# {{{1 Imports
from itertools import product
import pandas as pd

# {{{1 Params
# Default params
max_threads = 30

# {{{1 Project configuration
# Configure the pipeline
config_file = 'config.yaml'
configfile: config_file
# Metadata specific configurations
_library = pd.read_table(config['_meta_library'], index_col='library_id')
_asmbl_group = pd.read_table(config['_meta_asmbl_group'])
config['library'] = {}
for library_id, row in _library.iterrows():
    config['library'][library_id] = {}
    config['library'][library_id]['r1'] = row['file_r1']
    config['library'][library_id]['r2'] = row['file_r2']
config['asmbl_group'] = {}
for group, d in _asmbl_group.groupby('asmbl_group'):
    config['asmbl_group'][group] = list(d['library_id'])

# {{{1 Universal submission configuration
localrules: print_config, link_raw_reads, download_salask_reference,
    download_illumina_adapters, download_mouse_reference, alias_default_read_processing,

# {{{1 Utility rules
rule print_config:
    shell:
        '{config}'

# {{{1 Downloading and linking raw data
rule link_raw_reads:
    output:
        r1='seq/{library}.mgen.r1.fq.gz',
        r2='seq/{library}.mgen.r2.fq.gz'
    input:
        unpack(lambda wildcards: {'r1': 'raw/mgen/' + config['library'][wildcards.library]['r1'],
                                  'r2': 'raw/mgen/' + config['library'][wildcards.library]['r2']})
    shell:
        """
        ln -rs {input.r1} {output.r1}
        ln -rs {input.r2} {output.r2}
        """

rule download_salask_reference:
    output: 'raw/ref/salask.fn'
    params:
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000356.1&rettype=fasta&retmode=text'
    shell: "curl '{params.url}' > {output}"

rule download_illumina_adapters:
    output: 'raw/ref/illumina_adapters.fn'
    params:
        url='https://raw.githubusercontent.com/vsbuffalo/scythe/master/illumina_adapters.fa'
    shell: "curl '{params.url}' > {output}"

rule download_mouse_reference:
    output: 'raw/ref/mouse.fn'
    params:
        url='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCA_001632575.1_C3H_HeJ_v1/GCA_001632575.1_C3H_HeJ_v1_genomic.fna.gz'
    shell: "curl '{params.url}' | zcat > {output}"

rule download_sra_data:
    output: 'raw/sra/{sra_id}.fn'
    shell:
        """
        fastq-dump -Z {wildcards.sra_id} | seqtk seq -A > {output}
        """

# {{{1 Data pre-processing
rule deduplicate_reads:
    output:
        r1='seq/{stem}.mgen.r1.dedup.fq.gz',
        r2='seq/{stem}.mgen.r2.dedup.fq.gz'
    input:
        script='scripts/fastuniq_wrapper.sh',
        r1='seq/{stem}.mgen.r1.fq.gz',
        r2='seq/{stem}.mgen.r2.fq.gz'
    resources:
        mem_mb=24000
    shell: "{input.script} {input.r1} {input.r2} {output.r1} {output.r2}"

rule trim_adapters:
    output: 'seq/{stem}.deadapt.fq.gz'
    input:
        adapters='raw/ref/illumina_adapters.fn',
        seqs='seq/{stem}.fq.gz'
    log: 'log/{stem}.scythe.log'
    shell:
        """
        scythe -a {input.adapters} {input.seqs} > {output} 2>{log}
        ! grep -Fxq 'Blank FASTA header or sequence in adapters file.' {log}
        """

rule quality_trim:
    output:
        r1='seq/{stem}.r1.{proc}.qtrim.fq.gz',
        r2='seq/{stem}.r2.{proc}.qtrim.fq.gz',
        r3='seq/{stem}.r3.{proc}.qtrim.fq.gz'
    input:
        r1='seq/{stem}.r1.{proc}.fq.gz',
        r2='seq/{stem}.r2.{proc}.fq.gz'
    params:
        ends='pe',
        qual_type='sanger'
    shell:
        """
        sickle {params.ends} -t {params.qual_type} --gzip-output \
            --pe-file1 {input.r1} --pe-file2 {input.r2} \
            --output-pe1 {output.r1} --output-pe2 {output.r2} \
            --output-single {output.r3}
        """

# Alias assembly (drop line-noise)
rule alias_default_read_processing:
    output: 'seq/{library_id}.mgen.{r}.proc.fq.gz'
    input: 'seq/{library_id}.mgen.{r}.dedup.deadapt.qtrim.fq.gz'
    shell: 'ln -rs {input} {output}'


# {{{1 Assembly
rule assemble_mgen:
    output:
        contigs='seq/{group,[^.]+}.asmbl.{proc}.contigs.fn',
        outdir=temp('seq/{group}.asmbl.{proc}.megahit.d'),
        lengths='res/{group}.asmbl.{proc}.contigs.length.tsv'
    input:
        lambda wildcards: [f'seq/{library}.mgen.{read}.{wildcards.proc}.fq.gz'
                           for library, read
                           in product(config['asmbl_group'][wildcards.group],
                                      ['r1', 'r2', 'r3'])
                          ]
    log: 'log/{group}.asmbl.{proc}.log'
    threads: max_threads
    params:
        r1=lambda wildcards: ','.join([f'seq/{library}.mgen.r1.{wildcards.proc}.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
        r2=lambda wildcards: ','.join([f'seq/{library}.mgen.r2.{wildcards.proc}.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
    shell:
        r"""
        megahit \
            -1 {params.r1} \
            -2 {params.r2} \
            --k-min 21 --k-max 161 --k-step 20 \
            --out-dir {output.outdir} \
            --num-cpu-threads {threads} \
            --verbose
        sed 's:^>k:>{wildcards.group}-k:' {output.outdir}/final.contigs.fa > {output.contigs}
        cp {output.outdir}/log {log}
        printf 'contig_id\tlength\n' > {output.lengths}
        grep '^>' {output.contigs} \
            | sed 's:^>\(\S*\) flag=\([0123]\) multi=\(\S*\) len=\([0-9]*\)$:\\1\t\\4:' \
            >> {output.lengths}
        """

# rule quality_asses_assembly
# rule fragment_contigs  # To mitigate missassembly

# {{{1 Mapping
rule bowtie_index_build:
    output:
        'seq/{stem}.1.bt2',
        'seq/{stem}.2.bt2',
        'seq/{stem}.3.bt2',
        'seq/{stem}.4.bt2',
        'seq/{stem}.rev.1.bt2',
        'seq/{stem}.rev.2.bt2'
    input: 'seq/{stem}.fn'
    log: 'seq/{stem}.bowtie2-build.log'
    threads: max_threads
    shell:
        """
        bowtie2-build --threads {threads} {input} seq/{wildcards.stem} 2>&1 >{log}
        """

rule backmap_reads_to_assembly:
    output: 'res/{library}.mgen.{proc}.{group}-map.sort.bam'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        r1='seq/{library}.mgen.r1.{proc}.fq.gz',
        r2='seq/{library}.mgen.r2.{proc}.fq.gz',
        inx_1='seq/{group}.asmbl.{proc}.contigs.1.bt2',
        inx_2='seq/{group}.asmbl.{proc}.contigs.2.bt2',
        inx_3='seq/{group}.asmbl.{proc}.contigs.3.bt2',
        inx_4='seq/{group}.asmbl.{proc}.contigs.4.bt2',
        inx_rev1='seq/{group}.asmbl.{proc}.contigs.rev.1.bt2',
        inx_rev2='seq/{group}.asmbl.{proc}.contigs.rev.2.bt2'
    threads: max_threads
    shell:
        r"""
        bowtie2 --threads {threads} -x seq/{wildcards.group}.asmbl.{wildcards.proc}.contigs -1 {input.r1} -2 {input.r2} \
            | samtools sort --output-fmt=BAM -o {output}
        """

# # TODO: Drop this rule aand ruleorder fter all my BAMs are sorted.
# rule sort_bam:
#     output: 'res/{stem}.sort.bam'
#     input: 'res/{stem}.bam'
#     shadow: 'full'
#     shell:
#         """
#         samtools sort --output-fmt=BAM -o {output} {input}
#         """
# ruleorder: sort_bam > backmap_reads_to_assembly


# {{{1 Transform mapping into tabular data
# NOTE: The depth file format is lacking a header.
rule calculate_mapping_depth:
    output: temp('res/{stem}.depth.tsv')
    input: 'res/{stem}.sort.bam'
    shadow: 'full'
    shell:
        """
        samtools depth {input} >> {output}
        """

rule estimate_contig_cvrg:
    output: 'res/{library}.mgen.{proc}.{group}-map.cvrg.tsv'
    input:
        script='scripts/estimate_contig_coverage.py',
        depth='res/{library}.mgen.{proc}.{group}-map.depth.tsv',
        length='res/{group}.asmbl.{proc}.contigs.length.tsv'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    params:
        float_fmt='%.6g'
    shell:
        """
        {input.script} {input.depth} {input.length} {params.float_fmt} > {output}
        """

rule combine_cvrg:
    output: 'res/{group}.asmbl.{proc}.contigs.cvrg.tsv'
    input:
        script='scripts/concat_tables.py',
        tables=lambda wildcards: [f'res/{library}.mgen.{wildcards.proc}.{wildcards.group}-map.cvrg.tsv'
                                  for library
                                  in config['asmbl_group'][wildcards.group]
                                 ]
    wildcard_constraints:
        group='[^.]+'
    shell:
        r"""
        printf 'library_id\tcontig_id\tcoverage\n' > {output}
        for file in {input.tables}; do
            awk -v OFS='\t' -v library_id=$(basename ${{file%%.*}}) 'NR != 1 {{print library_id, $0}}' $file >> {output}
        done
        """

rule unstack_cvrg:
    output: 'res/{stem}.contigs.cvrg.unstack.tsv'
    input:
        script='scripts/unstack_cvrg.py',
        cvrg='res/{stem}.contigs.cvrg.tsv',
    params:
        float_format='%.6g'
    shell:
        """
        {input.script} {input.cvrg} '{params.float_format}' > {output}
        """

# {{{1 Binning
rule transform_contig_space:
    output:
        out='res/{stem}.contigs.pca.csv',
        dir=temp('res/{stem}.contigs.concoct.d')
    input:
        cvrg='res/{stem}.contigs.cvrg.unstack.tsv',
        seqs='seq/{stem}.contigs.fn'
    threads: max_threads
    params:
        length_threshold=1000
    shadow: 'full'
    shell:
        r"""
        concoct --coverage_file={input.cvrg} --composition_file={input.seqs} \
                --length_threshold={params.length_threshold} \
                --basename {output.dir}/ --iterations=1 \
                --epsilon=100 --converge_out
        mv {output.dir}/PCA_transformed_data_gt{params.length_threshold}.csv {output.out}
        """

# rule summarize_contig_coverage
# rule bin_contigs_metabat

rule bin_contigs_gmm:
    output:
        out='res/{stem}.contigs.bins.tsv',
        summary='res/{stem}.contigs.bins.summary.tsv'
    input:
        script='scripts/cluster_contigs.py',
        pca='res/{stem}.contigs.pca.csv',
        length='res/{stem}.contigs.length.tsv'
    params:
        method='gmm',
        min_total_length=100000,
        min_contig_length=1000
    shell:
        r"""
        {input.script} {input.pca} {input.length} \
                --min-length {params.min_contig_length} \
                --min-bin-size {params.min_total_length} \
                --summary {output.summary} \
                {params.method} \
                > {output.out}
        """

rule split_out_bins:
    output: 'seq/{stem}.bins.d'
    input: bins='res/{stem}.bins.tsv', contigs='seq/{stem}.fn'
    shell:
        r"""
        rm -rf {output}
        mkdir {output}
        for num in $(sed '1,1d' {input.bins} | cut -f 2 | sort | uniq); do
            outfile=$(printf 'bin_%04d.fn' $num)
            seqtk subseq {input.contigs} \
                         <(awk -v bin=$num '$2==bin' {input.bins} | cut -f1) \
                > {output}/$outfile
        done
        """

rule checkm_bins:
    output:
        outdir=temp('res/{stem}.bins.checkm.d'),
        summary='res/{stem}.bins.checkm.tsv'
    input: 'seq/{stem}.bins.d'
    threads: max_threads
    shell:
        r"""
        rm -rf {output}
        checkm lineage_wf -x fn \
                --threads {threads} --pplacer_threads {threads} \
                --file {output.summary} --tab_table \
                {input} {output.outdir}
        """

rule generate_checkm_markerset:
    output:
        'res/{level}_{taxon}.ms'
    shell:
        'checkm taxon_set {wildcards.level} {wildcards.taxon} {output}'

rule checkm_merge:
    output:
        checkm_work=temp('res/{stem}.bins.checkm_merge.d'),
        result='res/{stem}.bins.merge.tsv'
    input:
        bins='seq/{stem}.bins.d',
        markerset='res/kingdom_Bacteria.ms'
    threads: max_threads
    shell:
        """
        rm -rf {output.checkm_work}
        checkm merge --threads {threads} -x fn {input.markerset} {input.bins} {output.checkm_work}
        head -1 {output.checkm_work}/merger.tsv > {output.result}
        sed '1,1d' {output.checkm_work}/merger.tsv | sort -k9,9rn >> {output.result}
        """

# rule refine_bins:

# {{{1 Annotation
rule annotate_bin:
    output: 'res/{stem}.bins.prokka.d/{bin_id}.d'
    input: 'seq/{stem}.bins.d/{bin_id}.fn'
    log: 'log/{bin_id}.prokka.log'
    threads: max_threads
    shell:
        r"""
        prokka --force --cpus {threads} {input} \
                --outdir {output} --prefix prokka
                --locustag {wildcards.bin_id} --rawproduct \
                >{log} 2>&1
        """

rule extract_ec_numbers:
    output: 'res/{bin}.ec.list'
    input: 'res/{bin}.prokka.d'
    shell:
        """
        grep eC_number {input}/prokka.gff  | cut -f9 | cut -d';' -f2 | cut -d'=' -f2 | awk '{{print "ec:"$1}}' > {output}
        """

# {{{1 Compare to Ormerod2016
# {{{2 Metadata

config['ormerod'] = \
        {
        'GP1'                                         : 'DAAI01000000',
        'GP2'                                         : 'DAAJ01000000',
        'GP3'                                         : 'DAAK01000000',
        'GP4'                                         : 'DAAL01000000',
        'H1'                                          : 'DAAM01000000',
        'H10'                                         : 'DAAV01000000',
        'H2'                                          : 'DAAO01000000',
        'H3'                                          : 'DAAN01000000',
        'H4'                                          : 'DAAP01000000',
        'H5'                                          : 'DAAQ01000000',
        'H6'                                          : 'DAAR01000000',
        'H7'                                          : 'DAAS01000000',
        'H8'                                          : 'DAAT01000000',
        'H9'                                          : 'DAAU01000000',
        'K1'                                          : 'LUJZ01000000',
        'K10'                                         : 'LUKA01000000',
        'M1'                                          : 'LUJL01000000',
        'M10'                                         : 'LUJU01000000',
        'M11'                                         : 'LUJV01000000',
        'M12'                                         : 'LUJW01000000',
        'M13'                                         : 'LUJX01000000',
        'M14'                                         : 'LUJY01000000',
        'M2'                                          : 'LUJM01000000',
        'M3'                                          : 'LUJN01000000',
        'M5'                                          : 'LUJP01000000',
        'M6'                                          : 'LUJQ01000000',
        'M7'                                          : 'LUJR01000000',
        'M8'                                          : 'LUJS01000000',
        'M9'                                          : 'LUJT01000000',
        'Candidatus_Homeothermus_arabinoxylanisolvens': 'LUJO01000000',
        }

config['marker_gene'] = {}
config['marker_gene']['search_string'] = \
    {'rpoB': 'DNA-directed RNA polymerase subunit beta$',
     'gyrB': 'DNA gyrase subunit B$',
     'rrnS': '16S ribosomal RNA$'}

# }}}

rule link_ormerod_bin:
    output: 'seq/ormerod.bins.d/{bin_id}.fn'
    input: lambda wildcards: 'raw/sra/{accession}.fn'.format(accession=config['ormerod'][wildcards.bin_id])
    shell:
        """
        ln -rs {input} {output}
        """

rule process_ormerod_bins:
    input:
        expand('res/ormerod.bins.prokka.d/{bin_id}.d', bin_id=config['ormerod'].keys()),
        'res/core.asmbl.dedup.deadapt.qtrim.contigs.bins.prokka.d/bin_0090.d'

phylogenetic_marker_genes = []

rule pull_out_phylogenetic_marker_genes_nucl_ormerod:
    output: 'seq/{stem}.prokka.{gene_id}.fn',
    input:
        expand('res/{{stem}}.prokka.d/{bin_id}.d/prokka.ffn', bin_id=config['ormerod'].keys()),
        'res/core.asmbl.dedup.deadapt.qtrim.contigs.bins.prokka.d/bin_0090.d/prokka.ffn'
    params:
        search_string=lambda wildcards: config['marker_gene']['search_string'][wildcards.gene_id],
    shell:
        r"""
        rm -f {output}
        for file in {input}; do
            seqtk subseq $file <(grep '{params.search_string}' $file | sed 's:^>\([^ ]*\).*:\1:') >> {output}
        done
        """

rule pull_out_phylogenetic_marker_genes_amino_ormerod:
    output: 'seq/{stem}.prokka.{gene_id}.fa',
    input:
        expand('res/{{stem}}.prokka.d/{bin_id}.d/prokka.faa', bin_id=config['ormerod'].keys()),
        'res/core.asmbl.dedup.deadapt.qtrim.contigs.bins.prokka.d/bin_0090.d/prokka.faa'
    params:
        search_string=lambda wildcards: config['marker_gene']['search_string'][wildcards.gene_id],
    shell:
        r"""
        rm -f {output}
        for file in {input}; do
            seqtk subseq $file <(grep '{params.search_string}' $file | sed 's:^>\([^ ]*\).*:\1:') >> {output}
        done
        """
