# {{{1 Imports

from itertools import product
import pandas as pd

# {{{1 Configuration

# {{{2 Params

# Default params
max_threads = 24

# {{{2 Project configuration

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

# {{{3 Local includes
include: 'local.snake'

# {{{3 Sub-project includes
include: 'ormerod.snake'

# {{{1 Utility rules/recipes/templates

rule print_config:
    shell:
        '{config}'

localrules: print_config

rule drop_header:
    output: '{stem}.noheader.{ext,(tsv|csv)}'
    input: '{stem}.{ext}'
    shell: 'sed 1,1d {input} > {output}'

rule drop_header_meta:
    output: 'res/{stem}.noheader.{ext,(tsv|csv)}'
    input: 'meta/{stem}.{ext}'
    shell: 'sed 1,1d {input} > {output}'

# Here we have a template for aliasing
alias_recipe = "ln -rs {input} {output}"
alias_fmt = lambda input, output: alias_recipe.format(input=input, output=output)
curl_recipe = "curl '{params.url}' > {output}"
curl_unzip_recipe = "curl '{params.url}' | zcat > {output}"

# {{{1 Downloading and linking data

# {{{2 Reference data

rule download_salask_reference:
    output: 'raw/ref/salask.fn'
    params:
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000356.1&rettype=fasta&retmode=text'
    shell: curl_recipe

rule download_illumina_adapters:
    output: 'raw/ref/illumina_adapters.fn'
    params:
        url='https://raw.githubusercontent.com/vsbuffalo/scythe/master/illumina_adapters.fa'
    shell: curl_recipe

rule download_mouse_reference:
    output: 'raw/ref/mouse.fn'
    params:
        url='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCA_001632575.1_C3H_HeJ_v1/GCA_001632575.1_C3H_HeJ_v1_genomic.fna.gz'
    shell: curl_unzip_recipe

localrules: download_salask_reference, download_illumina_adapters, download_mouse_reference

rule download_sra_data:
    output: 'raw/sra/{sra_id}.fn'
    shell:
        """
        fastq-dump -Z {wildcards.sra_id} | seqtk seq -A > {output}
        """

# {{{2 Raw data
rule alias_raw_reads:
    output:
        r1='seq/{library}.m.r1.fq.gz',
        r2='seq/{library}.m.r2.fq.gz'
    input:
        unpack(lambda wildcards: {'r1': 'raw/mgen/' + config['library'][wildcards.library]['r1'],
                                  'r2': 'raw/mgen/' + config['library'][wildcards.library]['r2']})
    shell:
        alias_fmt("{input.r1}", "{output.r1}") +
        alias_fmt("{input.r2}", "{output.r2}")

# See scripts/query_longev_rrs_relative_abundance.sql
rule link_rabund_info:
    output: 'res/core.r.rabund.tsv'
    input: 'raw/longev_rrs_relative_abundance.tsv'
    shell: alias_recipe

localrules: alias_raw_reads, link_rabund_info

# {{{1 Metagenomics

# {{{2 Data pre-processing

rule deduplicate_reads:
    output:
        r1='seq/{stem}.m.r1.dedup.fq.gz',
        r2='seq/{stem}.m.r2.dedup.fq.gz'
    input:
        script='scripts/fastuniq_wrapper.sh',
        r1='seq/{stem}.m.r1.fq.gz',
        r2='seq/{stem}.m.r2.fq.gz'
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

rule quality_trim_reads:
    output:
        r1='seq/{stem}.r1.{proc}.qtrim.fq.gz',
        r2='seq/{stem}.r2.{proc}.qtrim.fq.gz',
        r3='seq/{stem}.r3.{proc}.qtrim.fq.gz'
    input:
        r1='seq/{stem}.r1.{proc}.fq.gz',
        r2='seq/{stem}.r2.{proc}.fq.gz'
    params:
        qual_type='sanger',
        qual_thresh=20
    shell:
        """
        sickle pe -t {params.qual_type} -q {params.qual_thresh} --gzip-output \
            --pe-file1 {input.r1} --pe-file2 {input.r2} \
            --output-pe1 {output.r1} --output-pe2 {output.r2} \
            --output-single {output.r3}
        """

# {{{3 Aliasing after processing

# This processing intended for mapping.
rule alias_read_processing:
    output: 'seq/{library_id}.m.{r}.proc.fq.gz'
    input: 'seq/{library_id}.m.{r}.dedup.deadapt.qtrim.fq.gz'
    shell: alias_recipe

localrules: alias_read_processing

# {{{2 Assembly

# {{{3 MEGAHIT

rule assemble_mgen:
    output:
        contigs='seq/{group,[^.]+}.a.{proc}.contigs.fn',
        outdir=temp('seq/{group}.a.{proc}.megahit.d'),
    input:
        lambda wildcards: [f'seq/{library}.m.{read}.{wildcards.proc}.fq.gz'
                           for library, read
                           in product(config['asmbl_group'][wildcards.group],
                                      ['r1', 'r2'])
                          ]
    log: 'log/{group}.a.{proc}.log'
    threads: max_threads
    params:
        r1=lambda wildcards: ','.join([f'seq/{library}.m.r1.{wildcards.proc}.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
        r2=lambda wildcards: ','.join([f'seq/{library}.m.r2.{wildcards.proc}.fq.gz'
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
        """

rule fragment_contigs:
    output: 'seq/{stem}.scontigs.fn'
    input: seqs='seq/{stem}.contigs.fn', script='scripts/fragment_seqs.py'
    params:
        length=10000, overlap=1000
    shell:
        """
        {input.script} {params.length} {params.overlap} {input.seqs} > {output}
        """

localrules: fragment_contigs

# {{{3 QC Assembly

rule quality_asses_assembly_with_spike:
    output: 'res/{stem}.{contigs,s?contigs}.metaquast.d'
    input: contigs='seq/{stem}.{contigs}.fn', ref='raw/ref/salask.fn'
    threads: max_threads
    params: min_contig_length=1000,
    shell:
        """
        metaquast.py --threads={threads} --min-contig {params.min_contig_length} -R {input.ref} --output-dir {output} {input.contigs}
        """

rule quality_asses_assembly:
    output: 'res/{stem}.{contigs,s?contigs}.metaquast2.d'
    input: contigs='seq/{stem}.{contigs}.fn'
    threads: max_threads
    params: min_contig_length=1000,
    shell:
        """
        metaquast.py --threads={threads} --min-contig {params.min_contig_length} --output-dir {output} {input.contigs}
        """

# {{{2 Mapping

# {{{3 Index building

rule bowtie_index_build:
    output:
        'seq/{stem}.1.bt2',
        'seq/{stem}.2.bt2',
        'seq/{stem}.3.bt2',
        'seq/{stem}.4.bt2',
        'seq/{stem}.rev.1.bt2',
        'seq/{stem}.rev.2.bt2'
    input: 'seq/{stem}.fn'
    log: 'log/{stem}.bowtie2-build.log'
    threads: max_threads
    shell:
        """
        bowtie2-build --threads {threads} {input} seq/{wildcards.stem} 2>&1 >{log}
        """

# {{{3 Backmap to an assembly

# Backmap with assembly from reads processed identically
rule backmap_reads_to_assembly:
    output: 'res/{library}.m.{proc}.{group}-map.sort.bam'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        r1='seq/{library}.m.r1.{proc}.fq.gz',
        r2='seq/{library}.m.r2.{proc}.fq.gz',
        inx_1='seq/{group}.a.{proc}.contigs.1.bt2',
        inx_2='seq/{group}.a.{proc}.contigs.2.bt2',
        inx_3='seq/{group}.a.{proc}.contigs.3.bt2',
        inx_4='seq/{group}.a.{proc}.contigs.4.bt2',
        inx_rev1='seq/{group}.a.{proc}.contigs.rev.1.bt2',
        inx_rev2='seq/{group}.a.{proc}.contigs.rev.2.bt2'
    threads: max_threads
    shell:
        r"""
        bowtie2 --threads {threads} \
                -x seq/{wildcards.group}.a.{wildcards.proc}.contigs \
                -1 {input.r1} -2 {input.r2} \
            | samtools sort --output-fmt=BAM -o {output}
        """

# Backmap with assembly from reads processed identically
rule backmap_reads_to_fragmented_assembly:
    output: 'res/{library}.m.{proc}.{group}-map.frag.sort.bam'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        r1='seq/{library}.m.r1.{proc}.fq.gz',
        r2='seq/{library}.m.r2.{proc}.fq.gz',
        inx_1='seq/{group}.a.{proc}.scontigs.1.bt2',
        inx_2='seq/{group}.a.{proc}.scontigs.2.bt2',
        inx_3='seq/{group}.a.{proc}.scontigs.3.bt2',
        inx_4='seq/{group}.a.{proc}.scontigs.4.bt2',
        inx_rev1='seq/{group}.a.{proc}.scontigs.rev.1.bt2',
        inx_rev2='seq/{group}.a.{proc}.scontigs.rev.2.bt2'
    threads: max_threads
    shell:
        r"""
        bowtie2 --threads {threads} \
                -x seq/{wildcards.group}.a.{wildcards.proc}.scontigs \
                -1 {input.r1} -2 {input.r2} \
            | samtools sort --output-fmt=BAM -o {output}
        """

# {{{2 Scaffolding


rule combine_bams:
    output: 'res/{group}.a.{proc}.{group}-map.sort.bam'
    input: lambda wildcards: [f'res/{library}.m.{wildcards.proc}.{wildcards.group}-map.sort.bam' for library in config['asmbl_group'][wildcards.group]]
    threads: 8
    shell:
        """
        tmpfile=$(mktemp -p $TMPDIR)
        samtools merge -f --threads {threads} $tmpfile {input}
        mv $tmpfile {output}
        """

rule tally_links:
    output: 'res/{stem}.linkage_tally.tsv'
    input: 'res/{stem}.bam'
    threads: 4
    params: min_hits=3, min_quality=40
    shell:
        r"""
        printf 'contig_id_A\tcontig_id_B\ttally\n' > {output}
        # For Flags explained see: https://broadinstitute.github.io/picard/explain-flags.html
        samtools view --threads={threads} -F3852 {input} | awk '($7 != "=") && ($7 != "*")' \
            | awk '{{print $3, $7, $1, $5}}' \
            | awk -v OFS='\t' \
                  '{{if ($1 > $2) {{print $3, "LEFT", $2, $1, $4}} \
                     else {{print $3, "RIGHT", $1, $2, $4}} \
                   }}' \
            | sort \
            | paste - - \
            | awk -v OFS='\t' -v q={params.min_quality} \
                  '($5 >= q) && ($10 >= q) {{print $3, $4}}' \
            | sort | uniq -c \
            | awk -v OFS='\t' '$1 >= {params.min_hits} {{print $2, $3, $1}}' \
            >> {output}
        """


# {{{2 Calculate statistics

# {{{3 Sequence lengths

rule count_seq_lengths_nucl:
    output: 'res/{stem}.nlength.tsv'
    input: script='scripts/count_seq_lengths.py', seqs='seq/{stem}.fn'
    shell:
        r"""
        printf "contig_id\tlength\n" > {output}
        {input.script} {input.seqs} >> {output}
        """

rule count_seq_lengths_aa:
    output: 'res/{stem}.alength.tsv'
    input: script='scripts/count_seq_lengths.py', seqs='seq/{stem}.fa'
    shell:
        """
        printf "seq_id\tlength\n" > {output}
        {input.script} {input.seqs} >> {output}
        """

localrules: count_seq_lengths_nucl, count_seq_lengths_aa

# {{{3 Coverage

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
    output: 'res/{library}.m.{proc}.{group}-map.cvrg.tsv'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        script='scripts/estimate_contig_coverage.py',
        depth='res/{library}.m.{proc}.{group}-map.depth.tsv',
        length='res/{group}.a.{proc}.contigs.nlength.tsv'
    params:
        float_fmt='%.6g'
    shell:
        """
        {input.script} {input.depth} {input.length} {params.float_fmt} > {output}
        """

rule estimate_fragmented_contig_cvrg:
    output: 'res/{library}.m.{proc}.{group}-map.frag.cvrg.tsv'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        script='scripts/estimate_contig_coverage.py',
        depth='res/{library}.m.{proc}.{group}-map.frag.depth.tsv',
        length='res/{group}.a.{proc}.scontigs.nlength.tsv'
    params:
        float_fmt='%.6g'
    shell:
        """
        {input.script} {input.depth} {input.length} {params.float_fmt} > {output}
        """

rule combine_cvrg:
    output: 'res/{group}.a.{proc}.contigs.cvrg.tsv'
    wildcard_constraints:
        group='[^.]+'
    input:
        script='scripts/concat_tables.py',
        tables=lambda wildcards: [f'res/{library}.m.{wildcards.proc}.{wildcards.group}-map.cvrg.tsv'
                                  for library
                                  in config['asmbl_group'][wildcards.group]
                                 ]
    shell:
        r"""
        printf 'library_id\tcontig_id\tcoverage\n' > {output}
        for file in {input.tables}; do
            awk -v OFS='\t' -v library_id=$(basename ${{file%%.*}}) 'FNR > 1 {{print library_id, $0}}' $file >> {output}
        done
        """

rule combine_fragmented_cvrg:
    output: 'res/{group}.a.{proc}.scontigs.cvrg.tsv'
    wildcard_constraints:
        group='[^.]+'
    input:
        script='scripts/concat_tables.py',
        tables=lambda wildcards: [f'res/{library}.m.{wildcards.proc}.{wildcards.group}-map.frag.cvrg.tsv'
                                  for library
                                  in config['asmbl_group'][wildcards.group]
                                 ]
    shell:
        r"""
        printf 'library_id\tcontig_id\tcoverage\n' > {output}
        for file in {input.tables}; do
            awk -v OFS='\t' -v library_id=$(basename ${{file%%.*}}) 'FNR > 1 {{print library_id, $0}}' $file >> {output}
        done
        """

rule unstack_cvrg:
    output: 'res/{stem}.{contigs,s?contigs}.cvrg.unstack.tsv'
    input:
        script='scripts/unstack_cvrg.py',
        cvrg='res/{stem}.{contigs}.cvrg.tsv',
    params:
        float_format='%.6g'
    shell:
        """
        {input.script} {input.cvrg} '{params.float_format}' > {output}
        """

# {{{2 Binning

# {{{3 Prepare input data

rule transform_contig_space:
    output:
        out='res/{stem}.{contigs}.pca.csv',
        dir=temp('res/{stem}.{contigs}.concoct.d')
    wildcard_constraints:
        contigs='s?contigs'
    input:
        cvrg='res/{stem}.{contigs}.cvrg.unstack.tsv',
        seqs='seq/{stem}.{contigs}.fn'
    threads: max_threads
    params:
        length_threshold=1000
    shadow: 'full'
    shell:
        r"""
        concoct --coverage_file={input.cvrg} --composition_file={input.seqs} \
                --length_threshold={params.length_threshold} \
                --basename={output.dir}/ \
                --cluster=10 --iterations=1 --num-threads={threads} --epsilon=1 --converge_out
        mv {output.dir}/PCA_transformed_data_gt{params.length_threshold}.csv {output.out}
        """

# {{{3 Clustering

rule bin_contigs:
    output:
        out='res/{stem}.{contigs,s?contigs}.cluster.tsv',
        summary='res/{stem}.{contigs}.cluster.summary.tsv',
    input:
        script='scripts/cluster_contigs.py',
        pca='res/{stem}.{contigs}.pca.csv',
        length='res/{stem}.{contigs}.nlength.tsv',
    params:
        min_total_length=100000,
        min_contig_length=1000,
        frac=0.10,
    shell:
        r"""
        {input.script} {input.pca} {input.length} \
                --min-length {params.min_contig_length} \
                --min-bin-size {params.min_total_length} \
                --summary {output.summary} \
                --frac {params.frac} \
                > {output.out}
        """

rule rename_clusters_to_bins:
    output: 'res/{stem}.{contigs,s?contigs}.bins.tsv'
    input: 'res/{stem}.{contigs}.cluster.tsv'
    shell:
        r"""
        awk -v OFS='\t' \
            'BEGIN   {{print "contig_id", "bin_id"}}
             FNR > 1 {{printf "%s\tbin%03d\n", $1, $2}}
            ' \
            < {input} \
            > {output}
        """

rule split_out_bins:
    output: 'seq/{stem}.{bins,m?bins}.d'
    input:
        script='scripts/fetch_bin.sh',
        bins='res/{stem}.{bins}.tsv',
        contigs='seq/{stem}.fn',
    threads: max_threads
    shell:
        r"""
        rm -rf {output}
        mkdir {output}
        bins=$(sed '1,1d' {input.bins} | cut -f 2 | sort | uniq)
        parallel --progress --jobs {threads} {input.script} {input.contigs} {input.bins} {{1}} {output}/{{1}}.fn ::: $bins
        """

# TODO: refine bins (scaffolds, contig extension, splitting/merging

# {{{3 QC bins

rule checkm_bins:
    output:
        outdir=temp('res/{stem}.{bins}.checkm.d'),
        summary='res/{stem}.{bins}.checkm.tsv'
    wildcard_constraints: bins='m?bins'
    input: 'seq/{stem}.{bins}.d'
    threads: max_threads
    shell:
        r"""
        rm -rf {output}
        checkm lineage_wf -x fn \
                --threads {threads} --pplacer_threads {threads} \
                --file {output.summary} --tab_table \
                {input} {output.outdir}
        """

rule reformat_checkm_output:
    output: 'res/{stem}.checkm_details.tsv'
    input: 'res/{stem}.checkm.tsv'
    shell:
        """
        cut -f1,4,12-15 {input} > {output}
        """

rule generate_checkm_markerset:
    output:
        'res/{level}_{taxon}.ms'
    shell:
        'checkm taxon_set {wildcards.level} {wildcards.taxon} {output}'

# {{{3 Refine bins
rule checkm_content_merge:
    output:
        checkm_work=temp('res/{stem}.bins.checkm_merge.d'),
        merge_stats='res/{stem}.bins.checkm_merge_stats.tsv',
    input:
        bins='seq/{stem}.bins.d',
        markerset='res/domain_Bacteria.ms',
    threads: max_threads
    shell:
        """
        rm -rf {output.checkm_work}
        checkm merge --threads {threads} -x fn {input.markerset} {input.bins} {output.checkm_work}
        head -1 {output.checkm_work}/merger.tsv > {output.result}
        printf 'bin_id_1\tbin_id_2\tscore\n' > {output.result}
        sed '1,1d' {output.checkm_work}/merger.tsv | cut -f1,2,9 >> {output.merge_stats}
        """

rule checkm_merge_bins:
    output: 'res/{stem}.mbins.tsv'
    input:
        script='scripts/merge_bins_by_checkm.py',
        bins='res/{stem}.bins.tsv',
        merge_list='res/{stem}.bins.checkm_merge_stats.tsv',
    params:
        min_complete_delta=30,
        max_contam_delta=2,
    shell:
        """
        {input.script} {input.bins} {input.merge_list} {params.min_complete_delta} {params.max_contam_delta} > {output}
        """


# {{{2 Annotation

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
    output: 'res/{stem}.{bin_id}.ec.list'
    input: 'res/{stem}.bins.prokka.d/{bin_id}.d'
    shell:
        """
        grep eC_number {input}/prokka.gff  | cut -f9 | cut -d';' -f2 | cut -d'=' -f2 | awk '{{print "ec:"$1}}' > {output}
        """


# {{{1 Compile all data

rule generate_database:
    output: 'res/{group}.results.db'
    input:
        schema='schema.sql',
        library='res/library.noheader.tsv',
        contig='res/{group}.a.proc.contigs.nlength.noheader.tsv',
        contig_bin='res/{group}.a.proc.contigs.bins.noheader.tsv',
        contig_coverage='res/{group}.a.proc.contigs.cvrg.noheader.tsv',
        bin_checkm='res/{group}.a.proc.contigs.bins.checkm_details.noheader.tsv',
        rrs_taxon_rabund='res/{group}.r.rabund.noheader.tsv',
        contig_linkage='res/{group}.a.proc.core-map.sort.linkage_tally.noheader.tsv',
        checkm_merge='res/{group}.a.proc.contigs.bins.checkm_merge_stats.noheader.tsv',
    shell:
        r"""
        rm -f {output}
        echo '
.bail ON
.read {input.schema}
PRAGMA cache_size = 1000000;
PRAGMA foreign_keys = TRUE;
.separator \t
.import {input.library} library
.import {input.checkm_merge} _bin_complementarity
.import {input.contig} contig
.import {input.contig_linkage} _contig_linkage
.import {input.contig_bin} contig_bin
.import {input.contig_coverage} contig_coverage
.import {input.bin_checkm} bin_checkm
.import {input.rrs_taxon_rabund} rrs_taxon_rabund
ANALYZE;
             ' \
        | sqlite3 {output}
        """
