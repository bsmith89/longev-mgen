# {{{1 Imports

from itertools import product
import pandas as pd

# {{{1 Configuration

# {{{2 Nomenclature

wildcard_constraints:
    group='[^.]+',
    library='[^.]+'

# {{{2 Params

# Default params
max_threads = 30


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
include: 'snake/local.snake'

# {{{3 Sub-project includes
include: 'snake/reference_genomes.snake'
include: 'snake/cazy.snake'

# {{{1 Utility rules/recipes/templates

rule all:
    shell:
        """
        echo "Figure out what you really want!"
        """

localrules: all

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

localrules: drop_header, drop_header_meta

ruleorder: drop_header_meta > drop_header

rule start_jupyter:
    shell: "jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/"

rule configure_git:
    shell:
        """
        # IPYNB Filter
        git config --local filter.dropoutput_ipynb.clean scripts/ipynb_output_filter.py
        git config --local filter.dropoutput_ipynb.smudge cat
        # Pager Config
        git config --local core.pager 'less -x4'
        # Daff config
        git config --local diff.daff-csv.command "daff.py diff --git"
        git config --local merge.daff-csv.name "daff.py tabular merge"
        git config --local merge.daff-csv.driver "daff.py merge --output %A %O %A %B"
        """

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

rule download_m_intestinale_genome:
    output: 'raw/ref/Muribaculum_intestinale_yl27.fn'
    params:
        url='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/201/515/GCA_002201515.1_ASM220151v1/GCA_002201515.1_ASM220151v1_genomic.fna.gz'
    shell: curl_unzip_recipe

rule download_tigrfam:
    output: "raw/ref/TIGRFAMs_14.0_HMM.tar.gz"
    params:
        url="ftp://ftp.jcvi.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz"
    shell:
        curl_recipe

rule extract_tigrfam_hmm:
    output: "ref/hmm/TIGR{num}.hmm"
    input: "raw/ref/TIGRFAMs_14.0_HMM.tar.gz"
    shell:
        """
        tar -xzf {input} TIGR{wildcards.num}.HMM -O > {output}
        """

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

rule download_sra_data:
    output: 'raw/sra/{sra_id}.fn'
    shell:
        """
        fastq-dump -Z {wildcards.sra_id} | seqtk seq -A > {output}
        """

rule download_cog_function_mapping:
    output: 'raw/ref/cognames2003-2014.tab'
    params:
        url="ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab"
    shell: curl_recipe

rule process_cog_function_mapping:
    output: 'ref/cog_function.tsv'
    input: 'raw/ref/cognames2003-2014.tab'
    shell: "iconv -f LATIN1 -t UTF-8 {input} | sed '1,1s:^# COG\tfunc\tname:cog_id\tfunction_categories\tfunction_name:' > {output}"

rule download_cog_to_ko_mapping:
    output: 'raw/ref/cog_from_string7_to_ko20080319_filtered_005.txt'
    params:
        url="http://pathways2.embl.de/data/cog_from_string7_to_ko20080319_filtered_005.txt.gz"
    shell: curl_unzip_recipe

rule alias_cog_to_ko_mapping:
    output: 'ref/cog_to_ko.tsv'
    input: 'raw/ref/cog_from_string7_to_ko20080319_filtered_005.txt'
    shell: alias_recipe


localrules: download_salask_reference, download_illumina_adapters,
            download_mouse_reference, download_sra_data, download_tigrfam,
            extract_tigrfam_hmm, download_cog_to_ko_mapping, alias_cog_to_ko_mapping


# {{{2 Raw data

rule alias_raw_read_r1:
    output: 'seq/{library}.m.r1.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r1'])
    shell: alias_recipe

rule alias_raw_read_r2:
    output: 'seq/{library}.m.r2.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r2'])
    shell: alias_recipe

localrules: alias_raw_read_r1, alias_raw_read_r2, link_rabund_info

# {{{3 Import results from 16S libraries

rule query_count_info:
    output: 'res/core.r.count.tsv'
    input: script='scripts/query_longev_rrs_count.sql', db='raw/longev_rrs_results.db'
    shell: "sqlite3 -header -separator '\t' {input.db} < {input.script} > {output}"

rule query_taxonomy_info:
    output: 'res/core.r.taxonomy.tsv'
    input: script='scripts/query_longev_rrs_taxonomy.sql', db='raw/longev_rrs_results.db'
    shell: "sqlite3 -header -separator '\t' {input.db} < {input.script} > {output}"

localrules: query_count_info, query_taxonomy_info

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
        scythe -a {input.adapters} {input.seqs} >{output} 2>{log}
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
        r"""
        sickle pe -t {params.qual_type} -q {params.qual_thresh} --gzip-output \
            --pe-file1 {input.r1} --pe-file2 {input.r2} \
            --output-pe1 {output.r1} --output-pe2 {output.r2} \
            --output-single {output.r3}
        """

# {{{3 Aliasing after processing

# This processing intended for mapping.
rule alias_read_processing_r1:
    output: 'seq/{library_id}.m.r1.proc.fq.gz'
    input: 'seq/{library_id}.m.r1.dedup.deadapt.qtrim.fq.gz'
    shell: alias_recipe

rule alias_read_processing_r2:
    output: 'seq/{library_id}.m.r2.proc.fq.gz'
    input: 'seq/{library_id}.m.r2.dedup.deadapt.qtrim.fq.gz'
    shell: alias_recipe

localrules: alias_read_processing_r1, alias_read_processing_r2

# {{{2 Assembly

# {{{3 MEGAHIT

rule assemble_mgen:
    output:
        fasta='seq/{group}.a.contigs.fn',
        outdir=temp('seq/{group}.a.megahit.d'),
        fastg='seq/{group}.a.contigs.fg',
    wildcard_constraints:
        group='[^.]+',
    input:
        lambda wildcards: [f'seq/{library}.m.{read}.proc.fq.gz'
                           for library, read
                           in product(config['asmbl_group'][wildcards.group],
                                      ['r1', 'r2'])
                          ]
    log: 'log/{group}.a.megahit.log'
    threads: max_threads
    params:
        r1=lambda wildcards: ','.join([f'seq/{library}.m.r1.proc.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
        r2=lambda wildcards: ','.join([f'seq/{library}.m.r2.proc.fq.gz'
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
        sed 's:^>k:>{wildcards.group}-k:' {output.outdir}/final.contigs.fa > {output.fasta}
        # TODO: Fix this hard-coding of k-parameters.
        megahit_toolkit contig2fastg 141 {output.outdir}/intermediate_contigs/k141.contigs.fa > {output.fastg}
        cp {output.outdir}/log {log}
        """


rule alias_mbin_contig_file_for_spades:
    output: 'seq/{group}.a.mbins.d/{mbin_id}.fasta'
    input: 'seq/{group}.a.mbins.d/{mbin_id}.fn'
    shell: alias_recipe

localrules: alias_mbin_contig_file_for_spades

rule assemble_mbin:
    output: 'seq/{group}.a.mbins.d/{mbin_id}.a.spades.d'
    input: reads='seq/{group}.a.mbins.d/{mbin_id}.m.pe.fq.gz', contigs='seq/{group}.a.mbins.d/{mbin_id}.fasta'
    threads: max_threads
    shell:
        """
        spades.py --threads {threads} --careful --12 {input.reads} --untrusted-contigs {input.contigs} -o {output}
        """

# {{{3 QC Assembly

rule quality_asses_assembly_with_spike:
    output: 'res/{group}.a.contigs.metaquast.d'
    wildcard_constraints:
        group='[^.]+',
    input: contigs='seq/{group}.a.contigs.fn', ref='raw/ref/salask.fn'
    threads: max_threads
    params: min_contig_length=1000,
    shell:
        """
        metaquast.py --threads={threads} --min-contig {params.min_contig_length} -R {input.ref} --output-dir {output} {input.contigs}
        """

# {{{2 Mapping

# {{{3 Index building

rule bowtie_index_build:
    output:
        'seq/{group}.a.contigs.1.bt2',
        'seq/{group}.a.contigs.2.bt2',
        'seq/{group}.a.contigs.3.bt2',
        'seq/{group}.a.contigs.4.bt2',
        'seq/{group}.a.contigs.rev.1.bt2',
        'seq/{group}.a.contigs.rev.2.bt2'
    wildcard_constraints:
        group='[^.]+',
    input: 'seq/{group}.a.contigs.fn'
    log: 'log/{group}.a.contigs.bowtie2-build.log'
    threads: min(max_threads, 14)
    shell:
        """
        bowtie2-build --threads {threads} {input} seq/{wildcards.group}.a.contigs >{log} 2>&1
        """

# {{{3 Backmap to an assembly

# Backmap with assembly from reads processed identically
rule map_reads_to_assembly:
    output: 'res/{library}.m.{group}-map.sort.bam'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        r1='seq/{library}.m.r1.proc.fq.gz',
        r2='seq/{library}.m.r2.proc.fq.gz',
        inx_1='seq/{group}.a.contigs.1.bt2',
        inx_2='seq/{group}.a.contigs.2.bt2',
        inx_3='seq/{group}.a.contigs.3.bt2',
        inx_4='seq/{group}.a.contigs.4.bt2',
        inx_rev1='seq/{group}.a.contigs.rev.1.bt2',
        inx_rev2='seq/{group}.a.contigs.rev.2.bt2'
    threads: max_threads
    shell:
        r"""
        bowtie2 --threads {threads} \
                -x seq/{wildcards.group}.a.contigs \
                -rg-id {wildcards.library} \
                -1 {input.r1} -2 {input.r2} \
            | samtools sort --output-fmt=BAM -o {output}
        """

rule combine_read_mappings:
    output: 'res/{group}.a.contigs.map.sort.bam'
    input:
        lambda wildcards: [f'res/{library}.m.{wildcards.group}-map.sort.bam'
                           for library
                           in config['asmbl_group'][wildcards.group]
                          ]
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output} {input}
        """

rule index_read_mappings:
    output: 'res/{stem}.sort.bam.bai'
    input: 'res/{stem}.sort.bam'
    shell: 'samtools index {input} {output}'


# {{{2 Scaffolding


rule tally_links:
    output: 'res/{library}.m.{group}-map.linkage_tally.tsv'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input: 'res/{library}.m.{group}-map.sort.bam'
    params: min_hits=1, min_quality=40
    shell:
        r"""
        printf 'contig_id_1\tcontig_id_2\tlibrary_id\tread_count\n' > {output}
        # For Flags explained see: https://broadinstitute.github.io/picard/explain-flags.html
        samtools view -F3852 {input} | awk '($7 != "=") && ($7 != "*")' \
            | awk '{{print $3, $7, $1, $5}}' \
            | awk -v OFS='\t' \
                  '{{if ($1 > $2) {{print $3, "LEFT", $2, $1, $4}} \
                     else {{print $3, "RIGHT", $1, $2, $4}} \
                   }}' \
            | sort \
            | paste - - \
            | awk -v OFS='\t' -v q={params.min_quality} \
                  '   ($5 >= q) \
                   && ($10 >= q) \
                   && ($2 == "LEFT") \
                   && ($7 == "RIGHT") \
                   {{print $3, $4}} \
                  ' \
            | sort | uniq -c \
            | awk -v OFS='\t' \
                  -v library_id='{wildcards.library}' \
                  -v min_hits='{params.min_hits}' \
                  '$1 >= min_hits {{print $2, $3, library_id, $1}}' \
            >> {output}
        """

rule combine_linkage_tallies:
    output: 'res/{group}.a.contigs.linkage_tally.tsv'
    input:
        lambda wildcards: [f'res/{library}.m.{wildcards.group}-map.linkage_tally.tsv'
                           for library in config['asmbl_group'][wildcards.group]
                          ]
    shell:
        r"""
        printf 'contig_id_1\tcontig_id_2\tlibrary_id\tread_count\n' > {output}
        for file in {input}; do
            sed 1,1d $file >> {output}
        done
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
    output: 'res/{library}.m.{group}-map.cvrg.tsv'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        script='scripts/estimate_contig_coverage.py',
        depth='res/{library}.m.{group}-map.depth.tsv',
        length='res/{group}.a.contigs.nlength.tsv'
    params:
        float_fmt='%.6g'
    shell:
        """
        {input.script} {input.depth} {input.length} {params.float_fmt} > {output}
        """

rule combine_cvrg:
    output: 'res/{group}.a.contigs.cvrg.tsv'
    wildcard_constraints:
        group='[^.]+'
    input:
        script='scripts/concat_tables.py',
        tables=lambda wildcards: [f'res/{library}.m.{wildcards.group}-map.cvrg.tsv'
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
    output: 'res/{group}.a.contigs.cvrg.unstack.tsv'
    wildcard_constraints:
        group='[^.]+'
    input:
        script='scripts/unstack_cvrg.py',
        cvrg='res/{group}.a.contigs.cvrg.tsv',
    params:
        float_format='%.6g'
    shell:
        """
        {input.script} {input.cvrg} '{params.float_format}' > {output}
        """

# {{{2 Binning

# {{{3 Prepare input data


# rule count_tetramers:
#     output: "res/{stem}.contigs.4mers.tsv"
#     input: script="scripts/count_kmers.py", seq="seq/{stem}.contigs.fn"
#     threads: 1
#     shell:
#         """
#         {input.script} 4 ACGT {input.seq} | sed '1,1s:seq_id:contig_id:' > {output}
#         """
#
# TODO: scripts/count_kmers.py
# TODO: scripts/transform_contig_space.py
rule transform_contig_space:
    output:
        pca='res/{group}.a.contigs.concoct.pca.tsv',
        raw='res/{group}.a.contigs.concoct.tsv',
        dir=temp('res/{group}.a.contigs.concoct.d')
    wildcard_constraints:
        group='[^.]+'
    input:
        cvrg='res/{group}.a.contigs.cvrg.unstack.tsv',
        seqs='seq/{group}.a.contigs.fn'
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
        sed 's:,:\t:g' {output.dir}/original_data_gt{params.length_threshold}.csv | sed '1,1s:^:contig_id:' > {output.raw}
        sed 's:,:\t:g' {output.dir}/PCA_transformed_data_gt{params.length_threshold}.csv > {output.pca}
        """

# {{{3 Clustering

# TODO: Don't rename as a separate stem.
# TODO: Refine this script.
rule cluster_contigs:
    output:
        out='res/{group}.a.contigs.cluster.tsv',
        summary='res/{group}.a.contigs.cluster.summary.tsv',
    wildcard_constraints:
        group='[^.]+'
    input:
        script='scripts/cluster_contigs.py',
        pca='res/{group}.a.contigs.concoct.pca.tsv',
        length='res/{group}.a.contigs.nlength.tsv',
    log: 'log/{group}.a.contigs.cluster.log'
    params:
        min_contig_length=1000,
        frac=0.10,
        alpha=1,
        max_clusters=2000,
        seed=1,
    shell:
        r"""
        {input.script} {input.pca} {input.length} \
                --min-length {params.min_contig_length} \
                --frac {params.frac} \
                --max-nbins {params.max_clusters} \
                --alpha {params.alpha} \
                --seed {params.seed} \
                --summary {output.summary} \
                > {output.out} 2> {log}
        """

rule rename_clusters_to_bins:
    output: 'res/{group}.a.contigs.bins.tsv'
    wildcard_constraints:
        group='[^.]+'
    input: 'res/{group}.a.contigs.cluster.tsv'
    params:
        padding=5
    shell:
        r"""
        awk -v OFS='\t' \
            'BEGIN   {{print "contig_id", "bin_id"}}
             FNR > 1 {{printf "%s\tbin%0{params.padding}d\n", $1, $2}}
            ' \
            < {input} \
            > {output}
        """

rule split_out_bins:
    output: 'seq/{group}.a.bins.d'
    input:
        script='scripts/fetch_bin.sh',
        bins='res/{group}.a.contigs.bins.tsv',
        contigs='seq/{group}.a.contigs.fn',
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

rule checkm_bins_or_mags:
    output:
        outdir=temp('res/{stem}.{bins_or_mbins}.checkm.d'),
        summary='res/{stem}.{bins_or_mbins}.checkm.tsv'
    input: 'seq/{stem}.{bins_or_mbins}.d'
    wildcard_constraints:
        bins_or_mbins = 'bins|mbins|mags'
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
    output: 'res/{group}.checkm_details.tsv'
    input: 'res/{group}.checkm.tsv'
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
# TODO: Understand what field 9 in checkM output file is.
# (I _think_ it's the difference between the increased completeness and the
# increased contamination.)
rule checkm_content_merge:
    output:
        checkm_work=temp('res/{group}.a.bins.checkm_merge.d'),
        merge_stats='res/{group}.a.bins.checkm_merge_stats.tsv',
    input:
        bins='seq/{group}.a.bins.d',
        markerset='res/domain_Bacteria.ms',
    threads: max_threads
    shadow: 'full'
    shell:
        """
        rm -rf {output.checkm_work}
        checkm merge --threads {threads} -x fn {input.markerset} {input.bins} {output.checkm_work}
        printf 'bin_id_1\tbin_id_2\tscore\n' > {output.merge_stats}
        sed '1,1d' {output.checkm_work}/merger.tsv | cut -f1,2,9 >> {output.merge_stats}
        """

rule query_merge_stats:
    output: 'res/{group}.a.bin_merge_stats.tsv'
    input: sql='scripts/query_bin_merger.sql', db='res/{group}.1.denorm.db'
    shell:
        """
        sqlite3 -header -separator '\t' {input.db} < {input.sql} > {output}
        """

localrules: query_merge_stats

rule construct_metabin:
    output: 'seq/{group}.a.mbins.d/{bin_id}.fn'
    input: 'res/{group}.a.bins.checkm_merge_stats.tsv'
    shell:
        """
        false  # {input} is new.  Create {output} or touch it to declare that it's up-to-date.
        """

localrules: construct_metabin

rule get_mbin_contig_ids:
    output: 'res/{group}.a.mbins.d/{bin_id}.list'
    input: 'seq/{group}.a.mbins.d/{bin_id}.fn'
    shell: 'ls_ids < {input} > {output}'

# TODO: Also extract paired-ends that didn't map to mbin directly.
# TODO: Filter out unmatched lines before the fastq step, so that we can keep
# all of the metadata.
# TODO: How to get at those reads?
# TODO: How to make -F3844 more expressive? https://broadinstitute.github.io/picard/explain-flags.html
rule extract_mbin_reads:
    output:
        fqgz='seq/{group}.a.mbins.d/{bin_id}.m.pe.fq.gz',
        sam=temp('res/{group}.a.mbins.d/{bin_id}.m.sam'),
        tmp1=temp('res/{group}.a.mbins.d/{bin_id}.m.sam.1.temp'),
        tmp2=temp('res/{group}.a.mbins.d/{bin_id}.m.sam.2.temp'),
        tmp3=temp('res/{group}.a.mbins.d/{bin_id}.m.sam.3.temp')
    input:
        script='scripts/match_paired_reads.py',
        bam='res/{group}.a.contigs.map.sort.bam',
        bai='res/{group}.a.contigs.map.sort.bam.bai',
        contig_ids='res/{group}.a.mbins.d/{bin_id}.list'
    threads: 2
    shell:
        r"""
        echo "Outputting header for {wildcards.bin_id}"
        samtools view -H {input.bam} > {output.sam}

        echo "Collecting intra-bin linking pairs for {wildcards.bin_id}"
        samtools view -@ {threads} -f 1 -F 3842 {input.bam} $(cat {input.contig_ids}) \
            | {input.script} {output.tmp1} >> {output.sam}

        # TODO: Figure out how to control memory usage for such a large 'samtool | grep'
        # echo "Collecting inter-bin linking pairs for {wildcards.bin_id}"
        # Find all of the linked out-of-bin contigs (output.tmp1).
        # Pull discordant reads that map to these contigs.
        # Then filter by the list of contigs that are in the bin and output
        # these lines to the sam-file.
        # samtools view -@ {threads} -f 1 -F 3854 {input.bam} \
        #         $(cut -f7 {output.tmp1} | sort | uniq) \
        #     | grep -wf {input.contig_ids} >> {output.tmp2}
        # cat {output.tmp1} {output.tmp2} | {input.script} \
        #     >> {output.sam}
        touch {output.tmp2}

        echo "Collecting singly mapped pairs for {wildcards.bin_id}"
        samtools view -@ {threads} -f 5 -F 3848 {input.bam} \
            | grep -wf {input.contig_ids} > {output.tmp3}
        samtools view -@ {threads} -f 9 -F 3844 {input.bam} $(cat {input.contig_ids}) \
            >> {output.tmp3}
        {input.script} < {output.tmp3} \
            >> {output.sam}

        echo "Collecting proper pairs for {wildcards.bin_id}"
        samtools view -@ {threads} -f 3 -F 3852 {input.bam} $(cat {input.contig_ids}) \
            | {input.script} >> {output.sam}

        echo "Converting to GZIPed FASTQ for {wildcards.bin_id}"
        samtools view -u {output.sam} | samtools fastq - | gzip > {output.fqgz}
        """


# {{{2 Annotation

# TODO: Be more explicit than {stem} in the below:

# TODO: Output the individual files rather than the directory
rule annotate_mag:
    output:
        fa="seq/{stem}.mags.annot.d/{mag_id}.prokka.fa",
        fn="seq/{stem}.mags.annot.d/{mag_id}.prokka.fn",
        gbk="seq/{stem}.mags.annot.d/{mag_id}.prokka.gbk",
        tbl="res/{stem}.mags.annot.d/{mag_id}.prokka.tbl",
        tsv="res/{stem}.mags.annot.d/{mag_id}.prokka.tsv",
        gff="res/{stem}.mags.annot.d/{mag_id}.prokka.gff",
    input: "seq/{stem}.mags.d/{mag_id}.fn"
    threads: max_threads
    shell:
        r"""
        prokka --force --cpus {threads} {input} \
                --outdir res/{wildcards.stem}.mags.prokka_temp.d --prefix {wildcards.mag_id} \
                --locustag {wildcards.mag_id}
        mv res/{wildcards.stem}.mags.prokka_temp.d/{wildcards.mag_id}.faa {output.fa}
        mv res/{wildcards.stem}.mags.prokka_temp.d/{wildcards.mag_id}.ffn {output.fn}
        mv res/{wildcards.stem}.mags.prokka_temp.d/{wildcards.mag_id}.gbk {output.gbk}
        mv res/{wildcards.stem}.mags.prokka_temp.d/{wildcards.mag_id}.tbl {output.tbl}
        mv res/{wildcards.stem}.mags.prokka_temp.d/{wildcards.mag_id}.tsv {output.tsv}
        mv res/{wildcards.stem}.mags.prokka_temp.d/{wildcards.mag_id}.gff {output.gff}
        """

rule extract_ec_numbers:
    output: 'res/{stem}.ec.tsv'
    input: 'res/{stem}.prokka.tsv'
    shell:
        """
        cut -f1,5 {input} | awk -v OFS='\t' 'NR > 1 && $2 != "" {{print $1,$2}}' > {output}
        """

rule extract_cogs:
    output: 'res/{stem}.cog.tsv'
    input: 'res/{stem}.prokka.tsv'
    shell:
        """
        cut -f1,6 {input} | awk -v OFS='\t' 'NR > 1 && $2 != "" {{print $1,$2}}' > {output}
        """

rule convert_cogs_to_ko:
    output: 'res/{stem}.ko.tsv'
    input: mapping='ref/cog_to_ko.tsv', cogs='res/{stem}.cog.tsv'
    shell:
        """
        join -t '\t' <(sort -k2 {input.cogs}) -1 2 <(sort -k1 {input.mapping}) -2 1 | cut -f2,3 > {output}
        """


localrules: extract_cogs, extract_ec_numbers

rule find_orfs:
    output: nucl="{stem}.orfs.fn", prot="{stem}.orfs.fa"
    input: "{stem}.fn"
    shell:
        "prodigal -q -p meta -i {input} -o /dev/null -d {output.nucl} -a {output.prot}"

rule press_hmm:
    output: "ref/hmm/{stem}.hmm.h3f",
            "ref/hmm/{stem}.hmm.h3i",
            "ref/hmm/{stem}.hmm.h3m",
            "ref/hmm/{stem}.hmm.h3p"
    input: "ref/hmm/{stem}.hmm"
    shell:
        "hmmpress {input}"

rule search_hmm:
    output: "res/{stem}.{hmm}-hits.hmmer-{cutoff}.tsv"
    wildcard_constraints:
        cutoff='ga|nc|tc'
    input:
        faa = "seq/{stem}.fa",
        hmm = "ref/hmm/{hmm}.hmm",
        h3f = "ref/hmm/{hmm}.hmm.h3f",
        h3i = "ref/hmm/{hmm}.hmm.h3i",
        h3m = "ref/hmm/{hmm}.hmm.h3m",
        h3p = "ref/hmm/{hmm}.hmm.h3p"
    threads: 10
    shell:
        """
        echo "orf_id\tmodel_id\tscore" > {output}
        hmmsearch --cut_{wildcards.cutoff} \\
                  --cpu {threads} \\
                  --tblout >(grep -v '^#' | sed 's:\s\+:\\t:g' | cut -f1,3,6 >> {output}) \\
                  {input.hmm} {input.faa} > /dev/null
        """

# {{{2 Sequences Analysis

# "Permissive" means that we include low scoring and partial hits.
# TODO: Are there cases where I want to pull this gene not permissively?
rule gather_hit_orfs_permissive:
    output:
        nucl="seq/{stem}.orfs.{hmm}-hits.fn",
        prot="seq/{stem}.orfs.{hmm}-hits.fa"
    input:
        hit_table="res/{stem}.orfs.{hmm}-hits.hmmer-nc.tsv",
        nucl="seq/{stem}.orfs.fn",
        prot="seq/{stem}.orfs.fa"
    shell:
        """
        seqtk subseq {input.nucl} <(sed 1,1d {input.hit_table} | cut -f 1) > {output.nucl}
        seqtk subseq {input.prot} <(sed 1,1d {input.hit_table} | cut -f 1) > {output.prot}
        """

# rule gather_hit_complete_orfs:
#     output:
#         nucl="seq/{stem}.orfs.{hmm}-hits.fn",
#         prot="seq/{stem}.orfs.{hmm}-hits.fa"
#     input:
#         hit_table="res/{stem}.orfs.{hmm}-hits.hmmer-ga.tsv",
#         nucl="seq/{stem}.orfs.fn",
#         prot="seq/{stem}.orfs.fa"
#     shell:
#         """
#         seqtk subseq {input.nucl} <(sed 1,1d {input.hit_table} | cut -f 1) | grep --no-group-separator -A1 'partial=00' > {output.nucl}
#         seqtk subseq {input.prot} <(sed 1,1d {input.hit_table} | cut -f 1) | grep --no-group-separator -A1 'partial=00' > {output.prot}
#         """
#

rule hmmalign_orfs:
    output: "seq/{stem}.orfs.{hmm}-hits.afa",
    input:
        prot="seq/{stem}.orfs.{hmm}-hits.fa",
        hmm="ref/hmm/{hmm}.hmm"
    shell:
        """
        hmmalign --informat FASTA {input.hmm} {input.prot} | convert -f stockholm -t fasta > {output}
        """

rule codonalign:
    output: "seq/{stem}.codonalign.afn"
    input:
        prot="seq/{stem}.afa",
        nucl="seq/{stem}.fn"
    shell:
        "codonalign {input.prot} {input.nucl} > {output}"

rule squeeze_codon_alignment:
    output: "seq/{stem}.codonalign.sqz.afn"
    input:
        script="scripts/squeeze_alignment.py",
        seq="seq/{stem}.codonalign.afn"
    shell: "{input.script} '-.acgtu' < {input.seq} > {output}"

rule squeeze_hmmalign_alignment:
    output: "seq/{stem}.orfs.{hmm}-hits.sqz.afa"
    wildcard_constraints:
        hmm="[^.]*"
    input:
        script="scripts/squeeze_alignment.py",
        seq="seq/{stem}.orfs.{hmm}-hits.afa"
    shell: "{input.script} '-.*abcdefghijklmnopqrstuvwxyz' < {input.seq} > {output}"

rule estimate_phylogeny_afn:
    output: "res/{stem}.sqz.nucl.nwk"
    input: "seq/{stem}.sqz.afn"
    shell: "FastTree -nt < {input} > {output}"

rule estimate_phylogeny_afa:
    output: "res/{stem}.sqz.prot.nwk"
    input: "seq/{stem}.sqz.afa"
    shell: "FastTree < {input} > {output}"

# TODO: Am I sure I want to use the nucleotide tree for sorting?
rule tree_sort_afn:
    output: "seq/{stem}.tree-sort.afn"
    input:
        script="scripts/get_ordered_leaves.py",
        tree="res/{stem}.nucl.nwk",
        seqs="seq/{stem}.afn"
    shell: "fetch_seqs --match-order <({input.script} {input.tree}) {input.seqs} > {output}"

rule tree_sort_afa:
    output: "seq/{stem}.tree-sort.afa"
    input:
        script="scripts/get_ordered_leaves.py",
        tree="res/{stem}.prot.nwk",
        seqs="seq/{stem}.afa"
    shell: "fetch_seqs --match-order <({input.script} {input.tree}) {input.seqs} > {output}"





# {{{1 Compile all data

# Base database, containing static metadata.
rule generate_database_0:
    output: 'res/{group}.0.db'
    input:
        schema='schema.sql',
        mouse='res/mouse.noheader.tsv',
        sample='res/sample.noheader.tsv',
        extraction='res/extraction.noheader.tsv',
        library='res/library.noheader.tsv',
        asmbl_group='res/asmbl_group.noheader.tsv',
        rrs_taxon_count='res/{group}.r.count.noheader.tsv',
        rrs_taxonomy='res/{group}.r.taxonomy.noheader.tsv',
    shell:
        r"""
        rm -f {output}
        echo '
.bail ON
PRAGMA cache_size = 1000000;
PRAGMA foreign_keys = TRUE;
.read {input.schema}
.separator \t
.import {input.mouse} mouse
.import {input.sample} sample
.import {input.extraction} extraction
.import {input.library} library
.import {input.asmbl_group} library_asmbl_group
.import {input.rrs_taxon_count} rrs_taxon_count
.import {input.rrs_taxonomy} rrs_taxonomy
ANALYZE;
             ' \
        | sqlite3 {output}
        """

localrules: generate_database_0

# First iteration of results db; this might be used to e.g. find scaffolds.
rule generate_database_1:
    output: 'res/{group}.1.db'
    input:
        db='res/{group}.0.db',
        contig='res/{group}.a.contigs.nlength.noheader.tsv',
        contig_bin='res/{group}.a.contigs.bins.noheader.tsv',
        contig_coverage='res/{group}.a.contigs.cvrg.noheader.tsv',
        bin_checkm='res/{group}.a.bins.checkm_details.noheader.tsv',
        contig_linkage='res/{group}.a.contigs.linkage_tally.noheader.tsv',
        checkm_merge='res/{group}.a.bins.checkm_merge_stats.noheader.tsv',
    shell:
        r"""
        tmpfile=$(mktemp -p $TMPDIR)
        cp {input.db} $tmpfile
        echo '
.bail ON
PRAGMA cache_size = 1000000;
PRAGMA foreign_keys = TRUE;
.separator \t
.import {input.checkm_merge} _bin_complementarity
.import {input.contig} contig
.import {input.contig_linkage} _contig_linkage
.import {input.contig_bin} contig_bin
.import {input.contig_coverage} contig_coverage
.import {input.bin_checkm} bin_checkm
ANALYZE;
             ' \
        | sqlite3 $tmpfile
        cp $tmpfile {output}
        """

rule denormalize_database:
    output: 'res/{stem}.denorm.db'
    input:
        db='res/{stem}.db',
        script='scripts/denormalize_db.sql',
    shell:
        """
        tmpfile=$(mktemp -p $TMPDIR)
        cp {input.db} $tmpfile
        cat <(echo "PRAGMA cache_size = 1000000;") {input.script} | sqlite3 $tmpfile
        sqlite3 $tmpfile "VACUUM; ANALYZE;"
        cp $tmpfile {output}
        """
