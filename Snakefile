# {{{1 Imports

from itertools import product
import pandas as pd
from snake.misc import alias_recipe, alias_fmt, curl_recipe, curl_unzip_recipe

# {{{1 Configuration

# {{{2 Nomenclature

one_word_wc_constraint = '[^./]+'
wildcard_constraints:
    group = one_word_wc_constraint,
    library = one_word_wc_constraint,
    mag = one_word_wc_constraint,
    bin = one_word_wc_constraint,
    genomes = one_word_wc_constraint

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
include: 'snake/genome_comparison.snake'

# {{{2 Params

# Default params
MAX_THREADS = 999
if 'MAX_THREADS' in config:
    MAX_THREADS = config['MAX_THREADS']



# {{{1 Utility rules/recipes/templates

# TODO: Figure out what I really want.
mag_proc_infix = 'rsmbl.scaffolds.pilon.ctrim'
rule all:
    input:
        [ f"data/core.a.mags.muri.g.{mag_proc_infix}.genome_stats.tsv"
        , f"data/core.a.mags.muri.g.{mag_proc_infix}.marker_genes.refine.gb.prot.nwk"
        , f"data/core.a.mags.muri.g.{mag_proc_infix}.ec-annot.count.tsv"
        , f"data/core.a.mags.muri.g.{mag_proc_infix}.ec-minpath.count.tsv"
        , f"data/core.a.mags.muri.g.{mag_proc_infix}.cog-annot.count.tsv"
        , f"data/core.a.mags.muri.g.{mag_proc_infix}.dbCAN-hits.annot_table.tsv"
        , f"data/core.a.mags.muri.g.{mag_proc_infix}.dbCAN-hits.denovo-clust.count.tsv"
        , f"data/core.a.mags.muri.g.{mag_proc_infix}.dbCAN-hits.domain-clust.count.tsv"
        ]

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
    output: 'data/{stem}.noheader.{ext,(tsv|csv)}'
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

rule display_rulegraph:
    output: "data/snake.rulegraph.pdf"
    input: "Snakefile", "snake/genome_comparison.snake"
    shell:
        """
        snakemake -n --rulegraph \
                data/core.a.mags.muri.g.rsmbl.scaffolds.pilon.ctrim.dbCAN-hits.denovo-clust.tsv \
                data/core.a.mags.muri.g.rsmbl.scaffolds.pilon.ctrim.genome_stats.tsv \
                data/core.a.mags.muri.g.rsmbl.scaffolds.pilon.ctrim.marker_genes.gb.prot.nwk \
            | dot -Tpdf > {output}
        """

rule display_dag:
    output:
        pdf="data/snake.dag.pdf",
        dot="data/snake.dag.dot",
    input: "Snakefile", "snake/genome_comparison.snake"
    shell:
        "snakemake -n --forceall --dag data/core.a.mags.muri.dbCAN-hits.denovo-clust.tsv | tee {output.dot} | dot -Tpdf > {output.pdf}"


# {{{1 Downloading and linking data

# {{{2 Reference data

# {{{3 General purpose reference genomes

rule download_salask_reference:
    output: 'raw/ref/salask.fn'
    params:
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000356.1&rettype=fasta&retmode=text'
    shell: curl_recipe

rule alias_salask_reference:
    output: 'ref/salask.fn'
    input: 'raw/ref/salask.fn'
    shell: alias_recipe

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

localrules: download_salask_reference, alias_salask_reference,
            download_mouse_reference, download_sra_data,

# {{{3 TIGRFAM

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

rule download_pfam:
    output: "raw/ref/Pfam-31.0.hmm"
    params:
        url="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz"
    shell: curl_unzip_recipe

rule link_pfam:
    output: "ref/hmm/Pfam.hmm"
    input: "raw/ref/Pfam-31.0.hmm"
    shell: alias_recipe

rule alias_hmm:
    output: "ref/hmm/{alias}.hmm"
    input: lambda wildcards: "ref/hmm/{}.hmm".format(config['gene_to_hmm'][wildcards.alias])
    shell: alias_recipe

localrules: download_tigrfam, extract_tigrfam_hmm, download_pfam, link_pfam, alias_hmm

# {{{3 Metadata for sequence processing

rule download_illumina_adapters:
    output: 'raw/ref/illumina_adapters.fn'
    params:
        url='https://raw.githubusercontent.com/vsbuffalo/scythe/master/illumina_adapters.fa'
    shell: curl_recipe

localrules: download_illumina_adapters

# {{{3 dbCAN

rule download_dBCAN_hmms:
    output: "raw/ref/dbCAN.hmm"
    params:
        url="http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt.v6"
    shell: curl_recipe

rule download_dbCAN_meta:
    output: "raw/ref/dbCAN.tsv"
    params:
        url="http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt"
    shell: curl_recipe

rule filter_dbCAN_hmms:
    output: "ref/hmm/dbCAN.hmm"
    input: "raw/ref/dbCAN.hmm"
    shell:
        r"""
        sed 's:^NAME  \(.*\).hmm$:NAME  \1\nDESC  hypothetical carbohydrate-active domain (\1) containing protein:' {input} > {output}
        """

rule filter_dbCAN_meta:
    output: "ref/dbCAN.tsv"
    input: "raw/ref/dbCAN.tsv"
    shell:
        r"""
        sed '1s:^#Family\t::' {input} > {output}
        """

# {{{3 COG

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

localrules: download_cog_to_ko_mapping, download_cog_function_mapping,
            alias_cog_to_ko_mapping, process_cog_function_mapping

# {{{3 Enzyme Commission

rule download_ec_mapping:
    output: 'raw/ref/expasy_enzyme.dat'
    params: url='ftp://ftp.expasy.org/databases/enzyme/enzyme.dat'
    shell: curl_recipe

rule transform_ec_mapping:
    output: 'ref/expasy.tsv'
    input: 'raw/ref/expasy_enzyme.dat'
    shell: """
    grep --no-group-separator -A1 '^ID\>' {input} \
        | paste - - \
        | sed -e 's:\<\(DE\|ID\)\>   ::g' \
        > {output}
    """


# {{{3 MetaCyc

rule download_metacyc_pathways_page:
    output: 'raw/ref/metacyc_pathway_page.html'
    params: url='https://biocyc.org/META/class-subs-instances?object=Pathways'
    shell:
        """
        wget -O {output} {params.url}
        """

# TODO: Improve this
rule scrape_metacyc_pathways_table:
    output: 'ref/metacyc_pathway_descriptions.tsv'
    input: 'raw/ref/metacyc_pathway_page.html'
    shell:
        r"""
        grep '^<li><a href=/META/NEW-IMAGE?type=PATHWAY&object=' {input} \
            | sed 's:.*object=\([^<>]*\)>\(.*\)</a>$:\1\t\2:' \
            | sed -e 's:\(<[iI]>\|</[iI]>\|<sup>\|</sup>\|<sub>\|</sub>\|<SUB>\|</SUB>\)::g' \
            | sort | uniq \
            > {output}
        """

rule download_picrust_metacyc_pathways:
    output: 'raw/ref/ec2path.picrust_prokaryotic.tsv'
    params: url='https://raw.githubusercontent.com/picrust/picrust2/master/MinPath/ec2metacyc_picrust_prokaryotic.txt'
    shell: curl_recipe

rule reformat_picrust_metacyc_pathways:
    output: 'ref/ec2path.picrust.tsv'
    input: 'raw/ref/ec2path.picrust_prokaryotic.tsv'
    shell:
        """
        printf '#Metacyc pathway and ec mapping file\n' > {output}
        printf '#Pathway\tEC\n' >> {output}
        sed 's/EC://' < {input} >> {output}
        """


localrules: scrape_metacyc_pathways_table,
            download_dBCAN_hmms, download_dbCAN_meta,
            filter_dbCAN_hmms, filter_dbCAN_meta

# {{{2 Raw data

rule alias_raw_read_r1:
    output: 'data/{library}.m.r1.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r1'])
    shell: alias_recipe

rule alias_raw_read_r2:
    output: 'data/{library}.m.r2.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r2'])
    shell: alias_recipe

localrules: alias_raw_read_r1, alias_raw_read_r2, link_rabund_info

# {{{3 Import results from 16S libraries

rule query_count_info:
    output: 'data/core.r.count.tsv'
    input: script='scripts/query_longev_rrs_count.sql', db='raw/longev_rrs_results.db'
    shell: "sqlite3 -header -separator '\t' {input.db} < {input.script} > {output}"

rule query_taxonomy_info:
    output: 'data/core.r.taxonomy.tsv'
    input: script='scripts/query_longev_rrs_taxonomy.sql', db='raw/longev_rrs_results.db'
    shell: "sqlite3 -header -separator '\t' {input.db} < {input.script} > {output}"

rule alias_rrs_reps:
    output: 'ref/core.r.reps.fn'
    input: 'raw/longev_rrs_reps.fn'
    shell: alias_recipe

localrules: query_count_info, query_taxonomy_info

# {{{1 Metagenomic Assembly

# {{{2 Data pre-processing

rule deduplicate_reads:
    output:
        r1='data/{stem}.m.r1.dedup.fq.gz',
        r2='data/{stem}.m.r2.dedup.fq.gz'
    input:
        script='scripts/fastuniq_wrapper.sh',
        r1='data/{stem}.m.r1.fq.gz',
        r2='data/{stem}.m.r2.fq.gz'
    resources:
        mem_mb=24000
    shell: "{input.script} {input.r1} {input.r2} {output.r1} {output.r2}"

rule trim_adapters:
    output: 'data/{stem}.deadapt.fq.gz'
    input:
        adapters='raw/ref/illumina_adapters.fn',
        seqs='data/{stem}.fq.gz'
    log: 'log/{stem}.scythe.log'
    shell:
        """
        scythe -a {input.adapters} {input.seqs} >{output} 2>{log}
        ! grep -Fxq 'Blank FASTA header or sequence in adapters file.' {log}
        """

rule quality_trim_reads:
    output:
        r1='data/{stem}.r1.{proc}.qtrim.fq.gz',
        r2='data/{stem}.r2.{proc}.qtrim.fq.gz',
        r3='data/{stem}.r3.{proc}.qtrim.fq.gz'
    input:
        r1='data/{stem}.r1.{proc}.fq.gz',
        r2='data/{stem}.r2.{proc}.fq.gz'
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
    output: 'data/{library}.m.r1.proc.fq.gz'
    input: 'data/{library}.m.r1.dedup.deadapt.qtrim.fq.gz'
    shell: alias_recipe

rule alias_read_processing_r2:
    output: 'data/{library}.m.r2.proc.fq.gz'
    input: 'data/{library}.m.r2.dedup.deadapt.qtrim.fq.gz'
    shell: alias_recipe

localrules: alias_read_processing_r1, alias_read_processing_r2

# {{{2 Metagenome Assembly

# {{{3 MEGAHIT

# TODO: Consider dropping the contig renaming since it's redundant with the filename.
rule assemble_mgen:
    output:
        fasta='data/{group}.a.contigs.fn',
        dir='data/{group}.a.megahit.d',
        fastg='data/{group}.a.contigs.fg',
    input:
        lambda wildcards: [f'data/{library}.m.{read}.proc.fq.gz'
                           for library, read
                           in product(config['asmbl_group'][wildcards.group],
                                      ['r1', 'r2'])
                          ]
    log: 'data/{group}.a.megahit.log'
    threads: MAX_THREADS
    params:
        r1=lambda wildcards: ','.join([f'data/{library}.m.r1.proc.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
        r2=lambda wildcards: ','.join([f'data/{library}.m.r2.proc.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
    shell:
        r"""
        megahit \
            -1 {params.r1} \
            -2 {params.r2} \
            --k-min 21 --k-max 161 --k-step 20 \
            --out-dir {output.dir} \
            --num-cpu-threads {threads} \
            --verbose \
            2>&1 > {log}
        sed 's:^>k161_\([0-9]\+\):>{wildcards.group}_\1:' {output.dir}/final.contigs.fa > {output.fasta}
        # TODO: Fix this hard-coding of k-parameters.
        megahit_toolkit contig2fastg 141 {output.dir}/intermediate_contigs/k141.contigs.fa > {output.fastg}
        """

# {{{3 QC

rule quality_asses_assembly_with_spike:
    output: 'data/{group}.a.quast.d'
    input: contigs='data/{group}.a.contigs.fn', ref='ref/salask.fn'
    threads: MAX_THREADS
    params: min_contig_length=1000,
    shell:
        """
        metaquast.py --threads={threads} --min-contig {params.min_contig_length} -R {input.ref} --output-dir {output} {input.contigs}
        """

# {{{2 Mapping

# {{{3 Index building

rule bowtie_index_build:
    output:
        'data/{stem}.1.bt2',
        'data/{stem}.2.bt2',
        'data/{stem}.3.bt2',
        'data/{stem}.4.bt2',
        'data/{stem}.rev.1.bt2',
        'data/{stem}.rev.2.bt2'
    input: 'data/{stem}.fn'
    log: 'log/{stem}.bowtie2-build.log'
    threads: min(5, MAX_THREADS)
    shell:
        """
        bowtie2-build --threads {threads} {input} data/{wildcards.stem} >{log} 2>&1
        """

# {{{3 Backmap

# TODO: Consider matching the MAG back-mapping recipe by including a {proc} wildcard after '-map.'
rule map_reads_to_metagenome_assembly:
    output: 'data/{library}.m.{group}-map.sort.bam'
    input:
        r1='data/{library}.m.r1.proc.fq.gz',
        r2='data/{library}.m.r2.proc.fq.gz',
        inx_1='data/{group}.a.contigs.1.bt2',
        inx_2='data/{group}.a.contigs.2.bt2',
        inx_3='data/{group}.a.contigs.3.bt2',
        inx_4='data/{group}.a.contigs.4.bt2',
        inx_rev1='data/{group}.a.contigs.rev.1.bt2',
        inx_rev2='data/{group}.a.contigs.rev.2.bt2'
    threads: MAX_THREADS
    shell:
        r"""
        tmp=$(mktemp)
        echo "Writing temporary bamfile to $tmp for {output}"
        bowtie2 --threads {threads} \
                -x data/{wildcards.group}.a.contigs \
                --rg-id {wildcards.library} \
                -1 {input.r1} -2 {input.r2} \
            | samtools view -@ {threads} -b - \
            > $tmp
        samtools sort -@ {threads} --output-fmt=BAM -o {output} $tmp
        rm $tmp
        """

rule combine_read_mappings:
    output: 'data/{group}.a.map-{group}.sort.bam'
    input:
        lambda wildcards: [f'data/{library}.m.{wildcards.group}-map.sort.bam'
                           for library
                           in config['asmbl_group'][wildcards.group]
                          ]
    threads: min(10, MAX_THREADS)
    shell:
        """
        samtools merge -@ {threads} {output} {input}
        """
rule index_read_mappings:
    output: 'data/{stem}.sort.bam.bai'
    input: 'data/{stem}.sort.bam'
    threads: min(6, MAX_THREADS)
    shell: 'samtools index -@ {threads} {input} {output}'


# {{{2 Scaffolding

rule tally_links:
    output: 'data/{library}.m.{group}-map.linkage_tally.tsv'
    input: 'data/{library}.m.{group}-map.sort.bam'
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
    output: 'data/{group}.a.contigs.linkage_tally.tsv'
    input:
        lambda wildcards: [f'data/{library}.m.{wildcards.group}-map.linkage_tally.tsv'
                           for library in config['asmbl_group'][wildcards.group]
                          ]
    shell:
        r"""
        printf 'contig_id_1\tcontig_id_2\tlibrary_id\tread_count\n' > {output}
        for file in {input}; do
            sed 1,1d $file
        done >> {output}
        """

# {{{2 Calculate statistics

# {{{3 Sequence lengths

rule count_seq_lengths_nucl:
    output: 'data/{stem}.nlength.tsv'
    input: script='scripts/count_seq_lengths.py', seqs='data/{stem}.fn'
    shell:
        r"""
        printf "contig_id\tlength\n" > {output}
        {input.script} {input.seqs} >> {output}
        """

rule count_seq_lengths_aa:
    output: 'data/{stem}.alength.tsv'
    input: script='scripts/count_seq_lengths.py', seqs='data/{stem}.fa'
    shell:
        """
        printf "seq_id\tlength\n" > {output}
        {input.script} {input.seqs} >> {output}
        """

localrules: count_seq_lengths_nucl, count_seq_lengths_aa

# {{{3 Coverage

# NOTE: The depth file format is lacking a header.
# TODO params: -Q flag (mapping_quality_thresh), -d 0 flag (no maximum mapping depth)
rule calculate_mapping_depth:
    output: 'data/{stem}.depth.tsv'
    input: 'data/{stem}.sort.bam'
    shell:
        """
        samtools depth -d 0 {input} > {output}
        """

rule estimate_contig_cvrg:
    output: 'data/{library}.m.{group}-map.cvrg.tsv'
    input:
        script='scripts/estimate_contig_coverage.py',
        depth='data/{library}.m.{group}-map.depth.tsv',
        length='data/{group}.a.contigs.nlength.tsv'
    params:
        float_fmt='%.6g'
    shell:
        """
        {input.script} {input.depth} {input.length} {params.float_fmt} > {output}
        """

rule combine_cvrg:
    output: 'data/{group}.a.contigs.cvrg.tsv'
    input:
        script='scripts/concat_tables.py',
        tables=lambda wildcards: [f'data/{library}.m.{wildcards.group}-map.cvrg.tsv'
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
    output: 'data/{group}.a.contigs.cvrg.unstack.tsv'
    input:
        script='scripts/unstack_cvrg.py',
        cvrg='data/{group}.a.contigs.cvrg.tsv',
    params:
        float_format='%.6g'
    shell:
        """
        {input.script} {input.cvrg} '{params.float_format}' > {output}
        """

# {{{1 MAGs
# {{{2 MAG Finding

# {{{3 Prepare input data

rule transform_contig_space:
    output:
        pca='data/{group}.a.contigs.concoct.pca.tsv',
        raw='data/{group}.a.contigs.concoct.tsv',
        dir=temp('data/{group}.a.contigs.concoct.d')
    input:
        cvrg='data/{group}.a.contigs.cvrg.unstack.tsv',
        seqs='data/{group}.a.contigs.fn'
    threads: MAX_THREADS
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
        out='data/{group}.a.contigs.cluster.tsv',
        summary='data/{group}.a.contigs.cluster.summary.tsv',
    input:
        script='scripts/cluster_contigs.py',
        pca='data/{group}.a.contigs.concoct.pca.tsv',
        length='data/{group}.a.contigs.nlength.tsv',
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
    output: 'data/{group}.a.contigs.bins.tsv'
    input: 'data/{group}.a.contigs.cluster.tsv'
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
    output: 'data/{group}.a.bins.d'
    input:
        script='scripts/fetch_bin.sh',
        bins='data/{group}.a.contigs.bins.tsv',
        contigs='data/{group}.a.contigs.fn',
    threads: MAX_THREADS
    shell:
        r"""
        rm -rf {output}
        mkdir {output}
        bins=$(sed '1,1d' {input.bins} | cut -f 2 | sort | uniq)
        parallel --progress --jobs {threads} {input.script} {input.contigs} {input.bins} {{1}} {output}/{{1}}.fn ::: $bins
        """

# TODO: refine bins (scaffolds, contig extension, splitting/merging

# {{{3 QC bins

rule checkm_seqs:
    output:
        dir=temp('data/{stem}.checkm.d'),
        summary='data/{stem}.checkm.tsv'
    input: 'data/{stem}.d'
    threads: MAX_THREADS
    shell:
        r"""
        rm -rf {output.dir}
        checkm lineage_wf -x fn \
                --threads {threads} --pplacer_threads {threads} \
                --file {output.summary} --tab_table \
                {input} {output.dir}
        """

rule checkm_ctrim_set:
    output:
        dir=temp('data/{group}.a.mags/{mag}.g.{proc}.ctrim_check.checkm.d'),
        summary='data/{group}.a.mags/{mag}.g.{proc}.ctrim_check.checkm.tsv'
    input:
        [ 'data/{group}.a.mags/{mag}.g.{proc}.fn'
        , 'data/{group}.a.mags/{mag}.g.{proc}.ctrim-50.fn'
        , 'data/{group}.a.mags/{mag}.g.{proc}.ctrim-60.fn'
        , 'data/{group}.a.mags/{mag}.g.{proc}.ctrim-70.fn'
        , 'data/{group}.a.mags/{mag}.g.{proc}.ctrim-80.fn'
        , 'data/{group}.a.mags/{mag}.g.{proc}.ctrim-90.fn'
        ]
    threads: 2
    shell:
        r"""
        tmpdir=$(mktemp -d)
        ln -rst $tmpdir {input}
        rm -rf {output.dir}
        checkm lineage_wf -x fn \
                --threads {threads} --pplacer_threads {threads} \
                --file {output.summary} --tab_table \
                $tmpdir {output.dir}
        """

ruleorder: checkm_ctrim_set > checkm_seqs

rule reformat_checkm_output:
    output: 'data/{stem}.checkm_details.tsv'
    input: 'data/{stem}.checkm.tsv'
    shell:
        """
        cut -f1,4,12-15 {input} > {output}
        """

rule generate_checkm_markerset:
    output:
        'data/{level}_{taxon}.ms'
    shell:
        'checkm taxon_set {wildcards.level} {wildcards.taxon} {output}'

# {{{3 Bins to MAGs

# TODO: Understand what field 9 in checkM output file is.
# (I _think_ it's the difference between the increased completeness and the
# increased contamination.)
rule checkm_content_merge:
    output:
        checkm_work=temp('data/{group}.a.bins.checkm_merge.d'),
        merge_stats='data/{group}.a.bins.checkm_merge_stats.tsv',
    input:
        bins='data/{group}.a.bins.d',
        markerset='data/domain_Bacteria.ms',
    threads: MAX_THREADS
    shadow: 'full'
    shell:
        """
        rm -rf {output.checkm_work}
        checkm merge --threads {threads} -x fn {input.markerset} {input.bins} {output.checkm_work}
        printf 'bin_id_1\tbin_id_2\tscore\n' > {output.merge_stats}
        sed '1,1d' {output.checkm_work}/merger.tsv | cut -f1,2,9 >> {output.merge_stats}
        """

rule query_merge_stats:
    output: 'data/{group}.a.bin_merge_stats.tsv'
    input: sql='scripts/query_bin_merger.sql', db='data/{group}.1.denorm.db'
    shell:
        """
        sqlite3 -header -separator '\t' {input.db} < {input.sql} > {output}
        """

localrules: query_merge_stats

# TODO: Automate this??
# TODO: Swap checkm_merge_stats input for bin_merge_stats? (this will introduce a dependencies on databases and therefore schema.sql
rule select_curated_mags:
    output:
        contigs='data/{group}.a.mags/{mag}.g.contigs.list',
        libraries='data/{group}.a.mags/{mag}.g.libraries.list',
    input: 'data/{group}.a.bins.checkm_merge_stats.tsv'
    shell:
        """
        cat <<EOF
        Select the contigs that belong to {wildcards.mag} and the libraries
        from which to collect reads for refinement.

        This can be accomplished by merely touching
        `{output.contigs}`
        and `{output.libraries}`
        if they already exist and you are confident that they reflects the
        current state of `data/{group}.a.bins.*` .

        However, since `{input}`
        is newer, this may not be the case, or you may
        want to refine your curation.  Instead you should manually identify all
        of contigs that belong to this MAG  and save their names to this list.

EOF
        false
        """

rule get_mag_contigs:
    output: 'data/{group}.a.mags/{bin}.g.contigs.fn'
    input:
        ids='data/{group}.a.mags/{bin}.g.contigs.list',
        seqs='data/{group}.a.contigs.fn'
    shell: 'seqtk subseq {input.seqs} {input.ids} > {output}'

localrules: select_curated_mags, get_mag_contigs

# {{{2 MAG Refinement
# {{{3 Mapping

rule map_reads_to_mag:
    output: 'data/{group}.a.mags/{library}.m.{mag}-map.{proc}.sort.bam'
    input:
        r1='data/{library}.m.r1.proc.fq.gz',
        r2='data/{library}.m.r2.proc.fq.gz',
        inx_1='data/{group}.a.mags/{mag}.g.{proc}.1.bt2',
        inx_2='data/{group}.a.mags/{mag}.g.{proc}.2.bt2',
        inx_3='data/{group}.a.mags/{mag}.g.{proc}.3.bt2',
        inx_4='data/{group}.a.mags/{mag}.g.{proc}.4.bt2',
        inx_rev1='data/{group}.a.mags/{mag}.g.{proc}.rev.1.bt2',
        inx_rev2='data/{group}.a.mags/{mag}.g.{proc}.rev.2.bt2',
    params:
        inx_name="data/{group}.a.mags/{mag}.g.{proc}"
    shell:
        """
        tmp1=$(mktemp)
        tmp2=$(mktemp)
        tmp3=$(mktemp)
        echo "$tmp1 $tmp2 $tmp3 ({output})"
        bowtie2 -x {params.inx_name} \
                --rg-id {wildcards.library} \
                -1 {input.r1} -2 {input.r2} \
            | samtools view -G 12 -b - > $tmp1
        # Header
        samtools view -H $tmp1 > $tmp2
        # Both mapped
        samtools view -F 3852 $tmp1 >> $tmp2
        # Only other read mapped
        samtools view -f 4 -F 3848 $tmp1 >> $tmp2
        # Only other read unmapped
        samtools view -f 8 -F 3844 $tmp1 >> $tmp2
        samtools view -b $tmp2 > $tmp3
        samtools sort --output-fmt=BAM -o {output} $tmp3
        rm $tmp1 $tmp2 $tmp3
        """

rule combine_mag_read_mappings:
    output: 'data/{group}.a.mags/{mag}.g.{proc}.map-{group}.sort.bam'
    input:
        lambda wildcards: [f'data/{wildcards.group}.a.mags/{library}.m.{wildcards.mag}-map.{wildcards.proc}.sort.bam'
                           for library
                           in config['asmbl_group'][wildcards.group]
                          ]
    threads: min(10, MAX_THREADS)
    shell:
        """
        samtools merge -@ {threads} {output} {input}
        """

rule extract_relevant_read_mappings:
    output: 'data/{group}.a.mags/{mag}.g.{proc}.map.sort.bam'
    input:
        bam='data/{group}.a.mags/{mag}.g.{proc}.map-{group}.sort.bam',
        libs='data/{group}.a.mags/{mag}.g.library.list',
    threads: min(10, MAX_THREADS)
    shell: "samtools view -@ {threads} -b -R {input.libs} {input.bam} > {output}"

rule extract_all_read_mappings:
    output: 'data/{group}.a.mags/{mag}_v0.g.{proc}.map.sort.bam'
    input: 'data/{group}.a.mags/{mag}_v0.g.{proc}.map-{group}.sort.bam',
    shell: alias_recipe

ruleorder: extract_all_read_mappings > extract_relevant_read_mappings


# {{{3 Prepare Data

rule bam_to_fastq:
    output:
        r1='data/{stem}.r1.fq.gz',
        r2='data/{stem}.r2.fq.gz',
    input:
        bam='data/{stem}.sort.bam',
        script='scripts/match_paired_reads.py',
    threads: min(10, MAX_THREADS)
    shell:
        r"""
        # TODO
        (
            samtools view -H {input.bam}
            samtools view -@ {threads} {input.bam} | {input.script}
        ) \
            | samtools view -@ {threads} -u - \
            | samtools fastq -@ {threads} -c 6 -1 {output.r1} -2 {output.r2} -
        """

rule diginorm_reads:
    output:
        r1='data/{stem}.r1.dnorm.fq.gz',
        r2='data/{stem}.r2.dnorm.fq.gz'
    input:
        r1='data/{stem}.r1.fq.gz',
        r2='data/{stem}.r2.fq.gz'
    params:
        ksize=32,
        tablecols=1024,
        tablerows=10
    threads: 2
    shell:
        """
        Bignorm -1 {input.r1} -2 {input.r2} -k {params.ksize} --both -m {params.tablecols} -t {params.tablerows} -z
        mv {input.r1}_keep.gz {output.r1}
        mv {input.r2}_keep.gz {output.r2}
        """

# {{{3 Reassembly

# TODO: Try larger minimum kmers to reduce missassembly using -k 21,33,55,77
# TODO: Filter contigs by length (some minimum) and by estimated coverage (remove outliers)
# TODO: Rename contigs with mag stem.
# TODO: Consider assembling from reads in all (not just trusted) libraries.
rule reassemble_mag:
    output:
        scaffolds='data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.fn',
        contigs='data/{group}.a.mags/{mag}.g.rsmbl.contigs.fn',
        dir='data/{group}.a.mags/{mag}.g.spades.d',
    input:
        r1='data/{group}.a.mags/{mag}.g.contigs.map.r1.dnorm.fq.gz',
        r2='data/{group}.a.mags/{mag}.g.contigs.map.r2.dnorm.fq.gz'
    threads: min(15, MAX_THREADS)
    shell:
        r"""
        spades.py --tmp-dir $TMPDIR --threads {threads} --careful -1 {input.r1} -2 {input.r2} -o {output.dir}
        sed 's:^>NODE_\([0-9]\+\)_length_[0-9]\+_cov_[0-9]\+.[0-9]\+$:>{wildcards.mag}_\1:' {output.dir}/scaffolds.fasta > {output.scaffolds}
        sed 's:^>NODE_\([0-9]\+\)_length_[0-9]\+_cov_[0-9]\+.[0-9]\+$:>{wildcards.mag}_\1:' {output.dir}/contigs.fasta > {output.contigs}
        """

# {{{3 Pilon Refinement

rule alias_seqs_for_pilon:
    output: 'data/{stem}.fasta'
    input: 'data/{stem}.fn'
    shell: alias_recipe

rule pilon_refine:
    output: "data/{group}.a.mags/{mag}.g.{proc}.pilon.fn",
    input:
        contigs="data/{group}.a.mags/{mag}.g.{proc}.fasta",
        bam="data/{group}.a.mags/{mag}.g.{proc}.map.sort.bam",
        bai="data/{group}.a.mags/{mag}.g.{proc}.map.sort.bam.bai",
    resources:
        mem_mb=100000
    params:
        prefix="data/{group}.a.mags/{mag}.g.{proc}.pilon"
    threads: min(12, MAX_THREADS)
    shell:
        r"""
        pilon -Xms1024m -Xmx{resources.mem_mb}m --threads {threads} \
                --fix snps,indels,gaps,local,breaks,amb \
                --genome {input.contigs} --frags {input.bam} \
                --output {params.prefix}
        mv {params.prefix}.fasta {output}
        """

# {{{3 Depth Trimming

# See rule: calculate_mapping_depth.


# {{{3 Correlation Trimming

# TODO: Use this combined naming scheme for all other "combine_*" recipes.
# i.e. don't combine from a bunch of library specific files into a group
# file with the same name; instead use the naming scheme for the sequences
# they're mapped to.
# TODO: I may want to move all such files to the new namespace so I don't need
# to rebuild them.

rule combine_depths:
    output: 'data/{group}.a.mags/{mag}.g.{proc}.library-depth.tsv.gz'
    input:
        lambda wildcards: [f'data/{wildcards.group}.a.mags/{library}.m.{wildcards.mag}-map.{wildcards.proc}.depth.tsv'
                           for library
                           in config['asmbl_group'][wildcards.group]
                          ]
    shell:
        """
        for file in {input}
        do
            library=$(basename --suffix=.m.{wildcards.mag}-map.{wildcards.proc}.depth.tsv $file)
            echo "Writing $library depths to {output}" >&2
            awk -v OFS='\t' -v library=$library '{{print library,$0}}' $file
        done | gzip > {output}
        """

rule estimate_mag_contig_cvrg:
    output: 'data/{stem}.cvrg.unstack.tsv'
    input:
        script='scripts/calculate_contig_cvrg_all_libs.py',
        depth='data/{stem}.library-depth.tsv.gz',
        length='data/{stem}.nlength.tsv'
    shell:
        """
        {input.script} {input.depth} {input.length} > {output}
        """

rule calculate_position_correlation_stats:
    output: "data/{group}.a.mags/{mag}.g.{proc}.pcorr.tsv"
    input:
        script="scripts/calculate_per_position_stats.py",
        trusted="data/{group}.a.mags/{mag}.g.trusted_depth.tsv",
        depth="data/{group}.a.mags/{mag}.g.{proc}.library-depth.tsv.gz",
        libs='data/{group}.a.mags/{mag}.g.library.list',
    shell:
        """
        {input.script} {input.depth} {input.trusted} {input.libs} > {output}
        """

rule calculate_position_correlation_stats_all_libs:
    output: "data/{group}.a.mags/{mag}_v0.g.{proc}.pcorr.tsv"
    input:
        script="scripts/calculate_per_position_stats.py",
        trusted="data/{group}.a.mags/{mag}.g.trusted_depth.tsv",
        depth="data/{group}.a.mags/{mag}_v0.g.{proc}.library-depth.tsv.gz",
    shell:
        """
        {input.script} {input.depth} {input.trusted} > {output}
        """

ruleorder: calculate_position_correlation_stats_all_libs > calculate_position_correlation_stats

rule plot_position_distribution_plots:
    output: "data/{group}.a.mags/{mag}.g.{proc}.pcorr.hist.pdf"
    input:
        script="scripts/plot_position_correlations_histogram.py",
        corrs="data/{group}.a.mags/{mag}.g.{proc}.pcorr.tsv"
    shell:
        """
        {input.script} {input.corrs} {output}
        """

rule correlation_trim_contigs:
    output:
        fn="data/{stem}.ctrim-{thresh}.fn",
    input:
        script="scripts/correlation_trim_contigs.py",
        scaffolds="data/{stem}.fn",
        corr="data/{stem}.pcorr.tsv",
    wildcard_constraints:
        thresh='[0-9][0-9]'
    params:
        window=100,
        flank=100,
        min_len=300,
    shell:
        """
        {input.script} --corr-thresh=0.{wildcards.thresh} \
                --window-size={params.window} \
                --flank-size={params.flank} \
                --min-length={params.min_len} \
                {input.scaffolds} {input.corr} \
                > {output.fn}
        """

rule select_correlation_trim_threshold:
    output:
        fn="data/{stem}.ctrim.fn",
    input:
        script="scripts/correlation_trim_contigs.py",
        plot="data/{stem}.pcorr.hist.pdf",
        seqs=[f"data/{{stem}}.ctrim-{cutoff}.fn"
              for cutoff in [50, 60, 70, 90]],
        checkm="data/{stem}.ctrim_check.checkm_details.tsv",
    params:
        window=100,
        flank=100,
        min_len=1000,
    shell:
        """
        cat <<EOF
        Pick a threshold for the canonical ctrimmed file.

        Already available:
        {input.seqs}

        If you want a different cutoff
        (replacing THRESHOLD with the chosen threshold):

        snakemake data/{wildcards.stem}.ctrim-THRESHOLD.fn

        You may want to use the following to guide your decision:
        -   Completeness/contamination ({input.checkm})
        -   Number of contigs/scaffolds
        -   Total sequence length
        -   Nucleotide correlation histogram ({input.plot})

        Once you have selected the final result:

        cp data/{wildcards.stem}.ctrim-THRESHOLD.fn {output}

EOF
        false
        """

# {{{3 QC

# TODO: Custom (non-16S) blast db for reference finding
# see http://quast.bioinf.spbau.ru/manual.html#faq_q12
rule quality_asses_mag:
    output: 'data/{group}.a.mags/{mag}.g.quast.d'
    input:
        asmbl=[
               'data/{group}.a.mags/{mag}.g.contigs.fn',
               'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.fn',
               'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.fn',
               'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-50.fn',
               'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-60.fn',
               'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-70.fn',
               'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-80.fn',
               'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-90.fn',
               'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-50.fn',
               'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-60.fn',
               'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-70.fn',
               'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-80.fn',
               'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-90.fn',
               ]
    threads: min(7, MAX_THREADS)
    shell:
        r"""
        quast.py --threads={threads} --min-contig 0 --output-dir {output} \
                -R {input.asmbl[4]} --fragmented \
                {input.asmbl}
        """

# TODO: Custom (non-16S) blast db for reference finding
# see http://quast.bioinf.spbau.ru/manual.html#faq_q12
rule quality_asses_spike_mag:
    output: 'data/{group}.a.mags/{mag}.g.quast-spike.d'
    input:
        ref='ref/salask.fn',
        asmbl=[
               'data/{group}.a.mags/{mag}.g.contigs.fn',
               'data/{group}.a.mags/{mag}.g.scaffolds.fn',
               'data/{group}.a.mags/{mag}.g.scaffolds.pilon.fn',
               'data/{group}.a.mags/{mag}.g.scaffolds.pilon.ctrim.fn',
               ]
    threads: min(6, MAX_THREADS)
    shell:
        r"""
        quast.py --threads={threads} --min-contig 0 --output-dir {output} \
                -R {input.ref} \
                --labels "MAG, Reassembled, Reassembled+Refined, Reassembled+Refined+CTrimmed" \
                {input.asmbl}
        """

localrules: quality_asses_reassembly, quality_asses_spike_reassembly

# {{{3 Strain Comparison
rule compare_strains:
    output:
        delta="data/{group_stem}/{magA}.g.{proc}.{magB}-align.delta",
    input:
        magA="data/{group_stem}/{magA}.g.{proc}.fn",
        magB="data/{group_stem}/{magB}.g.{proc}.fn",
    wildcard_constraints:
        magA=one_word_wc_constraint,
        magB=one_word_wc_constraint,
    shell:
        """
        nucmer --mum --delta {output.delta} {input.magA} {input.magB}
        """

# TODO
# FROM http://mummer.sourceforge.net/manual/#coords
# When run with the -B option, output format will consist of 21 tab-delimited
# columns. These are as follows: [1] query sequence ID [2] date of alignment
# [3] length of query sequence [4] alignment type [5] reference file [6]
# reference sequence ID [7] start of alignment in the query [8] end of
# alignment in the query [9] start of alignment in the reference [10] end of
# alignment in the reference [11] percent identity [12] percent similarity [13]
# length of alignment in the query [14] 0 for compatibility [15] 0 for
# compatibility [16] NULL for compatibility [17] 0 for compatibility [18]
# strand of the query [19] length of the reference sequence [20] 0 for
# compatibility [21] and 0 for compatibility.
rule process_strain_comparison_table:
    output:
        coords="data/{group_stem}/{magA}.g.{proc}.{magB}-align.coords",
    input:
        delta="data/{group_stem}/{magA}.g.{proc}.{magB}-align.delta",
    wildcard_constraints:
        magA=one_word_wc_constraint,
        magB=one_word_wc_constraint,
    shell:
        """
        show-coords -B {input.delta} > {output.coords}
        """

rule plot_strain_comparison:
    output:
        pdf="data/{group_stem}/{magA}.g.{proc}.{magB}-align.pdf",
    input:
        script="scripts/plot_nucmer_comparison.py",
        coords="data/{group_stem}/{magA}.g.{proc}.{magB}-align.coords",
        length1="data/{group_stem}/{magA}.g.{proc}.nlength.tsv",
        length2="data/{group_stem}/{magB}.g.{proc}.nlength.tsv",
        depth1="data/{group_stem}/{magA}.g.{proc}.cvrg.unstack.tsv",
        depth2="data/{group_stem}/{magB}.g.{proc}.cvrg.unstack.tsv",
        lib1="data/{group_stem}/{magA}.g.library.list",
        lib2="data/{group_stem}/{magB}.g.library.list",
    wildcard_constraints:
        magA=one_word_wc_constraint,
        magB=one_word_wc_constraint,
    params:
        alignment_length_thresh=150
    shell:
        """
        {input.script} {input.coords} \
                {input.length1} {input.length2} \
                {input.depth1} {input.depth2} \
                {input.lib1} {input.lib2} \
                {params.alignment_length_thresh} \
                {output}
        """


# {{{2 Annotation

# TODO: Be more explicit than {stem} in the below:

# {{{3 Untargetted
# TODO: Output the individual files rather than the directory
# TODO: Redo annotation now that --metagenome flag has been removed.
rule annotate_mag:
    output:
        fa="data/{stem}.mags.annot/{mag_stem}.cds.fa",
        fn="data/{stem}.mags.annot/{mag_stem}.cds.fn",
        gbk="data/{stem}.mags.annot/{mag_stem}.prokka-annot.gbk",
        tbl="data/{stem}.mags.annot/{mag_stem}.prokka-annot.tbl",
        tsv="data/{stem}.mags.annot/{mag_stem}.prokka-annot.tsv",
        gff="data/{stem}.mags.annot/{mag_stem}.prokka-annot.gff",
        dir=temp("data/{stem}.mags.annot/{mag_stem}.prokka-annot.d"),
    input: "data/{stem}.mags/{mag_stem}.fn"
    threads: min(10, MAX_THREADS)
    shell:
        r"""
        prokka --force --cpus {threads} {input} \
                --outdir {output.dir} --prefix {wildcards.mag_stem} \
                --locustag {wildcards.mag_stem} \
                --metagenome --cdsrnaolap
        cp {output.dir}/{wildcards.mag_stem}.faa {output.fa}
        cp {output.dir}/{wildcards.mag_stem}.ffn {output.fn}
        cp {output.dir}/{wildcards.mag_stem}.gbk {output.gbk}
        cp {output.dir}/{wildcards.mag_stem}.tbl {output.tbl}
        cp {output.dir}/{wildcards.mag_stem}.tsv {output.tsv}
        cp {output.dir}/{wildcards.mag_stem}.gff {output.gff}
        """

rule summarize_annotation:
    output: 'data/{group_stem}.annot/{mag_stem}.prokka.summary.tsv'
    input: 'data/{group_stem}.annot/{mag_stem}.prokka-annot.tsv'
    shell:
        r"""
        echo '
CREATE TABLE annotation
( locus_tag PRIMARY KEY
, ftype
, length_bp
, gene
, ec_number
, cog
, product
);
.separator \t
.import {input} annotation
-- Drop the header
--DELETE FROM annotation WHERE annotation.locus_tag = "locus_tag"

-- Number of loci
SELECT "n_loci", COUNT(*) FROM annotation;
-- Count different feature types
SELECT "n_features_CDS", COUNT(*) FROM annotation WHERE ftype = "CDS";
SELECT "n_features_tRNA", COUNT(*) FROM annotation WHERE ftype = "tRNA";
SELECT "n_features_rRNA", COUNT(*) FROM annotation WHERE ftype = "rRNA";
SELECT "n_features_tmRNA", COUNT(*) FROM annotation WHERE ftype = "tmRNA";
-- Count different annotations of interest
SELECT "n_annotated_hypothetical", COUNT(*) FROM annotation WHERE product = "hypothetical protein";
SELECT "n_annotated_16S", COUNT(*) FROM annotation WHERE product = "16S ribosomal RNA" OR product = "16S ribosomal RNA (partial)";
SELECT "n_with_gene_name", COUNT(*) FROM annotation WHERE gene != "";
SELECT "n_with_ec_number", COUNT(*) FROM annotation WHERE ec_number != "";
SELECT "n_with_cog", COUNT(*) FROM annotation WHERE cog != "";
SELECT "n_with_function", COUNT(*) FROM annotation WHERE gene != "" OR ec_number != "" OR cog != "";
SELECT "n_product_not_hypothetical", COUNT(*) FROM annotation WHERE product != "hypothetical protein";
        ' | sqlite3 > {output}


"""

rule extract_ec_numbers:
    output: 'data/{stem}.ec.tsv'
    input: 'data/{stem}.prokka-annot.tsv'
    shell:
        """
        cut -f1,5 {input} | awk -v OFS='\t' 'NR > 1 && $2 != "" {{print $1,$2}}' > {output}
        """

rule extract_cogs:
    output: 'data/{stem}.cog.tsv'
    input: 'data/{stem}.prokka-annot.tsv'
    shell:
        """
        cut -f1,6 {input} | awk -v OFS='\t' 'NR > 1 && $2 != "" {{print $1,$2}}' > {output}
        """

rule convert_cogs_to_ko:
    output: 'data/{stem}.ko.tsv'
    input: mapping='ref/cog_to_ko.tsv', cogs='data/{stem}.cog.tsv'
    shell:
        """
        join -t '\t' <(sort -k2 {input.cogs}) -1 2 <(sort -k1 {input.mapping}) -2 1 | cut -f2,3 > {output}
        """

rule filter_metacyc_pathways:
    output:
        keep='ref/ec2path.picrust.filt.tsv',
        drop='ref/metacyc_pathway_descriptions.drop.tsv',
    input:
        ecmap='ref/ec2path.picrust.tsv',
        meta='ref/metacyc_pathway_descriptions.tsv',
    params:
        dropregex="superpathway"
    shell:
        """
        grep "{params.dropregex}" {input.meta} > {output.drop}
        grep -v -f <(cut -f 1 {output.drop}) {input.ecmap} > {output.keep}
        """

rule infer_pathways:
    output:
        report='data/{stem}.ec-minpath.report.tsv',
        details='data/{stem}.ec-minpath.details.txt'
    input:
        ec_list='data/{stem}.ec.tsv',
        ec_map='ref/ec2path.picrust.filt.tsv',
    log:
        'data/{stem}.ec-minpath.log'
    threads: MAX_THREADS  # TODO: Figure out how to stop MinPath from overwriting its own processing files (which makes parallel jobs impossible).
    shell:
        """
        MinPath1.4.py -any {input.ec_list} -map {input.ec_map} -report {output.report} -details {output.details} >{log} 2>&1
        """

rule extract_metacyc_list:
    output: 'data/{stem}.ec-minpath.list'
    input: 'data/{stem}.ec-minpath.report.tsv'
    shell:
        """
        awk '$8==1{{print $14}}' {input} > {output}
        """

localrules: extract_cogs, extract_ec_numbers, infer_pathways, extract_metacyc_list

# {{{3 Targetted

# rule pull_annotated_genes:
#     output: fn='data/{stem}.mags.annot/{mag_id}.{gene_id}-hits.fn',
#             fa='data/{stem}.mags.annot/{mag_id}.{gene_id}-hits.fa'
#     input:
#         fn='data/{stem}.mags.annot/{mag_id}.prokka.fn',
#         fa='data/{stem}.mags.annot/{mag_id}.prokka.fa',
#         tsv='data/{stem}.mags.annot/{mag_id}.prokka.tsv'
#     params:
#         search_string=lambda wildcards: config['gene_to_search_string'][wildcards.gene_id]
#     shell:
#         """
#         # Nucleotide
#         seqtk subseq {input.fn} \
#             <(awk -v FS='\t' \
#                   -v search_string='{params.search_string}' \
#                   '$7~search_string{{print $1}}' {input.tsv} \
#              ) \
#             > {output.fn}
#         # Amino-acid
#         seqtk subseq {input.fa} \
#             <(awk -v FS='\t' \
#                   -v search_string='{params.search_string}' \
#                   '$7~search_string{{print $1}}' {input.tsv} \
#              ) \
#             > {output.fa}
#         """
#
# localrules: pull_annotated_genes

rule press_hmm:
    output: "ref/hmm/{stem}.hmm.h3f",
            "ref/hmm/{stem}.hmm.h3i",
            "ref/hmm/{stem}.hmm.h3m",
            "ref/hmm/{stem}.hmm.h3p"
    input: "ref/hmm/{stem}.hmm"
    shell:
        "hmmpress {input}"

rule search_hmm:
    output: "data/{stem}.{hmm}-hits.hmmer-{cutoff}.tsv"
    wildcard_constraints:
        cutoff='ga|nc|tc'
    input:
        faa = "data/{stem}.fa",
        hmm = "ref/hmm/{hmm}.hmm",
        h3f = "ref/hmm/{hmm}.hmm.h3f",
        h3i = "ref/hmm/{hmm}.hmm.h3i",
        h3m = "ref/hmm/{hmm}.hmm.h3m",
        h3p = "ref/hmm/{hmm}.hmm.h3p"
    threads: min(2, MAX_THREADS)
    shell:
        """
        echo "orf_id\tmodel_id\tscore" > {output}
        hmmsearch --cut_{wildcards.cutoff} \\
                  --cpu {threads} \\
                  --tblout >(grep -v '^#' | sed 's:\s\+:\\t:g' | cut -f1,3,6 >> {output}) \\
                  {input.hmm} {input.faa} > /dev/null
        """

rule build_blast_db:
    output:
        nhr='{stem}.fn.nhr',
        nin='{stem}.fn.nin',
        nsq='{stem}.fn.nsq',
    input:
        '{stem}.fn'
    shell:
        """
        makeblastdb -dbtype nucl -parse_seqids -in {input} -out {input}
        """

rule blast_rrs_reps:
    output:
        'data/{stem}.rrs-blastn.tsv'
    input:
        mag='data/{stem}.fn',
        ref='ref/core.r.reps.fn',
        nhr='ref/core.r.reps.fn.nhr',
        nin='ref/core.r.reps.fn.nin',
        nsq='ref/core.r.reps.fn.nsq',
    shell:
        """
        blastn -subject {input.ref} -query {input.mag} -max_target_seqs 1 -outfmt 6 > {output}
        """

rule gather_hit_cds_strict:
    output:
        nucl="data/{stem}.cds.{hmm}-hits.fn",
        prot="data/{stem}.cds.{hmm}-hits.fa"
    input:
        hit_table="data/{stem}.cds.{hmm}-hits.hmmer-tc.tsv",
        nucl="data/{stem}.cds.fn",
        prot="data/{stem}.cds.fa"
    shell:
        """
        seqtk subseq {input.nucl} <(sed 1,1d {input.hit_table} | cut -f 1 | sort -u) > {output.nucl}
        seqtk subseq {input.prot} <(sed 1,1d {input.hit_table} | cut -f 1 | sort -u) > {output.prot}
        """

rule identify_rrna_seqs:
    output: "data/{stem}.16S-blastn.gff"
    input: "data/{stem}.fn"
    shell:
        """
        barrnap {input} > {output}
        """

# {{{2 Sequences Analysis

# {{{3 Domain Analysis

rule hmmscan_domains:
    output: "data/{stem}.{hmm}-hits.domtblout"
    input: "data/{stem}.{hmm}-hits.fa"
    log: "data/{stem}.{hmm}-hits.domtblout.log"
    shell:
        """
        hmmscan --domtblout {output} ref/hmm/{wildcards.hmm}.hmm {input} > {log}
        """

rule compile_domain_info:
    output: "data/{stem}.domains.tsv"
    input:
        script="scripts/hmmscan-parser.sh",
        domtbl="data/{stem}.domtblout"
    shell:
        "{input.script} {input.domtbl} > {output}"

# {{{3 Alignment

rule hmmalign:
    output: "data/{stem}.{hmm}-hits.afa"
    input: fa="data/{stem}.{hmm}-hits.fa", hmm="ref/hmm/{hmm}.hmm"
    shell:
        """
        hmmalign --informat fasta {input.hmm} {input.fa} | convert -f stockholm -t fasta > {output}
        """

rule codonalign:
    output: "data/{stem}.codonalign.afn"
    input:
        prot="data/{stem}.afa",
        nucl="data/{stem}.fn"
    shell:
        "codonalign {input.prot} {input.nucl} > {output}"

# {{{3 Filter Alignment

rule squeeze_codon_alignment:
    output: "data/{stem}.codonalign.sqz.afn"
    input:
        script="scripts/squeeze_alignment.py",
        seq="data/{stem}.codonalign.afn"
    shell: "{input.script} '-.acgtu' < {input.seq} > {output}"

rule squeeze_hmmalign_alignment:
    output: "data/{stem}.{hmm}-hits.sqz.afa"
    wildcard_constraints:
        hmm="[^.]*"
    input:
        script="scripts/squeeze_alignment.py",
        seq="data/{stem}.{hmm}-hits.afa"
    shell: "{input.script} '-.*abcdefghijklmnopqrstuvwxyz' < {input.seq} > {output}"

rule gblocks_afa:
    output: '{stem}.gb.afa'
    input:
        script='scripts/Gblocks.py',
        seq='{stem}.afa'
    shell:
        """
        scripts/Gblocks.py < {input.seq} > {output}
        """

rule refine_alignment:
    output: '{stem}.refine.afa'
    input: '{stem}.afa'
    shell: 'muscle -refine < {input} > {output}'


rule gblocks_afn:
    output: '{stem}.gb.afn'
    input: '{stem}.afn'
    shell:
        """
        Gblocks {input} -t=c -p=y -v=150 || [ $? == 1 ]
        mv {input}-gb {output}
        rm {input}-gb.htm
        """

# {{{3 Phylogenetics
rule estimate_phylogeny_afn:
    output: "data/{stem}.nucl.nwk"
    input: "data/{stem}.afn"
    shell: "FastTree -nt < {input} > {output}"

rule estimate_phylogeny_afa:
    output: "data/{stem}.prot.nwk"
    input: "data/{stem}.afa"
    shell: "FastTree < {input} > {output}"

# {{{3 Sort alignment
# TODO: Am I sure I want to use the nucleotide tree for sorting?
rule tree_sort_afn:
    output: "data/{stem}.tree-sort.afn"
    input:
        script="scripts/get_ordered_leaves.py",
        tree="data/{stem}.nucl.nwk",
        seqs="data/{stem}.afn"
    shell: "fetch_seqs --match-order <({input.script} {input.tree}) {input.seqs} > {output}"

rule tree_sort_afa:
    output: "data/{stem}.tree-sort.afa"
    input:
        script="scripts/get_ordered_leaves.py",
        tree="data/{stem}.prot.nwk",
        seqs="data/{stem}.afa"
    shell: "fetch_seqs --match-order <({input.script} {input.tree}) {input.seqs} > {output}"

# {{{2 Protein Clustering

# {{{3 De novo
rule make_diamond_db:
    output: "{stem}.fa.dmnd"
    input: "{stem}.fa"
    shell: "diamond makedb --in {input} --db {input}"

rule all_by_all_blastp:
    output: "data/{stem}.self_blastp.tsv"
    input:
        fa='data/{stem}.fa',
        db='data/{stem}.fa.dmnd',
    threads: MAX_THREADS
    shell:
        "diamond blastp --threads {threads} --db {input.db} --max-target-seqs 1000 --outfmt 6 --query {input.fa} --out {output}"

# TODO: Check out https://micans.org/mcl/man/mcl.html
rule denovo_cluster_proteins:
    output:
        clust="data/{stem}.denovo-clust.tsv",
        diss="data/{stem}.blastp_diss.tsv",
    input:
        script='scripts/cluster_proteins.py',
        data='data/{stem}.self_blastp.tsv',
    params:
        n_clusters=2000
    shell:
        '{input.script} {input.data} {params.n_clusters} {output.diss} > {output}'

rule join_genome_by_cluster_table:
    output: "data/{stem}.{clust_type}-clust.count.tsv"
    input:
        clust='data/{stem}.{clust_type}-clust.tsv',
        meta='data/{stem}.gene_genome_map.tsv'
    wildcard_constraints:
        clust_type=one_word_wc_constraint
    shell:
        """
        join -t'\t' -1 2 -2 1 <(sort -k2,2 {input.meta}) <(sort -k1,1 {input.clust}) \
                | sort -k1,1 -k3,3 -u \
                | cut -f2,3 | sort \
                | uniq -c | awk -v OFS='\t' '{{print $2,$3,$1}}' \
                > {output}
        """

rule join_genome_by_annotation_table:
    output: "data/{stem}.{annot}-annot.count.tsv"
    input:
        clust='data/{stem}.{annot}-annot.tsv',
        meta='data/{stem}.gene_genome_map.tsv'
    wildcard_constraints:
        annot=one_word_wc_constraint
    shell:
        """
        join -t'\t' -1 2 -2 1 <(sort -k2,2 {input.meta}) <(sort -k1,1 {input.clust}) \
                | sort -k1,1 -k3,3 -u \
                | cut -f2,3 | sort \
                | uniq -c | awk -v OFS='\t' '{{print $2,$3,$1}}' \
                > {output}
        """

rule domain_structure_cluster_proteins:
    output:
        "data/{stem}.domain-clust.tsv"
    input:
        script="scripts/group_by_domain_structure.py",
        domains="data/{stem}.domains.tsv"
    shell:
        "{input.script} {input.domains} > {output}"


rule compile_all_dbCAN_hit_annotations:
    output: "data/{stem}.dbCAN-hits.annot_table.tsv"
    input:
        denovo="data/{stem}.dbCAN-hits.denovo-clust.tsv",
        domain="data/{stem}.dbCAN-hits.domain-clust.tsv",
        prokka="data/{stem}.prokka-annot.tsv",
    shell:
        r"""
        # join -t '\t' <(sort -k1,1 {input.denovo}) <(join -t '\t' <(sort -k1,1 {input.domain}) <(sort -k1,1 {input.prokka})) > {output}
        tmp=$(mktemp)
        echo '
.separator "\t"

CREATE TABLE prokka
( locus_tag PRIMARY KEY
, ftype
, length_bp
, gene
, ec_number
, cog
, product
);
.import {input.prokka} prokka

CREATE TABLE domain
( locus_tag
, domain_structure
);
.import {input.domain} domain

CREATE TABLE denovo
( locus_tag
, denovo_clust
);
.import {input.denovo} denovo

SELECT
    locus_tag
  , denovo_clust
  , domain_structure
  , length_bp
  , gene
  , ec_number
  , cog
  , product
FROM prokka
JOIN denovo USING (locus_tag)
LEFT JOIN domain USING (locus_tag)
;

        ' | sqlite3 $tmp -header > {output}
"""

# {{{1 Databases

# Base database, containing static metadata.
rule generate_database_0:
    output: 'data/{group}.0.db'
    input:
        schema='schema.sql',
        mouse='data/mouse.noheader.tsv',
        sample='data/sample.noheader.tsv',
        extraction='data/extraction.noheader.tsv',
        library='data/library.noheader.tsv',
        asmbl_group='data/asmbl_group.noheader.tsv',
        rrs_taxon_count='data/{group}.r.count.noheader.tsv',
        rrs_taxonomy='data/{group}.r.taxonomy.noheader.tsv',
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
    output: 'data/{group}.1.db'
    input:
        db='data/{group}.0.db',
        contig='data/{group}.a.contigs.nlength.noheader.tsv',
        contig_bin='data/{group}.a.contigs.bins.noheader.tsv',
        contig_coverage='data/{group}.a.contigs.cvrg.noheader.tsv',
        bin_checkm='data/{group}.a.bins.checkm_details.noheader.tsv',
        contig_linkage='data/{group}.a.contigs.linkage_tally.noheader.tsv',
        checkm_merge='data/{group}.a.bins.checkm_merge_stats.noheader.tsv',
    shell:
        r"""
        tmp=$(mktemp)
        cp {input.db} $tmp
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
        | sqlite3 $tmp
        mv $tmp {output}
        """

rule denormalize_database:
    output: 'data/{stem}.denorm.db'
    input:
        db='data/{stem}.db',
        script='scripts/denormalize_db.sql',
    shell:
        """
        tmp=$(mktemp)
        cp {input.db} $tmp
        cat <(echo "PRAGMA cache_size = 1000000;") {input.script} | sqlite3 $tmp
        sqlite3 $tmp "VACUUM; ANALYZE;"
        mv $tmp {output}
        """
