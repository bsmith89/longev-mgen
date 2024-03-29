# {{{1 Imports

from itertools import product
import pandas as pd
from snake.misc import alias_recipe, alias_fmt, curl_recipe, curl_unzip_recipe

# {{{1 Configuration

# {{{2 Nomenclature

one_word_wc_constraint = '[^./]+'
integer_wc_constraint = '[0-9]+'
wildcard_constraints:
    group = one_word_wc_constraint,
    library = one_word_wc_constraint,
    mag = one_word_wc_constraint,
    bin = one_word_wc_constraint,
    genomes = one_word_wc_constraint

limit_numpy_procs = \
        """
        export MKL_NUM_THREADS={threads}
        export OPENBLAS_NUM_THREADS={threads}
        export NUMEXPR_NUM_THREADS={threads}
        export OMP_NUM_THREADS={threads}
        export VECLIB_MAXIMUM_THREADS={threads}

        """

# {{{2 Project configuration

# Configure the pipeline
config_file = 'config.yaml'
configfile: config_file
# Metadata specific configurations
_library = pd.read_csv(config['_meta_library'], index_col='library_id', sep='\t')
_asmbl_group = pd.read_csv(config['_meta_asmbl_group'], sep='\t')
config['library'] = {}
for library_id, row in _library.iterrows():
    config['library'][library_id] = {}
    config['library'][library_id]['r1'] = row['file_r1']
    config['library'][library_id]['r2'] = row['file_r2']
config['asmbl_group'] = {}
for group, d in _asmbl_group.groupby('asmbl_group'):
    config['asmbl_group'][group] = list(d['library_id'])

# TODO: Build the other files from the database.
rule all:
    input:
        [ "data/core.a.mags.muri2.g.final.marker_genes.gb.prot.nwk"
        , "data/core.a.mags.muri2.g.final.cds.TIGR02013-hits.hmmer-tc.sqz.gb.prot.nwk"
        , "data/core.muri2.2.denorm.db"
        , "data/core.a.mags.muri2.g.final.genome_stats.tsv"
        , "data/core.a.mags.muri2.g.final.cds.fa"
        , "data/core.a.mags/B1.g.final.quast-align.tsv"
        , "data/core.muri2.2.query_B1_genomic_variants.tsv"
        , "data/core.muri2.2.query_putative_PUL_regions.tsv"
        ]
    shell:
        "# {input}"
localrules: all

rule submission:
    input:
        docx='build/submission.docx',
        fig1='build/figure1.pdf',
        fig2='build/figure2.pdf',
        fig3='build/figure3.pdf',
        figS1='build/figureS1.pdf',
        figS2='build/figureS2.pdf',
        tblS1='build/tableS1.tsv',
        tblS2='build/tableS2.tsv',
        gh13_tree='build/gh13_tree.nwk',
        gtdb_assign='build/gtdbtk_classify_summary.tsv',
        ipynbs=[
            'build/B1_inter_strain_comparison.ipynb.html',
            'build/annotation_ordination.ipynb.html',
            'build/genome_statistics_table.ipynb.html',
            'build/genomic_diagrams.ipynb.html',
            'build/gh13_families_to_opfs_mapping.ipynb.html',
            'build/match_2019_2020_otus.ipynb.html',
            'build/metabinning_and_strain_variation_muribaculaceae.ipynb.html',
            'build/otu_correlation_and_aca_response.ipynb.html',
            'build/phylogenetic_tree_figure.ipynb.html',
        ],

rule render_muri_comparison_figure_to_pdf:
    output:
        "fig/muri_comparison.pdf",
    input:
        "doc/static/muri_comparison.svg",
    shell:
        """
        inkscape {input} --export-filename {output}
        """

rule rename_figure1:
    output: 'build/figure1.pdf'
    input: 'fig/muri_comparison.pdf'
    shell: 'ln {input} {output}'

rule rename_figure2:
    output: 'build/figure2.pdf'
    input: 'fig/starch_related_PULs.pdf'
    shell: 'ln {input} {output}'

rule rename_figure3:
    output: 'build/figure3.pdf'
    input: 'fig/B1_strain_comparison_heatmap.pdf'
    shell: 'ln {input} {output}'

rule rename_tableS1:
    output: 'build/tableS1.tsv'
    input: 'data/core.muri2.2.query_all_genome_annotations.tsv'
    shell:
        'ln {input} {output}'
localrules: rename_tableS1

rule rename_tableS2:
    output: 'build/tableS2.tsv'
    input: 'data/core.muri2.2.query_B1_genomic_variants.tsv'
    shell:
        'ln {input} {output}'
localrules: rename_tableS2

rule rename_figureS1:
    output: 'build/figureS1.pdf'
    input: 'fig/phylogenetics_concatenated_with_accession.pdf'
    shell:
        'ln {input} {output}'
localrules: rename_figureS1
ruleorder: rename_figureS1 > build_pdf_documentation

rule rename_figureS2:
    output: 'build/figureS2.pdf'
    input: 'fig/phylogenetics_rpoB.pdf'
    shell:
        'ln {input} {output}'
localrules: rename_figureS2
ruleorder: rename_figureS2 > build_pdf_documentation

# Modify taxa names in a GH13 domain tree with 'Genome_FeatureNumber|OpfNNNNN'
rule rename_tree_modify_names:
    output: 'build/gh13_tree.nwk'
    input:
        tree='data/core.a.mags.muri2.g.final.cds.GH13-hits.hmmer-nothresh.sqz.prot.nwk',
        meta='data/core.muri2.2.query_gh13_features.tsv',
    shell:
        """
        sed --file <(awk -v OFS=':' '{{print "s",$1,$1"|"$2"|"$3,"g"}}' {input.meta}) {input.tree} > {output}
        """

rule rename_gtdb_classification_summary:
    output: 'build/gtdbtk_classify_summary.tsv'
    input: dir='data/core.a.mags.here.g.final.gtdbtk.d'
    params:
        basename='gtdbtk.bac120.summary.tsv',
    shell:
        "cp {input.dir}/{params.basename} {output}"
localrules: rename_gtdb_classification_summary


# {{{2 Tooling configuration

rule setup_git_locally:
    shell:
        """
        git update-index --skip-worktree snake/local.snake
        """
localrules: setup_git_locally

# {{{3 Local includes
include: 'snake/local.snake'

# {{{3 Sub-project includes
include: 'snake/genome_comparison.snake'
include: 'snake/cazy.snake'
include: 'snake/strain_variation.snake'
# include: 'snake/curation.snake'

# {{{2 Params

# Default params
MAX_THREADS = 999
if 'MAX_THREADS' in config:
    MAX_THREADS = config['MAX_THREADS']



# {{{1 Utility rules/recipes/templates


rule print_config:
    shell:
        '{config}'
localrules: print_config

rule drop_header:
    output: '{stem}.noheader.{ext,(tsv|csv)}'
    input: '{stem}.{ext}'
    shell: 'sed 1,1d {input} > {output}'
localrules: drop_header

rule drop_header_meta:
    output: 'data/{stem}.noheader.{ext,(tsv|csv)}'
    input: 'meta/{stem}.{ext}'
    shell: 'sed 1,1d {input} > {output}'
localrules: drop_header_meta
ruleorder: drop_header_meta > drop_header

rule start_ipython:
    threads: MAX_THREADS
    shell: limit_numpy_procs + 'ipython'

rule start_jupyter:
    threads: MAX_THREADS
    shell: limit_numpy_procs + "jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/"
localrules: start_jupyter

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
localrules: configure_git

rule display_rulegraph:
    output: "data/snake.rulegraph.pdf"
    input: "Snakefile", "snake/genome_comparison.snake"
    shell:
        """
        snakemake -n --rulegraph \
                data/core.a.mags.muri.g.rsmbl.scaffolds.pilon.ctrim.dbCAN-hits.denovo50-clust.tsv \
                data/core.a.mags.muri.g.rsmbl.scaffolds.pilon.ctrim.genome_stats.tsv \
                data/core.a.mags.muri.g.rsmbl.scaffolds.pilon.ctrim.marker_genes.gb.prot.nwk \
            | dot -Tpdf > {output}
        """
localrules: display_rulegraph

rule display_dag:
    output:
        pdf="data/snake.dag.pdf",
        dot="data/snake.dag.dot",
    input: "Snakefile", "snake/genome_comparison.snake"
    shell:
        "snakemake -n --forceall --dag data/core.a.mags.muri.dbCAN-hits.denovo50-clust.tsv | tee {output.dot} | dot -Tpdf > {output.pdf}"
localrules: display_dag

rule build_html_documentation:
    output: "build/{stem}.html"
    input:
        source="doc/{stem}.md",
        mathjax="doc/static/MathJax.js",
        bib="doc/bibliography.bib",
        csl="doc/static/style.csl",
    shell:
        """
        pandoc --from markdown --to html5 \
               --standalone --self-contained --mathjax={input.mathjax} \
               --filter pandoc-crossref --csl {input.csl} \
               --bibliography={input.bib} -s {input.source} -o {output}
        """
localrules: build_html_documentation

rule build_latex_documentation:
    output: "build/{stem}.tex"
    input:
        source="doc/{stem}.md",
        bib="doc/bibliography.bib",
        csl="doc/static/style.csl",
    shell:
        """
        pandoc --from markdown --to latex \
               --filter pandoc-crossref --csl {input.csl} \
               --bibliography={input.bib} -s {input.source} -o {output}
        """
localrules: build_pdf_documentation

rule build_pdf_documentation:
    output: "build/{stem}.pdf"
    input:
        source="doc/{stem}.md",
        bib="doc/bibliography.bib",
        csl="doc/static/style.csl",
    shell:
        """
        pandoc --from markdown --to pdf \
               --pdf-engine=pdflatex \
               --filter pandoc-crossref --csl {input.csl} \
               --bibliography={input.bib} -s {input.source} -o {output}
        """
localrules: build_pdf_documentation

rule build_docx_documentation:
    output: "build/{stem}.docx"
    input:
        source="doc/{stem}.md",
        bib="doc/bibliography.bib",
        template="doc/static/example_style.docx",
        csl="doc/static/style.csl",
    shell:
        """
        pandoc --from markdown --to docx \
               --standalone --self-contained --reference-doc {input.template} \
               --filter pandoc-crossref --csl {input.csl} \
               --bibliography={input.bib} -s {input.source} -o {output}
        """
localrules: build_docx_documentation


rule jupyter_notebook_to_html:
    output: 'build/{stem}.ipynb.html',
    input: 'nb/{stem}.ipynb',
    shell: "jupyter nbconvert --to html --stdout {input} > {output}"
localrules: jupyter_notebook_to_html


# {{{1 Downloading and linking data

# {{{2 Reference data

# {{{3 General purpose reference genomes

rule download_salask_reference:
    output: 'raw/ref/salask.fn'
    params:
        url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000356.1&rettype=fasta&retmode=text'
    shell: curl_recipe
localrules: download_salask_reference

rule alias_salask_reference:
    output: 'ref/salask.fn'
    input: 'raw/ref/salask.fn'
    shell: alias_recipe
localrules: alias_salask_reference

rule download_mouse_reference:
    output: 'raw/ref/mouse.fn'
    params:
        url='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCA_001632575.1_C3H_HeJ_v1/GCA_001632575.1_C3H_HeJ_v1_genomic.fna.gz'
    shell: curl_unzip_recipe
localrules: download_mouse_reference

rule download_sra_data:
    output: 'raw/sra/{sra_id}.fn'
    shell:
        """
        fastq-dump -Z {wildcards.sra_id} | seqtk seq -A > {output}
        """
localrules: download_sra_data

# {{{3 TIGRFAM

rule download_tigrfam:
    output: "raw/ref/TIGRFAMs_14.0_HMM.tar.gz"
    params:
        url="ftp://ftp.jcvi.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz"
    shell:
        curl_recipe
localrules: download_tigrfam

rule extract_tigrfam:
    output: "ref/hmm/TIGRFAM.hmm"
    input: "raw/ref/TIGRFAMs_14.0_HMM.tar.gz"
    shell:
        """
        tar -O -xzf {input} > {output}
        """

rule extract_single_tigrfam_model:
    output: "ref/hmm/TIGR{number}.hmm"
    input: "ref/hmm/TIGRFAM.hmm"
    shell:
        """
        hmmfetch {input} TIGR{wildcards.number} > {output}
        """

rule download_pfam:
    output: "raw/ref/Pfam-31.0.hmm"
    params:
        url="ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz"
    shell: curl_unzip_recipe
localrules: download_pfam

rule link_pfam:
    output: "ref/hmm/Pfam.hmm"
    input: "raw/ref/Pfam-31.0.hmm"
    shell: alias_recipe
localrules: link_pfam

# {{{3 Parse HMM DBs

rule format_hmm_database_descriptions:
    output: "ref/{db}.hmm.tsv"
    input: script="scripts/format_hmm_database_descriptions.py", hmm="ref/hmm/{db}.hmm"
    shell:
        """
        {input.script} {input.hmm} > {output}
        # grep '^NAME\|^DESC' {input} | sed 's:\(NAME\|DESC\)  ::' | paste - - > {output}
        """



# {{{3 Metadata for sequence processing

rule download_illumina_adapters:
    output: 'raw/ref/illumina_adapters.fn'
    params:
        url='https://raw.githubusercontent.com/vsbuffalo/scythe/master/illumina_adapters.fa'
    shell: curl_recipe
localrules: download_illumina_adapters

# {{{3 CheckM DB

rule download_checkm_reference:
    output: 'raw/ref/checkm_data_2015_01_16.tar.gz'
    params:
        url='https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz'
    shell: curl_recipe

rule extract_checkm_reference:
    output: "raw/ref/checkm_db"
    input: 'raw/ref/checkm_data_2015_01_16.tar.gz'
    shell:
        """
        tar -C {output} -xzf {input}
        """

# {{{3 GTDB

rule download_gtdbtk_db:
    output: 'raw/ref/gtdbtk_r202_data.tar.gz'
    params:
        url='https://data.gtdb.ecogenomic.org/releases/release202/202.0/auxillary_files/gtdbtk_r202_data.tar.gz'
    shell: curl_recipe

localrules: download_gtdbtk_db

rule unpack_gtdbtk_db:
    output: directory('raw/ref/release202')
    input: 'raw/ref/gtdbtk_r202_data.tar.gz'
    shell:
        """
        tar -xzvf {input}
        """

rule alias_gtdbtk_db:
    output: directory('ref/gtdbtk_db')
    input: 'raw/ref/release202'
    shell: alias_recipe

localrules: alias_gtdbtk_db


# {{{3 dbCAN

# FIXME: Update manuscript with new dbCAN version
# http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt
# and re-run downstream
# analyses.
rule download_dBCAN_hmms:
    output: "raw/ref/dbCAN.hmm"
    params:
        url="http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-old@UGA/dbCAN-fam-HMMs.txt.v6"
    shell: curl_recipe
localrules: download_dBCAN_hmms

# rule download_dbCAN_meta:
#     output: "raw/ref/dbCAN.tsv"
#     params:
#         url="http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt"
#     shell: curl_recipe
# localrules: download_dbCAN_meta

rule filter_dbCAN_hmms:
    output: "ref/hmm/dbCAN.hmm"
    input: "raw/ref/dbCAN.hmm"
    shell:
        r"""
        sed 's:^NAME  \(.*\).hmm$:NAME  \1\nDESC  hypothetical carbohydrate-active domain (\1) containing protein:' {input} > {output}
        """

rule download_dbCAN_seqs:
    output: "raw/ref/CAZyDB.07202017.fa"
    params:
        url="http://csbl.bmb.uga.edu/dbCAN/download/CAZyDB.07202017.fa"
    shell: curl_recipe
localrules: download_dbCAN_seqs

rule link_dbCAN_seqs:
    output: "ref/dbCAN.fa"
    input: "raw/ref/CAZyDB.07202017.fa"
    shell: alias_recipe
localrules: link_dbCAN_seqs

# {{{3 COG

rule download_cog_function_mapping:
    output: 'raw/ref/cognames2003-2014.tab'
    params:
        url="ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/cognames2003-2014.tab"
    shell: curl_recipe
localrules: download_cog_function_mapping

rule process_cog_function_mapping:
    output: 'ref/cog_function.tsv'
    input: 'raw/ref/cognames2003-2014.tab'
    shell: "iconv -f LATIN1 -t UTF-8 {input} | sed '1,1s:^# COG\tfunc\tname:cog_id\tfunction_categories\tfunction_name:' > {output}"

# TODO: Remove these rules; they're unecessary now that I have KEGG data.
rule download_cog_to_ko_mapping:
    output: 'raw/ref/cog_from_string7_to_ko20080319_filtered_005.txt'
    params:
        url="http://pathways2.embl.de/data/cog_from_string7_to_ko20080319_filtered_005.txt.gz"
    shell: curl_unzip_recipe
localrules: download_cog_to_ko_mapping

rule alias_cog_to_ko_mapping:
    output: 'ref/cog_to_ko.tsv'
    input: 'raw/ref/cog_from_string7_to_ko20080319_filtered_005.txt'
    shell: alias_recipe
localrules: alias_cog_to_ko_mapping

rule alias_kegg_fasta:
    output: 'ref/kegg.fa'
    input: 'raw/ref/kegg.fa'
    shell: alias_recipe
localrules: alias_kegg_fasta

# TODO: Where does raw/ref/kegg.tsv come from??
# raw/ref/kegg.fa comes from parse_kegg_files...
# The current version of raw/ref/kegg.tsv was read out of an archived copy of
# the output DB.
rule alias_kegg_tsv:
    output: 'ref/kegg.tsv'
    input: 'raw/ref/kegg.tsv'
    shell: alias_recipe
localrules: alias_kegg_tsv

# {{{3 Enzyme Commission

rule download_ec_mapping:
    output: 'raw/ref/expasy_enzyme.dat'
    params: url='ftp://ftp.expasy.org/databases/enzyme/enzyme.dat'
    shell: curl_recipe
localrules: download_ec_mapping

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
localrules: download_metacyc_pathways_page

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
localrules: download_picrust_metacyc_pathways

rule reformat_picrust_metacyc_pathways:
    output: 'ref/ec2path.picrust.tsv'
    input: 'raw/ref/ec2path.picrust_prokaryotic.tsv'
    shell:
        """
        printf '#Metacyc pathway and ec mapping file\n' > {output}
        printf '#Pathway\tEC\n' >> {output}
        sed 's/EC://' < {input} >> {output}
        """


# {{{2 Raw data

rule alias_raw_read_r1:
    output: 'data/{library}.m.r1.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r1'])
    shell: alias_recipe
localrules: alias_raw_read_r1

rule alias_raw_read_r2:
    output: 'data/{library}.m.r2.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r2'])
    shell: alias_recipe
localrules: alias_raw_read_r2

# {{{3 Import results from 16S libraries

rule query_count_info:
    output: 'data/core.r.count.tsv'
    input: script='scripts/query_longev_rrs_count.sql', db='raw/longev_rrs_results.db'
    shell: "sqlite3 -header -separator '\t' {input.db} < {input.script} > {output}"
localrules: query_count_info

rule query_taxonomy_info:
    output: 'data/core.r.taxonomy.tsv'
    input: script='scripts/query_longev_rrs_taxonomy.sql', db='raw/longev_rrs_results.db'
    shell: "sqlite3 -header -separator '\t' {input.db} < {input.script} > {output}"
localrules: query_taxonomy_info

rule alias_rrs_reps:
    output: 'ref/core.r.reps.fn'
    input: 'raw/longev_rrs_reps.fn'
    shell: alias_recipe
localrules: alias_rrs_reps

# {{{1 Metagenomic Assembly

# {{{2 Data pre-processing


rule count_library_size:
    output: 'data/{library}.m.{proc}.library_size.tsv'
    input: r1='data/{library}.m.r1.{proc}.fq.gz',
           r2='data/{library}.m.r2.{proc}.fq.gz',
    shell:
        """
        zcat {input.r1} {input.r2} \
                | awk -v OFS='\t' -v library='{wildcards.library}' \
                      'NR%4==2 {{total+=length($1)}}
                       END{{print library, total}}' \
                > {output}
        """

rule combine_library_sizes:
    output: 'data/{group}.a.{proc}.library_size.tsv'
    input:
        lambda wildcards: [f'data/{library}.m.{wildcards.proc}.library_size.tsv'
                           for library in config['asmbl_group'][wildcards.group]],
    shell:
        """
        cat {input} > {output}
        """

rule deduplicate_reads:
    output:
        r1='data/{stem}.m.r1.dedup.fq.gz',
        r2='data/{stem}.m.r2.dedup.fq.gz'
    input:
        script='scripts/fastuniq_wrapper.sh',
        r1='data/{stem}.m.r1.fq.gz',
        r2='data/{stem}.m.r2.fq.gz'
    resources:
        mem_mb=24000,
        nas_io_threads=1,
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

# Why does the existence of these aliased files result in a "ChildIOException" (or
# something like that)?
# This processing intended for mapping.
rule alias_read_processing_r1:
    output: 'data/{library}.m.r1.proc.fq.gz'
    input: 'data/{library}.m.r1.dedup.deadapt.qtrim.fq.gz'
    shell: alias_recipe
localrules: alias_read_processing_r1
rule alias_read_processing_r2:
    output: 'data/{library}.m.r2.proc.fq.gz'
    input: 'data/{library}.m.r2.dedup.deadapt.qtrim.fq.gz'
    shell: alias_recipe
localrules: alias_read_processing_r2

# {{{2 Metagenome Assembly

# {{{3 MEGAHIT

# TODO: Consider dropping the contig renaming since it's redundant with the filename.
rule assemble_mgen:
    output:
        fastg='data/{group}.a-k{k}.fg',
    input:
        r1=lambda w: [f'data/{library}.m.r1.proc.fq.gz'
                      for library in config['asmbl_group'][w.group]],
        r2=lambda w: [f'data/{library}.m.r2.proc.fq.gz'
                      for library in config['asmbl_group'][w.group]],
    wildcard_constraints:
        k=integer_wc_constraint,
    params:
        r1_list=lambda w, input: ','.join(input.r1),
        r2_list=lambda w, input: ','.join(input.r2),
        k_min=21,
        k_max=lambda w: int(w.k),
        k_step=20,
        dir='data/{group}.a-k{k}.megahit.d',
    threads: 18
    resources:
        mem_mb=int(400e3),
        mem_b=int(400e3 * 1e6),  # Should match the mem_mb resource * 1e6
        pmem=int(400e3 / 18)
    log: 'data/{group}.a-k{k}.megahit.d/log'
    shell:
        r'''
        rm -rf {params.dir}
        du -L -hsc {input.r1} {input.r2}
        megahit \
            -1 {params.r1_list} \
            -2 {params.r2_list} \
            --k-min {params.k_min} --k-max {params.k_max} --k-step {params.k_step} \
            --min-contig-len {params.k_max} \
            --out-dir {params.dir} \
            --num-cpu-threads {threads} \
            --memory {resources.mem_b} --mem-flag 2
        megahit_toolkit contig2fastg {params.k_max} {params.dir}/intermediate_contigs/k{params.k_max}.contigs.fa > {output.fastg}
        '''

rule fastg_to_gfa:
    output: '{stem}.a-k{k}.gfa'
    input: '{stem}.a-k{k}.fg'
    params:
        k=lambda w: int(w.k)
    shell: "fastg2gfa {input} | sed 's:\<NODE_\([0-9]\+\)_[^\\t]*\>:\\1:g' | sed 's:\\<0M\\>:{params.k}M:' > {output}"

# TODO: Ensure this hardcoded k is appropriate.
rule gfa_to_fn:
    output: '{stem}.a.contigs.fn'
    input: '{stem}.a-k161.gfa'
    shell:
        """
        awk '/^S/{{print ">"$2"\\n"$3}}' < {input} > {output}
        """

# {{{3 QC

rule quality_asses_assembly_with_spike:
    output: directory('data/{group}.a.quast.d')
    input: contigs='data/{group}.a.contigs.fn', ref='ref/salask.fn'
    threads: min(20, MAX_THREADS)
    params: min_contig_length=1000,
    conda: 'conda/quast.yaml'
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
    threads: 3
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
    shell: 'samtools index -b {input} {output}'


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

rule filter_sequences_by_length:
    output: 'data/{stem}.gt{length}.fn'
    input: nlength='data/{stem}.nlength.tsv', fn='data/{stem}.fn'
    wildcard_constraints:
        length=integer_wc_constraint
    params: length=lambda w: int(w.length)
    shell:
        """
        seqtk subseq {input.fn} <(awk -v thresh={params.length} '$2>=thresh {{print $1}}' {input.nlength}) > {output}
        """

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

# {{{3 Binning

rule cut_up_contigs_for_concoct_binning:
    output:
        fn='{stem}.shard-{maxshard}.fn',
        bed='{stem}.shard-{maxshard}.bed',
    input: '{stem}.fn'
    params:
        maxshard=lambda w: int(float(w.maxshard)),
        overlap=0,
    conda: 'conda/concoct.yaml'
    shell:
        """
        cut_up_fasta.py {input} -c {params.maxshard} -o {params.overlap} --merge_last -b {output.bed} > {output.fn}
        """

rule construct_concoct_coverage_table:
    output:
        'data/{group}.a.contigs.shard-{maxshard}.concoct_cvrg.tsv'
    input:
        bam=lambda wildcards: [f'data/{library}.m.{wildcards.group}-map.sort.bam'
                               for library
                               in config['asmbl_group'][wildcards.group]
                              ],
        bai=lambda wildcards: [f'data/{library}.m.{wildcards.group}-map.sort.bam.bai'
                               for library
                               in config['asmbl_group'][wildcards.group]
                              ],
        bed='data/{group}.a.contigs.shard-{maxshard}.bed',
    conda: 'conda/concoct.yaml'
    resources:
        mem_mb=int(100e3)
    shell:
        """
        concoct_coverage_table.py {input.bed} {input.bam} > {output}
        """

# TODO: Be explicit about parameters to this.
# FIXME: Why does '--converge_out' not give me any debug info in the log file?
# TODO: Consider if the concoct default of 400 clusters is appropriate.  Other defaults?
# TODO: Fix the output file names. Currently `-b X.d` produces outputs of `X.d_*`
rule concoct_cluster_contigs:
    output: directory('{stem}.concoct.d')
    input:
        fn='{stem}.fn',
        cvrg='{stem}.concoct_cvrg.tsv',
    params:
        clusters=1000,
        max_iter=1000,
    conda: 'conda/concoct.yaml'
    threads: 18
    shell:
        """
        concoct --threads {threads} --clusters {params.clusters} --iterations {params.max_iter} --converge_out --composition_file {input.fn} --coverage_file {input.cvrg} -b {output}
        """

# FIXME: 'gt1000' in recipe means that I've hard-coded the minimum contig length.
rule merge_concoct_sharded_clusters:
    output: 'data/{group}.a.contigs.shard-{maxshard}.concoct.unshard.csv'
    input:  'data/{group}.a.contigs.shard-{maxshard}.concoct.d'
    conda: 'conda/concoct.yaml'
    shell:
        """
        merge_cutup_clustering.py {input}/clustering_gt1000.csv > {output}
        """

# This reformatting saves me from having to rework any of my rules downstream
# of the clustering.
rule reformat_concoct_clusters:
    output: '{stem}.concoct.rfmt.tsv'
    input: '{stem}.concoct.csv'
    shell:
        """
        sed '1,1d' {input} | (echo 'contig_id\tcluster'; tr ',' '\t') > {output}
        """

rule alias_concoct_binning:
    output: 'data/{group}.a.contigs.cluster.tsv'
    input: 'data/{group}.a.contigs.shard-1e4.unshard.concoct.rfmt.tsv'
    shell: alias_recipe


# {{{3 Prepare input data

# # {{{3 Clustering

# TODO: Don't rename as a separate stem.
# TODO: Refine this script.
rule cluster_contigs:
    output:
        out='data/{group}.a.contigs.shard-{maxshard}.cluster.csv',
        summary='data/{group}.a.contigs.shard-{maxshard}.cluster.summary.tsv',
    input:
        script='scripts/cluster_contigs.py',
        concoct_data='data/{group}.a.contigs.shard-{maxshard}.concoct.d',
        length='data/{group}.a.contigs.shard-{maxshard}.nlength.tsv',
    log: 'log/{group}.a.contigs.shard-{maxshard}.cluster.log'
    params:
        frac=0.10,
        alpha=1,
        max_clusters=2000,
        seed=1,
        min_contig_length=1000,
    threads: 12
    resources:
        memory_mb=100000,
        pmem=100000 // 12
    shell:
        limit_numpy_procs + r"""
        {input.script} {input.concoct_data}/PCA_transformed_data_gt1000.csv \
                {input.length} \
                --min-length {params.min_contig_length} \
                --frac {params.frac} \
                --max-nbins {params.max_clusters} \
                --alpha {params.alpha} \
                --seed {params.seed} \
                --summary {output.summary} \
                > {output.out} 2> {log}
        """

# FIXME: 'gt1000' in recipe means that I've hard-coded the minimum contig length.
rule merge_sharded_clusters:
    output: 'data/{group}.a.contigs.shard-{maxshard}.cluster.unshard.csv'
    input:  'data/{group}.a.contigs.shard-{maxshard}.cluster.csv'
    conda: 'conda/concoct.yaml'
    shell:
        """
        merge_cutup_clustering.py {input} > {output}
        """

rule clusters_to_bins:
    output: 'data/{group}.a.contigs.bins.tsv'
    input: 'data/{group}.a.contigs.shard-1e4.cluster.unshard.csv'
    params:
        padding=5
    shell:
        r"""
        awk -v FS=',' -v OFS='\t' \
            'BEGIN   {{print "contig_id", "bin_id"}}
             FNR > 1 {{printf "%s\tbin%0{params.padding}d\n", $1, $2}}
            ' \
            < {input} \
            > {output}
        """

rule split_out_bins:
    output: directory('data/{group}.a.bins.d')
    input:
        script='scripts/fetch_bin.sh',
        bins='data/{group}.a.contigs.bins.tsv',
        contigs='data/{group}.a.contigs.fn',
    threads: MAX_THREADS
    shell:
        r"""
        njobs=`expr {threads} \* 4`
        rm -rf {output}
        mkdir {output}
        bins=$(sed '1,1d' {input.bins} | cut -f 2 | sort | uniq)
        parallel --progress --jobs $njobs {input.script} {input.contigs} {input.bins} {{1}} {output}/{{1}}.fn ::: $bins
        """

# TODO: refine bins (scaffolds, contig extension, splitting/merging

# {{{3 QC bins

rule run_gtdbtk_classify_on_mags:
    output:
        dir=directory('data/{stem}.gtdbtk.d')
    input:
        genomes_dir='data/{stem}.for_checkm.d',
        db_dir='ref/gtdbtk_db',
    threads: 8
    resources:
        mem_mb=250000,
    conda: 'conda/gtdbtk.yaml'
    shell:
        """
        GTDBTK_DATA_PATH={input.db_dir} \
                gtdbtk classify_wf \
                --genome_dir {input.genomes_dir} \
                --out_dir {output} -x .fn \
                --cpus {threads}
        """


rule checkm_seqs:
    output:
        dir=directory('data/{stem}.checkm.d'),
        summary='data/{stem}.checkm.tsv'
    input: 'data/{stem}.for_checkm.d'
    threads: min(20, MAX_THREADS)
    conda: 'conda/checkm.yaml'
    shell:
        r"""
        rm -rf {output.dir}
        checkm lineage_wf -x fn \
                --threads {threads} --pplacer_threads {threads} \
                --file {output.summary} --tab_table \
                {input} {output.dir}
        """

rule checkm_simple_bins:
    output:
        dir=directory('data/{stem}.a.bins.checkm.d'),
        summary='data/{stem}.a.bins.checkm.tsv'
    input: 'data/{stem}.a.bins.d'
    threads: min(20, MAX_THREADS)
    conda: 'conda/checkm.yaml'
    shell:
        r"""
        rm -rf {output.dir}
        checkm lineage_wf -x fn \
                --threads {threads} --pplacer_threads {threads} \
                --file {output.summary} --tab_table \
                {input} {output.dir}
        """

rule checkm_refinements:
    output:
        dir=directory('data/{group}.a.mags/{mag}.g.rfn_check.checkm.d'),
        summary='data/{group}.a.mags/{mag}.g.rfn_check.checkm.tsv'
    input:
        seqs=[
               'data/{group}.a.mags/{mag}.g.contigs.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-50.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-60.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-70.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-80.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-90.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-50.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-60.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-70.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-80.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-90.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-50.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-60.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-70.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-80.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-90.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-50.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-60.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-70.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-80.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-90.fn'
             ],
    threads: 6
    conda: 'conda/checkm.yaml'
    shell:
        r"""
        tmpdir=$(mktemp -d)
        ln -rst $tmpdir {input.seqs}
        rm -rf {output.dir}
        checkm lineage_wf -x fn \
                --threads {threads} --pplacer_threads {threads} \
                --file {output.summary} --tab_table \
                $tmpdir {output.dir}
        """

ruleorder: checkm_refinements > checkm_seqs

rule reformat_checkm_output:
    output: 'data/{stem}.checkm_details.tsv'
    input: 'data/{stem}.checkm.tsv'
    shell:
        """
        cut -f1,4,12-15 {input} > {output}
        """
localrules: reformat_checkm_output

rule generate_checkm_markerset:
    output:
        'data/{level}_{taxon}.ms'
    conda: 'conda/checkm.yaml'
    shell:
        'checkm taxon_set {wildcards.level} {wildcards.taxon} {output}'

# {{{3 Bins to MAGs

# TODO: Understand what field 9 in checkM output file is.
# (I _think_ it's the difference between the increased completeness and the
# increased contamination.)
rule checkm_content_merge:
    output:
        checkm_work=directory('data/{group}.a.bins.checkm_merge.d'),
        merge_stats='data/{group}.a.bins.checkm_merge_stats.tsv',
    input:
        bins='data/{group}.a.bins.d',
        markerset='data/domain_Bacteria.ms',
    threads: min(20, MAX_THREADS)
    conda: 'conda/checkm.yaml'
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

# TODO: Automate this??
# TODO: Swap checkm_merge_stats input for bin_merge_stats? (this will introduce
# a dependencies on databases and therefore schema.sql
# TODO: Write curated data out to chkpt/ dir, a location for intermeidate data
# files that SHOULD be version controlled (e.g. manual curation)
rule manual_choose_curated_mags:
    output:
        contigs='chkpt/{group}.a.mags/{mag}.g.contigs.list',
        libraries='chkpt/{group}.a.mags/{mag}.g.library.list',
        trusted='chkpt/{group}.a.mags/{mag}.g.trusted_depth.tsv'
    input:
        contig='data/{group}.a.contigs.nlength.tsv',
        contig_bin='data/{group}.a.contigs.bins.tsv',
        contig_coverage='data/{group}.a.contigs.cvrg.tsv',
        bin_checkm='data/{group}.a.bins.checkm_details.tsv',
        checkm_merge='data/{group}.a.bins.checkm_merge_stats.tsv',
    shell:
        """
        cat <<EOF
        Select the contigs that belong to {wildcards.mag} ({output.contigs})
        and the libraries from which to collect reads for refinement
        ({output.libraries}), and provide an estimate of the expected depth of
        single-copy genes ({output.trusted}).

        If these already exist and are up-to-date, this can be accomplished by
        `touch {output}`
        if they already exist and you are confident that they reflects the
        current state of `data/{group}.a.bins.*` .

        However, since one of {input} is newer, this may not be the case, or
        you may want to refine your curation.  Instead you should manually
        identify all of contigs that belong to this MAG and save their names to
        this list.

EOF
        false
        """
localrules: manual_choose_curated_mags

rule link_curated_mag:
    output:
        contigs='data/{group}.a.mags/{mag}.g.contigs.list',
        libraries='data/{group}.a.mags/{mag}.g.library.list',
        depth='data/{group}.a.mags/{mag}.g.trusted_depth.tsv'
    input:
        contigs='chkpt/{group}.a.mags/{mag}.g.contigs.list',
        libraries='chkpt/{group}.a.mags/{mag}.g.library.list',
        depth='chkpt/{group}.a.mags/{mag}.g.trusted_depth.tsv'
    shell:
        """
        ln -rs {input.contigs} {output.contigs}
        ln -rs {input.libraries} {output.libraries}
        ln -rs {input.depth} {output.depth}
        """
localrules: link_curated_mag

rule get_mag_contigs:
    output: 'data/{group}.a.mags/{bin}.g.contigs.fn'
    input:
        ids='data/{group}.a.mags/{bin}.g.contigs.list',
        seqs='data/{group}.a.contigs.fn'
    shell: 'seqtk subseq {input.seqs} {input.ids} > {output}'

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
    threads: 2
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


rule extract_library_bam_for_strain_contigs:
    output: 'data/{group}.a.mags/{library}.m.{mag}-map.contigs.sort.bam'
    input:
        bam='data/{library}.m.{group}-map.sort.bam',
        contigs='data/{group}.a.mags/{mag}.g.contigs.list'
    params:
        max_contig_length=9999999999
    shell:
        """
        tmp=$(mktemp)
        awk -v dummy_len={params.max_contig_length} \
                '{{print $1, 0, dummy_len}}' \
                {input.contigs} \
                > $tmp
        samtools view -b \
                -L $tmp \
                {input.bam} \
            > {output}
        """
ruleorder: extract_library_bam_for_strain_contigs > map_reads_to_mag


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
        r1='data/{stem}.map.r1.fq.gz',
        r2='data/{stem}.map.r2.fq.gz',
    input:
        bam='data/{stem}.map.sort.bam',
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
            | samtools fastq -c 6 -1 {output.r1} -2 {output.r2} -
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
        dir=directory('data/{group}.a.mags/{mag}.g.spades.d'),
    input:
        r1='data/{group}.a.mags/{mag}.g.contigs.map.r1.fq.gz',  # TODO: Check if dnorm helped.
        r2='data/{group}.a.mags/{mag}.g.contigs.map.r2.fq.gz'
    threads: min(15, MAX_THREADS)
    resources:
        mem_mb=200000,
        pmem=200000 // 15,
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
localrules: alias_seqs_for_pilon

rule pilon_refine:
    output: "data/{group}.a.mags/{mag}.g.{proc}.pilon.fn",
    input:
        contigs="data/{group}.a.mags/{mag}.g.{proc}.fasta",
        bam="data/{group}.a.mags/{mag}.g.{proc}.map.sort.bam",
        bai="data/{group}.a.mags/{mag}.g.{proc}.map.sort.bam.bai",
    resources:
        mem_mb=int(1e5),
        pmem=int(1e5 // 12),
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

rule estimate_mag_feature_cvrg:
    output: 'data/{group}.a.mags.annot/{mag_stem}.feature_cvrg.tsv'
    input:
        depth='data/{group}.a.mags/{mag_stem}.library-depth.tsv.gz',
        feature='data/{group}.a.mags.annot/{mag_stem}.features.tsv'
    threads: 8
    resources:
        memory_mb=256000,
        pmem=256000 // 8,
    shell:
        """
script=$(mktemp)
db=$(mktemp)
echo $db

cat <<EOF > $script
.bail on
.separator '\t'

CREATE TABLE _feature
( feature_id
, sequence_id
, left INT
, right INT
);
.import {input.feature} _feature

CREATE TABLE feature
( feature_id PRIMARY KEY
, sequence_id
, start_coord INT
, stop_coord INT
);
CREATE INDEX feature__sequence_id_idx ON feature(sequence_id);
CREATE INDEX feature__start_coord_idx ON feature(start_coord);
CREATE INDEX feature__stop_coord_idx ON feature(stop_coord);

INSERT INTO feature
SELECT
    feature_id
  , sequence_id
  , MIN(left, right) AS start_coord
  , MAX(left, right) AS stop_coord
FROM _feature
;

CREATE TABLE depth
( library_id
, sequence_id
, position INT
, tally INT
);
.import /dev/stdin depth

SELECT
    library_id
  , feature_id
  , SUM(tally) * 1.0 / (stop_coord - start_coord) AS depth
FROM depth AS d
JOIN feature AS f
    ON f.sequence_id = d.sequence_id
    AND position >= start_coord
    AND position <= stop_coord
GROUP BY library_id, feature_id
;
EOF

zcat {input.depth} | tqdm | sqlite3 -init $script $db | tqdm > {output}
rm $script
rm $db
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

# {{{3 Refinement Selection

# NOTE: Checkpoint rule. Comment out to speed up dependency analysis.
rule manual_choose_mag_refinement:
    output:
        fn="data/{group}.a.mags/{mag}.g.rfn.fn",
    input:
        seqs=[
               'data/{group}.a.mags/{mag}.g.contigs.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-50.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-60.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-70.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-80.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-90.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-50.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-60.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-70.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-80.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-90.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-50.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-60.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-70.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-80.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-90.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-50.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-60.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-70.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-80.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-90.fn'
             ],
        checkm="data/{group}.a.mags/{mag}.g.rfn_check.checkm_details.tsv",
        quast='data/{group}.a.mags/{mag}.g.quast.d',
    shell:
        """
        cat <<EOF
        Pick a canonical MAG genome.

        Already available:
        {input.seqs}

        You may want to use the following to guide your decision:
        -   Completeness/contamination ({input.checkm})
        -   Number and size of contigs/scaffolds ({input.quast})
        -   Total sequence length

        Once you have selected the final MAG:

        ln -rsf data/{wildcards.group}.a.mags/{wildcards.mag}.g.REFINEMENT.fn {output}

EOF
        false
        """
localrules: manual_choose_mag_refinement

rule rename_mag_sequences:
    output: "data/{group}.a.mags/{mag}.g.final.fn"
    input: "data/{group}.a.mags/{mag}.g.rfn.fn"
    shell:
        """
        awk -v name={wildcards.mag} -v i=1 \
                '/^>/{{print ">" name "." i; i+=1}}
                !/^>/{{print $0}}' \
                < {input} > {output}
        """
localrules: rename_mag_sequences

# {{{2 QC

# TODO: Custom (non-16S) blast db for reference finding
# see http://quast.bioinf.spbau.ru/manual.html#faq_q12
rule quality_asses_mag:
    output: directory('data/{group}.a.mags/{mag}.g.quast.d')
    input:
        asmbl=[
               'data/{group}.a.mags/{mag}.g.contigs.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-50.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-60.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-70.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-80.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.ctrim-90.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-50.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-60.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-70.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-80.fn'
             , 'data/{group}.a.mags/{mag}.g.contigs.pilon.ctrim-90.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-50.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-60.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-70.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-80.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.ctrim-90.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-50.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-60.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-70.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-80.fn'
             # , 'data/{group}.a.mags/{mag}.g.rsmbl.scaffolds.pilon.ctrim-90.fn'
             ],
    threads: min(7, MAX_THREADS)
    conda: 'conda/quast.yaml'
    shell:
        r"""
        quast.py --threads={threads} --min-contig 0 --output-dir {output} \
                -R {input.asmbl[4]} --fragmented \
                {input.asmbl}
        """

# TODO: Custom (non-16S) blast db for reference finding
# see http://quast.bioinf.spbau.ru/manual.html#faq_q12
rule quality_asses_spike_mag:
    output: directory('data/{group}.a.mags/{mag}.g.quast-spike.d')
    input:
        ref='ref/salask.fn',
        asmbl=[
               'data/{group}.a.mags/{mag}.g.contigs.fn',
               'data/{group}.a.mags/{mag}.g.scaffolds.fn',
               'data/{group}.a.mags/{mag}.g.scaffolds.pilon.fn',
               'data/{group}.a.mags/{mag}.g.scaffolds.pilon.ctrim.fn',
               ]
    threads: min(6, MAX_THREADS)
    conda: 'conda/quast.yaml'
    shell:
        r"""
        quast.py --threads={threads} --min-contig 0 --output-dir {output} \
                -R {input.ref} \
                --labels "MAG, Reassembled, Reassembled+Refined, Reassembled+Refined+CTrimmed" \
                {input.asmbl}
        """


# {{{2 Annotation

# TODO: Be more explicit than {stem} in the below:

# {{{3 Untargetted

# {{{4 Prokka

rule annotate_mag:
    output:
        fa="data/{stem}.mags.annot/{mag}.g.{proc}.cds.fa",
        fn="data/{stem}.mags.annot/{mag}.g.{proc}.cds.fn",
        gbk="data/{stem}.mags.annot/{mag}.g.{proc}.prokka-annot.gbk",
        tbl="data/{stem}.mags.annot/{mag}.g.{proc}.prokka-annot.tbl",
        tsv="data/{stem}.mags.annot/{mag}.g.{proc}.prokka-annot.tsv",
        gff="data/{stem}.mags.annot/{mag}.g.{proc}.prokka-annot.gff",
        dir=directory("data/{stem}.mags.annot/{mag}.g.{proc}.prokka-annot.d"),
    input: "data/{stem}.mags/{mag}.g.{proc}.fn"
    threads: min(10, MAX_THREADS)
    conda: 'conda/prokka.yaml'
    shell:
        r"""
        prokka --force --cpus {threads} {input} \
                --outdir {output.dir} --prefix {wildcards.mag} \
                --locustag {wildcards.mag} \
                --cdsrnaolap
        cp {output.dir}/{wildcards.mag}.faa {output.fa}
        cp {output.dir}/{wildcards.mag}.ffn {output.fn}
        cp {output.dir}/{wildcards.mag}.gbk {output.gbk}
        cp {output.dir}/{wildcards.mag}.tbl {output.tbl}
        sed '/repeat_region/d' {output.dir}/{wildcards.mag}.tsv > {output.tsv}
        cp {output.dir}/{wildcards.mag}.gff {output.gff}
        """

# TODO: Can I drop this?
rule annotate_reference_mag:
    output:
        fa="data/ref.mags.annot/{mag}.g.cds.fa",
        fn="data/ref.mags.annot/{mag}.g.cds.fn",
        gbk="data/ref.mags.annot/{mag}.g.prokka-annot.gbk",
        tbl="data/ref.mags.annot/{mag}.g.prokka-annot.tbl",
        tsv="data/ref.mags.annot/{mag}.g.prokka-annot.tsv",
        gff="data/ref.mags.annot/{mag}.g.prokka-annot.gff",
        dir=directory("data/ref.mags.annot/{mag}.g.prokka-annot.d"),
    input: "data/ref.mags/{mag}.g.fn"
    threads: min(10, MAX_THREADS)
    conda: 'conda/prokka.yaml'
    shell:
        r"""
        prokka --force --cpus {threads} {input} \
                --outdir {output.dir} --prefix {wildcards.mag} \
                --locustag {wildcards.mag} \
                --metagenome --cdsrnaolap
        cp {output.dir}/{wildcards.mag}.faa {output.fa}
        cp {output.dir}/{wildcards.mag}.ffn {output.fn}
        cp {output.dir}/{wildcards.mag}.gbk {output.gbk}
        cp {output.dir}/{wildcards.mag}.tbl {output.tbl}
        sed '/repeat_region/d' {output.dir}/{wildcards.mag}.tsv > {output.tsv}
        cp {output.dir}/{wildcards.mag}.gff {output.gff}
        """

rule extract_feature_details:
    output: 'data/{group_stem}/{strain_stem}.feature_details.tsv'
    input: 'data/{group_stem}/{strain_stem}.prokka-annot.tsv'
    shell:
        """
        cut -f1,2,3,7 {input} | sed 1,1d > {output}
        """

rule parse_feature_table:
    output: 'data/{stem}.features.tsv'
    input:
        script='scripts/parse_prokka_features.py',
        feature_table='data/{stem}.prokka-annot.tbl'
    shell:
        """
        {input.script} {input.feature_table} > {output}
        """

rule list_sequence_ids:
    output: 'data/{stem}.sequence.list'
    input: 'data/{stem}.fn'
    shell:
        """
        ls_ids < {input} > {output}
        """


rule summarize_annotation:
    output: 'data/{group_stem}.annot/{mag_stem}.prokka.summary.tsv'
    input:
        prokka='data/{group_stem}.annot/{mag_stem}.prokka-annot.tsv',
        ko="data/{group_stem}.annot/{mag_stem}.ec.tsv",
    shell:
        r"""
        echo '
CREATE TABLE prokka
( locus_tag PRIMARY KEY
, ftype
, length_bp
, gene
, ec_number
, cog
, product
);
.separator \t
.import {input.prokka} prokka
-- Drop the header
--DELETE FROM annotation WHERE annotation.locus_tag = "locus_tag"

CREATE TABLE ko
( locus_tag PRIMARY KEY
, ko_id TEXT
);
.import {input.ko} ko

CREATE VIEW annotation AS
SELECT * FROM prokka LEFT JOIN ko USING (locus_tag)
;

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
SELECT "n_with_ko", COUNT(*) FROM annotation WHERE ko_id != "";
SELECT "n_with_function", COUNT(*) FROM annotation WHERE gene != "" OR ec_number != "" OR cog != "" OR ko_id != "";
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

# {{{4 Not Prokka

rule kegg_blastp:
    output: 'data/{stem}.kegg-blastp.tsv'
    input:
        db='ref/kegg.fa.dmnd',
        seq='data/{stem}.cds.fa'
    threads: min(6, MAX_THREADS)
    shell:
        """
        diamond blastp --threads {threads} --tmpdir $TMPDIR \
                --db {input.db} --query {input.seq} \
                --max-target-seqs 1 \
                --outfmt 6 --out {output}
        """

rule parse_kegg_blastp_ko:
    output: 'data/{stem}.ko.tsv'
    input:
        script='scripts/parse_kegg_blastp.py',
        table='data/{stem}.kegg-blastp.tsv'
    params:
        max_evalue=1e-10
    shell:
        """
        awk -v OFS='\t' '$11 < {params.max_evalue} {{print $1,$2}}' {input.table} | {input.script} > {output}
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

# {{{3 Targetted

rule press_hmm:
    output: "ref/hmm/{stem}.hmm.h3f",
            "ref/hmm/{stem}.hmm.h3i",
            "ref/hmm/{stem}.hmm.h3m",
            "ref/hmm/{stem}.hmm.h3p"
    input: "ref/hmm/{stem}.hmm"
    shell:
        "hmmpress {input}"

rule search_hmm:
    output:
        tbl="data/{stem}.{hmm}-hits.hmmer-{cutoff}.tblout",
        domtbl="data/{stem}.{hmm}-hits.hmmer-{cutoff}.domtblout"
    wildcard_constraints:
        cutoff='ga|nc|tc'
    input:
        faa = "data/{stem}.cds.fa",
        hmm = "ref/hmm/{hmm}.hmm",
        h3f = "ref/hmm/{hmm}.hmm.h3f",
        h3i = "ref/hmm/{hmm}.hmm.h3i",
        h3m = "ref/hmm/{hmm}.hmm.h3m",
        h3p = "ref/hmm/{hmm}.hmm.h3p"
    threads: 2
    shell:
        r"""
        printf "orf_id\tmodel_id\tscore" > {output.tbl}
        hmmsearch --cut_{wildcards.cutoff} \
                  --cpu {threads} \
                  --tblout {output.tbl} \
                  --domtblout {output.domtbl} \
                  {input.hmm} {input.faa} > /dev/null
        """

rule search_hmm_no_threshold:
    output:
        tbl="data/{stem}.{hmm}-hits.hmmer-nothresh.tblout",
        domtbl="data/{stem}.{hmm}-hits.hmmer-nothresh.domtblout"
    input:
        faa = "data/{stem}.cds.fa",
        hmm = "ref/hmm/{hmm}.hmm",
        h3f = "ref/hmm/{hmm}.hmm.h3f",
        h3i = "ref/hmm/{hmm}.hmm.h3i",
        h3m = "ref/hmm/{hmm}.hmm.h3m",
        h3p = "ref/hmm/{hmm}.hmm.h3p"
    threads: 2
    shell:
        r"""
        printf "orf_id\tmodel_id\tscore" > {output.tbl}
        hmmsearch --cpu {threads} \
                  --tblout {output.tbl} \
                  --domtblout {output.domtbl} \
                  {input.hmm} {input.faa} > /dev/null
        """

rule parse_hmmsearch_tblout:
    output: "data/{stem}.hmmer-{cutoff}.tsv"
    input:  "data/{stem}.hmmer-{cutoff}.tblout"
    wildcard_constraints:
        cutoff='ga|nc|tc|nothresh'
    shell:
        r"""
        grep -v '^#' {input} | sed 's:\s\+:\t:g' | cut -f1,3,6 >> {output}
        """

rule parse_hmmsearch_domtblout:
    output: "data/{stem}.{hmm}-domain.tsv"
    input:
        script="scripts/parse_domtbl_to_domains.py",
        domtbl="data/{stem}.{hmm}-hits.hmmer-nothresh.domtblout",
    shell:
        """
        {input.script} {input.domtbl} > {output}
        """


rule find_minimal_domains:
    output: "data/{stem}.{hmm}-domain-best.tsv"
    input:
        script="scripts/pick_minimal_domain_set.py",
        domains="data/{stem}.{hmm}-domain.tsv",
    wildcard_constraints:
        hmm=one_word_wc_constraint
    params:
        max_overlap_frac = 0.4
    shell:
        """
        {input.script} {input.domains} {params.max_overlap_frac} > {output}
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

rule gather_hit_cds_fa:
    output:
        prot="data/{stem}.cds.{hmm}-hits.hmmer-{cutoff}.fa"
    input:
        hit_table="data/{stem}.{hmm}-hits.hmmer-{cutoff}.tsv",
        prot="data/{stem}.cds.fa"
    shell:
        """
        seqtk subseq {input.prot} <(sed 1,1d {input.hit_table} | cut -f 1 | sort -u) > {output.prot}
        """


rule alias_hit_cds_fa_strict:
    output:
        prot="data/{stem}.cds.{hmm}-hits.fa"
    input:
        hit_table="data/{stem}.{hmm}-hits.hmmer-tc.tsv",
    shell: alias_recipe
localrules: alias_hit_cds_fa_strict


rule gather_hit_cds_fn_strict:
    output:
        nucl="data/{stem}.cds.{hmm}-hits.fn",
    input:
        hit_table="data/{stem}.{hmm}-hits.hmmer-tc.tsv",
        nucl="data/{stem}.cds.fn",
    shell:
        """
        seqtk subseq {input.nucl} <(sed 1,1d {input.hit_table} | cut -f 1 | sort -u) > {output.nucl}
        """


rule identify_rrna_seqs:
    output: "data/{stem}.16S-blastn.gff"
    input: "data/{stem}.fn"
    shell:
        """
        barrnap {input} > {output}
        """

rule identify_signal_peptides:
    output: "data/{stem}.signalp-raw.txt"
    input: "data/{stem}.cds.fa"
    params:
        gram='gram-',
        cutoffs="-U 0.42 -u 0.42",  # More sensitive cutoffs because false-positives are bad.
    threads: 2
    shell:
        """
        signalp -t '{params.gram}' {params.cutoffs} {input} > {output}
        """

rule format_signalp_data:
    output:
        "data/{stem}.signalp.tsv"
    input:
        script="scripts/format_signalp_data.py",
        tsv="data/{stem}.signalp-raw.txt",
        fasta="data/{stem}.cds.fa"
    shell:
        """
        {input.script} {input.tsv} {input.fasta} > {output}
        """

rule identify_transmembrane_proteins:
    output: "data/{stem}.tmhmm-raw.txt"
    input: "data/{stem}.cds.fa"
    shell:
        """
        tmhmm {input} > {output}
        """

rule format_tmhmm_data:
    output: "data/{stem}.tmhmm.tsv"
    input: script="scripts/format_tmhmm_data.py", data="data/{stem}.tmhmm-raw.txt"
    shell: "{input.script} {input.data} > {output}"


rule identify_other_secreted_proteins:
    output: "data/{stem}.secretomep-raw.txt"
    input: "data/{stem}.cds.fa"
    shell:
        """
        echo "SecretomeP 2.0 must be run from the webserver."
        echo "Version 1.0 is not appropriate for bacterial proteins."
        false
        """

rule identify_lipoproteins:
    output: "data/{stem}.lipop-raw.txt"
    input: "data/{stem}.cds.fa"
    shell:
        """
        LipoP < {input} > {output}
        """

# Take the summary line and format into a table.
# feature_id, lipop_type, score, margin, cleavage_site, aa_at_position_plus_2
rule format_lipop_data:
    output: "data/{stem}.lipop.tsv"
    input: "data/{stem}.lipop-raw.txt"
    shell:
        """
        grep 'score=' {input} \
            | awk -v OFS='\t' \
                      '$3~/SpII/{{sub("-[0-9]+$", "", $6);
                                 print $2, $3, substr($4,7), substr($5,8), substr($6,10), substr($7,7)
                               }}' \
            > {output}
        """


# {{{2 Sequences Analysis

# {{{3 Domain Analysis

rule hmmsearch_domains:
    output: "data/{stem}.{hmm}.domtblout"
    input:
        fa="data/{stem}.cds.fa",
        hmm="ref/hmm/{hmm}.hmm"
    log: "data/{stem}.{hmm}.domtblout.log"
    threads: min(6, MAX_THREADS)
    shell:
        """
        hmmsearch --cpu {threads} --cut_nc --domtblout {output} {input.hmm} {input.fa} > {log}
        """

# {{{3 Alignment

rule hmmalign:
    output: "data/{stem}.{hmm}-hits.afa"
    input: fa="data/{stem}.{hmm}-hits.fa", hmm="ref/hmm/{hmm}.hmm"
    shell:
        """
        hmmalign --informat fasta {input.hmm} {input.fa} | convert -f stockholm -t fasta > {output}
        """

rule hmmalign_with_other_cutoff:
    output: "data/{stem}.{hmm}-hits.hmmer-{cutoff}.afa"
    input: fa="data/{stem}.{hmm}-hits.hmmer-{cutoff}.fa", hmm="ref/hmm/{hmm}.hmm"
    shell:
        """
        hmmalign --informat fasta {input.hmm} {input.fa} | seqret -sfortmat1 stockholm -osformat2 fasta -filter > {output}
        """

rule codonalign:
    output: "data/{stem}.codonalign.afn"
    input:
        prot="data/{stem}.afa",
        nucl="data/{stem}.fn"
    shell:
        "codonalign {input.prot} {input.nucl} > {output}"

rule muscle_align:
    output: "data/{stem}.muscle.afa"
    input: "data/{stem}.fa"
    shell:
        "muscle < {input} > {output}"

# {{{3 Filter Alignment

rule squeeze_codon_alignment:
    output: "data/{stem}.codonalign.sqz.afn"
    input:
        script="scripts/squeeze_alignment.py",
        seq="data/{stem}.codonalign.afn"
    shell: "{input.script} '-.acgtu' < {input.seq} > {output}"

rule squeeze_hmmalign_alignment:
    output: "data/{stem}.{hmm}-hits.hmmer-{cutoff}.sqz.afa"
    wildcard_constraints:
        hmm="[^.]*"
    input:
        script="scripts/squeeze_alignment.py",
        seq="data/{stem}.{hmm}-hits.hmmer-{cutoff}.afa"
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

rule cdhit_dereplicate_proteins:
    output:
        fasta="{stem}.cds.derep.fa",
        mapping="{stem}.cds.derep.fa.clstr",
    input: "{stem}.cds.fa"
    params:
        cutoff=0.99,
    threads: 4
    shell:
        """
        cd-hit -T {threads} -i {input} -o {output.fasta} -c {params.cutoff} -d 0
        """

rule cdhit_to_map:
    output: '{stem}.cds.derep.tsv'
    input:
        script='scripts/cdhit_to_cluster_mapping.py',
        clstr='{stem}.cds.derep.fa.clstr',
    shell:
        """
        {input.script} {input.clstr} > {output}
        """

rule make_diamond_db:
    output: "{stem}.fa.dmnd"
    input: "{stem}.fa"
    shell: "diamond makedb --tmpdir $TMPDIR --in {input} --db {input}"

rule all_by_all_blastp:
    output: "data/{stem}.self-blastp.tsv"
    input:
        fa='data/{stem}.fa',
        db='data/{stem}.fa.dmnd',
    threads: 12
    params:
        # Should probably be substantially greater than the effective number of
        # genomes being compared.
        # This is the maximum degree of any node in the blast-network.
        max_target_seqs=1000
    shell:
        """
        diamond blastp --threads {threads} --tmpdir $TMPDIR \
                --db {input.db} --query {input.fa} \
                --max-target-seqs {params.max_target_seqs} \
                --outfmt 6 --out {output}
        """

# TODO: Fix this to deal with overlapping or multiple hits and normalize to self-blast.
# Current version just takes the best hit for any pair of sequences.
rule transform_blastp_to_similarity:
    output: '{stem}.protsim.tsv'
    input:
        script='scripts/transform_blastp_to_similarity.py',
        blastp='{stem}.self-blastp.tsv'
    threads: 4
    shell:
        """
        {input.script} {threads} {input.blastp} > {output}
        """

rule clip_protein_similarity_graph:
    output: '{stem}.protsim.clip.tsv'
    input: '{stem}.protsim.tsv'
    params:
        min_score = 0.2
    shell:
        """
        awk '$3>{params.min_score}' {input} > {output}
        """

# TODO: Check out https://micans.org/mcl/man/mcl.html
rule denovo_cluster_proteins:
    output: "data/{stem}.protclusts-{infl}.mcl",
    input: 'data/{stem}.protsim.clip.tsv',
    params:
        inflation=lambda wildcards: int(wildcards.infl) / 10
    shell:
        """
        mcl {input} --abc -I {params.inflation} -o {output}
        """

rule mcl_to_opf_mapping:
    output: 'data/{stem}.denovo{infl}-clust.tsv',
    input:
        script='scripts/mcl_output_to_map.py',
        mcl='data/{stem}.protclusts-{infl}.mcl',
    shell:
        """
        {input.script} {input.mcl} > {output}
        """

rule underep_clust_mapping:
    output: 'data/{stemA}.derep.{stemB}.underep.tsv'
    input:
        derep='data/{stemA}.derep.tsv',
        clust='data/{stemA}.derep.{stemB}.tsv',
    shell:
        """
        join -1 2 -2 1 \
                <(sort -k2,2 -k1,1 {input.derep}) \
                <(sort -k1,1 -k2,2 {input.clust}) \
            | awk -v OFS='\t' '{{print $2,$3}}' \
            > {output}
        """

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

# TODO: Rename script to format_domain_structure.py
rule architecture_annotate_proteins:
    output:
        "data/{stem}.Pfam-architecture.tsv"
    input:
        script="scripts/group_by_domain_structure.py",
        domains="data/{stem}.Pfam-domain-best.tsv"
    shell:
        "{input.script} {input.domains} > {output}"


# {{{1 DBs

# Base database, containing static metadata.
rule generate_database_0:
    output: 'data/{group}.0.db'
    input:
        schema='schema.0.sql',
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

# First iteration of results db; this might be used to e.g. find scaffolds.
rule generate_database_1:
    output: 'data/{group}.1.db'
    input:
        db='data/{group}.0.db',
        schema='schema.1.sql',
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
.read {input.schema}
.separator \t
.import {input.contig} contig
.import {input.contig_coverage} contig_coverage
.import {input.contig_bin} contig_bin
.import {input.bin_checkm} bin_checkm
.import {input.contig_linkage} _contig_linkage
.import {input.checkm_merge} _bin_complementarity
ANALYZE;
             ' \
        | sqlite3 $tmp
        mv $tmp {output}
        """

# TODO: Do I need to denormalize more tables?  Fewer?
# TODO: Do I need to add any indices?
rule denormalize_database_1:
    output: 'data/{stem}.1.denorm.db'
    input:
        db='data/{stem}.1.db',
    shell:
        """
        tmp=$(mktemp --dry-run)
        echo $tmp
        cp {input.db} $tmp
        ls -l $tmp
        echo '
PRAGMA legacy_alter_table = 1;
PRAGMA cache_size = 1000000;

CREATE TABLE __bin_coverage AS SELECT * FROM bin_coverage;
DROP VIEW bin_coverage;
ALTER TABLE __bin_coverage RENAME TO bin_coverage;

CREATE TABLE __library_total_nucleotides_mapping AS SELECT * FROM library_total_nucleotides_mapping;
DROP VIEW library_total_nucleotides_mapping;
ALTER TABLE __library_total_nucleotides_mapping RENAME TO library_total_nucleotides_mapping;

VACUUM; ANALYZE;
        ' | sqlite3 -bail -echo $tmp
        mv $tmp {output}
        """

rule generate_database_2:
    output: 'data/{group}.{genomes}.2.db'
    input:
        db='data/{group}.0.db',
        schema='schema.2.sql',
        genome='data/genome.noheader.tsv',
        genome_group='data/genome_group.noheader.tsv',
        library_size='data/{group}.a.proc.library_size.tsv',
        checkm='data/{group}.a.mags.{genomes}.g.final.checkm_details.noheader.tsv',
        quast='data/{group}.a.mags.{genomes}.g.final.quast.noheader.tsv',
        sequence='data/{group}.a.mags.{genomes}.g.final.sequence_to_genome.tsv',
        sequence_length='data/{group}.a.mags.{genomes}.g.final.nlength.noheader.tsv',
        ko='ref/kegg.noheader.tsv',
        cog='ref/cog_function.noheader.tsv',
        pfam_domain='ref/Pfam.hmm.tsv',
        cazy_domain='ref/dbCAN.hmm.tsv',
        tigr_domain='ref/TIGRFAM.hmm.tsv',
        feature='data/{group}.a.mags.{genomes}.g.final.features.tsv',
        feature_details='data/{group}.a.mags.{genomes}.g.final.feature_details.tsv',
        feature_x_ko='data/{group}.a.mags.{genomes}.g.final.ko-annot.tsv',
        feature_to_cog='data/{group}.a.mags.{genomes}.g.final.cog-annot.tsv',
        feature_to_opf='data/{group}.a.mags.{genomes}.g.final.cds.derep.denovo50-clust.underep.tsv',
        feature_pfam_domain='data/{group}.a.mags.{genomes}.g.final.Pfam-domain-annot.tsv',
        feature_cazy_domain='data/{group}.a.mags.{genomes}.g.final.dbCAN-domain-annot.tsv',
        feature_cazy_min_domain='data/{group}.a.mags.{genomes}.g.final.dbCAN-domain-best-annot.tsv',
        feature_tigr_domain='data/{group}.a.mags.{genomes}.g.final.TIGRFAM-domain-annot.tsv',
        feature_to_architecture='data/{group}.a.mags.{genomes}.g.final.Pfam-architecture-annot.tsv',
        signal_peptide='data/{group}.a.mags.{genomes}.g.final.signalp-annot.tsv',
        feature_tmhmm='data/{group}.a.mags.{genomes}.g.final.tmhmm-annot.tsv',
        feature_lipop='data/{group}.a.mags.{genomes}.g.final.lipop-annot.tsv',
        feature_library_cvrg='data/{group}.a.mags.{genomes}.g.final.feature_cvrg.tsv',
        variant_cross_cvrg='data/core.a.mags.annot/B1.g.final.cvrg-ratio.tsv',
    shell:
        r"""
        tmp=$(mktemp -u)
        cp {input.db} $tmp
        echo '
.bail ON
PRAGMA cache_size = 1000000;
PRAGMA foreign_keys = TRUE;
.read {input.schema}
.separator \t
.import {input.genome} genome
.import {input.genome_group} genome_group
.import {input.library_size} library_size
.import {input.checkm} checkm
.import {input.quast} quast
.import {input.sequence} _sequence
.import {input.sequence_length} _sequence_length
.import {input.ko} ko
.import {input.cog} cog
.import {input.pfam_domain} pfam_domain
.import {input.cazy_domain} cazy_domain
.import {input.tigr_domain} tigr_domain
.import {input.feature} feature
.import {input.feature_details} _feature_details
.import {input.feature_x_ko} feature_x_ko
.import {input.feature_to_cog} feature_to_cog
.import {input.feature_to_opf} feature_to_opf
.import {input.feature_pfam_domain} feature_x_pfam_domain
.import {input.feature_cazy_domain} feature_x_cazy_domain
.import {input.feature_cazy_min_domain} feature_x_cazy_minimal_domain
.import {input.feature_tigr_domain} feature_x_tigr_domain
.import {input.feature_to_architecture} feature_to_architecture
.import {input.signal_peptide} feature_signal_peptide
.import {input.feature_tmhmm} feature_tmh
.import {input.feature_lipop} feature_lipop
.import {input.feature_library_cvrg} feature_library_coverage
.import {input.variant_cross_cvrg} variant_cross_coverage
ANALYZE;
             ' \
        | sqlite3 $tmp
        mv $tmp {output}
        """

rule denormalize_database_2:
    output: 'data/{stem}.2.denorm.db'
    input:
        db='data/{stem}.2.db',
    shell:
        """
        tmp=$(mktemp)
        cp {input.db} $tmp
        echo '
.bail on
PRAGMA legacy_alter_table = 1;
PRAGMA cache_size = 1000000;

CREATE TABLE __feature_details AS SELECT * FROM feature_details;
DROP VIEW feature_details;
ALTER TABLE __feature_details RENAME TO feature_details;

DROP TABLE feature;
CREATE VIEW feature AS
SELECT feature_id, sequence_id, feature_start, feature_stop FROM feature_details
;
DROP TABLE _feature_details;
CREATE VIEW _feature_details AS
SELECT feature_id, ftype, nlength, product_description FROM feature_details
;
DROP VIEW feature_localization;
CREATE VIEW feature_localization AS
SELECT
    feature_id, localization, lp_score, lp_margin
  , sp_score, closest_cysteine, tmhelix_count
FROM feature_details
;
DROP TABLE feature_to_opf;
CREATE VIEW feature_to_opf AS
SELECT feature_id, opf_id FROM feature_details
;
DROP TABLE feature_to_architecture;
CREATE VIEW feature_to_architecture AS
SELECT feature_id, architecture FROM feature_details
;
DROP TABLE feature_to_cog;
CREATE VIEW feature_to_cog AS
SELECT feature_id, cog_id FROM feature_details
;

CREATE TABLE __feature_possible_ko AS SELECT * FROM feature_possible_ko;
DROP VIEW feature_possible_ko;
ALTER TABLE __feature_possible_ko RENAME TO feature_possible_ko;

CREATE TABLE __pul_susC AS SELECT * FROM pul_susC;
DROP VIEW pul_susC;
ALTER TABLE __pul_susC RENAME TO pul_susC;

CREATE TABLE __closest_PUL_susC AS SELECT * FROM closest_PUL_susC;
DROP VIEW closest_PUL_susC;
ALTER TABLE __closest_PUL_susC RENAME TO closest_PUL_susC;

-- TODO: Replace the lost indices, where necessary

VACUUM; ANALYZE;
        ' | sqlite3 -echo $tmp
        mv $tmp {output}
        """

rule run_db0_query:
    output: 'data/{db}.0.query_{query}.tsv'
    input: db='data/{db}.0.db', query='scripts/queries/0/{query}.sql'
    shell:
        """
        sqlite3 -header -separator '	' {input.db} < {input.query} > {output}
        """
localrules: run_db0_query

rule run_db1_query:
    output: 'data/{db}.1.query_{query}.tsv'
    input: db='data/{db}.1.denorm.db', query='scripts/queries/1/{query}.sql'
    shell:
        """
        sqlite3 -header -separator '	' {input.db} < {input.query} > {output}
        """
localrules: run_db1_query

rule run_db2_query:
    output: 'data/{db}.2.query_{query}.tsv'
    input: db='data/{db}.2.denorm.db', query='scripts/queries/2/{query}.sql'
    shell:
        """
        sqlite3 -header -separator '	' {input.db} < {input.query} > {output}
        """
