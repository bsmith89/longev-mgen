from itertools import product
import pandas as pd

# Default params
max_threads = 40

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

localrules: print_config, link_raw_reads, download_salask_reference,
    download_illumina_adapters, download_mouse_reference

rule print_config:
    shell:
        '{config}'

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

rule assemble_mgen:
    output:
        contigs='seq/{group,[^.]+}.{proc}.asmbl.fn',
        outdir='seq/{group}.{proc}.asmbl.d'
    input:
        lambda wildcards: [f'seq/{library}.mgen.{read}.{wildcards.proc}.fq.gz'
                           for library, read
                           in product(config['asmbl_group'][wildcards.group],
                                      ['r1', 'r2', 'r3'])
                          ]
    log: 'log/{group}.{proc}.asmbl.log'
    threads: max_threads
    params:
        r1=lambda wildcards: ','.join([f'seq/{library}.mgen.r1.{wildcards.proc}.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
        r2=lambda wildcards: ','.join([f'seq/{library}.mgen.r2.{wildcards.proc}.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]]),
        r3=lambda wildcards: ','.join([f'seq/{library}.mgen.r3.{wildcards.proc}.fq.gz'
                                      for library in config['asmbl_group'][wildcards.group]])
    shell:
        """
        megahit \
            -1 {params.r1} \
            -2 {params.r2} \
            -r {params.r3} \
            --k-min 21 --k-max 161 --k-step 20 \
            --out-dir {output.outdir} \
            --num-cpu-threads {threads} \
            --verbose
        sed 's:^>k:>{wildcards.group}-k:' {output.outdir}/final.contigs.fa > {output.contigs}
        cp {output.outdir}/log {log}
        """

rule make_contig_table:
    output: 'res/{stem}.asmbl.contigs.tsv'
    input: 'seq/{stem}.asmbl.fn'
    shell:
        """
        printf 'contig_id\tflag\tmulti\tlength\n' > {output}
        grep '^>' {input} \
        | sed 's:^>\(\S*\) flag=\([0123]\) multi=\(\S*\) len=\([0-9]*\)$:\\1\t\\2\t\\3\t\\4:' \
        >> {output}
        """


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
    output: 'res/{library}.mgen.{proc}.{group}.asmbl.map.bam'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    input:
        inx_1='seq/{group}.{proc}.asmbl.1.bt2',
        inx_2='seq/{group}.{proc}.asmbl.2.bt2',
        inx_3='seq/{group}.{proc}.asmbl.3.bt2',
        inx_4='seq/{group}.{proc}.asmbl.4.bt2',
        inx_rev1='seq/{group}.{proc}.asmbl.rev.1.bt2',
        inx_rev2='seq/{group}.{proc}.asmbl.rev.2.bt2',
        r1='seq/{library}.mgen.r1.{proc}.fq.gz',
        r2='seq/{library}.mgen.r2.{proc}.fq.gz'
    threads: max_threads
    shell:
        """
        bowtie2 --threads {threads} -x seq/{wildcards.group}.{wildcards.proc}.asmbl -1 {input.r1} -2 {input.r2} \
            | samtools view -bS - \
            > {output}
        """

rule calculate_mapping_depth:
    output: 'res/{stem}.depth.tsv'
    input: 'res/{stem}.map.bam'
    shadow: 'full'
    shell:
        """
        printf 'contig_id\tposition\tdepth\n' > {output}
        samtools sort {input} | samtools depth - >> {output}
        """

rule estimate_contig_cvrg:
    output: 'res/{library}.mgen.{proc}.{group}.asmbl.cvrg.tsv'
    input:
        script='scripts/estimate_contig_coverage.py',
        depth='res/{library}.mgen.{proc}.{group}.asmbl.depth.tsv',
        contigs='res/{group}.{proc}.asmbl.contigs.tsv'
    wildcard_constraints:
        library='[^.]+',
        group='[^.]+'
    shell:
        """
        {input.script} {input.depth} {input.contigs} {wildcards.library} > {output}
        """

rule combine_cvrg_tables:
    output: 'res/{group}.{proc}.cvrg.tsv'
    input:
        script='scripts/concat_tables.py',
        all_cvrg_tables=lambda wildcards: [f'res/{library}.mgen.{wildcards.proc}.cvrg.tsv'
                                           for library
                                           in config['asmbl_group'][wildcards.group]
                                          ]
    wildcard_constraints:
        group='[^.]+'
    shell:
        """
        {input.script} {input.all_cvrg_tables} > {output}
        """

