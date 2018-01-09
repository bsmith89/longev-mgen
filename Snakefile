max_threads = 99
just_symlink_cmd = 'ln -rs {input} {output}'
config_file = 'config.yml'
configfile: config_file

rule make_config:
    run:
        import pandas as pd
        import yaml
        library = pd.read_table('meta/library.tsv', index_col='mgen_library_id')
        library_dict = {}
        for lib, data in library.iterrows():
            library_dict[lib] = {'r1': 'raw/mgen/' + data['file_r1'],
                                'r2': 'raw/mgen/' + data['file_r2']}
        with open(config_file, 'w') as handle:
            handle.write(yaml.dump({'library': library_dict}))

rule link_raw_reads:
    output: r1='seq/{library}.mgen.r1.fq.gz', r2='seq/{library}.mgen.r2.fq.gz'
    input:
        unpack(lambda wildcards: config['library'][wildcards.library])
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

rule bowtie_index_build:
    output:
        'seq/{stem}.1.bt2',
        'seq/{stem}.2.bt2',
        'seq/{stem}.3.bt2',
        'seq/{stem}.4.bt2',
        'seq/{stem}.rev.1.bt2',
        'seq/{stem}.rev.2.bt2'
    input: 'seq/{stem}.fn'
    threads: max_threads
    shell: "bowtie2-build --threads {threads} {input} seq/{stem}"

rule deduplicate_reads:
    output:
        r1='seq/{stem}.mgen.r1.dedup.fq.gz',
        r2='seq/{stem}.mgen.r2.dedup.fq.gz'
    input:
        script='scripts/fastuniq_wrapper.sh',
        r1='seq/{stem}.mgen.r1.fq.gz',
        r2='seq/{stem}.mgen.r2.fq.gz'
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
