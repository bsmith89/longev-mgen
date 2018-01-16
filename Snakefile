from itertools import product
import pandas as pd

max_threads = 99
just_symlink_cmd = 'ln -rs {input} {output}'
# Default params
max_threads = 40

# Configure the pipeline
config_file = 'config.yml'
qsub_resources="nodes=1:ppn=1,pmem=24gb,walltime=00:10:00:00,qos=flux"
configfile: config_file
# Metadata specific configurations
_library = pd.read_table('meta/library.tsv', index_col='library_id')
_asmbl_group = pd.read_table('meta/asmbl_group.tsv')
config['library'] = {}
for library_id, row in _library.iterrows():
    config['library'][library_id] = {}
    config['library'][library_id]['r1'] = row['file_r1']
    config['library'][library_id]['r2'] = row['file_r2']
config['asmbl_group'] = {}
for group, d in _asmbl_group.groupby('asmbl_group'):
    config['asmbl_group'][group] = list(d['library_id'])

rule print_config:
    shell:
        '{config}'

# rule make_config:
#     run:
#         import pandas as pd
#         import yaml
#         library = pd.read_table('meta/library.tsv', index_col='mgen_library_id')
#         library_dict = {}
#         for lib, data in library.iterrows():
#             library_dict[lib] = {'r1': 'raw/mgen/' + data['file_r1'],
#                                 'r2': 'raw/mgen/' + data['file_r2']}
#         with open(config_file, 'w') as handle:
#             handle.write(yaml.dump({'library': library_dict}))

rule link_raw_reads:
    output: r1='seq/{library}.mgen.r1.fq.gz', r2='seq/{library}.mgen.r2.fq.gz'
    input:
        unpack(lambda wildcards: config['library'][wildcards.library])
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
    params: qsub_resources="nodes=1:ppn=1,pmem=24gb,walltime=00:01:00:00,qos=flux"
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
        'seq/{group}.{proc}.asmbl.fn'
        'seq/{group,[^.]+}.{proc}.asmbl.fn'
    input:
        lambda wildcards: ['seq/{library}.{read}.{{proc}}.fq.gz'
        lambda wildcards: [f'seq/{library}.mgen.{read}.{wildcards.proc}.fq.gz'
                           for library, read
                           in product(config['group'][wildcards.group],
                           in product(config['asmbl_group'][wildcards.group],
                                      ['r1', 'r2', 'r3'])
                          ]
    shell:
        """
        echo {input} > {output}
        """

+# ruleorder: deduplicate_reads > trim_adapters >  quality_trim > assemble_mgen
+
+
+# rule qsub_test:
+#     output: touch('qsub_test.{n}')
+#     log: 'log/qsub_test.{n}.log'
+#     threads: 2
+#     shell:
+#         """
+#         qstat -f $PBS_JOBID >>{log}
+#         """
+#
+# rule qsub_test_dep:
+#     output: touch('qsub_test_dep')
+#     input: 'qsub_test.1', 'qsub_test.2'
+#     log: 'log/qsub_test_dep.log'
+#     shell:
+#         """
+#         qstat -f $PBS_JOBID >>{log}
+#         """
