snakemake --cluster-config cluster.yml \
          --cluster 'qsub -v PATH="$PATH" \
                          -A {cluster.alloc} -q {cluster.queue} \
                          -l nodes=1:ppn={cluster.ppn} -l pmem={cluster.pmem}mb \
                          -l walltime={cluster.walltime} \
                          -m {cluster.mail_options} \
                          -o {cluster.stdout_file} -e {cluster.stderr_file} \
                     ' $@
