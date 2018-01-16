snakemake --cluster-config cluster.yml \
          --cluster 'qsub -v PATH="$PATH" \
                          -A {cluster.alloc} -q {cluster.queue} \
                          -l nodes=1:ppn={threads} -l pmem={cluster.pmem}mb \
                          -l walltime={cluster.walltime} \
                          -m {cluster.mail_options} \
                          -o /dev/null -e /dev/null \
                     ' $@
