snakemake --jobs 1 --cluster 'qsub -A schmidti_fluxm -q fluxm -l nodes=1:ppn=1,pmem=24gb,walltime=00:10:00:00,qos=flux -v PATH="$PATH"' $@
