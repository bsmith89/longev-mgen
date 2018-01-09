snakemake --jobs 1 --cluster 'qsub -A schmidti_fluxm -q fluxm -l {qsub_resources} -v PATH="$PATH"' $@
