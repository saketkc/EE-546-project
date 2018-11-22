#!/bin/bash
snakemake --snakefile Snakefile.kallisto.se\
    --config config_path=configs/$1.py\
    --js $PWD/jobscript.sh\
    --printshellcmds\
    --cluster-config $PWD/cluster.yaml\
    --jobname "{rulename}.{jobid}.$1"\
    --keep-going\
    --stats $PWD/stats/$1.riboraptor.stats\
    --rerun-incomplete\
    -j 400\
    --cluster 'sbatch --nodes=1 --partition={cluster.partition} --ntasks={cluster.cores} --mem={cluster.mem} --time={cluster.time} -o {cluster.logout} -e {cluster.logerror}'
