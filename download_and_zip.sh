#!/bin/bash

# Inputs:
# arg1 is the name of the "to_download.txt" file which contains single SRA ID on each line
# arg2 is the name of the ".ngc" file
# arg3 is the name of the user (e.g. eowen)

# get number of jobs to run
num_jobs=$(sed -n '$=' $1)

# Get first .sra file from $1
sra_file=$(sed '1!d' $1)

# Downloads the first fastq file(s) and
# Returns the "Job ID" from .out message
job_id=$(sbatch -p 20 --mail-user $3 --mail-type fail -J ${sra_file} --mem=32G -c 12 --wrap "echo ${sra_file} | fasterq-dump --ngc $2 --split-files ${sra_file} -p -e 12" | cut -d ' ' -f4)
job_id=$(sbatch --dependency=afterok:${job_id} -p 20 --mail-user $3 --mail-type fail -J zip --mem=8G -c 24 --wrap "pigz -9 '${sra_file}_1.fastq' | pigz -9 '${sra_file}_2.fastq'" | cut -d ' ' -f4)

# Sequentially run the remainder of the jobs
for k in $(seq 2 ${num_jobs});
    do temp="${k}"
        sra_file=$(sed "${temp}q;d" $1)
        job_id=$(sbatch --dependency=afterok:${job_id} -p 20 --mail-user $3 --mail-type fail -J ${sra_file} --mem=32G -c 12 --wrap "echo ${sra_file} | fasterq-dump --ngc $2 --split-files ${sra_file} -p -e 12" | cut -d ' ' -f4)
        job_id=$(sbatch --dependency=afterok:${job_id} -p 20 --mail-user $3 --mail-type fail -J zip --mem=8G -c 24 --wrap "pigz -9 '${sra_file}_1.fastq' | pigz -9 '${sra_file}_2.fastq'" | cut -d ' ' -f4)
    done