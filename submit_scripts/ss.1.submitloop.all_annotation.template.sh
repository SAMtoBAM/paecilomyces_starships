#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=shared       # default "shared", if not specified
#SBATCH --time=3-00:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1               # require 1 nodes
#SBATCH --ntasks-per-node=4    # cpus per node (by default, "ntasks"="cpus")
#SBATCH --mem=32000            # RAM per node in megabytes
#SBATCH --error=submit_script.1.all_annotation.XXXXX.job.%J.err
#SBATCH --output=submit_script.1.all_annotation.XXXXX.job.%J.out
# Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

##### STEP 0: SETUP #####

##source the conda env ability to run in bash
source ~/miniconda3/etc/profile.d/conda.sh

dataset="paecilomyces"



cd /scratch/saodonnell/projects/${dataset}/annotation/
###launching a single submit script do to all jobs for each assembly
ls 1.curation/1.raw/ | grep fa$ | while read genome
do
genome2=$( echo $genome | awk -F "." '{print $1}'  )
##move into submit scripts for submission so the report files are automatically placed here
cd submit_scripts
sed "s/XXXXX/${genome2}/g"  ../../ss.1.all_annotation.template.sh > ss.1.all_annotation.${genome2}.sh
sbatch ss.1.all_annotation.${genome2}.sh
cd ..
##there are issues with running so many simultaneously
sleep 30m
done
