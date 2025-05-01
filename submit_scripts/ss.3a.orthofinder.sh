#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=shared       # default "shared", if not specified
#SBATCH --time=3-00:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1               # require 1 nodes
#SBATCH --ntasks-per-node=64    # cpus per node (by default, "ntasks"="cpus")
#SBATCH --mem=256000            # RAM per node in megabytes
#SBATCH --error=ss.3b.orthofinder.job.%J.err
#SBATCH --output=ss.3b.orthofinder.job.%J.out
# Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output


#######################################################
#################### STEP 0: SETUP ####################
#######################################################

##source the conda env ability to run in bash
source ~/miniconda3/etc/profile.d/conda.sh

dataset="paecilomyces"
wrkdir="PATH/${dataset}"
cd ${wrkdir}/


##################################################################
#################### STEP 1: Running orthofinder ####################
##################################################################
threads=64
cd starfish
mkdir orthofinder
cd orthofinder

conda activate orthofinder
orthofinder -f ../proteome/ -t ${threads} 
mv ../proteome/OrthoFinder/ orthofinder/orthofinder_results/
