#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=shared       # default "shared", if not specified
#SBATCH --time=3-00:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1               # require 1 nodes
#SBATCH --ntasks-per-node=16    # cpus per node (by default, "ntasks"="cpus")
#SBATCH --mem=128000            # RAM per node in megabytes
#SBATCH --error=ss.2.busco_check.XXXXX.job.%J.err
#SBATCH --output=ss.2.busco_check.XXXXX.job.%J.out
# Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

#######################################################
#################### STEP 0: SETUP ####################
#######################################################

##source the conda env ability to run in bash
source ~/miniconda3/etc/profile.d/conda.sh

dataset="paecilomyces"
cd ${wrkdir}/annotation/

##the genome accession/strain name given to the raw assembly and kept throughout as the name
##this will be automatically modified and saved as another submittable script
genome2=XXXXX

##############################################################################################################
#################### STEP 1: Checking BUSCO completeness of the resulting predicted genes ####################
##############################################################################################################
threads=16
cd ${wrkdir}/annotation/

##because we have no real truth dataset it is hard to evaluate the annotations
##one way is we can just see if the BUSCOs, found in the genome are also found in the annotations
##below BUSCO is run on the predicted protein dataset (with all transcripts)
##NOTE because some genes have multiple transcripts, this increases the duplicate copy number
##then take these results and the results of the previous BUSCO analyses and generate a table so it is easy to compare all genomes
conda activate busco

##run BUSCO; using the eurotiales lineage database
busco -i 3.annotate/3.gff3/${genome2}/${genome2}.proteins.fa \
-o 3.annotate/3.gff3/${genome2}/${genome2}.busco_proteins \
-l eurotiales -m prot -c ${threads}
rm -r 3.annotate/3.gff3/${genome2}/${genome2}.busco_proteins/run_eurotiales_odb10/*_output
rm -r 3.annotate/3.gff3/${genome2}/${genome2}.busco_proteins/logs
rm -r 3.annotate/3.gff3/${genome2}/${genome2}.busco_proteins/run_eurotiales_odb10/busco_sequences
##extract total, complete, fragmented and missing buscos from each dataset and calculate the proportion of each from the total
cat 1.curation/1.raw/${genome2}.busco/short_summary.specific.eurotiales_odb10.${genome2}.busco.txt | grep "Complete BUSCOs\|Fragmented BUSCOs\|Missing BUSCOs\|Total BUSCO" | awk '{line=line";"$1} END{print line}' | sed 's/^;//' | tr ';' '\t' | awk -v genome="$genome2" '{print genome"\traw\t"$4"\t"$1"\t"$1/$4"\t"$2"\t"$2/$4"\t"$3"\t"$3/$4}' >> 3.annotate/3.gff3/BUSCO_results.tsv
cat 1.curation/2.sorted/${genome2}.busco/short_summary.specific.eurotiales_odb10.${genome2}.busco.txt | grep "Complete BUSCOs\|Fragmented BUSCOs\|Missing BUSCOs\|Total BUSCO" | awk '{line=line";"$1} END{print line}' | sed 's/^;//' | tr ';' '\t' | awk -v genome="$genome2" '{print genome"\tsorted\t"$4"\t"$1"\t"$1/$4"\t"$2"\t"$2/$4"\t"$3"\t"$3/$4}' >> 3.annotate/3.gff3/BUSCO_results.tsv
cat 3.annotate/3.gff3/${genome2}/${genome2}.busco_proteins/short_summary.specific.eurotiales_odb10.${genome2}.busco_proteins.txt | grep "Complete BUSCOs\|Fragmented BUSCOs\|Missing BUSCOs\|Total BUSCO" | awk '{line=line";"$1} END{print line}' | sed 's/^;//' | tr ';' '\t' | awk -v genome="$genome2" '{print genome"\tproteins\t"$4"\t"$1"\t"$1/$4"\t"$2"\t"$2/$4"\t"$3"\t"$3/$4}' >> 3.annotate/3.gff3/BUSCO_results.tsv
conda deactivate

