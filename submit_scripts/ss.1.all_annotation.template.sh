#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=shared       # default "shared", if not specified
#SBATCH --time=3-00:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1               # require 1 nodes
#SBATCH --ntasks-per-node=32    # cpus per node (by default, "ntasks"="cpus")
#SBATCH --mem=256000            # RAM per node in megabytes
#SBATCH --error=ss.1.all_annotation.XXXXX.job.%J.err
#SBATCH --output=ss.1.all_annotation.XXXXX.job.%J.out
# Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output


#######################################################
#################### STEP 0: SETUP ####################
#######################################################

##source the conda env ability to run in bash
source ~/miniconda3/etc/profile.d/conda.sh

dataset="paecilomyces"
wrkdir="PATH/${dataset}"
cd ${wrkdir}/annotation/

##the genome accession/strain name given to the raw assembly and kept throughout as the name
##this will be automatically modified and saved as another submittable script
genome2=XXXXX


#######################################################################
#################### STEP 1: PREPPING GENOMES STEP ####################
#######################################################################
threads=32
cd ${wrkdir}/annotation/

conda activate funannotate
funannotate sort -i 1.curation/1.raw/${genome2}.fa -o 1.curation/2.sorted/${genome2}.sorted.fa --minlen 1000
conda deactivate

conda activate busco
busco -i 1.curation/1.raw/${genome2}.fa -o 1.curation/1.raw/${genome2}.busco -l eurotiales --mode genome -c ${threads}
rm -r 1.curation/1.raw/${genome2}.busco/run_eurotiales_odb10/*_output
rm -r 1.curation/1.raw/${genome2}.busco/run_eurotiales_odb10/busco_sequences
busco -i 1.curation/2.sorted/${genome2}.sorted.fa -o 1.curation/2.sorted/${genome2}.busco -l eurotiales --mode genome -c ${threads}
#rm -r 1.curation/2.sorted/${genome2}.busco/run_eurotiales_odb10
rm -r 1.curation/2.sorted/${genome2}.busco/run_eurotiales_odb10/hmmer_output
rm -r 1.curation/2.sorted/${genome2}.busco/run_eurotiales_odb10/miniprot_output
rm -r 1.curation/2.sorted/${genome2}.busco/run_eurotiales_odb10/busco_sequences/single_copy_busco_sequences/*.gff
rm -r 1.curation/2.sorted/${genome2}.busco/logs
conda deactivate

gzip 1.curation/1.raw/${genome2}.fa

########################################################################################
#################### STEP 2: EARLGREY REPEAT ANNOTATION AND MASKING ####################
########################################################################################
threads=32
cd ${wrkdir}/annotation/


conda activate earlgrey
earlGrey -g 1.curation/2.sorted/${genome2}.sorted.fa -s ${genome2} -o 1.curation/3.earlgrey/ -t ${threads} -d yes
mv 1.curation/3.earlgrey/${genome2}_EarlGrey 1.curation/3.earlgrey/${genome2}.EarlGrey
cp 1.curation/3.earlgrey/${genome2}.EarlGrey/${genome2}_summaryFiles/${genome2}.softmasked.fasta 1.curation/4.softmasked/${genome2}.softmasked.fa
conda deactivate

##earlgrey generates a low of index files and blast database files etc so we will just remove these now, leaving just the genomes
rm 1.curation/2.sorted/${genome2}.sorted.fa.*

gzip 1.curation/2.sorted/${genome2}.sorted.fa

################################################################################
#################### STEP 3: BRAKER GENE/PROTEIN ANNOTATION ####################
################################################################################
##set the thread number to be used (on the HPC we can use 64 easily)
threads=32
cd ${wrkdir}/annotation/

##add path to compleasm kit
export PATH=$PATH:/home/saodonnell/compleasm_kit

##activate env
conda activate braker3

##now we have bams for the RNAseq data, protein data and tools that were all downloaded and set up above
##now run braker3
braker.pl --fungus --gff3 \
--genome=1.curation/4.softmasked/${genome2}.softmasked.fa \
--workingdir=2.braker3/${genome2} \
--threads=16 \
--busco_lineage=eurotiales_odb10 \
--COMPLEASM_PATH=/home/saodonnell/compleasm_kit/ \
--prot_seq=/home/saodonnell/fungi_odb11/Fungi.fa \
--AUGUSTUS_CONFIG_PATH=/home/saodonnell/Augustus/config/ \
--AUGUSTUS_SCRIPTS_PATH=/home/saodonnell/Augustus/bin/ \
--GENEMARK_PATH=/home/saodonnell/ETP/bin/ \
--PROTHINT_PATH=/home/saodonnell/gmes_linux_64_4/ProtHint/bin/ \
--SRATOOLS_PATH=/home/saodonnell/sratoolkit.3.1.0-ubuntu64/bin/
##the output is stored in annotation/2.braker3/${genome2}
##the important information is essentially the predicted proteins and gff3 which are here : 2.braker3/${genome2}/braker.aa  2.braker3/${genome2}/braker.gff3
##however the output of braker has these stars at the end of the protein sequences and we just need to remove them and use this 'star_mod' file now
cat 2.braker3/${genome2}/braker.aa | sed 's/\*//g' > 2.braker3/${genome2}/braker.star_mod.faa
conda deactivate

##remove some larger folders from the output to save space
rm -r 2.braker3/${genome2}/Augustus
rm -r 2.braker3/${genome2}/bbc


###################################################################################
#################### STEP 4: FUNANNOTATE FUNCTIONAL ANNOTATION ####################
###################################################################################
threads=32
cd ${wrkdir}/annotation/
##reduce thread count for interproscan
threads2=$( echo $threads | awk '{print $0/8}' | awk -F "." '{print $1}' )
##export funannotate db path
export FUNANNOTATE_DB=/home/saodonnell/funannotate_db_ascomycota

##activate env
conda activate funannotate

##feed iprscan the braker3 proteins (with the path to the iprscan tool)
funannotate iprscan -i 2.braker3/${genome2}/braker.star_mod.faa \
-m local \
-o 3.annotate/1.iprscan/${genome2}.iprscan \
-c ${threads2} \
--iprscan_path /home/saodonnell/my_interproscan/interproscan-5.67-99.0/interproscan.sh

##feed eggnog mapper the braker3 proteins
emapper.py -i 2.braker3/${genome2}/braker.star_mod.faa \
-o 3.annotate/2.eggnog/${genome2} \
--data_dir /home/saodonnell/eggnog_db/ \
--cpu $threads

##compile and generate functional annotations using funannotate
funannotate annotate --gff 2.braker3/${genome2}/braker.gff3 \
--fasta 1.curation/4.softmasked/${genome2}.softmasked.fa \
--species "Paecilomyces" \
-o 3.annotate/3.gff3/${genome2} \
--strain ${genome2} \
--busco_db ascomycota \
--cpus ${threads2} \
--iprscan 3.annotate/1.iprscan/${genome2}.iprscan \
--eggnog 3.annotate/2.eggnog/${genome2}.emapper.annotations \
--force

##because of the terrible naming system we can just move around some files and change the names in order to make them more easy to collate etc afterwards
cp 3.annotate/3.gff3/${genome2}/annotate_results/*_${genome2}.proteins.fa 3.annotate/3.gff3/${genome2}/${genome2}.proteins.fa
cp 3.annotate/3.gff3/${genome2}/annotate_results/*_${genome2}.mrna-transcripts.fa 3.annotate/3.gff3/${genome2}/${genome2}.mrna-transcripts.fa
cp 3.annotate/3.gff3/${genome2}/annotate_results/*_${genome2}.cds-transcripts.fa 3.annotate/3.gff3/${genome2}/${genome2}.cds-transcripts.fa
cp 3.annotate/3.gff3/${genome2}/annotate_results/*_${genome2}.gff3 3.annotate/3.gff3/${genome2}/${genome2}.gff3
cp 3.annotate/3.gff3/${genome2}/annotate_results/*_${genome2}.gbk 3.annotate/3.gff3/${genome2}/${genome2}.gbk
cp 3.annotate/3.gff3/${genome2}/annotate_results/*_${genome2}.agp 3.annotate/3.gff3/${genome2}/${genome2}.agp
conda deactivate

## remove iprscan and annotation results as they take up a lot of space and are duplicated in the funannotate results folder
rm 3.annotate/1.iprscan/${genome2}.iprscan
rm 3.annotate/2.eggnog/${genome2}.*
rm -r 3.annotate/3.gff3/${genome2}/annotate_misc
