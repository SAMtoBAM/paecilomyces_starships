#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=shared       # default "shared", if not specified
#SBATCH --time=3-00:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1               # require 1 nodes
#SBATCH --ntasks-per-node=8    # cpus per node (by default, "ntasks"="cpus")
#SBATCH --mem=384000            # RAM per node in megabytes
#SBATCH --error=ss.X.pezizomycotina_ncbi_db.job.%J.err
#SBATCH --output=ss.X.pezizomycotina_ncbi_db.job.%J.out
# Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output


#######################################################
#################### STEP 0: SETUP ####################
#######################################################

##source the conda env ability to run in bash
source ~/miniconda3/etc/profile.d/conda.sh

dataset="paecilomyces"
cd /scratch/saodonnell/projects/${dataset}/


#############################################################################################
#################### STEP 1: DOWNLOADING ALL NCBI PEZIZOMYCOTINA GENOMES ####################
#############################################################################################
threads=8
cd /home/saodonnell/

conda activate ncbi_datasets
datasets download genome taxon pezizomycotina --filename pezizomycotina_ncbi.zip --exclude-atypical --dehydrated --assembly-source genbank --mag exclude
unzip pezizomycotina_ncbi.zip -d pezizomycotina_ncbi
datasets rehydrate --gzip --directory pezizomycotina_ncbi/

## go in database folder
cd pezizomycotina_ncbi/

## get the assembly metadata
## convert the assembly jsonl file to a tsv (selecting only for specific fields, the jsonl has so much and many tools get confused with fields etc so have to reduce it anyway)
dataformat tsv genome --fields accession,current-accession,organism-name,organism-infraspecific-isolate,organism-infraspecific-strain,assmstats-contig-l50,assmstats-contig-n50,assmstats-gc-count,assmstats-gc-percent,assmstats-genome-coverage,assmstats-number-of-component-sequences,assmstats-number-of-contigs,assmstats-number-of-organelles,assmstats-number-of-scaffolds,assmstats-scaffold-l50,assmstats-scaffold-n50,assmstats-total-sequence-len,assmstats-total-ungapped-len,organism-tax-id,source_database,type_material-display_text,type_material-label,wgs-contigs-url,wgs-project-accession,wgs-url --inputfile ncbi_dataset/data/assembly_data_report.jsonl > pezizomycotina_ncbi.assembly_data_report.tsv
## save the jsonl in case
cp ncbi_dataset/data/assembly_data_report.jsonl pezizomycotina_ncbi.jsonl
conda deactivate


####################################################################
#################### STEP 2: FORMATTING DATASET ####################
####################################################################
threads=8
cd /home/saodonnell/

## go in database folder
cd pezizomycotina_ncbi/

## rename the headers for each assembly as '>genomeaccession_contigaccession' so with the accession we will be able to easily extract more info from the BLAST results
## combine all the assemblies into a single fasta file (now easily recognisable by the accession) so it can be used to build the blast/diamond database
rm pezizomycotina_ncbi.fa.gz
ls ncbi_dataset/data/ | grep -v json | while read accession
do
genome=$( ls ncbi_dataset/data/${accession}/*.fna.gz )
echo $genome
zcat ${genome} | awk -F " " -v accession="$accession" '{if($0 ~ ">") {print ">"accession"_"$1} else {print}}' | sed 's/_>/_/g' | gzip >> pezizomycotina_ncbi.fa.gz
rm -r ncbi_dataset/data/${accession}
done

############################################################################
#################### STEP 3: BLAST SEARCH FOR STARSHIPS ####################
############################################################################
threads=8
cd /home/saodonnell/

## go in database folder
cd pezizomycotina_ncbi/

## now use BLAST to generate a database of all the genomes and search against it using Starship sequences
## generate a blastdb with compressed files (decompress into stdin first)
## use starfish env
conda activate starfish
gunzip -c pezizomycotina_ncbi.fa.gz | makeblastdb -in - -dbtype nucl -title pezizomycotina_ncbi -out pezizomycotina_ncbi

## now search the database using starships
## first make a header for the output then search
echo "query reference pident length mismatch gapopen qstart qend sstart send evalue bitscore" | tr ' ' '\t' > starships.pezizomycotina_ncbi.blastn.raw.tsv
blastn -db pezizomycotina_ncbi -query /scratch/saodonnell/projects/${dataset}/starfish/elementFinder/paecilomyces_seq.elements.fna -outfmt 6 -num_threads ${threads} >> starships.pezizomycotina_ncbi.blastn.raw.tsv 

##basic filters to remove sequences and no aggregation
minlength="10000"
minidentity="80"

##now use the filters
minlength2=$( echo $minlength | awk '{print $0/1000}' )
tail -n+2 starships.pezizomycotina_ncbi.blastn.raw.tsv  | awk -F "\t" -v minlength="$minlength" -v minidentity="$minidentity" '{if($3 > minidentity && $4 > minlength) print}'  >> starships.pezizomycotina_ncbi.blastn.${minidentity}pid_${minlength2}kb_filt.tsv


##we can now use the raw output to stitch together alignments and averaging out the identity
##we will print the query_genome_contig query_start query_end ref_starship ref_start ref_end pident length pident*length
##FIRST: select only starships with at least one region mapping > the minimum identity and > minimum length
##SECOND: get all alignments for that starship after applying a light filter for short fragments of > the minimum identity and 1kb
##THIRD: use bedtools merge to merge any overlapping regions in the query, taking the min ref_start, max ref_end, median pident, sum length and sum pident*length
##NOTE: bedtools merge will be allowed to merge regions with a maximum of 500bp gap (-d 500)
##FORTH: recalculate the length and then divide the sum pident*length for the new pident noramlised by the length of each match as opposed to the median pident
echo "query reference pident length qstart qend sstart send" | tr ' ' '\t' > starships.pezizomycotina_ncbi.blastn.filt_agg.tsv
tail -n+2 starships.pezizomycotina_ncbi.blastn.raw.tsv | awk -F "\t" -v minlength="$minlength" -v minidentity="$minidentity" '{if($3 > minidentity && $4 > minlength) print}' | cut -f1 | sort -u | while read starship
do
cat starships.pezizomycotina_ncbi.blastn.raw.tsv  | awk -F "\t" -v minidentity="$minidentity" '{if($3 > minidentity && $4 > 1000) print}'  | awk -v starship="$starship" '{if($1 == starship) print}' | awk -F "\t" '{print $2"\t"$9"\t"$10"\t"$1"\t"$7"\t"$8"\t"$3"\t"$4"\t"($3*$4)}' | awk '{if($2 > $3) {print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} else {print}}'  | sort -k1,1 -k2,2n | bedtools merge -d 500 -c 4,5,6,7,8,9 -o distinct,min,max,median,sum,sum | awk -F "\t" '{print $0"\t"$9/($3-$2)}'
done | awk '{print $4"\t"$1"\t"($7+$10)/2"\t"$3-$2"\t"$5"\t"$6"\t"$2"\t"$3}' | awk '{if($3 > 100){print $1"\t"$2"\t100\t"$4"\t"$5"\t"$6"\t"$7"\t"$8} else {print}}' >> starships.pezizomycotina_ncbi.blastn.filt_agg.tsv 

##ISSUES WITH SHORT EMBEDDED SEQUENCES


##apply a harsher filter to the aggregated alignments
cat starships.pezizomycotina_ncbi.blastn.filt_agg.tsv  | awk -F "\t" -v minlength="$minlength" -v minidentity="$minidentity" '{if($3 > minidentity && $4 > minlength) print}'  >> starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.tsv


##now we need to add metadata to the file
##grab all the matches and add the species on the match, the length of the starship and then how much of the starship covered 
echo "query reference pident length qstart qend sstart send reference_species query_size percent_covered" | tr ' ' '\t' > starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv
tail -n+2 starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.tsv  | while read line
do
##extract the accession then use this to extract the species name
accession=$( echo $line | awk -F " " '{print $2}' | awk -F "_" '{print $1"_"$2}' )
species=$( cat pezizomycotina_ncbi.assembly_data_report.tsv | awk -v acc="$accession" '{if($1 == acc) print $3"_"$4}' )
##get the starship and its size
starship=$( echo $line | awk -F "|" '{print $1}' )
size=$( cat /scratch/saodonnell/projects/paecilomyces/starfish/paecilomyces.starships.summary.manual_mod.tsv | awk -v starship="$starship" '{if($1 == starship) print $11}'  )
##remove any with the genus and Paecilomyces of Byssochlamys, keep the species name and add the size then divide by the length of the alignment
echo "${line};${species};${size}" | tr ';' '\t'  | awk '{print $0"\t"$4/$10}'
done >> starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv

## summarise the total alignments per starship for each matching genome
echo "starship_query starship_species starship_genome ref_species ref_genome average_identity total_length starship_length proportion_covered" | tr ' ' '\t' > starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.sum_per_starship_genome.tsv
tail -n+2 starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv | cut -f1 | sort -u | while read starship
do
## set variables for the genome and species of the starship
ssgenome=$( echo $starship | awk -F "_" '{print $1}' )
ssspecies=$( cat ../../metadata.tsv | awk -F "\t" -v ssgenome="$ssgenome" '{if($6 == ssgenome) print $7}'   )
##now get all matched genomes and cycle through them
cat starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv | grep "${starship}" | cut -f2 | awk -F "_" '{print $1"_"$2}' | sort -u | while read genome
do

species=$( cat starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv | grep "${starship}" | grep "${genome}" | cut -f9 | sort -u)
## sum the match size
size=$( cat starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv | grep "${starship}" | grep "${genome}" | cut -f4 | awk '{sum=sum+$0} END{print sum}' )
## average the identity
identity=$( cat starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv | grep "${starship}" | grep "${genome}" | cut -f3 | awk '{sum=sum+$0; count++} END{print sum/count}' )

cat starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.tsv | grep "${starship}" | grep "${genome}"  | cut -f1,9,10 | sort -u | awk -v ssspecies="$ssspecies" -v ssgenome="$ssgenome" -v genome="$genome" -v size="$size" -v identity="$identity" '{print $1"\t"ssspecies"\t"ssgenome"\t"$2"\t"genome"\t"identity"\t"size"\t"$3"\t"size/$3}'
done
done >> starships.pezizomycotina_ncbi.blastn.filt_agg.${minidentity}pid_${minlength2}kb_filt.metadata_plus.sum_per_starship_genome.tsv



mkdir /scratch/saodonnell/projects/paecilomyces/starfish/pezizomycotina_BLAST
cp starships.pezizomycotina_ncbi.* /scratch/saodonnell/projects/paecilomyces/starfish/pezizomycotina_BLAST/
cp pezizomycotina_ncbi.assembly_data_report.tsv /scratch/saodonnell/projects/paecilomyces/starfish/pezizomycotina_BLAST/