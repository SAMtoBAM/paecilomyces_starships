# paecilomyces_starships
## Describing the Starship repertoire of the Paecilomyces genus using 81 (mostly public in 2024) genomes

This reposistory is for the data analysis done for XXXXX et al. 20XX in order to describe all Paecilomyces Starships <br/>
The dataset uses mostly public available genomes (at the moment) <br/>
There are only a handful of other required inputs available in this repository <br/>
  1. A metadata file containing information on the set of genomes used AND those that were not used (and why) <br/>
  2. Genomes currently unavailable through a public repository (just using ncbi here but they are not available anywhere)

This analysis could techinically be reproduced on any genus; just need to modify the clade downloaded from ncbi <br/>
However in this analysis we did not directly download all Paecilomyces genomes (although the steps are indicated below) as we had further information as to other genomes not named Paecilomyces in ncbi <br/>
Hence here we use the metadata file with a curated list

Most computationally expensive steps were run on an HPC with the SLURM workload manager <br/>
The scripts are provided in XXXXXX and indicated as to what they will perform

#### Step 0: Set up and download genomes
Variables to manually set: <br/>
  $wrkdir = where all files and directories etc will be created; an absolute path is best especially when dealing with many HPC configurations
  $dataset = a prefix for the output (here set as 'paecilomyces')
```
######################################################################
######################### STEP 0: SETUP ##############################
######################################################################
#### Need to set up directories and assemblies for the dataset
#### Also assigning some variables

## move to the working directory already created
dataset="paecilomyces"
wrkdir="/scratch/saodonnell/projects/${dataset}/"
cd $wrkdir

## create the directory structure for analysis and annotation
mkdir annotation
cd annotation
mkdir 1.curation
mkdir 1.curation/1.raw
mkdir 1.curation/2.sorted
mkdir 1.curation/3.earlgrey
mkdir 1.curation/4.softmasked
mkdir 2.braker3
mkdir 3.annotate
mkdir 3.annotate/1.iprscan
mkdir 3.annotate/2.eggnog
mkdir 3.annotate/3.gff3

## create directory for the submit scripts
mkdir submit_scripts

######################## ASSEMBLIES ########################

###### DID NOT PERFORM BUT COULD BE DONE (note above) ######
##will download genomes using the ncbi toolkit; create a conda env for it
#mamba create -n ncbi_datasets
#conda activate ncbi_datasets
#mamba install -c conda-forge ncbi-datasets-cli
## download genomes on ncbi assigned the Paecilomyces genus, removing those that look strange according to ncbi or were assembled from metagenomic data
conda activate ncbi_datasets
datasets download genome taxon paecilomyces --filename paecilomyces_ncbi.zip --exclude-atypical --dehydrated --assembly-source genbank --mag exclude
unzip paecilomyces_ncbi.zip -d paecilomyces_ncbi
datasets rehydrate --gzip --directory paecilomyces_ncbi/
##extract metadata for the genomes just in case
dataformat tsv genome --fields accession,current-accession,organism-name,organism-infraspecific-isolate,organism-infraspecific-strain,assmstats-contig-l50,assmstats-contig-n50,assmstats-gc-count,assmstats-gc-percent,assmstats-genome-coverage,assmstats-number-of-component-sequences,assmstats-number-of-contigs,assmstats-number-of-organelles,assmstats-number-of-scaffolds,assmstats-scaffold-l50,assmstats-scaffold-n50,assmstats-total-sequence-len,assmstats-total-ungapped-len,organism-tax-id,source_database,type_material-display_text,type_material-label,wgs-contigs-url,wgs-project-accession,wgs-url --inputfile paecilomyces_ncbi/ncbi_dataset/data/assembly_data_report.jsonl > paecilomyces_ncbi.assembly_data_report.tsv
conda deactivate
##move to the directory genomes and rename the files with just the public genome accession (without the underscore)
mkdir genomes
ls paecilomyces_ncbi/ncbi_dataset/data/ | grep -v json | while read genome
do
genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
cat paecilomyces_ncbi/ncbi_dataset/data/$genome/$genome*.fna | sed "s/>/>${genome2}\_/g" > genomes/$genome2.fa
done

############################################################

## the assemblies have been given manually and therefore needed to be placed in the 1.curation/1.raw/ folder
## the naming system should be simply 'accession.fa' or 'strain.fa'
## all previously softmasked nucleotides have been made uppercase also ( awk '{if($0 ~ ">") {print $0} else {print toupper($0)}}' )

## currently have 82 Paecilomyces sp. genomes (removed 1 Monascus floridanus; incorrectly identified as Paecilomyces previously but too distant to be useful)
## 29 from Viergie et al. and are not yet public
## 1 long-read assembly from Andrew and not yet public (replacing the short read assembly for the same strain)
## only 4 long-read assemblies in total (2 variotti, 1 dactylethromorphus, 1 formosus)

## UPDATE Dec 2024; all Viergie assemblies are now public. Accessions have been added to the metadata file but the names have no been changed yet...

```
