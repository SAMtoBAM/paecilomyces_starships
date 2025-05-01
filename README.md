# paecilomyces_starships
## Describing the <i>Starship</i> repertoire of the Paecilomyces genus using 79 (public as of 2025) genomes

This reposistory is for the data analysis done for [Urquhart et al. 2025](https://doi.org/10.1101/2025.02.28.640899
) in order to describe all <i>Paecilomyces Starships</i> <br/>
The dataset uses mostly public available genomes (at the moment) <br/>
There are only a handful of other required inputs available in this repository <br/>
  1. A metadata file containing information on the set of genomes used AND those that were not used (and why) <br/>
  2. Genomes currently unavailable through a public repository (just using ncbi here but they are not available anywhere)

This analysis could techinically be reproduced on any genus; just need to modify the clade downloaded from ncbi <br/>
However in this analysis we did not directly download all <i>Paecilomyces</i> genomes (although the steps are indicated below) as we had further information as to other genomes not named <i>Paecilomyces</i> on NCBI <br/>
Hence here we use the metadata file with a curated list

Most computationally expensive steps were run on an HPC with the SLURM workload manager <br/>
The scripts are provided in submit_scripts and indicated as to what they will perform

#### Step 0: Set up and download genomes
Variables to manually set: <br/>
  $wrkdir = where all files and directories etc will be created; an absolute path is best especially when dealing with many HPC configurations (change in all submit files) <br/>
  $dataset = a prefix for the output (here set as 'paecilomyces') (change in all submit files)
```
######################################################################
######################### STEP 0: SETUP ##############################
######################################################################
#### Need to set up directories and assemblies for the dataset
#### Also assigning some variables

## move to the working directory already created
dataset="paecilomyces"
wrkdir="/PATH/${dataset}/"
mkdir $wrkdir
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

###### DID NOT PERFORM BELOW STEPS BUT YOU COULD INSTEAD OF USING A METADATA FILE AS I DID (note above) ######
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
##move to the directory for raw assemblies and rename the files with just the public genome accession (without the underscore)
ls paecilomyces_ncbi/ncbi_dataset/data/ | grep -v json | while read genome
do
genome2=$( echo $genome | sed 's/_//' | awk -F "." '{print $1}')
zcat paecilomyces_ncbi/ncbi_dataset/data/$genome/$genome*.fna.gz | sed "s/>/>${genome2}\_/g" > 1.curation/1.raw/$genome2.fa
done
## can remove all the ncbi data now
rm paecilomyces_ncbi.zip
rm -r paecilomyces_ncbi/

###### DID NOT PERFORM ABOVE STEPS BUT YOU  COULD INSTEAD OF USING A METADATA FILE AS I DID (note above) ######

## the assemblies were given manually and therefore needed to be placed in the 1.curation/1.raw/ folder
## the naming system should be simply 'accession.fa' or 'strain.fa'
## all previously softmasked nucleotides have been made uppercase also ( awk '{if($0 ~ ">") {print $0} else {print toupper($0)}}' )

## currently have 79 Paecilomyces sp. genomes (removed 1 Monascus floridanus; incorrectly identified as Paecilomyces previously but too distant to be useful)
## 1 long-read assembly from Andrew and not yet public (replacing the short read assembly for the same strain)
## only 4 long-read assemblies in total (2 variotti, 1 dactylethromorphus, 1 formosus)

## UPDATE Dec 2024; all Viergie assemblies are now public. Accessions have been added to the metadata file but the names have no been changed yet...
```


#### Step 1: Functionally annotate whole genomes

```
###########################################################################
######################### STEP 1: ANNOTATION ##############################
###########################################################################
#### In order to easily analyse Starship cargo we can annotate the entire genome of all assemblies present
#### A submission script template has been generated to launch all annotation steps for each assembly individually within the 1.curation/1.raw/ folder


## move into annotation folder
cd ${wrkdir}/annotation/
## launching a single submit script do to all jobs for each assembly
ls 1.curation/1.raw/ | grep fa$ | while read genome
do
genome2=$( echo $genome | awk -F "." '{print $1}'  )
## move into submit scripts for submission so the report files are automatically placed here
cd submit_scripts
sed "s/XXXXX/${genome2}/g"  ../../ss.1.all_annotation.template.sh > ss.1.all_annotation.${genome2}.sh
sbatch ss.1.all_annotation.${genome2}.sh
cd ..
## there are issues with running so many simultaneously
sleep 30m
done


## sometimes there are issues with the cluster etc
## check if they all finished and relunch those that didn't
ls 1.curation/1.raw/ | grep fa.gz$ | while read genome
do
genome2=$( echo $genome | awk -F "." '{print $1}'  )
if [ -f 3.annotate/3.gff3/${genome2}/${genome2}.gff3  ]
then
echo "${genome2} contains an annotated gff3, no worries"
else
echo "${genome2} needs to be re done (if it was done at all)"
gunzip 1.curation/1.raw/${genome}
echo "removing old filtered genome"
rm 1.curation/2.sorted/${genome2}.sorted.fa.gz
echo "removing old earlgrey output if present"
if [ -d 1.curation/3.earlgrey/${genome2}.EarlGrey/  ]
then
rm -r 1.curation/3.earlgrey/${genome2}.EarlGrey/
fi
echo "removing old braker output if present"
if [ -d 2.braker3/${genome2}/  ]
then
rm -r 2.braker3/${genome2}/
fi
echo "removing old funannotate output if present"
if [ -d 3.annotate/3.gff3/${genome2}/  ]
then
rm -r 3.annotate/3.gff3/${genome2}/
fi
echo "rerunning ${genome2}"
cd submit_scripts
sed "s/XXXXX/${genome2}/g"  ../../ss.1.all_annotation.template.sh > ss.1.all_annotation.${genome2}.sh
sbatch ss.1.all_annotation.${genome2}.sh
cd ..
## there are issues with running so many simultaneously
sleep 30m
fi
done

## compress assemblies
gzip 1.curation/1.raw/*.fa
gzip 1.curation/2.sorted/*.sorted.fa


#### RENAMING THE PUBLIC ASSEMBLIES USING THEIR PUBLIC ACCESSION FOR THE CONTIGS ####

##rename all the contigs to their original public accession naming using the ragtag scaffolding tool
##only applied to the public assemblies
##the old/ not-renamed are put into a newly named file whilst the renamed are left with the normal name for the file (like the others that were not renamed)
##install seqkit in conda env
conda activate general
cd ${wrkdir}/annotation/
mkdir 4.renaming
##first identify the genomes of interest for performing this
ls 1.curation/1.raw/ | grep .fa.gz$ | grep ^GCA | while read genome
do
genome2=$( echo $genome | awk -F "." '{print $1}')
gunzip 1.curation/1.raw/${genome}
## create new files which will be progressively renamed
cp 3.annotate/3.gff3/${genome2}/${genome2}.gff3 3.annotate/3.gff3/${genome2}/${genome2}.accession_not_renamed.gff3
cp 1.curation/4.softmasked/${genome2}.softmasked.fa 3.annotate/3.gff3/${genome2}/${genome2}.softmasked.fa
cp 3.annotate/3.gff3/${genome2}/${genome2}.softmasked.fa 3.annotate/3.gff3/${genome2}/${genome2}.softmasked.accession_not_renamed.fa
## now sort these genomes similarly, get the list of contig names in that order then lines up the two lists and make the swap
mkdir ${genome2}.naming
seqkit sort 1.curation/1.raw/${genome2}.fa -s -i | seqkit sort - -b | seqkit sort - -l -r | grep '>' > ${genome2}.naming/list1.txt
seqkit sort 1.curation/4.softmasked/${genome2}.softmasked.fa -s -i | seqkit sort - -b | seqkit sort - -l -r | grep '>' > ${genome2}.naming/list2.txt
paste ${genome2}.naming/list*.txt -d "\t" > ${genome2}.naming/list.comb.txt
cat ${genome2}.naming/list.comb.txt | while read list
do
accession=$( echo "${list}" | awk -F "\t" '{print $1}' | awk -F " " '{print $1}' | sed 's/>//g' )
current=$( echo "${list}" | awk -F "\t" '{print $2}' | sed 's/>//g' )
sed -i "s/^${current}\t/${accession}\t/g" 3.annotate/3.gff3/${genome2}/${genome2}.gff3
sed -i "s/^>${current}$/>${accession}/g" 3.annotate/3.gff3/${genome2}/${genome2}.softmasked.fa
done
#rm -r ${genome2}.naming
gzip 1.curation/1.raw/${genome2}.fa
done

##place a copy of the softmasked genome for the other assemblies in the annotation folder too
ls 1.curation/1.raw/ | grep .fa.gz$ | grep -v ^GCA | while read genome
do
genome2=$( echo $genome | awk -F "." '{print $1}')
cp 1.curation/4.softmasked/${genome2}.softmasked.fa 3.annotate/3.gff3/${genome2}/${genome2}.softmasked.fa
done

mv *.naming 4.renaming

## now there a few output folders of interest
## 	filtered and masked assembly = ${wrkdir}/annotation/1.curation/4.softmasked/*.softmasked.fa
## 	gff3 = ${wrkdir}/annotation/3.annotate/3.gff3/${genome2}/${genome2}.gff3
## 	proteins = ${wrkdir}/annotation/3.annotate/3.gff3/${genome2}/${genome2}.proteins.fa

```

#### Step 2: Validate genome completeness, annotation accuracy using BUSCOs and generate a phylogeny with single copy BUSCOs to validate species ID

```
############################################################################
######################### STEP 2: BUSCO CHECK ##############################
############################################################################
#### We can use BUSCO to check the completeness of the assemblies and the annotated proteins side by side
#### A submission script template has been generated to launch each assembly and annotated protein dataset individually

## move into annotation folder
cd ${wrkdir}/annotation/
## first create a header for the summary file, then launch the BUSCO analyses
echo "genome;dataset;total_BUSCOs;complete;complete_proportion;fragmented;fragmented_proportion;missing;missing_proportion" | tr ';' '\t' > 3.annotate/3.gff3/BUSCO_results.tsv

## launching a single submit script do to all jobs for each assembly
ls 1.curation/1.raw/ | grep fa.gz$ | while read genome
do
genome2=$( echo $genome | awk -F "." '{print $1}'  )
## move into submit scripts for submission so the report files are automatically placed here
cd submit_scripts
sed "s/XXXXX/${genome2}/g"  ../../ss.2.BUSCO_check.template.sh > ss.2.BUSCO_check.${genome2}.sh
sbatch ss.2.BUSCO_check.${genome2}.sh
cd ..
done

## the output file 3.annotate/3.gff3/BUSCO_results.tsv contains a combination of the BUSCO results from the raw assembly, the filtered assembly and the proteins
## there should be no large loss of BUSCO score after the filtering and in the annotated proteins
## in fact in most cases the annotated proteins usually have equal or better BUSCO scores


```

#### Step 3: Generate <i>Starship</i> database using Starfish

```
###################################################################################
######################### STEP 3: STARSHIP DISCOVERY ##############################
###################################################################################
#### Using Starfish we can take our newly annotated genomes and look for Starships using all genomes simultneously
#### A submission script is already generate for running all the steps required for Starfish

## install starfish first
## downloaded starfish, provides databases and some auxillary tools (all directed in the github)
#cd ~
#git clone https://github.com/egluckthaler/starfish.git

## created as below using the conda env (as on the wiki)
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda config --set channel_priority strict
#mamba create -c egluckthaler -n starfish python=3.8 starfish
## create an R and gggenome env (needed for later visualiations)
#mamba create -c conda-forge -n R_env r-base r-gggenomes r-svglite

## move into the project folder
cd ${wrkdir}/
## make a directory for the Starfish analysis output
mkdir starfish

##some paths for directing scripts to tools
#export PATH=$PATH:~/starfish/bin/
#export PATH=$PATH:~/starfish/CNEFinder

############################# STEP 3a: generating orthogroups ###################################

cd ${wrkdir}/starfish
### get the orthogroups that will then be used to look at alignments between genomes and detail multiple insertion sites of a single starship
##first create an environment
#mamba create -n orthofinder -c bioconda -c conda-forge orthofinder
#conda activate orthofinder
##copy the protein files and just rename them ${assembly}.faa (orthofinder will take everything before the suffix tag as the name so can't leave 'proteins' in the file name)
##add the genome name to the beginning of the protein name tags e.g. >GCAXXXXXXXX_g1234.t1 (used downstream to know which genome it has come from)
mkdir proteome
ls ${wrkdir}/annotation/3.annotate/3.gff3/*/*.proteins.fa | while read file
do
genome2=$( echo $file | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
cat $file | sed "s/>/>${genome2}_/g" > proteome/${genome2}.faa
done

mkdir transcriptome
ls ${wrkdir}/annotation/3.annotate/3.gff3/*/*.cds-transcripts.fa | while read file
do
genome2=$( echo $file | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
cat $file | sed "s/>/>${genome2}_/g" > transcriptome/${genome2}.fa
done

##run orthofinder
#mkdir orthofinder
#cd orthofinder
#conda activate orthofinder
#orthofinder -f ../proteome/ -t ${threads} 
#mv ../proteome/OrthoFinder/ orthofinder/orthofinder_results/

## use the submit file as this takes awhile
sbatch ss.3a.orthofinder.sh

###use a species tree or not????

## the orthofinder results are now in the directory ${wrkdir}/starfish/orthofinder/orthofinder_results/Results_*/
## the wildcard being the date it was run (hopefully there is just one run)


############################ STEP: 3b: finding the Starships ############################

## starfish will be run on all genomes simultaneously therefore we just need to submit a single script
## this script, as opposed to others. sequentially identifies starship insertions using closely related genomes then more distant genomes
## this is done by reducing the percent identity minimum for alignments in a second round, running on the output of the first
## by doing this all candidate starships from the first round are kept, and those without a candidate insertion site are searched again using the reduced identity min
sbatch ss.3b.starfish.sh

######### manual annotation ############

## once this script has finished there will be candidate starships to manually label as good, mid or poor candidates, based on the visulisation
##first produce a summary file containing a few features of interest
##this summary file will contain a column to be filled out that will then be used to automate the filtering for the final dataset
##columns of interest for a candidate starship are those in the elementFinder/*_seq.elements.named.stats and elementFinder/*_seq.elements.feat; plus the species names
## starship_id = name given to the starship candidate
## starship_family = family classification of the starship based on the captain
## starship_navisHap = the navis classification of the starship based on the captain similarity and its haplotype based on kmer similarity of the entire region
## species = as previously identified
## strain = all strain names known for the strain combined, split by ";"
## genome = the public accession if present otherwise NA
## assemblyid = the name used for this assembly in the analysis (can be strain or public accession if not yet public)
## contig:start-end = coordinates for the starship
## step_identified = whether 1st or sequential/2nd step (gives an idea of how it was found)
## nested_inside = the ID of a starship that this starship is said to be nested inside
## nested_contained = the ID of a starship said to be nested within this starship
## manual_classification = classify as either good; mid or poor by looking at the alignments
## manual_notes = notes about the classfication if necessary, perhaps why one that looks iffy by alignment is classified as good eventually
echo "starship_id;starship_family;starship_navisHap;species;strain;genome;assemblyid;contig;start;end;length;strand;boundary;emptysite_species;emptysite_strain;emptysite_contig;emptysite_start;emptysite_end;captain_id;captain_start;captain_end;nested_inside_id;nested_contained_id;step_identified;manual_classification;manual_notes" | tr ';' '\t' > ${dataset}.starships.summary.manual_mod.tsv
tail -n+2 elementFinder/${dataset}_seq.elements.ann.feat | awk '{print $1"\t"$2"\t"$3"\tXX\tXX\tXX\tXX\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\tXX\tXX\t"$12"\t"$13"\t"$14"\t"$5"\tXX\tXX\t"$22"\t"$23"\tXX\tXX\tXX"}' | while read line
do

## data on the strain/species
starship=$( echo "${line}" | awk -F "\t" '{print $1}' )
assemblyid=$( echo $starship | awk -F "_" '{print $1}' )
species=$( cat ../metadata.tsv | awk -F "\t" -v assemblyid="$assemblyid" '{if($6 == assemblyid) print $7}' | awk -F "_" '{print $2}'  )
strain=$( cat ../metadata.tsv | awk -F "\t" -v assemblyid="$assemblyid" '{if($6 == assemblyid) print $3";"$4";"$5"XXXXX"}' | sed 's/;;XXXXX//g' | sed 's/;XXXXX//g' | sed 's/;;/;/g' | sed 's/XXXXX//g'   )
genome=$( cat ../metadata.tsv | awk -F "\t" -v assemblyid="$assemblyid" '{if($6 == assemblyid) print $1}' )

##the same data for the genome compared against to identify the insertion (if present)
altassemblyid=$( echo "${line}" | awk -F "\t" '{print $16}' | awk -F "_" '{print $1}' )
altspecies=$( cat ../metadata.tsv | awk -F "\t" -v altassemblyid="$altassemblyid" '{if($6 == altassemblyid) print $7}' | awk -F "_" '{print $2}'  )
if [ "${altspecies}" == "" ]
then
altspecies="."
altstrain="."
else
altstrain=$( cat ../metadata.tsv | awk -F "\t" -v altassemblyid="$altassemblyid" '{if($6 == altassemblyid) print $3";"$4";"$5"XXXXX"}' | sed 's/;;XXXXX//g' | sed 's/;XXXXX//g' | sed 's/;;/;/g' | sed 's/XXXXX//g'  )
fi

##captain position
captainid=$( echo "${line}" | awk -F "\t" '{print $19}' )
captainstart=$( grep ${captainid} ${dataset}.filt_intersect.consolidated.gff | awk -F "\t" '{print $4}' )
captainend=$( grep ${captainid} ${dataset}.filt_intersect.consolidated.gff | awk -F "\t" '{print $5}' )

##which step this starship was found (can just see if it is present in just the initial ID files or both the initial and sequentially added)
step=$( tail -n+2 elementFinder/${dataset}.elements.feat elementFinder/${dataset}_seq.elements.feat | awk -v starship="$starship" '{if($2 == starship) print}' | wc -l | awk '{if($0 == "1") {print "step2"} else if($0 == "2"){print "step1"}}' )

echo "${line}" | awk -v species="$species" -v strain="$strain" -v genome="$genome" -v assemblyid="$assemblyid" -v altspecies="$altspecies" -v altstrain="$altstrain" -v captainstart="$captainstart" -v captainend="$captainend" -v step="$step" -F "\t" '{print $1"\t"$2"\t"$3"\t"species"\t"strain"\t"genome"\t"assemblyid"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"altspecies"\t"altstrain"\t"$16"\t"$17"\t"$18"\t"$19"\t"captainstart"\t"captainend"\t"$22"\t"$23"\t"step"\t\t"}' >> ${dataset}.starships.summary.manual_mod.tsv
done 


###### now the candidates need to be manually validated by inspecting the visualisation
## starships will be sorted into three groups
## easy candidates, when the insertion in clear
## on the fence, when it is not absolitely clear
## poor candidates, when the alignment is bad/ clearly a bad candidate (label as "poor")
## first use the pairViz_seq


##after certain starships have been labelled as poor (i.e. keep all others) then we can generate a matrix of presence/absence PER HAPLOTYPE per genome
cat ${dataset}.starships.summary.manual_mod.tsv | awk -F "\t" '{if($25 !~ "poor") print}'  > ${dataset}.starships.summary.manual_mod.good_only.tsv

##now get a list of all genomes then add sequentially each haplotype labelling presence if so
##the order of the haplotyes should be by family-navis-hap so they are clustered in the plot (if the order keeps)
#list of haplotypes/starships tail -n+2 ${dataset}.starships.summary.manual_mod.good_only.tsv | cut -f2,3 | sort -u
##create header to add to the matrix
tail -n+2 ${dataset}.starships.summary.manual_mod.good_only.tsv | cut -f2,3 | sort -u | awk '{print $1"-"$2}' | tr '\n' '\t' | tr ' ' '-' | awk '{print "genome\t"$0}' | sed 's/	$//g'  > ${dataset}.starships.summary.manual_mod.good_only.per_genome_matrix.tsv
ls assemblies/ | awk -F "." '{print $1}' | while read genome
do
##get species (need it to add the the row in order to connect it with the phylogeny tips) (needs to be 'P' instead of Paecilomyces and no gap e.g. Pvariotii)
species=$( cat ../metadata.tsv | awk -F "\t" -v genome="$genome" '{if($6 == genome) print $7}' | sed 's/aecilomyces_//g' )
strip=$( echo "${species}_${genome}" )
count=$( tail -n+2 ${dataset}.starships.summary.manual_mod.good_only.tsv | cut -f2,3 | sort -u | awk '{print $1"-"$2}' | wc -l )
tail -n+2 ${dataset}.starships.summary.manual_mod.good_only.tsv | cut -f2,3 | sort -u | awk '{print $1"-"$2}' | while read haplotype
do
((counter++))
presence=$( cat ${dataset}.starships.summary.manual_mod.good_only.tsv | awk -F "\t" -v haplotype="$haplotype" -v genome="$genome" 'BEGIN{hold="not_found"} {if(($2"-"$3) == haplotype && $7 == genome) {hold="found"}} END{if(hold=="found") {print "present"} else {print "NA"}}' )
strip=$( echo "${strip};${presence}" )
if [[ ${count} == ${counter} ]]
then
	echo ${strip} | tr ';' '\t' >> ${dataset}.starships.summary.manual_mod.good_only.per_genome_matrix.tsv
fi
done
done


##a simple summary file giving the number of filtered captains identified per strain and the number of good elements
echo "genome;captains;starships;ratio" | tr ';' '\t' > ${dataset}.starships.summary.manual_mod.good_only.starships_vs_captains.summary.tsv
cat ../metadata.tsv | awk -F "\t" '{if($9 == "yes" )print $6}'  | while read genome
do
species=$( cat ../metadata.tsv | awk -F "\t" -v genome="$genome" '{if($6 == genome) print $7}' | sed 's/Paecilomyces_/P/g' )
captains=$( cat geneFinder/${dataset}.filt.gff | grep "${genome}_" | wc -l )
starships=$( cat ${dataset}.starships.summary.manual_mod.good_only.tsv | grep ^"${genome}_" | wc -l )
echo "${species}_${genome};${captains};${starships}" | tr ';' '\t' | awk '{if($2 > 0) {print $0"\t"$3/$2} else {print $0"\t0"}}'
done >> ${dataset}.starships.summary.manual_mod.good_only.starships_vs_captains.summary.tsv


##extract only the good sequences from the set of all given by starfish
##can just use grep -A1 as the elements are on a single line and the naming is non-redundant in text so there is no overlap (so useful)
cat ${dataset}.starships.summary.manual_mod.good_only.tsv | cut -f1 | while read starship; do grep -A1 "${starship}" elementFinder/${dataset}_seq.elements.fna; done | sed 's/|[+-]//g' >> ${dataset}.starships.summary.manual_mod.good_only.fa

```

#### Step 4: Look for evidence of horizontal transfer within the entire Ascomycota subdivision Pezizomycotina

```

```

#### Step 5: Generate phylogeny and alignments for all candidates of horizontal transfer

```

```

