#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=shared       # default "shared", if not specified
#SBATCH --time=3-00:00:00       # run time in days-hh:mm:ss
#SBATCH --nodes=1               # require 1 nodes
#SBATCH --ntasks-per-node=64    # cpus per node (by default, "ntasks"="cpus")
#SBATCH --mem=256000            # RAM per node in megabytes
#SBATCH --error=ss.3.starfish.job.%J.err
#SBATCH --output=ss.3.starfish.job.%J.out
# Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output


#######################################################
#################### STEP 0: SETUP ####################
#######################################################

##source the conda env ability to run in bash
source ~/miniconda3/etc/profile.d/conda.sh

dataset="paecilomyces"
cd /scratch/saodonnell/projects/${dataset}/

##set up some paths
#export PATH=$PATH:~/starfish/bin/
#export PATH=$PATH:~/starfish/CNEFinder/


##################################################################
#################### STEP 1: Running starfish ####################
##################################################################
threads=64
cd starfish/

## this step is to take the genomes already annotated and run starfish on them
## activate the env
conda activate starfish

##need to recreate the softmasked files with the headers having the assembly name and removing the underscore between the 'scaffold' amnd the scaffold number
mkdir assemblies
ls ../annotation/3.annotate/3.gff3/*/*.softmasked.fa | while read file
do
assembly=$( echo $file | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'  )
cat $file | awk -v assembly="$assembly" '{if($0 ~ ">") {print ">"assembly"_"$0} else {print}}' | sed 's/_>/_/g' | sed 's/scaffold_/scaffold/g' > assemblies/${assembly}.fa
done

##do the same with the gff3 files
##here we also add another tag to the end of the file with the name of the gene being accession"_"gene
mkdir gff3
ls ../annotation/3.annotate/3.gff3/*/*.gff3 | while read file
do
assembly=$( echo $file | awk -F "/" '{print $NF}' | awk -F "." '{print $1}'  )
cat $file | while read line
do
gene=$( echo "$line" | awk -F "\t" '{print $9}' | awk -F ";" '{print $1}' | sed 's/ID=//g'   )
echo "$line" | awk -v assembly="$assembly" -v gene="$gene" '{if($0 ~ "#") {print $0} else {print assembly"_"$0"Name="assembly"_"gene";"}}' | sed 's/scaffold_/scaffold/g' 
done > gff3/${assembly}.gff3
done


##need a list of the assemblies and paths (using the softmasked genomes previously annotated)
realpath assemblies/*.fa | while read line
do
name=$( echo $line | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
echo ${name}";"${line} | tr ';' '\t'
done > ome2assembly.txt

##list of annotation files and paths also
realpath gff3/*.gff3 | while read line
do
name=$( echo $line | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
echo ${name}";"${line} | tr ';' '\t'
done > ome2gff.txt

## a consolidated file with all annotations to compare with Starfish
cat gff3/*.gff3 > ${dataset}.assemblies.gff3

## now run the genefinder step with the captains
mkdir geneFinder
starfish annotate -T ${threads} -x ${dataset} -a ome2assembly.txt -g ome2gff.txt -p ~/starfish/db/YRsuperfams.p1-512.hmm -P ~/starfish/db/YRsuperfamRefs.faa -i tyr -o geneFinder/
##get unique regions for each captain (cluster them based on distance)
##first need to generate a tsv with the path to the gff generated in the previous step (stupid step)
#realpath geneFinder/${dataset}.filt.gff | awk -v dataset="$dataset" '{print dataset"\t"$0}' > ome2consolidatedGFF.txt
##consolidate starfish captain positions with annotations (should spit out a file called XXX_tyr.filt_intersect.consolidated.gff)
starfish consolidate -o ./ -g ${dataset}.assemblies.gff3 -G geneFinder/${dataset}.filt_intersect.gff
##now generate a tsv with the path to the gff generated in the previous step
realpath ${dataset}.filt_intersect.consolidated.gff | awk -v dataset="$dataset" '{print dataset"\t"$0}' > ome2consolidatedGFF.txt
starfish sketch -m 10000 -q geneFinder/${dataset}.filt_intersect.ids -g ome2consolidatedGFF.txt -i s -x ${dataset} -o geneFinder/
##grab just the candidate captain coords
grep -P '\ttyr\t' geneFinder/${dataset}.bed > geneFinder/${dataset}.tyr.bed 


##now use the tyrosine elements to find starships
mkdir elementFinder
##generate a blastdb with the assemblies
mkdir blastdb
cut -f2 ome2assembly.txt | xargs cat > blastdb/${dataset}.assemblies.fa
makeblastdb -in blastdb/${dataset}.assemblies.fa -out blastdb/${dataset}.assemblies -parse_seqids -dbtype nucl
##look for the actual insertions
starfish insert -T ${threads} -a ome2assembly.txt -d blastdb/${dataset}.assemblies -b geneFinder/${dataset}.tyr.bed -i tyr -x ${dataset} -o elementFinder/
##sequentially running starfish insert, first with 90% identity and match size of 1kb; now 80% identity qnd match size of 750bp with those not previously found
starfish insert -T ${threads} -a ome2assembly.txt -d blastdb/${dataset}.assemblies -b elementFinder/${dataset}.insert.bed -i tyr -x ${dataset}_seq -o elementFinder/ --pid 80 --hsp 750

##find flanking elements
starfish flank -a ome2assembly.txt -b elementFinder/${dataset}.insert.bed -x ${dataset} -o elementFinder/
starfish flank -a ome2assembly.txt -b elementFinder/${dataset}_seq.insert.bed -x ${dataset}_seq -o elementFinder/

##need to concatenate the output stats files (which only output the stats for the new starships detected during the sequential run, not all)
mv elementFinder/${dataset}_seq.insert.stats elementFinder/${dataset}_temp.insert.stats
cat elementFinder/${dataset}.insert.stats elementFinder/${dataset}_temp.insert.stats > elementFinder/${dataset}_seq.insert.stats
#rm elementFinder/${dataset}_temp.insert.stats
mv elementFinder/${dataset}_seq.flank.singleDR.stats  elementFinder/${dataset}_temp.flank.singleDR.stats 
cat elementFinder/${dataset}.flank.singleDR.stats elementFinder/${dataset}_temp.flank.singleDR.stats  > elementFinder/${dataset}_seq.flank.singleDR.stats 
#rm elementFinder/${dataset}_seq.flank.singleDR.stats
mv elementFinder/${dataset}_seq.elements.named.stats  elementFinder/${dataset}_temp.elements.named.stats 
cat elementFinder/${dataset}.elements.named.stats elementFinder/${dataset}_temp.elements.named.stats  > elementFinder/${dataset}_seq.elements.named.stats 
#rm elementFinder/${dataset}_seq.elements.named.stats

##summarise the info (not using the flank data)
starfish summarize -a ome2assembly.txt -b elementFinder/${dataset}.insert.bed -x ${dataset} -o elementFinder/ -S elementFinder/${dataset}.insert.stats -g ome2consolidatedGFF.txt -t geneFinder/${dataset}.filt_intersect.ids 
starfish summarize -a ome2assembly.txt -b elementFinder/${dataset}_seq.insert.bed -x ${dataset}_seq -o elementFinder/ -S elementFinder/${dataset}_seq.insert.stats -g ome2consolidatedGFF.txt -t geneFinder/${dataset}.filt_intersect.ids 
##add colours to the features
awk '{
    if ($4 ~ /DR/) print $0 "\t255,255,0";        # Yellow for 'DR' in column 4
    else if ($4 ~ /TIR/) print $0 "\t255,165,0";  # Orange for 'TIR' in column 4
    else if ($5 ~ /tyr|cap/) print $0 "\t255,0,0"; # Red for 'tyr' or 'cap' in column 5
    else if ($5 != ".") print $0 "\t128,0,128";   # Purple if column 5 is not '.'
    else print $0 "\t169,169,169";                # Dark gray otherwise
}' elementFinder/${dataset}.elements.bed > elementFinder/${dataset}.elements.color.bed
awk '{
    if ($4 ~ /DR/) print $0 "\t255,255,0";        # Yellow for 'DR' in column 4
    else if ($4 ~ /TIR/) print $0 "\t255,165,0";  # Orange for 'TIR' in column 4
    else if ($5 ~ /tyr|cap/) print $0 "\t255,0,0"; # Red for 'tyr' or 'cap' in column 5
    else if ($5 != ".") print $0 "\t128,0,128";   # Purple if column 5 is not '.'
    else print $0 "\t169,169,169";                # Dark gray otherwise
}' elementFinder/${dataset}_seq.elements.bed > elementFinder/${dataset}_seq.elements.color.bed

##visualise 
mkdir pairViz
starfish pair-viz -m all -t empty -T ${threads} -A nucmer -a ome2assembly.txt -b elementFinder/${dataset}.elements.color.bed -S elementFinder/${dataset}.elements.named.stats -o pairViz/
##visualise the lenient version, modifying the nucmer parameters also (due to the leniency allowed, therefore allowing for more distant pairing)
mkdir pairViz_seq
starfish pair-viz -m all -t empty -T ${threads} -A nucmer --deltafilteropts \'-i 75\' -a ome2assembly.txt -b elementFinder/${dataset}_seq.elements.color.bed -S elementFinder/${dataset}_seq.elements.named.stats -o pairViz_seq/


## first assign the family using reference captain sequences HMM
hmmsearch --noali --notextw -E 0.001 --max --cpu 12 --tblout elementFinder/${dataset}_tyr_vs_YRsuperfams.out ~/starfish/db/YRsuperfams.p1-512.hmm geneFinder/${dataset}.filt_intersect.fas
perl -p -e 's/ +/\t/g' elementFinder/${dataset}_tyr_vs_YRsuperfams.out | cut -f1,3,5 | grep -v '#' | sort -k3,3g | awk '!x[$1]++' > elementFinder/${dataset}_tyr_vs_YRsuperfams_besthits.txt

## rename the captains by their starships
grep -P '\tcap\t' elementFinder/${dataset}_seq.elements.bed | cut -f4,7 > elementFinder/${dataset}_seq.cap2ship.txt
~/starfish/aux/searchReplace.pl --strict -i elementFinder/${dataset}_tyr_vs_YRsuperfams_besthits.txt -r elementFinder/${dataset}_seq.cap2ship.txt > elementFinder/${dataset}_elements_vs_YRsuperfams_besthits.txt

## group starships into naves using coverage threshold of 25% and identity of 50%
## make sure to use the filt_intersect.fas file for the putative captains as otherwise the name will not correspond downstream and naves information will be lost
mmseqs easy-cluster geneFinder/${dataset}.filt_intersect.fas elementFinder/${dataset}_tyr elementFinder/ --threads 2 --min-seq-id 0.5 -c 0.25 --alignment-mode 3 --cov-mode 0 --cluster-reassign
~/starfish/aux/mmseqs2mclFormat.pl -i elementFinder/${dataset}_tyr_cluster.tsv -g navis -o elementFinder/

## group starships into haplotypes using kmer similarity (95% threshold?)
~/starfish/bin/starfish sim -m element -t nucl -b elementFinder/${dataset}_seq.elements.bed -x ${dataset}_seq -o elementFinder/ -a ome2assembly.txt
#~/starfish/bin/starfish group -m mcl -s elementFinder/${dataset}_seq.element.nucl.sim -i hap -o elementFinder/ -t 0.05
~/starfish/bin/starfish group -s elementFinder/${dataset}_seq.element.nucl.sim -i hap -o elementFinder/ -t 0.05


## replace captain IDs with starship IDs in the naves file
~/starfish/aux/searchReplace.pl -i elementFinder/${dataset}_tyr_cluster.mcl -r elementFinder/${dataset}_seq.cap2ship.txt > elementFinder/${dataset}_seq.element_cluster.mcl
## merge navis with haplotype to create a navis-haplotype label for each Starship
~/starfish/aux/mergeGroupfiles.pl -t elementFinder/${dataset}_seq.element_cluster.mcl -q elementFinder/${dataset}_seq.element.nucl.I1.5.mcl > elementFinder/${dataset}_seq.element.navis-hap.mcl
## convert mcl to gene2og format to simplify downstream parsing:
awk '{ for (i = 2; i <= NF; i++) print $i"\t"$1 }' elementFinder/${dataset}_seq.element.navis-hap.mcl > elementFinder/${dataset}_seq.element.navis-hap.txt
## now add the family and navis-haplotype into to the element.feat file to consolidate metadata:
join -t$'\t' -1 1 -2 2 <(sort -t$'\t' -k1,1 elementFinder/${dataset}_seq.element.navis-hap.txt | grep -P '_e|_s') <(sort -t$'\t' -k2,2 elementFinder/${dataset}_seq.elements.feat) | awk -F'\t' '{print}' > elementFinder/${dataset}_seq.elements.temp.feat
echo -e "#elementID\tfamilyID\tnavisHapID\tcontigID\tcaptainID\telementBegin\telementEnd\telementLength\tstrand\tboundaryType\temptySiteID\temptyContig\temptyBegin\temptyEnd\temptySeq\tupDR\tdownDR\tDRedit\tupTIR\tdownTIR\tTIRedit\tnestedInside\tcontainNested" > elementFinder/${dataset}_seq.elements.ann.feat
join -t$'\t' -1 1 -2 1 <(sort -t$'\t' -k1,1 elementFinder/${dataset}_elements_vs_YRsuperfams_besthits.txt | grep -P '_e|_s' | cut -f1,2) <(sort -t$'\t' -k1,1 elementFinder/${dataset}_seq.elements.temp.feat) | awk -F'\t' '{print}' >> elementFinder/${dataset}_seq.elements.ann.feat


##now consider where the starships are in the genome
mkdir regionFinder

## create a file with tyrs that are not found in any elements (this will let us assign them to fragmented haplotypes in the dereplicate analysis, which can be helpful if annotating 'dead' or 'derelict' element copies):
grep -f <(comm -23 <(cut -f1 geneFinder/${dataset}.filt_intersect.ids | sort) <(grep -P '\tcap\t|\ttyr\t' elementFinder/${dataset}_seq.elements.bed | cut -f4| sort)) geneFinder/${dataset}.tyr.bed > regionFinder/unaffiliated_tyrs.bed
## filter for confident orthogroups
~/starfish/aux/filterOG.pl -O orthofinder/orthofinder_results/Results_*/Orthogroups/Orthogroups.txt -a 1 -c 5 -o ./

## now, dereplicate your data to identify independently segregating element insertions
starfish dereplicate -e elementFinder/${dataset}_seq.element.navis-hap.mcl -t regionFinder/unaffiliated_tyrs.bed -F elementFinder/${dataset}_seq.elements.feat -S elementFinder/${dataset}_seq.elements.named.stats -O Orthogroups.a1.c6.txt -g ome2gff.txt -x ${dataset} -o regionFinder/ --flanking 6 --mismatching 1

## count the number of regions each navis-haplotype is found in (i.e. transposed element)
grep -v '#' regionFinder/${dataset}.fog3.d600000.m1.dereplicated.txt | cut -f2 | sort | uniq -c | perl -pe 's/ +//' | sort -k1,1nr

## now visualise the loci
mkdir locusViz
## add gc content information to be used by the locus visualisation
~/starfish/aux/seq-gc.sh -Nbw 1000 blastdb/${dataset}.assemblies.fna > ${dataset}.assemblies.gcContent_w1000.bed
rm blastdb/${dataset}.assemblies.fna

starfish locus-viz -T 2 -m region -a ome2assembly.txt -b elementFinder/${dataset}_seq.elements.bed -x ${dataset} -o locusViz/ -A nucmer -r regionFinder/${dataset}.fog3.d600000.m1.regions.txt -d regionFinder/${dataset}.fog3.d600000.m1.dereplicated.txt -j regionFinder/${dataset}.fog3.d600000.m1.haplotype_jaccard.sim  -g ome2consolidatedGFF.txt --tags geneFinder/${dataset}.filt_intersect.ids --gc ${dataset}.assemblies.gcContent_w1000.bed

##now use the R scripts that were made and the R environment with gggenomes to generate the plots
conda activate R_env
ls locusViz/*.R | while read script
do
region=$( echo $script | awk -F "/" '{print $NF}' | awk -F "." '{print $2}' )
##modify the Rscript so that the legend is printed at the top and adjusting the width by adding the 25kb on each flank used for alignment
sed 's/geom_gene_tag/theme(legend.position = "top")+geom_gene_tag/g' $script | sed 's/regionSeqs$length/regionSeqs$length+50000/g' > locusViz/${dataset}.${element}.mod.R
Rscript locusViz/${dataset}.${element}.mod.R
done
conda deactivate

##now we have all the locations and their starships, but what about visualising the elements and all their regions
##pick the flanking size in kb
flank=100
mkdir elementViz_${flank}kbflank
##run the element visualisation by navis grouping (to see if haplotypes are truely different etc)
##to do this the list of elements given is the order of the plot, therefore sort by haplotype and locusViz region
cat elementFinder/${dataset}_seq.elements.ann.feat | cut -f2,3 | awk -F "-" '{print $1}' | sort -u | while read set
do
##new name for the set
set2=$( echo ${set} | awk '{print $1"-"$2}'  )
##get the list of element names per family-navis clade
grep "${set}" elementFinder/${dataset}_seq.elements.ann.feat | cut -f1 | while read element
do
##now get the regions associated with each element and sort by the region (here we lose haplotype information in the sorting...)
region=$( grep ${element} regionFinder/${dataset}.fog3.d600000.m1.dereplicated.txt | cut -f1 )
echo ${element}";"${region} | tr ';' '\t' | sort -k2 
done > elementViz_${flank}kbflank/${set2}.list
starfish locus-viz -T 2 -m element -a ome2assembly.txt -b elementFinder/${dataset}_seq.elements.bed -x ${set2} -U ${flank}000 -D ${flank}000 -l elementViz_${flank}kbflank/${set2}.list  -o elementViz_${flank}kbflank/ -A nucmer -r regionFinder/${dataset}.fog3.d600000.m1.regions.txt -d regionFinder/${dataset}.fog3.d600000.m1.dereplicated.txt -j regionFinder/${dataset}.fog3.d600000.m1.haplotype_jaccard.sim  -g ome2consolidatedGFF.txt --tags geneFinder/${dataset}.filt_intersect.ids --gc ${dataset}.assemblies.gcContent_w1000.bed
done

##now use the R scripts that were made and the R environment with gggenomes to generate the plots
conda activate R_env
ls elementViz_${flank}kbflank/*.R | while read script
do
element=$( echo $script | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' )
##modify the Rscript so that the legend is printed at the top and adjusting the width by adding the 25kb on each flank used for alignment
sed 's/geom_gene_tag/theme(legend.position = "top")+geom_gene_tag/g' $script | sed 's/regionSeqs$length/regionSeqs$length+50000/g' > elementViz_${flank}kbflank/${element}.mod.R
Rscript elementViz_${flank}kbflank/${element}.mod.R
done
conda deactivate


echo "################################################################"
echo "Need to now manually validate the candidate starships; have fun"
echo "################################################################"


