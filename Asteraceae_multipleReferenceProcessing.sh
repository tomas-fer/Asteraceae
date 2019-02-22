#!/bin/bash
#----------------MetaCentrum----------------
#PBS -l walltime=2:0:0
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=1gb
#PBS -j oe
#PBS -N combination
#PBS -m abe

#Combines results for mapping to the three Asteraceae references (lettuce, sunflower, safflowerLettuce)
#Works on results of HybPhyloMaker (https://github.com/tomas-fer/HybPhyloMaker)
#Tomas Fer, 2017
#v.0.1.0

#Requires three mapping analysis results in exons/71selected${MISSINGPERCENT}_${genome} where genome is lettuce, sunflower, and safflowerLettuce
#Requires three files for names translation in HybSeqSource (AT_lettuce.txt, AT_sunflower.txt, and AT_safflower.txt)

#Produces:
#- gene list (Arabidopsis codes) for each mapping and sample
#- a combined gene list (Arabidopsis codes) for each sample - ${genome}_list.txt
#- list of genes shared among samples - shared_genes.txt, shared_genes_modif.txt
#- table with number of recovered genes (combined) per samples - combined_genes.txt
#- table with missing data per mapping and gene - all_missdata_sorted.txt
#- table with least amount of missing data per gene and indication of mapping produced it - genome_with_least_missing_data.txt
#- list of genes producing best mapping (for each genome), both Arabidopsis and Asteraceae genome specific codes - ${genome}_best.txt
#- list of selected genes (genes selected by combination of the three mappings), with Asteraceae genome specific codes - selected_genes.txt

#Morever, copies selected genes to exons/71selected${MISSINGPERCENT} and calculates alignment properties using AMAS


#Complete path and set configuration for selected location
#Copy file with settings from home and set variables from settings.cfg
cp $PBS_O_WORKDIR/settings.cfg .
. settings.cfg
. /packages/run/modules-2.0/init/bash
path=/storage/$server/home/$LOGNAME/$data
source=/storage/$server/home/$LOGNAME/HybSeqSource
#Move to scratch
cd $SCRATCHDIR
#Add necessary modules
module add python-3.4.1-gcc

#Make a dir for results
mkdir $path/exons/combination

#Copy files
cp $source/AT_*.txt .
for genome in lettuce sunflower safflowerLettuce; do
	cp $path/exons/71selected${MISSINGPERCENT}_${genome}/MissingDataOverview.txt ${genome}_MissingDataOverview.txt
	cp $path/exons/71selected${MISSINGPERCENT}_${genome}/MissingDataOverview_${MISSINGPERCENT}.txt ${genome}_MissingDataOverview_${MISSINGPERCENT}.txt
	cp $path/exons/71selected${MISSINGPERCENT}_${genome}/selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt ${genome}_selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
done
mv safflowerLettuce_MissingDataOverview.txt safflower_MissingDataOverview.txt
mv safflowerLettuce_MissingDataOverview_${MISSINGPERCENT}.txt safflower_MissingDataOverview_${MISSINGPERCENT}.txt
mv safflowerLettuce_selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt safflower_selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt

#Rename lettuce, safflower and sunflower specific gene codes to Arabidopsis names
#Loop over all three genomes
echo "Renaming..."
for genome in lettuce sunflower safflower; do
	awk '{ print $1 }' ${genome}_MissingDataOverview_${MISSINGPERCENT}.txt > ${genome}_first.txt
	cat ${genome}_MissingDataOverview_${MISSINGPERCENT}.txt | cut -d' ' -f 2- > ${genome}_rest.txt
	#Make changes for all genome lists
	for geneset in lettuce sunflower safflower; do
		echo -e "$genome by $geneset"
		cat AT_${geneset}.txt | while read -r a b; do
			sed -i "s/Assembly_$b$/$a/" ${genome}_first.txt
		done
	done
	paste ${genome}_first.txt ${genome}_rest.txt | tr '\t' ' ' > ${genome}_MissingDataOverview_${MISSINGPERCENT}_modified.txt
	rm ${genome}_first.txt ${genome}_rest.txt
	cp ${genome}_MissingDataOverview_${MISSINGPERCENT}_modified.txt $path/exons/combination
done

#Make gene lists per species from MissingDataOverview
echo "Making species-specific gene lists..."
for genome in lettuce sunflower safflower; do
	#extract species names (take 1st line, change spaces to EOLs, delete last two lines, delete first line)
	head -1 ${genome}_MissingDataOverview_${MISSINGPERCENT}_modified.txt | tr ' ' '\n' | sed '$d' | sed '$d' | sed '1d' > ${genome}_list.txt
	#Loop over species
	count=1
	for i in $(cat ${genome}_list.txt); do
		count=$((count + 1))
		#print gene name (1st column) if value in count's column is not 'N/A'
		awk -F' ' -v val="$count" '{ if($val != "N/A") { print $1 } }' ${genome}_MissingDataOverview_${MISSINGPERCENT}_modified.txt | sed '$d' | sed '$d' | sed '1d' > ${i}_${genome}.txt
		cp ${i}_${genome}.txt $path/exons/combination
	done
done

#Combine gene lists
echo "Combining gene lists..."
for i in $(cat lettuce_list.txt); do
	cat ${i}* | sort | uniq > ${i}_combined_genes.txt
	cp ${i}_combined_genes.txt $path/exons/combination
	nrcombgenes=$(cat ${i}_combined_genes.txt | wc -l)
	echo -e "$i\t$nrcombgenes" >> combined_genes.txt
done
cp combined_genes.txt $path/exons/combination

#Make list of shared genes
awk '(++c[$0])==(ARGC-1)' *_combined* > shared_genes.txt
cp shared_genes.txt $path/exons/combination
cp shared_genes.txt shared_genes_modif.txt

#Make renaming files for genome-specific list of recovered genes (select only recovered genes for renaming back)
for genome in lettuce sunflower safflower; do
	for i in $(cat ${genome}_selected_genes_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt | cut -d'_' -f 2); do
		grep "$i$" AT_${genome}.txt >> AT_${genome}_back.txt
	done
	cp AT_${genome}_back.txt $path/exons/combination
done

#Rename Arabidopsis names back to genome specific codes
#Make changes for all genome lists
echo "Renaming gene names back..."
for geneset in lettuce sunflower safflower; do
	echo -e "$geneset"
	cat AT_${geneset}_back.txt | while read -r a b; do
		sed -i "s/$a$/Assembly_$b/" shared_genes_modif.txt
	done
done
cp shared_genes_modif.txt $path/exons/combination

#Rename lettuce, safflower and sunflower specific gene codes to Arabidopsis names
#Loop over all three genomes
echo "Renaming..."
for genome in lettuce sunflower safflower; do
	awk '{ print $1 }' ${genome}_MissingDataOverview.txt > ${genome}_first.txt
	cat ${genome}_MissingDataOverview.txt | cut -d' ' -f 2- > ${genome}_rest.txt
	#Make changes for all genome lists
	for geneset in lettuce sunflower safflower; do
		echo -e "$genome by $geneset"
		cat AT_${geneset}.txt | while read -r a b; do
			sed -i "s/Assembly_$b$/$a/" ${genome}_first.txt
		done
	done
	paste ${genome}_first.txt ${genome}_rest.txt | tr '\t' ' ' > ${genome}_MissingDataOverview_modified.txt
	rm ${genome}_first.txt ${genome}_rest.txt
	cp ${genome}_MissingDataOverview_modified.txt $path/exons/combination
done

echo "Preparing table of missing data..."
for i in lettuce sunflower safflower; do
	#print first and second last ($--NF), i.e. columns with gene name and amount of missing data, and delete last line in the file (contains unnecessary header)
	awk '{print $1 "\t" $--NF}' ${i}_MissingDataOverview_modified.txt | sed '$d'> ${i}_missdata.txt
	sort ${i}_missdata.txt > ${i}_missdata_sorted.txt
	cp ${i}_missdata_sorted.txt $path/exons/combination
done

echo "Making gene lists..."
#join lettuce and sunflower data (report also unpaired lines by '-a' and adds missing fields (na) by '-e' (this works only together with '-o' command !!!)
join -e'NA' -a1 -a2 -o 1.1,1.2,2.2 lettuce_missdata_sorted.txt sunflower_missdata_sorted.txt | sort > lettuce_sunflower_missdata_sorted.txt
#join and add safflower
join -e'NA' -a1 -a2 -o 1.1,1.2,1.3,2.2 lettuce_sunflower_missdata_sorted.txt safflower_missdata_sorted.txt | sort > all_missdata_sorted.txt
cp all_missdata_sorted.txt $path/exons/combination

#select minimum value of missing data on each line (looping from 2 to NF, i.e. number of fields), also reports the column in which this minimum value appears
awk '{min=$2;genome=2;for(i=2;i<=NF;i++){if ($i < min) genome = i;if($i < min) min = $i}print $1 "\t" min "\t" genome}' all_missdata_sorted.txt > genome_with_least_missing_data.txt
cp genome_with_least_missing_data.txt $path/exons/combination
#select only genes for which mapping to the lettuce reference gave least missing data (lines with '2' at the end of the line)
grep '2$' genome_with_least_missing_data.txt | cut -f1 > lettuce_best.txt
#select only genes for which mapping to the sunflower reference gave least missing data (lines with '3' at the end of the line)
grep '3$' genome_with_least_missing_data.txt | cut -f1 > sunflower_best.txt
#select only genes for which mapping to the safflower reference gave least missing data (lines with '4' at the end of the line)
grep '4$' genome_with_least_missing_data.txt | cut -f1 > safflower_best.txt
cp *_best.txt $path/exons/combination

#Make lists of best matching genes that are shared by all samples
for genome in lettuce sunflower safflower; do
	for i in $(cat shared_genes.txt); do
		grep ${i} ${genome}_best.txt >> ${genome}_best_shared.txt
	done
done
cp *_best_shared.txt $path/exons/combination

#Rename Arabidopsis names back to genome specific codes in all best_shared lists
echo "Making list of shared best genes..."

echo -e "...lettuce"
for genome in lettuce sunflower safflower; do
	cat AT_${genome}_back.txt | while read -r a b; do
		sed -i "s/$a$/Assembly_$b/" lettuce_best_shared.txt
	done
done
cp lettuce_best_shared.txt $path/exons/combination/lettuce_best_shared_modified.txt

echo -e "...sunflower"
for genome in sunflower lettuce safflower; do
	cat AT_${genome}_back.txt | while read -r a b; do
		sed -i "s/$a$/Assembly_$b/" sunflower_best_shared.txt
	done
done
cp sunflower_best_shared.txt $path/exons/combination/sunflower_best_shared_modified.txt

echo -e "...safflower"
for genome in safflower lettuce sunflower; do
	cat AT_${genome}_back.txt | while read -r a b; do
		sed -i "s/$a$/Assembly_$b/" safflower_best_shared.txt
	done
done
cp safflower_best_shared.txt $path/exons/combination/safflower_best_shared_modified.txt

#Make final lists
cat *_best_shared.txt > selected_genes.txt
cp selected_genes.txt $path/exons/combination

#Preparing data for concatenation
echo "Copying selected gene assemblies..."
mkdir $path/exons/71selected${MISSINGPERCENT}
mkdir $path/exons/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}

for i in $(cat lettuce_best_shared.txt); do
	cp $path/exons/71selected${MISSINGPERCENT}_lettuce/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas $path/exons/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
done
for i in $(cat sunflower_best_shared.txt); do
	cp $path/exons/71selected${MISSINGPERCENT}_sunflower/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas $path/exons/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
done
for i in $(cat safflower_best_shared.txt); do
	cp $path/exons/71selected${MISSINGPERCENT}_safflowerLettuce/deleted_above${MISSINGPERCENT}/${i}_modif${MISSINGPERCENT}.fas $path/exons/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
done

#Calculating alignment properties using AMAS
cd $path/exons/71selected${MISSINGPERCENT}/deleted_above${MISSINGPERCENT}
cp $source/AMAS.py .
python3 AMAS.py summary -f fasta -d dna -i *.fas
cp summary.txt ../summarySELECTED_${MISSINGPERCENT}_${SPECIESPRESENCE}.txt
rm AMAS.py

#Clean scratch/work directory
if [[ $PBS_O_HOST == *".cz" ]]; then
	#delete scratch
	if [[ ! $SCRATCHDIR == "" ]]; then
		rm -rf $SCRATCHDIR/*
	fi
fi

echo "Finished..."
