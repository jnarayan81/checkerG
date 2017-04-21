#!/bin/bash

#A pipeline to check the two genome assembled by different methods - Jitendra

#General location settings for checker
scriptLoc=/home/jitendra/ETC/TESTing/checkerG/scriptBase
rGenome=/home/jitendra/ETC/TESTing/checkerG/contigs.fa
tGenome=/home/jitendra/ETC/TESTing/checkerG/cleanMIRA_100.fa
genomeFileName=contigs.fa
config=/home/jitendra/ETC/TESTing/checkerG/lastZconfig

#General thresholds and folders
DIR=OutData

#/home/urbe/Tools/Alienomics_v1.1/Adineta_vaga_v2.0.scaffolds.fa

if [ -d "$DIR" ]; then
    printf '%s\n' "Removing Directory ($DIR)"
    rm -rf "$DIR"
fi

mkdir $DIR; cd $DIR

echo "Align to genome"
perl $scriptLoc/parallelLastZ.pl -q $rGenome -t $tGenome -c $config -s 4 -l 1000

echo "Move all alignment to ALN folder"
mkdir 'ALN'
mv *.lz  ALN

cd ALN

for f in *.lz
do
	echo "Processing $f"
	echo "Filtering direct overlaps/alignments"
	perl $scriptLoc/filterOverlaps.pl $f > $f.bed

done

#find . -name "*.lz" -size 0 -delete

cat *.bed > allALN.bed
cp allALN.bed ..

echo "Come out of the dir"
cd ..

echo "Merge overlapping blocks of interest +/-"
perl $scriptLoc/mergeOverlaps.pl allALN.bed general > finalALN.bed

echo "Create the faidx"
samtools faidx $rGenome

echo "extract all breaks location"
perl $scriptLoc/findBRK.pl finalALN.bed final_breakspoints.txt $rGenome.fai 1

echo "Looking for TRF"
$scriptLoc/trf/trf409.linux64 $rGenome 2 7 7 80 10 50 500 -f -d -m
perl $scriptLoc/trf/trfparser_v1.pl $genomeFileName.2.7.7.80.10.50.500.dat 1
rm -rf *.tmp *.html *.mask *.txt.parse *.dat.parse *.dat

echo "Remove the small contigs"
perl $scriptLoc/removeSmall.pl 1000 $rGenome > newGenome.fa

exit;

echo "Map the contigs"
bwa index newGenome.fa
#bwa mem -x pacbio newGenome.fa pacbio.fq > aln.sam

echo "Sam to Bam conversion"
samtools view -Sb  aln.sam  >  aln.bam

#Create genomecov file
echo "checking for genome coverage"
bedtools genomecov -ibam aln.bam -bga > allCoverage.bed

#Look into depth
#samtools depth $bamFile  |  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}'

#Create final GFF result file for visualization
echo "Creating final GFF file in detail"
perl $scriptLoc/createFinalGFF.pl final_breakspoints.txt $genomeFileName.2.7.7.80.10.50.500.final.parse allCoverage.bed 1000 newGenome.fa > final.gff


