#!/bin/bash

# finds group of lanes in a given directory using standard illumina filenames
# and merges BAM files to a specified directory
# param 1: directory of lane BAMs
# param 2: suffix to remove (with wildcards as necessary)
# param 3: output directory
# param 4: number of threads to use

cd $1
module load samtools/1.6

suffix=$2
outDir=$3
threads=$4

all=*.bam
sampleNames=''
for bam in $all
do
	#first remove suffix to get base name
	sample="${bam%${suffix}}"
	#dynamically build new list of base names
	sampleNames="$sampleNames $sample"
done

#remove duplicate base names
sampleNames=($sampleNames)
uniqueNames=($(echo "${sampleNames[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

#now iterate through each base name, find all lanes, and merge to new directory
#keeping RG tags
for s in ${uniqueNames[@]}
do
	laneFiles=${s}*bam
	echo files: $laneFiles
	outName="${outDir}/${s}_merged.bam"
	echo output: $outName
	samtools merge -r $outName $laneFiles --threads $threads
done