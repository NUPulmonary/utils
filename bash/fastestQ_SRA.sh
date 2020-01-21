#!/bin/bash
#Downloads a set of SRA files listed in  a return-separated text file
#parameter 1: filepath of the text file
#parameter 2: number of threads to use
#parameter 3: scratch directory (created on the fly)
#parameter 4: output directory (fastq files will be overwritten!)

mkdir -p "$3"

#note: automatically splits paired reads
module load sratoolkit/2.9.2
ids=$(<"$1") #put SRA/ERA IDs into single string to iterate through
for i in $ids
do
	fasterq-dump "$i" -e $2 -t "$3" -O "$4" -f
done

#gzip all files
cd $4
parallel -j $2 gzip ::: *.fastq #run in parallel

#now remove tmp files
rm -r "$3"
