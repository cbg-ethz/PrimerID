#!/bin/bash
FASTQC=/Applications/FastQC.app/Contents/Resources/Java/fastqc

for i in 3223 3236
do
	for j in a b c
	do
		echo ${i}${j}
		rm -rf ${i}${j}
		mkdir ${i}${j}
		
		for k in R1 R2
		do
			rm -rf tmp
			mkdir tmp
			
			cd tmp/
			$FASTQC -o . ../${i}${j}_S?_L001_${k}_001.fastq
			unzip *.zip
			cd ..
			
			cp tmp/${i}${j}_S?_L001_${k}_001_fastqc/Images/per_tile_quality.png ${i}${j}/${i}${j}_${k}.png
			rm -rf tmp
		done
	done
done