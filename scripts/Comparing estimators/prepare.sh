#!/bin/bash -l
AMPLICONCLIPPER=~/Applications/Tools/AmpliconClipper
PICARDTOOLS=~/Applications/Tools/picard-tools/picard.jar
BWA=~/Applications/root/bin/bwa

sam2bam() {
	STEM=${1%.sam}
	samtools view -bS "${STEM}.sam" > temp.bam 2> /dev/null
	samtools sort temp.bam "${STEM}" 2> /dev/null
	samtools index "${STEM}.bam" 2> /dev/null
	rm temp.bam;
}

# 1) Trim reads to fit and remove primerID simultaneously
for i in 3223 3236
do
	for j in a b c
	do
		rm -rf ${i}${j}
		mkdir ${i}${j}
		
		echo ${i}${j}
		cd ${i}${j}
		
		$AMPLICONCLIPPER -i ../../Alignments/${i}${j}/${i}${j}_QC_1.sam -o ${i}${j}_clipped.sam -a ../${i}_amplicon.txt
		
		cd ..
	done
done
echo "========================================================="


# 2) Extract reads with picard
for i in 3223 3236
do
	for j in a b c
	do
		echo ${i}${j}
		cd ${i}${j}
		
		java -jar $PICARDTOOLS SamToFastq "I=${i}${j}_clipped.sam" "F=${i}${j}_R1.fastq" "F2=${i}${j}_R2.fastq" FU=single.fastq VALIDATION_STRINGENCY=LENIENT
		rm -f single.fastq ${i}${j}_clipped.sam
		
		cd ..
	done
done
echo "========================================================="


# 3) perform alignment with BWA
for i in 3223 3236
do
	cp ../References/5VM_${i}_cons.fasta .
	$BWA index 5VM_${i}_cons.fasta
	
	for j in a b c
	do
		echo ${i}${j}
		cd ${i}${j}
		
		$BWA mem -t 4 -L 1000 -I 450,100000000 ../5VM_${i}_cons.fasta ${i}${j}_R1.fastq ${i}${j}_R2.fastq | samtools view -h -F 4 - | samtools view -h -F 8 - > ${i}${j}.sam
		sam2bam ${i}${j}.sam
		rm -f ${i}${j}.sam
		rm -f ${i}${j}_R?.fastq
		
		cd ..
	done
	
	rm -f 5VM_32??_cons.fasta.*
done
echo "========================================================="
