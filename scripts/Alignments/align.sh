#!/bin/bash
PIDALIGN=../../../pidalign

# 1) Perform Needleman-Wunsch alignments
for i in 3223 3236
do
	for j in a b c
	do
		echo ${i}${j}
		rm -rf ${i}${j}
		mkdir ${i}${j}
		cd ${i}${j}
		
		echo "Aligning after first QC step"
		echo -e "Aligning after first QC step\n============================" > Alignment_Statistics.txt
		$PIDALIGN \
			-r1 ../../PreprocessedData/${i}${j}/${i}${j}_QC_1_R_1.fastq \
			-r2 ../../PreprocessedData/${i}${j}/${i}${j}_QC_1_R_2.fastq \
			-o ${i}${j}_QC_1.sam \
			--ref ../../References/5VM_${i}_N_withPrimer_cons.fasta \
			-v 1 \
			-s 210 \
			-d 4 \
			-Q >> Alignment_Statistics.txt
		
		echo "Aligning after second nucleotide masking step"
		echo -e "\n\nAligning after second nucleotide masking step\n=============================================" >> Alignment_Statistics.txt
		$PIDALIGN \
			-r1 ../../PreprocessedData/${i}${j}/${i}${j}_nucMask_2_R_1.fastq \
			-r2 ../../PreprocessedData/${i}${j}/${i}${j}_nucMask_2_R_2.fastq \
			-o ${i}${j}_nucMask_2.sam \
			--ref ../../References/5VM_${i}_N_withPrimer_cons.fasta \
			-v 1 \
			-s 210 \
			-d 20 >> Alignment_Statistics.txt
		cd ..
	done
done
echo "========================================================="
