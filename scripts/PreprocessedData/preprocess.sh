#!/bin/bash
PRINSEQ=~/Applications/Tools/prinseq/prinseq-lite.pl
FASTQ_MASKER=~/Applications/root/bin/fastq_masker

MIN_PHRED_SCORE=30

# 1) Perform full prinseq quality filtering
for i in 3223 3236
do
	for j in a b c
	do
		echo ${i}${j}
		rm -rf ${i}${j}
		mkdir ${i}${j}
		cd ${i}${j}
		
		echo -e "PRINSEQ filtering\n=================" > Preprocess_Statistics.txt
		$PRINSEQ \
			-fastq ../../RawData/${i}${j}*R1*.fastq \
			-fastq2 ../../RawData/${i}${j}*R2*.fastq \
			-out_format 3 \
			-out_good ${i}${j}_QC_1_R \
			-out_bad null \
			-min_len 230 \
			-trim_qual_left ${MIN_PHRED_SCORE} \
			-trim_qual_right ${MIN_PHRED_SCORE} \
			-trim_qual_window 5 2>> Preprocess_Statistics.txt #-ns_max_n 6
		cd ..
	done
done
echo "========================================================="


# 2) Replace ALL low-quality nucleotides
for i in 3223 3236
do
	for j in a b c
	do
		echo ${i}${j}
		cd ${i}${j}
		
		echo -e "\n\nNucleotide masking\n==================" >> Preprocess_Statistics.txt
		for r in R_1 R_2
		do
			$FASTQ_MASKER \
				-q ${MIN_PHRED_SCORE} \
				-i ${i}${j}_QC_1_${r}.fastq \
				-o ${i}${j}_nucMask_2_${r}.fastq \
				-v >> Preprocess_Statistics.txt
		done
		
		# Discard reads with missing mate
		rm -rf *singletons*

		cd ..
	done
done
echo "========================================================="
