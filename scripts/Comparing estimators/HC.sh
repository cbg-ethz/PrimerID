#!/bin/bash -l
source /etc/profile
module load repo/beerenwinkel
module load HIV/haploclique

i=$1
j=$2

cd ${i}${j}
rm -rf results
mkdir results
cd results

for k in {1..30}
do
	touch ../../${i}${j}_cycle_${k}
	haploclique-assembly -r ../../5VM_${i}_cons.fasta -i ../${i}${j}.bam
	rm ../../${i}${j}_cycle_${k}
done
