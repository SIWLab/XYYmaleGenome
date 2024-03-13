#!/bin/bash

## redoing blast with the whoel genome not just the transcriptome
#fasta="NChap1.out.fasta"
#genelist="genelisthap1.out.txt"
#while read i; do 
#	samtools faidx $fasta $i >> hap1_sequences.fa
#done < $genelist

echo "creating hap2 genome db"

makeblastdb -in /ohta2/Rumex/Dovetail_xyy_pacbio-male/final_scaffolded_assemblies/NChap2_final.fa \
	-input_type fasta -dbtype nucl -out hap2genome

wait

echo "blasting hap1 sequences to hap2 geno"

blastn -db hap2genome -query hap1_sequences.fa \
	-max_hsps 1 \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore score" \
	-out blast_h2genome.tsv



