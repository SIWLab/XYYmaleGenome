#!/bin/bash
awk -F"-" '/^>/{print $1; next}1' NChap1_final.all.maker.transcripts.fasta > NChap1.out.fasta

awk '{print $2}' genelisthap1.txt > genelisthap1.out.txt 
echo "Filter out candidate loss genes from nchap1 transcript fasta"
genelist="genelisthap1.out.txt"
fasta="NChap1.out.fasta"
while read i; do
    samtools faidx $fasta $i  >> hap1_sequences.fa
done < $genelist

wait

echo "starting blastn"
makeblastdb -in NChap2_final.all.maker.transcripts.fasta -dbtype nucl \
  -input_type fasta -out hap2transcripts
  
wait



blastn -db hap2transcripts -query hap1_sequences.fa \
  -out blastout.txt
