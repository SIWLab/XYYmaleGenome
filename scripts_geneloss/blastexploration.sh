# the following are one-liners I used to parse and explore blast output

## blast results - blasting hap1 missing genes against hap2 transcript fasta ##
## first attempt:
blastn -db hap2transcripts -query hap1_sequences.fa   -out blastout.txt
blastn -db hap2transcripts -query hap1_sequences.fa  -outfmt 6 -out blastout_tab.txt
## in parsing the results I realized that there were multiple instances of blast hits for the same sequence pairs. Not very useful for me, because at least one strong blast his is enough for me to discount a gene a properly missing.
## setting max-hsps to 1 (default is 0, whihc means no limit)
blastn -db hap2transcripts -query hap1_sequences.fa -max_hsps 1  -outfmt 6 -out blast_max1.tsv
blastn -db hap2transcripts -query hap1_sequences.fa -max_hsps 1 -out blast_max1.xml.txt
## only the xml format reveals "no hits found" genes
grep -B 5 "No hits found" blast_max1.xml.txt |grep "Query" > nohits.max1.txt
## counting
less nohits.max1.txt |wc -l
#output: 400
awk '{print $2}' nohits.max1.txt | grep -v "|" >nohits.max1.genes.txt
# 397 - 3 of those were dups
##----------------------------------------------------------------------
## unique hap1 genes in blast output
awk '{print $1}' blast_max1.tsv | sort| uniq |wc -l
# 449, 122 entries are duplicate

## total genes in putative missing list 
#849

## numbers add up - out of the 849 genes, 400 have no hits, 449 have at least one hit, 122 of those have multiple hits in hap2 (on different genes)
##----------------------------------------------------------------------
awk '{ if ($3 < 90) { print } }' blast_max1.tsv > lowpercenthits_max1.txt
# 151 genes have low % ID hits
awk '{ if ($11 < 0.01) { print } }' blast_max1.tsv > goodEhits_max1.txt
# 663 (all hits) have good "e" lol
awk '{ if ($3 < 97) { print } }' blast_max1.tsv > reallylowpercenthits_max1.txt
# 240 genes have low ID % hits when we make it more stringent

