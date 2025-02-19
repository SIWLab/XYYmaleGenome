# Scripts and analyses for *R. hastatulus* XYY male genome paper

[Phased assembly of neo-sex chromosomes reveals extensive Y degeneration and rapid genome evolution in Rumex hastatulus](https://academic.oup.com/mbe/article/41/4/msae074/7644656) . Sacchi, Humphries et al. 2024

### Ks Plots
- input: anchorwave dot plots, Ks output from coge
- comparing haplotype A and B assemblies of a rumex hastatulus xyy male genome
  
### Genespace
- synteny plots and pangene annotations
  
### Gene loss
- tabulating gametologs present and absent between X and Y

### NCmaleHaplotypeTEs
Purpose: understand where TEs are in the genome

Specific goals for this file:
- Windows MAKER and EDTA data (summarizes along chromosomes)
- For the X and Y chromosomes, categorizes areas based on when they lost recombination: Pseudoautosomal Regions (PAR), old sex-linked, and new sex-linked
- Visualises the chromosomal "landscape": proportion of 1Mb windows occupied by genes/TEs
- Visualises the fraction of sex chromosome regions occupied by TEs

### RelaxedSelection
Purpose: look for signs of TE insertion/retention near genes

Specific goals for this file:
- Combine output files from =bedtools coverage= (1kb upstream, 1kb downstream, introns, exons)
- Visualise the chromosomal landscape of TEs near 1:1:1 orthologous genes
- Visualise the fraction of genes with nearby TE insertions, categorized by chromosomal location
- Per gene, understand whether it's more likely to be near TEs if on the X or Ys
