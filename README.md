# Scripts and analyses for *R. hastatulus* XYY male genome paper

Initial proposed formatting for plots

```r
pubTheme <-
  theme(title = element_text(size=10), #usually 14
        text = element_text(size=8),
        plot.background = element_rect(fill="lightgrey"),
        strip.background = element_rect(linetype=0,linewidth=8,
                                        fill="grey"), #facet boxes
        strip.text = element_text(size=10), #usually 11
        legend.background = element_rect(fill="grey"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        panel.background = element_rect(fill="#FFFFFF"), #plot body
        panel.grid = element_blank() #plot body pt2
  )
pubColours <- c("#2d9da6","#99b700","#8a53b6","#bf4824","#152944")
```

## Input data

### Gene annotation with MAKER

### TE annotation with EDTA

## Analysis

### anchorwaveplots :: Figure 1
Make the Hap1:Hap2 alignment look nice

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

# License
The content of this project itself is licensed under the Creative Commons Attribution 4.0 license, and the underlying source code used to format and display that content is licensed under the MIT license.
