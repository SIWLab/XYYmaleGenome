#R script
#library(devtools)
# detach("package:GENESPACE", unload = TRUE) # will return an error if GENESPACE isn't loaded
# if (!requireNamespace("devtools", quietly = TRUE))
#     install.packages("devtools")
# devtools::install_github("jtlovell/GENESPACE",force=TRUE)
# #install.packages()
#install.packages("tidyverse")
#install.packages("BiocManager")
#BiocManager::install(c("Biostrings", "rtracklayer"))
#install.packages(c("igraph","dbscan"))

##---- load dependencies! ----
library(tidyverse)
library(GENESPACE)
library(data.table)
## --- setwd ---
wd <- "/ohta2/bianca.sacchi/genespacev1_rumex/xyymaleall_final"
setwd(wd)
## --- run GENESPACE ---
## parsedPaths
#genomeRepo <-"/ohta2/bianca.sacchi/genespacev1_rumex/genomeRepo"
#parsedPaths <- parse_annotations(
    # rawGenomeRepo = genomeRepo, 
    #  genomeDirs = c("sal_genes","hap1_genes","hap2_genes","texas"),
    #  genomeIDs = c("salicifolius","hap1","hap2","texas"),
    #  gffString = "gff",
    #  faString = "fasta",
    #  presets = "none",
    #  genespaceWd = wd,
    #  gffIdColumn = "ID",
    #  gffStripText = "ID=",
    #  headerEntryIndex = 1,
    #  headerSep = " ",
    #  headerStripText = "-RA")
## init genespace
#gpar <- init_genespace(
#  genomeIDs = c("hap1","salicifolius","hap2","texas"),
#  onewayBlast = TRUE,
#  ploidy=1,
#  wd = wd, 
#  path2mcscanx ="~/bin/MCScanX-master")
#out <- run_genespace(gpar)
#saveRDS(out,"output_01_06_2023.rds")

## ---reload output from GENESPACE run (June 1, 2023) ---
#out<-readRDS("output_01_06_2023.rds")
load("/ohta2/bianca.sacchi/genespacev1_rumex/xyymaleall_final/results/gsParams.rda")
## original, unmodified plots are in the wd -> riparian folder

## --- plot for figure S3 ---

## renaming dictionary for Sal chr and Texas
## sal
salkey<-readr::read_delim("/ohta2/Rumex/Dovetail_R_salicifolius/scaffolded_assembly/input_to_output.txt",col_names =c("old","new"),show_col_types = FALSE)
salkey 
salkey<-salkey %>% mutate(fix=str_replace(old, "\\=", "\\_")) %>%
  mutate(.,fix=str_replace(fix, "\\;", "\\_")) %>% select(fix,new) %>%
  mutate(.,new =str_replace(new,"Scaffold\\_",""))
salkey<-as.data.frame(salkey)
## tx
fix<-c("LG1","LG3","LG5","LG2","LG4")
new<-c("A1","A2","A3","A4","XY")
texkey<-data.frame(fix,new)
#head(texkey)
#head(salkey)
mainkey<-rbind(salkey,texkey)
#print(mainkey)
## --- sal, hap1, hap2, tx - colored by tx ---

## invert chr
invchr <- data.frame(
  genome = c("hap1", "hap1","texas","texas","hap2","hap2"), 
  chr = c("A4","A2","LG5","LG4","Y1","Y2"))
# nice background
ggthemes<-ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"),
  axis.text.y = element_text(size = 20))

ripDat_2A<-plot_riparian(
  gsParam = gsParam,
  refGenome = "texas",
  genomeIDs = c("hap2","texas","hap1","salicifolius"),
  forceRecalcBlocks = FALSE, 
  invertTheseChrs = invchr,
  addThemes = ggthemes,
  chrLabFontSize = 14,
  chrExpand	=0.6,
  chrFill = "lightgrey",
  braidAlpha = 0.75,
  chrLabFun = function(x) matchmaker::match_vec(x,dictionary = mainkey),
  customRefChrOrder = c("LG1","LG3","LG2","LG4","LG5"))
ggsave(filename="Plot2a_04oct2023.pdf",height=5,width=10,units="in",device="pdf")

  #scalePlotWidth=10,
  #scalePlotHeight=5)
# ripDat_2A<-plot_riparian(
#   gsParam = out,
#   refGenome = "hap1",
#   genomeIDs = c("hap1","hap2","texas","salicifolius"),
#   forceRecalcBlocks = TRUE, 
#   invertTheseChrs = invchr,
#   addThemes = ggthemes,
#   chrLabFontSize = 6,
#   chrFill = "lightgrey",
#   braidAlpha = 0.75,
#   chrLabFun = function(x) matchmaker::match_vec(x,dictionary = mainkey),
#   #customRefChrOrder = c("LG1","LG3","LG2","LG4","LG5"),
#   scalePlotWidth=10,
#   scalePlotHeight=5)
#p1 <- ripDat_2A$plotData$ggplotObj
#p2 <- p1 + 
#  scale_y_discrete(NULL,limits = c("R. hastatulus XYY male \nPaternal haplotype", "R.hastatulus XY male", "R. hastatulus XYY male \nMaternal haplotype", "Rumex salicifolius"))
cowplot::plot_grid(ripDat_2A,ripDat_2B)
?grob
library("ggplotify")

#ggsave(filename="Plot2a_14Aug2023.pdf",height=5,width=10,units="in",device="pdf")
#ggsave(filename="Plot2a_18July2023.pdf",height=5,width=10,units="in",device="pdf")
#ggsave(filename="Plot2a_07June2023.tiff",device="tiff")
#ggsave(filename="Plot2a_07June2023.jpg",device="jpg")


## --- sex chr plot only ----
invchr2 <- data.frame(
  genome = c("hap2","hap2"), 
  chr = c("Y1","Y2"))
PARs <- data.table(
   genome = c("hap1","hap1","hap1","hap1"),
   chr = c("X","X","X","X"),
   start = c(0,71e6,261e6,366167432),
   end = c(71e6,261e6,366167433,483154172),
   color = c("#202020","darkblue","darkorange","#202020"))
## hap1,hap2,sal
## issue with PARs - one large block b/w X and Y2 is colored all orange when 
## it should be grey
## slight overlap in positions

### what do the blocks look like in the prev plot?
# df2a<-ripDat_2A$blks
# xyblks<-df2a %>% filter(genome1 == "hap1" & genome2 == "hap2" & chr1 == "X" & chr2 =="Y2")
# View(xyblks)

### trying to plot with useregions = false
### creating a plot function
plot2b<-function(region){plot_riparian(
  gsParam = gsParam,
  syntenyWeight = 1,
  refGenome = "hap1",
  genomeIDs = c("hap1","hap2","salicifolius"),
  forceRecalcBlocks = TRUE, 
  useRegions = FALSE,
  invertTheseChrs = invchr2,
  addThemes = ggthemes,
  chrFill = "lightgrey",
  chrLabFontSize = 14,
  chrExpand	=0.6,
  braidAlpha = 0.75,
  highlightBed = region,
  backgroundColor = NULL, 
  chrLabFun = function(x) matchmaker::match_vec(x,dictionary = mainkey))
  }

#ripDat_2B<-plot2b(PARs)

#ggsave("2b_aug30.pdf")
#fs::file_show(fs::path(wd, "2b_aug30.pdf"))



#362 blk 77 ends, which is still in region 42
#### which blocks in PAR boundary?
# df<-ripDat_2B$blks
# df1<-subset(
#   df, genome1=="hap1" & genome2 == "hap2" &
#   chr1 == "X" & chr2 == "Y2")
# df1
# write.table(df1,"XY2blks.txt")
## 42_! n;pcls - one ends at 667, the other ends at 366

### ok how does the regions file work
# myhits<-read_refGenomeSynHits(gsParam=gsParam, refGenome = "hap1")
# myhits2<-filter(myhits, chr1 == "X",chr2 =="Y2") %>% filter(start1 >300e6) %>% select(start1,end1,id1,regID,blkID)
# View(myhits2)
# overlap<-myhits2%>%filter(grepl("42",regID))
# View(overlap)
# write.csv(overlap,"geneparoverlaps.csv")
PARs2 <- data.table(
   genome = c("hap1","hap1","hap1","hap1"),
   chr = c("X","X","X","X"),
   start = c(0,71e6,261e6,364e6),
   end = c(71e6,261e6,364e6,483154172),
   color = c("#202020","darkblue","darkorange","#202020"))
plot2b(PARs2) # yas
ggsave("2b_final_04oct.pdf",height = 5,width=10)
fs::file_show(fs::path(wd, "2b_final_aug30.pdf"))
##367MB boundary - block49_1/51_1 on either side. block 367MB very slightly over - colored dark orange maybe shouldn't be
## adjusting in parfix below
# df2<-subset(
#   df, genome1=="hap2" & genome2=="salicifolius" &
#   chr1 == "Y2" & 
#   chr2 =="ScaOOZI_5_HRSCAF_7" &
#   (color=="#202020"|color=="darkorange")) 
# View(df2)
## PAR from fixed pos -  121 - first orange block (170_1) encompasses that (120.8Mb-121.8MB)
## all good
# <-function(region){plot_riparian(
#   gsParam = gsParam,
#   syntenyWeight = 1,
#   refGenome = "hap1",
#   genomeIDs = c("hap1","hap2","salicifolius"),
#   forceRecalcBlocks = TRUE, 
#   useRegions = TRUE,
#   invertTheseChrs = invchr2,
#   addThemes = ggthemes,
#   chrFill = "lightgrey",
#   braidAlpha = 0.75,
#   highlightBed = region,
#   backgroundColor = NULL, 
#   chrLabFun = function(x) matchmaker::match_vec(x,dictionary = mainkey))
#   }
# plot2b_reg(PARs2)
## plotting larger regions bad :(


## PAR2??
# df3<-subset(df, 
# genome1 == "hap1" & genome2 == "salicifolius" &
# chr1 == "X" & chr2 == "ScaOOZI_6_HRSCAF_9" & 
# (color == "darkblue" | color == "#202020"))
# View(df3)
## seems to coord with PAR1 boundary 71MB
## Y?
# df4<-subset(df, 
# genome1 == "hap2" & genome2 == "salicifolius" &
# chr1 == "Y1" & chr2 == "ScaOOZI_6_HRSCAF_9" & 
# (color == "darkblue" | color == "#202020"))
# View(df4)
# View(subset(df, blkID == "hap2_vs_salicifolius: 68_1"))
## checks out!
## still thin thread of scaffold_4/old y1 homology not on X
# df5<-subset(df, 
# genome1 == "salicifolius" &
# chr1 == "ScaOOZI_6_HRSCAF_9" & 
# (color == "darkblue"))
# View(df5)

# PARfix <- data.table(
#   genome = c("hap2","hap2","hap2",
#         "hap2","hap2","hap2"),
#   chr = c("Y1","Y1","Y1",
#           "Y2","Y2","Y2"),
#   start = c(0,64e6,272e6,
#           0,121e6,159e6),
#   end = c(64e6,272e6,345e6,
#           121e6,159e6,350e6),
#   color = c("darkorange","darkblue","#202020",
#             "#202020","darkorange","darkblue"))

# ripDat_2B_parfix<-plot_riparian(
#   gsParam = out,
#   syntenyWeight = 1,
#   refGenome = "hap1",
#   genomeIDs = c("hap1","salicifolius","hap2"),
#   forceRecalcBlocks = TRUE, 
#   invertTheseChrs = invchr2,
#   addThemes = ggthemes,
#   chrFill = "lightgrey",
#   braidAlpha = 0.75,
#   highlightBed = PARfix,
#   backgroundColor = NULL, 
#   chrLabFun = function(x) matchmaker::match_vec(x,dictionary = mainkey))
# #ggsave("2b_test_fix.pdf")

# df<-ripDat_2B_parfix$blks
# df1<-subset(
#   df, genome1=="hap1" & genome2=="salicifolius" &
#   chr1 == "X" & 
#   chr2 =="ScaOOZI_5_HRSCAF_7" &
#   (color=="#202020"|color=="darkorange")) 
# View(df1)

### add all positions to plot?????
# head(PARfix)
# head(PARs)
# allPos<-rbind(PARs,PARfix)
# allPos2<-rbind(PARfix,PARs)


# ripDat_2B_parfix<-plot_riparian(
#   gsParam = out,
#   syntenyWeight = 1,
#   refGenome = "hap1",
#   genomeIDs = c("hap1","salicifolius","hap2"),
#   forceRecalcBlocks = FALSE, 
#   invertTheseChrs = invchr2,
#   addThemes = ggthemes,
#   chrFill = "lightgrey",
#   braidAlpha = 0.75,
#   highlightBed = allPos2,
#   backgroundColor = NULL, 
#   chrLabFun = function(x) matchmaker::match_vec(x,dictionary = mainkey))

## order in table doesn't really matter
#this will be the final pdf!

#ggsave("2B_21Aug2023.pdf",height =5,width=10)

###----- adding legend -------
# library(cowplot)
