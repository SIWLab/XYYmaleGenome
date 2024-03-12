setwd("~/Library/CloudStorage/Dropbox/XYYmaleGenome/SC_coverage/tinycov")
scmalefile<-"male_tinycov.txt"

scfemalefile<-"female_tinycov.txt"

library(tidyverse)

colnames<- c("chr","start","end","coverage")

scMale<-read_table(scmalefile, col_names = colnames) 
scFemale<-read_table(scfemalefile, col_names = colnames) 

#boundaries<-
## skip all the scaffolds we ignored inthe bamqc run

chrname<-c("X","A4","A1","A2","Y2:121000001-347874410","Y1:1-273000000")
chrlength<-c(483481970,168693024,466657150,320470792,226874410,273000000)


# need file with chromosmes and posistons wrt 

ggplot(scFemale, aes(start, coverage))+geom_line()+ facet_grid(~chr)
