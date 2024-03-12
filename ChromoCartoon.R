### chromosome diagram
### use ggplot to make bars and regions
library(ggplot2)
library(tidyverse)
df<-read.csv("~/Dropbox/Positions_long.csv") %>%
  select(Genome:SLR) %>%
  mutate(regionsize = stop - start) 
ggplot(data =df, aes(x = regionsize,y = Chromosome,fill = SLR)) + geom_col()

ggplot() + geom_boxplot(df, mapping = aes(x = Chromosome, y = regionsize))
