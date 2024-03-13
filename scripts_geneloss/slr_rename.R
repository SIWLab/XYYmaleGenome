slr_rename<- function(x){mutate(x,
   region=if_else(genome == "hap1" & chr =="X" & end<=71000000,"PAR1",
          if_else(genome == "hap1" & chr =="X" & (end>71000000 & end<=261000000),"OldSLR",
          if_else(genome == "hap1" & chr=="X" & (end>261000000 & end<=367000000),"NewSLR", 
          if_else(genome == "hap1" & chr=="X" & end>367000000,"PAR2",
          if_else(genome == "hap1" & chr=="A1","A1",
          if_else(genome == "hap1" & chr=="A2","A2",
          if_else(genome == "hap1" & chr=="A4","A4",
          if_else(genome == "hap2" & chr=="Y1" & end<=64000000,"NewSLR",
          if_else(genome == "hap2" & chr =="Y1" & 
                    (end>64000000 & end<=272000000),"OldSLR",
          if_else(genome == "hap2" & chr =="Y1" & 
                    (end > 272000000 & end <= 344000000),"PAR1",
          if_else(genome == "hap2" & chr=="Y2" & end<=121000000,"PAR2",
          if_else(genome == "hap2" & chr=="Y2"& 
                    (end>121000000 & end<=159000000),"NewSLR",
          if_else(genome == "hap2" & chr=="Y2" & end>159000000, "OldSLR",
          if_else(genome == "hap2" & chr=="A1","A1",
          if_else(genome == "hap2" & chr=="A2","A2",
          if_else(genome == "hap2" & chr=="A4","A4",NA)))))))))))))))))}

slr_rename_PARs<-function(x){mutate(x,
   region=if_else(genome == "hap1" & chr =="X" & end<=71000000,"PAR",
          if_else(genome == "hap1" & chr =="X" & (end>71000000 & end<=261000000),"OldSLR",
          if_else(genome == "hap1" & chr=="X" & (end>261000000 & end<=367000000),"NewSLR", 
          if_else(genome == "hap1" & chr=="X" & end>367000000,"PAR",
          if_else(genome == "hap1" & chr=="A1","A1",
          if_else(genome == "hap1" & chr=="A2","A2",
          if_else(genome == "hap1" & chr=="A4","A4",
          if_else(genome == "hap2" & chr=="Y1" & end<=64000000,"NewSLR",
          if_else(genome == "hap2" & chr =="Y1" & 
                    (end>64000000 & end<=272000000),"OldSLR",
          if_else(genome == "hap2" & chr =="Y1" & 
                    (end > 272000000 & end <= 344000000),"PAR",
          if_else(genome == "hap2" & chr=="Y2" & end<=121000000,"PAR",
          if_else(genome == "hap2" & chr=="Y2"& 
                    (end>121000000 & end<=159000000),"NewSLR",
          if_else(genome == "hap2" & chr=="Y2" & end>159000000, "OldSLR",
          if_else(genome == "hap2" & chr=="A1","A1",
          if_else(genome == "hap2" & chr=="A2","A2",
          if_else(genome == "hap2" & chr=="A4","A4",NA)))))))))))))))))}
summarisebyregion<-function(x,y){ x %>% tibble() %>%
    group_by(region) %>% summarise({{y}}:=n()) %>% 
    filter(!is.na(region))}