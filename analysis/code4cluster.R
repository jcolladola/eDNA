library(tidyverse)
#library(here)

library(dada2)

base_dir <- "/home/meg/jcollado/data"

list.samples <- tibble (sample = list.dirs(base_dir, full.names = T)) |> 
  filter (str_detect(sample, "Nicaragua")) |>
  filter (!str_detect(sample, "filtered")) |> 
  mutate (R1 = file.path(sample, "with_primer.R1.fastq"),
          R2 = file.path(sample, "with_primer.R2.fastq")) 



 list.samples |> 
  mutate (filtf1 = file.path(sample, "filtered/goodQ.R1.fastq"),
          filtf2 = file.path(sample, "filtered/goodQ.R2.fastq")) |> 
  mutate (outFs = pmap_dfr(.l= list (R1, filtf1, R2, filtf2),
                       .f = function(a, b, c, d) {
                         filterAndTrim(a,b,c,d,
                                       
                                        truncLen=c(100,175),
                                       maxN=0, maxEE=c(2,2), rm.phix=TRUE, multithread=T) %>% 
                           as.data.frame()
                       } )) -> list.samples 
 print("finished filtering")

error_r1 <- learnErrors(list.samples$filtf1[1])
error_r2 <- learnErrors(list.samples$filtf2[1])
list.samples |> 
  mutate(derepF1 = map(filtF1, derepFastq),                   # dereplicate seqs
          derepR1 = map(filtR1, derepFastq),
          dadaF1  = map(derepF1, ~ dada(.x, err = error_r1, multithread = T)),  # dada2
          dadaR1  = map(derepR1, ~ dada(.x, err = error_r2, multithread = T)),
          mergers = pmap(.l = list(dadaF1,derepF1, dadaR1,derepR1),                 # merge things
                         .f = mergePairs,
                         minOverlap = 10)) -> output.dada2

output.dada2 |> write_rds(file.path(base_dir, "dada2.rds"))
