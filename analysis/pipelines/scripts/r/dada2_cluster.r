#!/usr/bin/env Rscript

msg <-"the script gets executed"

library (tidyverse)
library (dada2)
library (Biostrings)
library (digest)


print (msg)

version$version.string

arguments <- commandArgs(TRUE)

output_folder <- arguments[1] 

demult_folder <- arguments[2]



hashing <- arguments[4]



### Starting from here 

sample.map <- read_delim(file.path(output_folder,"sample_trans.tmp"),
                         col_names = c("Full_Id", "fastq_header","Sample"),
                         delim = "\t") |> 
  distinct()

head (sample.map)




## ----listing files

F1s <- sort(list.files(demult_folder, pattern="_FWD.1.fastq", full.names = TRUE))
F2s <- sort(list.files(demult_folder, pattern="_FWD.2.fastq", full.names = TRUE))
R1s <- sort(list.files(demult_folder, pattern="_REV.1.fastq", full.names = TRUE))
R2s <- sort(list.files(demult_folder, pattern="_REV.2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- str_replace(basename(F1s), "_FWD.1.fastq","")

#Now use only those that reflect a real sample

good.sample.names<-sample.names[sample.names %in% sample.map$fastq_header]
# Introduce here the biological counterpart of the fastq file

real.sample.name <- sample.map[,3][match(good.sample.names,sample.map$fastq_header),1]

# Now subset the number of files to be looked at

F1sgood<-F1s[sample.names %in% sample.map$fastq_header]
F2sgood<-F2s[sample.names %in% sample.map$fastq_header]
R1sgood<-R1s[sample.names %in% sample.map$fastq_header]
R2sgood<-R2s[sample.names %in% sample.map$fastq_header]

## ----filter and trim-----------------------------------------------------
filt_path <- file.path(output_folder, "/filtered") # Place filtered files in filtered/ subdirectory
dir.create(filt_path)
filtF1s <- file.path(filt_path, paste0(good.sample.names, "_F1_filt.fastq.gz"))
filtF2s <- file.path(filt_path, paste0(good.sample.names, "_F2_filt.fastq.gz"))
out_Fs <- filterAndTrim(F1sgood, filtF1s, F2sgood, filtF2s, truncLen=c(200,130),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=F) 

filtR1s <- file.path(filt_path, paste0(good.sample.names, "_R1_filt.fastq.gz"))
filtR2s <- file.path(filt_path, paste0(good.sample.names, "_R2_filt.fastq.gz"))
out_Rs <- filterAndTrim(R1sgood, filtR1s, R2sgood, filtR2s, truncLen=c(200,130),
                      maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=F)

## ----as tibble-----------------------------------------------------------
out_Fs_tibble<-tibble(file=dimnames(out_Fs)[[1]],
                      reads.in = out_Fs[,1],
                      reads.out=out_Fs[,2],
                      direction = "Fwd")
out_Rs_tibble<-tibble(file=dimnames(out_Rs)[[1]],
                      reads.in = out_Rs[,1],
                      reads.out=out_Rs[,2],
                      direction = "Rev")


## Reduce the dataset to those samples with more than 20 reads, in all pairs



filtF1s[out_Fs_tibble$reads.out>10] -> filtF1s_good
filtF2s[out_Fs_tibble$reads.out>10] -> filtF2s_good

filtR1s[out_Rs_tibble$reads.out>10] -> filtR1s_good
filtR2s[out_Rs_tibble$reads.out>10] -> filtR2s_good
## ----learning errors, echo=F---------------------------------------------
errF1 <- learnErrors(filtF1s_good, multithread=TRUE)
errF2 <- learnErrors(filtF2s_good, multithread=TRUE)
errR1 <- learnErrors(filtR1s_good, multithread=TRUE)
errR2 <- learnErrors(filtR2s_good, multithread=TRUE)

## ----dereplication, echo=F,message=FALSE---------------------------------
derepF1s <- derepFastq(filtF1s_good, verbose=TRUE)
derepF2s <- derepFastq(filtF2s_good, verbose=TRUE)
derepR1s <- derepFastq(filtR1s_good, verbose=TRUE)
derepR2s <- derepFastq(filtR2s_good, verbose=TRUE)

# Name the derep-class objects by the fastq header (that way we skip wrong naming practices)
names(derepF2s)<- names(derepF1s) <- str_remove_all(names(derepF1s),"_F1_filt.fastq.gz" )
names(derepR2s)<- names(derepR1s) <- str_remove_all(names(derepR1s),"_R1_filt.fastq.gz" )

# rownames(out_Fs) <- names(derepF1s) <- names(derepF2s) <- real.sample.name$Sample
# 
# rownames(out_Rs) <- names(derepR1s) <- names(derepR2s) <- real.sample.name$Sample


## ----dadaing, message=FALSE----------------------------------------------
dadaF1s <- dada(derepF1s, err = errF1, multithread = TRUE, verbose = T)
dadaF2s <- dada(derepF2s, err = errF2, multithread = TRUE)
dadaR1s <- dada(derepR1s, err = errR1, multithread = TRUE)
dadaR2s <- dada(derepR2s, err = errR2, multithread = TRUE)


## ----merging pairs-------------------------------------------------------

mergersF <- mergePairs(dadaF1s, derepF1s, dadaF2s, derepF2s, verbose=T)

#Run a for loop that adds the number of unique reads that went into each ASV

for (j in 1:length(mergersF)){

  dadaF1s[[j]]@.Data[[2]] %>% rownames_to_column(var="forward") %>% dplyr::select("forward", "nunq") ->Fwd
  Fwd$forward<-as.integer(Fwd$forward)
  dadaF2s[[j]]@.Data[[2]] %>% rownames_to_column(var="reverse") %>% dplyr::select("reverse", "nunq") ->Rev
  Rev$reverse<-as.integer(Rev$reverse)

  mergersF[[j]] <- left_join(mergersF[[j]],Fwd, by="forward") %>% left_join(Rev, by="reverse") %>% mutate(nunq=pmin(nunq.x,nunq.y)) %>% dplyr::select(-nunq.x,-nunq.y)


}

mergersR <- mergePairs(dadaR1s, derepR1s, dadaR2s, derepR2s, verbose = T)

for (j in 1:length(mergersR)){

  dadaR1s[[j]]@.Data[[2]] %>% rownames_to_column(var="forward") %>% dplyr::select("forward", "nunq") ->Fwd
  Fwd$forward<-as.integer(Fwd$forward)
  dadaR2s[[j]]@.Data[[2]] %>% rownames_to_column(var="reverse") %>% dplyr::select("reverse", "nunq") ->Rev
  Rev$reverse<-as.integer(Rev$reverse)

  mergersR[[j]] <- left_join(mergersR[[j]],Fwd, by="forward") %>% left_join(Rev, by="reverse") %>% mutate(nunq=pmin(nunq.x,nunq.y)) %>% dplyr::select(-nunq.x,-nunq.y)

}


## ----merging F and R (1)-------------------------------------------------

seqtabF <- makeSequenceTable(mergersF)

as.data.frame(seqtabF) |> 
  rownames_to_column("sample") |>
  as_tibble() |> 
  pivot_longer(-sample, names_to = "sequence", values_to= "nr") -> seqtab.F.df



## ----merging F and R (2)-------------------------------------------------
seqtabR <- makeSequenceTable(mergersR)

reversed_sequences<-as.character(reverseComplement(DNAStringSet(colnames(seqtabR))))

summary (colnames(seqtabF) %in% reversed_sequences)

summary (reversed_sequences %in% colnames(seqtabF))

colnames(seqtabR)<-reversed_sequences


## ----merging F and R (3)-------------------------------------------------

as.data.frame(seqtabR) |> 
  rownames_to_column("sample") |>
  as_tibble() |> 
  pivot_longer(-sample, names_to = "sequence", values_to= "nr") -> seqtab.R.df



final.seqtab <- full_join(seqtab.F.df, seqtab.R.df, by = c("sample", "sequence"), suffix = c(".F", ".R")) |> 
  group_by(sample,sequence) |> 
  summarise(nr =sum(nr.F, nr.R,  na.rm = T)) |> 
  pivot_wider(names_from = sequence,
              values_from = nr) |> 
  as.data.frame() |> 
  column_to_rownames("sample") |> 
  as.matrix()


## ----RemovingChimeras, message=F-----------------------------------------

seqtab.nochim <- removeBimeraDenovo(final.seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)  



## ----tidying and writing-------------------------------------------------

as.data.frame(seqtab.nochim) |> 
  rownames_to_column("sample") |>
  as_tibble() |> 
  pivot_longer(-sample, names_to = "sequence", values_to= "nReads") |> 
  filter (nReads > 0) -> seqtab.nochim.df

ASV_file <- file.path(output_folder, "ASV_table.csv")





if (grepl ("yes", hashing, ignore.case = TRUE)) {
  
  conv_file <-  file.path(output_folder,"hash_key.csv")
  conv_fasta <- file.path(output_folder, "hash_key.fasta")
  
  seqtab.nochim.df |>
    distinct(sequence) |> 
    mutate (Hash = map_chr(sequence, digest::sha1)) -> hash_key
  
  write_csv(hash_key, conv_file)
  
  eDNAfuns::fasta.writer(hash_key,
                         sequence= sequence,
                         header = Hash,
                         file.out = conv_fasta)
  
  seqtab.nochim.df |>
    inner_join(hash_key)   |> 
    inner_join(sample.map |> 
                 select(sample = fastq_header,
                        sample_name = Sample)) |> 
    select(sample_name, Hash, nReads) |> 
    write_csv(ASV_file)

  
}else{ # what happens if you don´t want hashes 

seqtab.nochim.df |> 
  inner_join(sample.map |> 
               select(sample = fastq_header,
                      sample_name = Sample)) |> 
    select(sample_name, sequence, nReads)-> ASV_table

 

write_csv(ASV_table, ASV_file)
}




## ----output_summary------------------------------------------------------

getN <- function(x) sum(getUniques(x))

list(sapply(dadaF1s, getN), sapply(dadaR1s, getN),
     sapply(mergersF, getN),sapply(mergersR, getN),
     rowSums(seqtabF),rowSums(seqtabR),
     rowSums(final.seqtab), rowSums(seqtab.nochim)) |> 
  set_names("2a.denoised_F", "2b.denoised_R",
            "3a.merged_F", "3b.merged_R",
            "4a.tabled_F", "4b.tabled_R",
            "5.tabled_together", "6.nonchim") |> 
  bind_rows(.id= "step") |> 
  pivot_longer(-step, names_to = "sample_name", values_to = "nReads")-> tracka

out_Fs_tibble |> 
  mutate (sample_name= str_remove(file,"_FWD.1.fastq" )) |> 
  select(sample_name,`0a.Input_F`=reads.in, `1a.Filtered_F`=reads.out) |> 
  pivot_longer(-sample_name,
               names_to="step" ,
               values_to= "nReads") -> outFsum

out_Rs_tibble |> 
  mutate (sample_name= str_remove(file,"_REV.1.fastq" )) |> 
  select(sample_name,`0b.Input_R`=reads.in, `1b.Filtered_R`=reads.out) |> 
  pivot_longer(-sample_name,
               names_to="step",
               values_to= "nReads") -> outRsum

bind_rows(outFsum, outRsum, tracka) |> 
  write_csv( file.path(output_folder,"dada2_summary.csv"))


