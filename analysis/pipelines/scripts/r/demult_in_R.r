# An Rscript to take the info files from cutadapt and the sample map and return the 
# list of fastqs /  or seq IDs to subset

# Rscript demultiplex_from_info.r <info_file.1.txt> <info.cutadapt.demult.txt> <info_file.2.txt>  <info_file.3.txt> <sample_map>
# <info.cutadapt.demult.R2f> <info.cutadapt.demult.R2r> <LIB> <DEMULT_DIR>

arguments <- commandArgs(TRUE)

library(tidyverse)

parsing_colnames <- c("Seq.id",
                      "n_errors",
                      "start_adap",
                      "end_adap",
                      "seq_before_adap",
                      "matching_seq",
                      "seq_after_adap",
                      "adap_name",
                      "QScores_seq_before",
                      "QScores_matching",
                      "QScores_after")

sample.map <- read_table(arguments[5],
                         col_names= c("Ids", "combo", "name") )



cutadapt.r1 <- vroom::vroom(arguments[1], 
                          col_names = parsing_colnames[c(1,2,5,7,8,11)],
                          delim = "\t", col_types ="ci--c-cc--c") %>% 
  filter (n_errors!= -1)

cutadapt.r1$key <- str_extract(cutadapt.r1$Seq.id, "\\S+")
cutadapt.r1 |> select(key, adap_name, seq_before_adap,seq_after_adap,
                      QScores_after) -> cutadapt.r1

# cutadapt_no_adap_r1 <- cutadapt.r1 %>% 
#   filter(n_errors!= -1 & start_adap == 0) %>% select(Seq.id) %>% 
#   separate(Seq.id, into = "key", sep = " ") 

cutadapt_demultR1 <- vroom::vroom(file = arguments[2], 
                                # file ="pipeline_output/last_demult/Lib01_demultF.txt", 
                                col_names =parsing_colnames[c(1,8)],
                                delim = "\t",
                                col_types = "c------c") 
cutadapt_demultR1$key <- str_extract(cutadapt_demultR1$Seq.id, "\\S+")

cutadapt_demultR1 |>   
  select(key, adap_name) %>%
  right_join(cutadapt.r1 %>% select(key, Primer =  adap_name, seq_before_adap)) -> cutadapt_demultR1
 


cutadapt.r2 <- vroom::vroom(#file ="pipeline_output/last_demult/Lib01_demultF.txt",
                          arguments[3],
                          col_names = parsing_colnames[c(1,2,5,7,11)],
                          delim = "\t", col_types ="ci--c-cc--c")  %>% 
  filter (n_errors!= -1)

cutadapt.r2$key <- str_extract(cutadapt.r2$Seq.id, "\\S+")


cutadapt.r3 <- vroom::vroom(arguments[4], 
                          col_names = parsing_colnames[c(1,2,5,7,11)],
                          delim = "\t", col_types ="ci--c-cc--c") %>% 
  filter (n_errors!= -1)

cutadapt.r3$key <- str_extract(cutadapt.r3$Seq.id, "\\S+")
# cutadapt_no_adap_R2 <- cutadapt.r2 %>%
#  
#   filter(n_errors!= -1 & start_adap == 0) %>% 
#   separate(Seq.id, into = "key", sep = " ")
# 
# cutadapt_no_adap_R3 <- cutadapt.r3 %>% 
#   filter(n_errors!= -1 & start_adap == 0) %>% 
#   separate(Seq.id, into = "key", sep = " ") 


cutadapt_demultR2 <- vroom::vroom(file = arguments[6], 
  # file ="pipeline_output/last_demult/Lib01_demultR2.txt", 
  col_names =parsing_colnames[c(1,8)],
  delim = "\t",
  col_types = "c------c")  

cutadapt_demultR2$key <-  str_extract(cutadapt_demultR2$Seq.id, "\\S+")

cutadapt_demultR2 %>% 
  select(key, adap_name) |> 
  mutate(Primer = "FWD") %>% 
  right_join(cutadapt.r2 %>% select(key, seq_before_adap)) -> cutadapt_demultR2
 


cutadapt_demultR3 <- vroom::vroom(file = arguments[7], 
  # file ="pipeline_output/last_demult/Lib01_demultR3.txt", 
  col_names =parsing_colnames[c(1,8)],
  delim = "\t",
  col_types = "c------c") 

cutadapt_demultR3$key <-  str_extract(cutadapt_demultR3$Seq.id, "\\S+")

cutadapt_demultR3 %>% 
  select(key, adap_name) |> 
  mutate(Primer = "REV") %>% 
  right_join(cutadapt.r3 %>% select(key,  seq_before_adap), by = "key") -> cutadapt_demultR3
 # 

bind_rows(cutadapt_demultR2, cutadapt_demultR3) %>% 
  inner_join(cutadapt_demultR1, ., by = "key", suffix = c(".1", ".2")) %>% 
  mutate(Primer =case_when(is.na(Primer.1) ~ Primer.2,
                           is.na(Primer.2) ~ Primer.1,
                           Primer.1 == Primer.2 ~ Primer.1) ,
        
         Final_sample = case_when( adap_name.1 == adap_name.2 ~ adap_name.1,
                                   is.na(adap_name.1) & is.na(seq_before_adap.1)& !is.na(adap_name.2) ~ adap_name.2,
                                   !is.na(adap_name.1) & is.na(seq_before_adap.2) & is.na(adap_name.2) ~ adap_name.1)) %>% 
  
  inner_join(cutadapt.r1 %>% 
               select(key,
                      seq_after_adap.1 = seq_after_adap,
                      QScores_after.1 = QScores_after), by = "key") %>% 
  inner_join(bind_rows(cutadapt.r2, cutadapt.r3) %>% 
               select(key,
                      seq_after_adap.2 = seq_after_adap,
                      QScores_after.2 = QScores_after), by = "key") %>%
  
group_by(Final_sample, Primer) %>%
  
  nest() %>%
#   ungroup() %>% 
#   slice(1) %>% 
#   unnest() -> temp
#   
# temp %>% 
#   select(Header = key,
#          Sequence = seq_after_adap.1,
#          Quality = QScores_after.1) %>% 
#   microseq::writeFastq(fdta = .,
#                        out.file = "test.fastq.gz")


mutate(Writing.1 = pmap(.l = list(Primer, Final_sample, data),
                            .f= function(a,b,c){
                             c %>%
                                select(Header = key,
                                       Sequence = seq_after_adap.1,
                                       Quality = QScores_after.1) %>%

                              microseq::writeFastq(fdta = .,
                                                   out.file = file.path(arguments[9],
                                                                        paste0(arguments[8],"_", b,"_",a,".1.fastq.gz")))
                              
                              
                            }),
       Writing.2 = pmap(.l = list(Primer, Final_sample, data),
                         .f= function(a,b,c){
                           c %>%
                             select(Header = key,
                                    Sequence = seq_after_adap.2,
                                    Quality = QScores_after.2) %>%
                             
                             microseq::writeFastq(fdta = .,
                                                  out.file = file.path(arguments[9],
                                                                       paste0(arguments[8],"_", b,"_",a,".2.fastq.gz")))
                           
                           
                         }))


# We are keeping the ones that are shared, and the ones that didn't have an ID on the opposite Read

# cutadapt_demultR1 %>%
#   inner_join(cutadapt_demultR2, by = "adap_name") %>% 
#   inner_join(cutadapt_demultR3, by = "adap_name") %>% 
#   mutate(AgreedFwd = map2 (data.x, data.y, ~ inner_join(.x, .y)),
#          AgreedRev = map2 (data.x, data  ,  ~ inner_join(.x, .y)),
#          Missing.oneF = map2(data.x, data.y, ~ anti_join(.x, .y) %>% 
#                                inner_join(cutadapt_no_adap_r1 ) %>% 
#                                bind_rows(anti_join(.y,.x) %>% 
#                                            inner_join(cutadapt_no_adap_R2))),
#          Missing.oneR = map2(data.x, data, ~ anti_join(.x, .y) %>% 
#                                inner_join(cutadapt_no_adap_r1 ) %>% 
#                                bind_rows(anti_join(.y,.x) %>% 
#                                            inner_join(cutadapt_no_adap_R3))),
#          Fs = map2(AgreedFwd, Missing.oneF, bind_rows),
#          Rs = map2(AgreedRev, Missing.oneR, bind_rows)) -> demult_object
  


