library(insect)
library(eDNAfuns)
library(tidyverse)
library(here)


hashes <- fasta_reader(here("data/hash_key_nica.fasta"))

hashes |> 
  mutate(seqlen = str_length(seq)) -> hashes
  
hashes |> 
ggplot(aes(seqlen))+
  geom_density()


hashes |> 
  filter(seqlen>250) -> hashes_good_length
  
hashes_good_length |> 
group_by(seqlen) |> 
  tally(sort = T)


hashes_good_length |> 
  dplyr::slice(1:100) |> 
  mutate(translation_2 = map_int(seq, ~count_stop_codons(sequence = .x, codon = 2,format = "character", dictionary = 5))) |> 
  group_by(translation_2, seqlen) |> 
  tally()

hashes_good_length |> 
  # dplyr::slice(1:100) |> 
  mutate(translation_2 = map_int(seq, ~count_stop_codons(sequence = .x, codon = 2,format = "character", dictionary = 5))) |> 
  filter(translation_2 == 0) -> good_trans


good_trans |> 
  mutate(translation_2 = map_chr(seq, ~count_stop_codons(sequence = .x, 
                                                         codon = 2,
                                                         format = "character",
                                                         dictionary = 5,
                                                         return = "sequence"))) -> good_trans 
  ## Traer la información de abundancia para ver si las cosas son reales

ASV <- read_csv(here("data/ASV_nicaragua_after_vsearch.csv"))

ASV |> 
  group_by(Hash) |> 
  summarise(n_occ = n(),
            nR    = sum(nReads)) -> ASV

###

hashes |> 
  anti_join(hashes_good_length) |> 
  mutate (status = "short") -> summary_hashes

hashes |> 
  semi_join(hashes_good_length) |> 
  anti_join(good_trans) |> 
  mutate (status = "notranslation") |> 
  bind_rows(summary_hashes)-> summary_hashes


hashes |> 
  semi_join(good_trans) |> 
  mutate(status = "translation") |> 
  bind_rows(summary_hashes) -> summary_hashes

summary_hashes |> 
  dplyr::rename(Hash = header) |> 
  inner_join(ASV)-> summary_hashes

summary_hashes |> 
  ggplot(aes(x= n_occ, y = nR))+
  geom_boxplot(aes(group = n_occ))+
  geom_smooth()+
  # scale_y_log10()+
  # scale_x_log10()+
  facet_wrap(~status)

## preparar un fasta de AA para blast

good_trans |> 
  distinct(translation_2) |> 
  rownames_to_column("ID") |> 
  dplyr::slice(1:10) |> 
  fasta_writer(sequence = translation_2,
               header = ID,
               file.out = here("data/ten_AAS.fasta"))

## Clasificador insect

classifier <- read_rds(here("data/classifier_Leray.rds"))


summary_hashes |> 
  filter(n_occ > 20) |> 
  arrange(desc(nR)) |> 
  dplyr::slice(1:10) |> 
  pull(seq) |> 
  insect::classify( tree = classifier) -> test
  