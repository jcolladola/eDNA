################### BLAST 27/10/2024 ######################

library(tidyverse)
library(Biostrings)
library(writexl)
library(here)

blast_results <- read_tsv(here("data/resultados_run2_271024.txt"), 
                          col_names = c("qseqid", "sseqid", "staxids", "pident", 
                                        "length", "mismatch", "evalue", "bitscore"))
View(blast_results)

##### Modificación del archivo rankedlineage.
path_taxdump <- "~/Desktop/Datosparaconservarciclidos/Ultima_asignacion_general/taxdump/new_taxdump/"
lineage_tbl <- read_delim(paste0(path_taxdump,"rankedlineage.dmp"), delim = "|", col_names = c("TaxID", "Tax_name" , "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Superkingdom","col11"))

# Eliminar columna 11 si todos son NAs.
if(nrow(filter(lineage_tbl,!is.na(col11)))==0){
  lineage_tbl <- select(lineage_tbl, -col11)
}


# lineage_tbl <- lineage_tbl %>%
  mutate_all(stringr::str_trim) %>%
  mutate_if(is.character, list(~na_if(., ""))) %>%
  select(TaxID, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species, Tax_name)

#Specify strings that should be given a lower score (escape special characters with '\\' for regex)
bad_results <- c("uncultured","sp\\.","fungal", "environmental", "unidentified")

# Ponemos en filas diferentes las secuencias asociadas a más de un taxid, usando el separador ;

clean_blast <- blast_results %>%
  separate_rows(staxids, sep=";")
View(clean_blast)

# Juntamos el Blast con la taxonomía.

clean_blast <- left_join(clean_blast, lineage_tbl, by = c("staxids" = "TaxID"))
View(clean_blast)

# Eliminamos los resultados que contengan los bad_results.

clean_blast <- clean_blast %>%
  filter(!str_detect(Tax_name, paste(bad_results, collapse = "|")))

# Guardamos el archivo como csv y xls.

write.csv(clean_blast, "clean_blast_all_hits.csv", row.names = FALSE)
write_xlsx(clean_blast, "clean_blast_all_hits.xlsx")

### Filtrar por un evalue < 10^-5

blast_filtered_by_evalue <- clean_blast %>%
  filter(evalue < 1e-5)

View(blast_filtered_by_evalue)

# Nos quedamos con hits ordenados, dentro de cada qseqid, con menor evalue, mayor bitscore, mayor porcentaje de identidad y por último menor número de mismatch.

best_hits <- blast_filtered_by_evalue %>%
  group_by(qseqid) %>%
  arrange(evalue, desc(bitscore), desc(pident), mismatch, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

View(best_hits)

### Colapsar a distintos niveles taxonómicos los resultados con hits pobres (>10^-5)
filtered_qseqid <- blast_filtered_by_evalue %>%
  distinct(qseqid) %>%
  pull(qseqid) 
filtered_qseqid

hits_high_evalue <- clean_blast %>%
  filter(!qseqid %in% filtered_qseqid)
View(hits_high_evalue)

write_xlsx(hits_high_evalue, "blast_poor_hits.xlsx")

best_hits_high_evalue <- hits_high_evalue %>%
  group_by(qseqid) %>%
  arrange(evalue, desc(bitscore), desc(pident), mismatch, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

View(best_hits_high_evalue)

## Juntamos ambos blasts.

blasts_todos_resultados <- rbind(best_hits, best_hits_high_evalue)
View(blasts_todos_resultados)

write_xlsx(blasts_todos_resultados, "blast_todos_resultados.xlsx")


##### Juntamos con resultados de Midori.

library(readxl)
asignacion_midori <- read_excel("asignacion_midori_simplificada.xlsx")

# Eliminamos resultados con NA a nivel de phylum.

asignacion_midori_filtered <- asignacion_midori[!is.na(asignacion_midori$Phylum), ]

# Nos quedamos con los qseqids no presentes en los resultados de blasts.

asignacion_midori_nuevas_coincidencias_sin_blasts <- asignacion_midori_filtered[!asignacion_midori_filtered$query_id %in% best_hits$qseqid, ]
View(asignacion_midori_nuevas_coincidencias_sin_blasts)

# Juntamos la asignación de BLAST con los resultados de Midori para los que no hay coincidencias en blast.

midori_blast_merged <- merge(best_hits, asignacion_midori_nuevas_coincidencias_sin_blasts, by.x = "qseqid", by.y = "query_id", all = TRUE)

# Cambiamos el nombre de Origen a BLAST para los resultados procedentes de BLAST y así distinguirlos de Midori.

midori_blast_merged <- midori_blast_merged %>%
  mutate(Origen = ifelse(is.na(Origen), "BLAST", Origen))

View(midori_blast_merged)


### Procesamos resultados de MEGAN.

library(purrr)
library(dplyr)

megan_asignation <- read.csv("megan_taxonomy_with_metadata.csv", header = TRUE)

megan_asignation <- megan_asignation %>%
  select(c(Hash, Phylum, Class, Order, Family, Genus))

View(megan_asignation)

# Eliminamos las filas con todos NAs.

megan_asignation_filtered <- megan_asignation %>%
  filter(if_any(c(Phylum, Class, Order, Family, Genus), ~ !is.na(.))) %>%
  distinct(Hash, .keep_all = TRUE)

View(megan_asignation_filtered)

# Nos quedamos con las filas con Hashes nuevos (no identificados por BLAST sin colapsar ni por Megan)

megan_asignation_new_hashes <- megan_asignation_filtered[!megan_asignation_filtered$Hash %in% midori_blast_merged$qseqid, ]
View(megan_asignation_new_hashes)

nrow(midori_blast_merged)
n_distinct(midori_blast_merged$qseqid)

# Comprobamos si hay algún hash más compartido entre una y otra tabla (midori-blast y megan)

hashes_comunes <- megan_asignation_new_hashes %>%
  inner_join(midori_blast_merged, by = c("Hash" = "qseqid")) %>%
  select(Hash)

if (nrow(hashes_comunes) > 0) {
  print("Existen hashes comunes entre ambos objetos:")
  print(hashes_comunes)
} else {
  print("No hay hashes comunes entre ambos objetos.")
}

# Juntamos ambas tablas por hash=qseqid.

asignacion_completa_blast_midori_megan <- merge(midori_blast_merged, megan_asignation_new_hashes, by.x = "qseqid", by.y = "Hash", all = TRUE)

# Cambiamos el origen a MEGAN para los nuevos resultados.

asignacion_completa_blast_midori_megan <- asignacion_completa_blast_midori_megan %>%
  mutate(Origen = ifelse(is.na(Origen), "MEGAN", Origen))
View(asignacion_completa_blast_midori_megan)

nrow(asignacion_completa_blast_midori_megan)
n_distinct(asignacion_completa_blast_midori_megan$qseqid)

write_xlsx(asignacion_completa_blast_midori_megan, "asignacion_completa_blast_midori_megan_1.xlsx")

# Eliminamos los números de referencia de cada grupo.

asignacion_completa_blast_midori_megan_limpio <- data.frame(lapply(asignacion_completa_blast_midori_megan, function(col) gsub("_\\d+$", "", col)))
View(asignacion_completa_blast_midori_megan_limpio)
write_xlsx(asignacion_completa_blast_midori_megan_limpio, "asignacion_completa_blast_midori_megan_limpio.xlsx")

# Leemos el excel nuevamente, ya que hemos realizado unos cambios para limpiar el df (eliminar class_, order_...)

library(readxl)
asignacion_completa_final <- read_excel("asignacion_completa_blast_midori_megan_limpio.xlsx")
View(asignacion_completa_final)

# Juntamos el archivo con otro que contiene los metadatos.
metadatos <- read.csv("metadatos.csv")
View(metadatos)
metadatos <- metadatos[,c(2:7)]

colnames(asignacion_completa_final)[colnames(asignacion_completa_final) == "qseqid"] <- "Hash"

asignacion_completa_final_con_metadatos <- metadatos %>%
  left_join(asignacion_completa_final, by = "Hash")

View(asignacion_completa_final_con_metadatos)
write_xlsx(asignacion_completa_final_con_metadatos, "asignacion_completa_final_con_metadatos.xlsx")
write.csv(asignacion_completa_final_con_metadatos, "asignacion_completa_final_con_metadatos.csv")

# Con esta opción nos quedamos tanto con los hashes asignados taxonómicamente como con los que no.

# Para quedarnos solo con los hashes asignados utilizamos merge.

asignacion_completa_final_con_metadatos_hashes_encontrados <- merge(metadatos, asignacion_completa_final, by = "Hash")
View(asignacion_completa_final_con_metadatos_hashes_encontrados)
write_xlsx(asignacion_completa_final_con_metadatos_hashes_encontrados, "asignacion_completa_final_con_metadatos_hashes_encontrados.xlsx")
write.csv(asignacion_completa_final_con_metadatos_hashes_encontrados, "asignacion_completa_final_con_metadatos_hashes_encontrados.csv")

# Comprobamos el número de hashes únicos.
length(unique(asignacion_completa_final_con_metadatos$Hash))
length(unique(asignacion_completa_final_con_metadatos_hashes_encontrados$Hash))
