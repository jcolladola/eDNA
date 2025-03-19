
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggVennDiagram")

BiocManager::install("Biostrings")

install.packages("vegan")

library(dplyr)
library(ggplot2)
library(vegan)
library(Biostrings)
library(ggVennDiagram)

# Nos quedamos con el mejor alineamiento de results.out (output de BLAST).

# Leer el archivo de resultados de BLAST sin nombres de columnas
blast_results <- read.table("results.out", header = FALSE, sep = "\t")

# Asignar nombres a las columnas según el formato de salida estándar de BLAST
colnames(blast_results) <- c("query_id", "subject_id", "percent_identity", "alignment_length",
                             "mismatches", "gap_opens", "query_start", "query_end",
                             "subject_start", "subject_end", "e_value", "bit_score")

# Verificar las primeras filas para asegurar que las columnas están correctamente asignadas
head(blast_results)

# Filtrar los mejores alineamientos por cada query_id
best_alignments <- blast_results %>%
  group_by(query_id) %>%
  arrange(e_value, desc(percent_identity)) %>%
  slice(1) %>%
  ungroup()  # Desagrupar después de la operación

# Verificar las primeras filas de los mejores alineamientos
head(best_alignments)

write.csv(best_alignments, "best_alignments.csv", row.names = FALSE)

## Asignación taxonómica

# Leer el archivo FASTA
fasta_file <- readDNAStringSet("fish_sequences.fasta")

# Extraer los nombres de las secuencias (headers)
headers <- names(fasta_file)

# Extraer la información taxonómica de los headers
taxonomy_data <- data.frame(
  subject_id = sapply(strsplit(headers, " "), `[`, 1),
  genus = sapply(strsplit(headers, " "), `[`, 2),
  species = sapply(strsplit(headers, " "), `[`, 3),
  info = sapply(strsplit(headers, " "), `[`, 4),
  stringsAsFactors = FALSE
)
head(taxonomy_data)



# Unir los resultados de BLAST con la información taxonómica
annotated_alignments <- best_alignments %>%
  left_join(taxonomy_data, by = "subject_id")

# Verificar las primeras filas de los alineamientos anotados
head(annotated_alignments)

# Guardar los alineamientos anotados en un archivo CSV
write.csv(annotated_alignments, "annotated_alignments.csv", row.names = FALSE)

# Juntamos los alineamientos anotados con los números de lecturas.
annotated_alignments_fishes_by_sample <- annotated_alignments %>%
  left_join(asv_data_agg, by = c("query_id" = "Hash"))
head(annotated_alignments_fishes_by_sample)

write.csv(annotated_alignments_fishes_by_sample, file="annotated_alignments_fishes_by_sample.csv")

# Añadimos la info del lago.

lagos_info <- data.frame(
  ID_tubo = c("21F110A", "21F110B", "21F109A", "21F109B", "21F111A", "21F111B", 
              "21F114A", "21F114B", "21F104A", "21F104B", "21F112A", "21F112B", 
              "21F107A", "21F107B", "21F103A", "21F103B", "2F101A", "21F101B", 
              "21F105A", "21F105B", "21F108A", "21F108B", "19F101A", "19F101B", 
              "19F101C", "19F102A", "19F102B", "19F102C", "19F103A", "19F103B", 
              "19F103C", "19F104A", "19F104B", "19F104C", "19F105A", "19F105B", 
              "19F105C", "V016", "V019", "V020", "V037", "V047", "V073", 
              "V089", "V095", "Xi", "Ma2s", "A2", "Ni5", "AL6"),
  Lago = c("Apoyo", "Apoyo", "Apoyo", "Apoyo", "Masaya", "Masaya", 
           "Masaya", "Masaya", "Nicaragua", "Nicaragua", "Managua", "Managua", 
           "Managua", "Managua", "Asososca León", "Asososca León", "Asososca León", 
           "Asososca León", "Xiloá", "Xiloá", "Xiloá", "Xiloá", "Apoyo", "Apoyo", 
           "Apoyo", "Nicaragua", "Nicaragua", "Nicaragua", "Asososca León", 
           "Asososca León", "Asososca León", "Xiloá", "Xiloá", 
           "Xiloá", "Xiloá", "Xiloá", "Xiloá", 
           "Asososca León", "Asososca León", "Nicaragua", "Xiloá", "Nicaragua", "Xiloá", 
           "Apoyo", "Apoyo", "Xiloá", "Masaya", "Apoyo", "Nicaragua", "Asososca León")
)

annotated_alignments_fishes_by_sample_with_lake_info <- annotated_alignments_fishes_by_sample %>%
  left_join(lagos_info, by = c("sample_id" = "ID_tubo"))
head(annotated_alignments_fishes_by_sample_with_lake_info)
write.csv(annotated_alignments_fishes_by_sample_with_lake_info, "annotated_alignments_fishes_by_sample_with_lake_info.csv")

# Comprobamos que a cada muestra le corresponde un lago.

lagos_sample_id <- annotated_alignments_fishes_by_sample_with_lake_info %>%
  select(Lago, sample_id) %>%
  distinct() %>%
  arrange(Lago)
print(lagos_sample_id, n=100)

# Filtramos por un porcentaje de identidad >=85

annotated_alignments_fishes_by_sample_with_lake_info_high_identity <- annotated_alignments_fishes_by_sample_with_lake_info %>%
  filter(percent_identity >= 85) 

# Juntamos la columna que contiene género y especie.

annotated_alignments_fishes_by_sample_with_lake_info_high_identity_species <- annotated_alignments_fishes_by_sample_with_lake_info_high_identity %>%
  mutate(genus_species = paste(genus, species, sep = " ")) %>%
  select(-genus, -species)

# Obtenemos la lista de géneros y especies por lago.

generos_especies_por_lago <- annotated_alignments_fishes_by_sample_with_lake_info_high_identity %>%
  select(Lago, genus, species) %>%
  distinct() %>%
  arrange(Lago, genus, species)

generos_especies_por_lago <- generos_especies_por_lago %>%
  mutate(genus_species = paste(genus, species, sep = " ")) %>%
  select(-genus, -species)


write.csv(generos_especies_por_lago, "generos_especies_por_lago.csv")


# Creamos un diagrama de venn.

generos_especies_por_lago <- generos_especies_por_lago %>%
  filter(!is.na(genus_species) & !is.na(Lago))

especies_por_lago <- generos_especies_por_lago %>% group_by(Lago) %>%
  summarize(especies = list(genus_species))

especies_list <- setNames(especies_por_lago$especies, especies_por_lago$Lago)


venn <- ggVennDiagram(especies_list, label="count")

ggsave(filename = "venn_diagram.png", plot = venn)


# Calculamos diversidad alpha general para cada lago.


diversity_data <- annotated_alignments_fishes_by_sample_with_lake_info_high_identity_species %>%
  group_by(Lago, genus_species) %>%
  filter(!is.na(genus_species) & !is.na(Lago)) %>%
  summarize(nReads = sum(nReads)) %>%
  group_by(Lago) %>%
  summarize(shannon_index = diversity(nReads, index = "shannon"))

ggplot(diversity_data, aes(x = Lago, y = shannon_index, fill = Lago)) +
  geom_bar(stat = "identity") +
  labs(title = "Diversidad Alfa (Índice de Shannon) por Lago", x = "Lago", y = "Índice de Shannon") +
  theme_minimal()

ggsave(filename = "shannon_diversity_per_lake.jpg")


# Diversidad alpha por muestra, boxplots.

annotated_alignments_fishes_by_sample_with_lake_info_high_identity_species <- annotated_alignments_fishes_by_sample_with_lake_info_high_identity_species %>% 
  filter(!is.na(Lago))

ggsave(filename = "shannon_diversity_boxplots.jpg")
diversity_data_2 <- annotated_alignments_fishes_by_sample_with_lake_info_high_identity_species %>%
  group_by(sample_id, genus_species, Lago) %>%
  summarize(nReads = sum(nReads)) %>%
  group_by(sample_id, Lago) %>%
  summarize(shannon_index = diversity(nReads, index = "shannon"))

ggplot(diversity_data_2, aes(x = Lago, y = shannon_index, fill = Lago)) +
  geom_boxplot() +
  labs(title = "Diversidad Alfa (Índice de Shannon) por Lago", x = "Lago", y = "Índice de Shannon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
