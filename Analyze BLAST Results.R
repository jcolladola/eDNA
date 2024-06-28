library(dplyr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings)

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

# Asumamos que los headers tienen la información en el formato ">subject_id species genus family"
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
