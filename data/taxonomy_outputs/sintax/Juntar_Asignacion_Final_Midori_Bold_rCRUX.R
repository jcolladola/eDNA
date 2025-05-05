#### Asignación final Midori, BOLD, rCRUX.

library(dplyr)
library(stringr)
library(readr)

# Cargar los archivos y añadir la columna DB
bold <- read_tsv("final_ids_BOLD.tsv") %>%
  mutate(DB = "BOLD")

midori <- read_tsv("final_ids_MIDORI.tsv") %>%
  mutate(DB = "MIDORI")

rcrux <- read_tsv("final_ids_rCRUX.tsv") %>%
  mutate(DB = "rCRUX")

# Unir los 3
all_assignments <- bind_rows(bold, midori, rcrux)
all_assignments <- as.data.frame(all_assignments)

# Contar cuántos niveles no son "NA" en cada final_id
all_assignments <- all_assignments %>%
  mutate(
    non_NA_levels = stringr::str_count(final_id, "%") + 1 - stringr::str_count(final_id, "NA"), # Sin generar lista
    DB_priority = case_when(
      DB == "MIDORI" ~ 1,
      DB == "BOLD" ~ 2,
      DB == "rCRUX" ~ 3
    )
  )

# Ahora elegir la mejor asignación
best_assignments <- all_assignments %>%
  dplyr::group_by(qseqid) %>%
  dplyr::arrange(desc(non_NA_levels), DB_priority) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

write_tsv(best_assignments, "best_taxonomic_assignments.tsv")


# Unimos con la tabla de ASVs con reads.

ASV_table <- read_csv("new_ASV_after_swarm.csv")

joined <- ASV_table %>%
  dplyr::inner_join(best_assignments, by = c("Hash" = "qseqid"))


#Agrupar por muestra y taxonomía, sumar nReads

grouped_table <- joined %>%
  dplyr::group_by(Sample, Hash, final_id, reason, DB, non_NA_levels, DB_priority) %>%
  dplyr::summarise(nReads = sum(nReads, na.rm = TRUE), .groups = "drop")

write_tsv(grouped_table, "taxonomic_assignments_with_hash_reads.tsv")


# ¿Cuántos hashes hay diferentes por muestra?

hashes_por_muestra <- joined %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(hashes_unicos = n_distinct(Hash))

print(hashes_por_muestra, n= 100)

# PROBLEMA CON ALGUNAS MUESTRAS. ¿ELIMINAMOS?

