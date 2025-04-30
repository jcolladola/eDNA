library(tidyverse)
library(insect)
library(taxonomizr)
library(here)

# Usage: Rscript parse_blast.r <blast_output>

params <- commandArgs(TRUE)
blast_output <- params[1]
output_folder <- dirname(blast_output)

# Extract file ID to suffix outputs (e.g., split_part01)
file_id <- tools::file_path_sans_ext(basename(blast_output))

message("Processing: ", blast_output)

# Read BLAST results
BLAST_results <- read_table(blast_output, 
    col_names = c("qseqid", "sseqid",  "pident", "length", "mismatch", "gapopen",
                  "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxid", "qlen"),
    col_types = cols())

# Expand multi-staxid hits
BLAST_results |> separate_rows(staxid, sep = ";") -> BLAST_results

# Load and clean lineage file
taxonomy_cols <-  c("phylum", "class", "order", "family", "genus", "species")
lineage_output <- file.path(output_folder, "lineages.txt")

all_lineages <- read_delim(lineage_output,
  col_names = c("taxID", taxonomy_cols),
  delim = ";",
  col_types = cols()) |>
  mutate(taxID = as.numeric(str_remove(taxID, "\t$"))) |>
  unite(full_tax, all_of(taxonomy_cols), sep = "%") |>
  filter (!str_detect(full_tax, "NA%NA%NA")) |> # Instead of removing all environemtal taxa, we remove those with no taxonomic information
  # filter(!str_detect(full_tax, regex("environmental|uncultured|unidentified|unknown|unclassified|unassigned|uncharacterized", ignore_case = TRUE))) |>
  separate(full_tax, into = taxonomy_cols, sep = "%") -> all_lineages_tibble

# Prep for LCA
BLAST_results |> 
  filter(length > 200) |>
  select(qseqid, pident, taxID = staxid) |>
  inner_join(all_lineages_tibble, by = "taxID") |>
  filter(!is.na(genus) | !is.na(family)) -> ready_to_roll

# Save max identities
BLAST_results |>
  filter(length > 200) |>
  group_by(qseqid) |> 
  semi_join(all_lineages_tibble, by = c("staxid" = "taxID")) |>
  summarise(
    max_pident = max(pident),
    max_length = max(qlen),
    .groups = "drop"
  ) -> max_pidents

write_tsv(max_pidents, file.path(output_folder, paste0("max_pidents_", file_id, ".tsv")))

# Custom LCA function
custom.lca <- function(df, cutoff = 90) {
  df %>%
    group_by(qseqid) %>%
    select(pident, all_of(taxonomy_cols)) %>%
    nest() %>%
    mutate(consensus = map_chr(data, function(.x) {
      if (any(.x$pident == 100)) {
        .x <- filter(.x, pident == 100)
      } else if (any(.x$pident > cutoff)) {
        .x <- filter(.x, pident > cutoff)
      }
      .x %>%
        select(-pident) %>%
        condenseTaxa() %>%
        paste(collapse = "%")
    })) %>%
    select(qseqid, consensus) %>%
    unnest(cols = c(consensus))
}

# Run LCA for each threshold
thresholds <- list("99" = 99, "97" = 97, "95" = 95, "90" = 90)

ids_thresholds <- map(thresholds, ~custom.lca(ready_to_roll, .x))

walk2(ids_thresholds, names(ids_thresholds), function(df, threshold) {
  write_tsv(df, file.path(output_folder, paste0("ids_threshold_", threshold, "_", file_id, ".tsv")))
})

# Output summary function
output.summary <- function(tibble) {
  tibble %>%
    rename(Hash = 1) %>%
    separate(consensus, into = taxonomy_cols, sep = "%", fill = "right") %>%
    transmute(final_rank = case_when(
      !is.na(species) ~ "species",
      !is.na(genus) ~ "genus",
      !is.na(family) ~ "family",
      TRUE ~ "Worse"
    ))
}

# Summarise all thresholds
map(ids_thresholds, output.summary) |>
  bind_rows(.id = "Threshold") |>
  pivot_wider(names_from = Threshold, values_from = final_rank) -> performance_id

write_tsv(performance_id, file.path(output_folder, paste0("performance_", file_id, ".tsv")))
