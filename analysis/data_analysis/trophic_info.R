# Get all trophic information at the class level

taxonomia |> 
  separate(taxa, into = taxonomy_cols, sep= "%") |> 
  distinct(phylum, class, order) |> 
  arrange(phylum, class, order) |>
  mutate (trophic_level = case_when(str_detect(phylum, "Rhodophyta") ~ "Producer")) |> 
  write_csv(here("data", "final_dataset", "trophic_info.csv"))


## Collapsing the sdataset at the taxa level

prune_taxa("species",x =  filter_samples)
rank_names(filter_samples)
## No sirve

# Importa los objetos: taxonomia, ASV y metadata_NA desde el Rmd

taxonomia |> 
  distinct(taxa) |>
  arrange(taxa) |> 
  rownames_to_column("taxon") |>
  mutate(taxon = glue::glue("Taxa_{str_pad(taxon,4,pad=0)}")) -> unique_tax

ASV |> 
  left_join(taxonomia) |> 
  left_join(unique_tax) |> 
  group_by(taxon, Sample) |> 
  summarise(nReads= sum(nReads)) -> unique_ASV

eDNAfuns::tidy2phyloseq(ASV_table = unique_ASV,
                        OTU_taxonomy = unique_tax,
                        metadata = metadata_no_NAs,
                        Taxa = "taxon", Sample = "Sample", Reads = "nReads") -> taxa_nica

