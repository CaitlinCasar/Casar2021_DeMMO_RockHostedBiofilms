pacman::p_load(qiime2R, tidyverse, lubridate, randomcoloR, vegan)
metadata <-read_csv("../data/metadata.csv")
ncbi_errors <- read_delim("github/data/asvs/ncbi_submission_error", delim = "\t", col_names = F) %>%
  rename(Sequence_ID = "X1")

unrarefied_data <- read_qza("../DeMMO_NativeRock/data/asvs/Osburn2020-asv-table.qza")
unrarefied_table <-  data.frame(unrarefied_data$data, check.names = F) %>%
  rownames_to_column("Feature.ID") %>%
  mutate_at(vars(-Feature.ID), funs(./sum(.)*100)) %>%
  pivot_longer(-Feature.ID, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  mutate(sample_id = str_replace_all(sample_id, "_", ".")) %>%
  inner_join(metadata %>%filter(coupon_id != "D3T15rep")) %>%
  left_join(as_tibble(taxonomy$data)) %>%
  mutate(taxonomy = str_remove_all(Taxon, "d__| p__| c__|o__| f__| g__| s__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species"))

#get unique features for filtering fasta 
unique_features <- unrarefied_table %>%
  dplyr::select(Feature.ID) %>%
  distinct() %>%
  rename(`#asv` = Feature.ID)

write_delim(unique_features, "github/data/asvs/unique_features.txt", delim = "\t")

biosample_assignment <- read_delim("github/data/asvs/biosample_assignment.tsv", delim = "\t") %>%
  left_join(unrarefied_table %>% 
              dplyr::select(Feature.ID, sample_id) %>% 
              rename(Sequence_ID = "Feature.ID")) %>%
  mutate(biosample_accession = if_else(str_detect(sample_id, "DeMMO1"), "SAMN12684770", "SAMN12684819")) %>%
  select(-sample_id) %>%
  distinct() %>%
  group_by(Sequence_ID) %>%
  mutate(biosample_accession = paste0(biosample_accession, collapse = ",")) %>%
  distinct() %>%
  mutate(biosample_accession = str_remove(biosample_accession, ",.*")) %>%
  anti_join(ncbi_errors)#keep only first BioSample ID for submission purposes, this will be fixed later by NCBI staff

write_delim(biosample_assignment, "github/data/asvs/biosample_assignment.tsv", delim = "\t")
