pacman::p_load(qiime2R, tidyverse, lubridate, randomcoloR, vegan)

asv_table <-read_qza("../data/asvs/Casar2021-asv-table-rarefied-47468.qza")
metadata <-read_csv("../data/metadata.csv")
taxonomy <- read_qza("../data/asvs/Osburn2020-taxonomy-Silva138.qza")

asv_seqs <- data.frame(asvs$data, check.names = F) %>%
  rownames_to_column("Feature.ID")

asv_abundance <- data.frame(asv_table$data, check.names = F) %>%
  rownames_to_column("Feature.ID") %>%
  mutate_at(vars(-Feature.ID), funs(./sum(.)*100)) %>%
  pivot_longer(-Feature.ID, names_to = "sample_id", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  mutate(sample_id = str_replace_all(sample_id, "_", ".")) %>%
  inner_join(metadata %>%filter(coupon_id != "D3T15rep")) %>%
  left_join(as_tibble(taxonomy$data)) %>%
  mutate(taxonomy = str_remove_all(Taxon, "d__| p__| c__|o__| f__| g__| s__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species"))

write_delim(asv_abundance, "../data/asvs/rarefied_asv_abundance.txt", delim = "\t")



n_asvs <- length(unique(asv_abundance$Feature.ID))


family_abundance <- asv_abundance %>%
  mutate(family = if_else(is.na(family) | str_detect(family, "uncultured"), order, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), class, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), phylum, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), domain, family),
         family = if_else(family == "GWF2-40-263", order, family),
         type = factor(experiment_type, levels = c("rep", "exp")),
         substrate = factor(substrate, levels = rev(c("fluid", "sand", "Yates", "Homestake", "Poorman", "Ellison")))) %>%
  #family = factor(family, levels = family_color$family)) %>% 
  group_by(site, `Date Sampled`, substrate, experiment_type, family) %>%
  summarise(abundance = sum(abundance))

less_abundant_taxa <- family_abundance %>%
  group_by(family) %>%
  filter(max(abundance) < 5) %>% #filter for families that represent less than 5% 
  group_by(site, `Date Sampled`, substrate, experiment_type) %>%
  summarise(abundance = sum(abundance)) %>%
  mutate(family = 'Less abundant taxa')


dominant_fams <- family_abundance %>%
  group_by(family) %>%
  filter(max(abundance) >= 5) %>%
  bind_rows(less_abundant_taxa)

write_delim(asv_abundance, "../data/asvs/dominant_fams.txt", delim = "\t")

#palette
palette <- distinctColorPalette(k = length(unique(dominant_fams$family)), altCol = FALSE, runTsne = FALSE)
names(palette) <- unique(dominant_fams$family)

tax <- asv_abundance %>%
  dplyr::select(domain, phylum, class, order, family) %>%
  distinct() %>%
  mutate(tax = paste(domain, phylum, class, order, family, sep = ";"),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), order, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), class, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), phylum, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), domain, family),
         family = if_else(family == "GWF2-40-263", order, family)) %>%
  dplyr::select(tax, family)
         
bar_plot <- dominant_fams %>%
  left_join(tax) %>%
  ggplot(aes(fill=family, y=abundance, x=experiment_type, label = tax)) +
  geom_bar(stat='identity', position='fill') +
  scale_fill_manual(values=palette) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() + 
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.25, "cm")) +
  facet_grid(cols = vars(site), rows = vars(substrate), switch = "y") + 
  guides(fill = guide_legend(ncol = 1))

#plotly::ggplotly(bar_plot)

phyla_plot <- asv_abundance %>%
  mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum)) %>%
  group_by(site, experiment_type, `Date Sampled`, substrate, phylum) %>%
  summarise(abundance = sum(abundance)) %>%
  ggplot(aes(fill=phylum, y=abundance, x=experiment_type, label = `Date Sampled`)) +
  geom_bar(stat='identity', position='fill') +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() + 
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.25, "cm")) +
  facet_grid(cols = vars(site), rows = vars(substrate), switch = "y") + 
  guides(fill = guide_legend(ncol = 1))

#plotly::ggplotly(phyla_plot)


# Get dominant asv’s for comparing to blast database ------------------------

dominant_asvs <- asv_abundance %>% 
  dplyr::filter(str_detect(Taxon, paste(dominant_fams %>% dplyr::select(family) %>% distinct() %>% dplyr::filter(!family == "Bacteria") %>% pull(), collapse = "|"))) %>%
  inner_join(asv_seqs)

unassigned_bacteria <- asv_abundance %>% 
  dplyr::filter(domain == "Bacteria" & is.na(phylum)) %>%
  inner_join(asv_seqs)

# rarefied alpha diversity ------------------------------------------------

rarefied_alpha_plot <- asv_abundance %>%
  group_by(sample_id, site, substrate) %>%
  summarise(value = n()) %>%
  ggplot(aes(reorder(substrate, value), value, color=substrate, group = substrate)) + 
  geom_line(size = 4) +
  stat_summary(fun.y=mean, geom="point", size=1, color="#262626") +
  coord_flip() +
  theme(legend.position = "none", 
        axis.title = element_blank()) +
  facet_wrap(~site)


# NMDS --------------------------------------------------------------------

NMDS <- asv_abundance %>%
  select(Feature.ID, sample_id, abundance) %>%
  pivot_wider(id_cols = sample_id, names_from = Feature.ID, values_from = abundance, values_fill = list(abundance = 0)) %>%
  column_to_rownames("sample_id") %>%
  metaMDS(k=2)

NMDS_coords <- NMDS[["points"]] %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(metadata %>% filter(coupon_id != "D3T15rep"))

#NMDS plot without controls 
NMDS_plot <- NMDS_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(color=site, shape = substrate, label = sample_id)) +
  #geom_text(aes(label=paste0(site, '.', date),hjust = 1, vjust = 1),  size=2, color="black") +
  #stat_ellipse(data=NMDS_coords, aes(color=interaction(site, substrate))) +
  theme(legend.key.size = unit(.5, "cm"))

#visualize interactive plot
#plotly::ggplotly(NMDS_plot)


#plot community data
top_panel <- gridExtra::grid.arrange(rarefied_alpha_plot + 
                                       theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                                       ggplot2::ggtitle("A. Alpha Diversity"), NMDS_plot + ggplot2::ggtitle("B. Beta Diversity"), nrow = 1)

gridExtra::grid.arrange(top_panel, bar_plot + 
                          ggplot2::ggtitle("C. Community Composition") +
                          theme(legend.key.size = unit(.5, "cm")), nrow = 2, heights = c(1, 1.5))

#simper analysis


asv_wide <- asv_abundance %>%
  select(Feature.ID, sample_id, abundance) %>%
  pivot_wider(id_cols = sample_id, names_from = Feature.ID, values_from = abundance, values_fill = list(abundance = 0)) %>%
  arrange(sample_id) 

asv_taxonomy <- asv_abundance %>% 
  dplyr::select(Feature.ID, domain, phylum, class, order, family, genus, species) %>%
  distinct()
  
D1_sim_meta <- metadata %>% 
  dplyr::filter(site == "D1") %>% 
  arrange(sample_id) %>% 
  column_to_rownames("sample_id") %>% 
  mutate(substrate = if_else(substrate %in% c("fluid", "sand"), substrate, "rock"))

(D1_sim <- with(D1_sim_meta, 
             simper(asv_wide %>% 
                      inner_join(metadata %>% 
                                   dplyr::filter(site == "D1") %>% 
                                   dplyr::select(sample_id)) %>%
                      column_to_rownames("sample_id"), substrate)))

D1_simper <- bind_rows(unclass(summary(D1_sim)), .id = "id") %>%
  rownames_to_column("Feature.ID") %>%
  separate(Feature.ID, into = c("Feature.ID", "row_id"), sep = "[.][.][.]") %>%
  left_join(asv_taxonomy) %>%
  mutate(tax = if_else(is.na(species) | str_detect(species,"uncultured") , genus, species),
         tax = if_else(is.na(tax) | str_detect(tax,"uncultured") | str_detect(tax, ".*[0-9].*"), family, tax),
         tax = if_else(is.na(tax) | str_detect(tax,"uncultured")  | str_detect(tax, ".*[0-9].*"), order, tax),
         tax = if_else(is.na(tax) | str_detect(tax,"uncultured")  | str_detect(tax, ".*[0-9].*"), class, tax),
         tax = if_else(is.na(tax) | str_detect(tax,"uncultured") | str_detect(tax, ".*[0-9].*"), phylum, tax),
         tax = if_else(is.na(tax) | str_detect(tax,"uncultured") , domain, tax))

D1_simper_plot <- D1_simper %>%
  filter(id == "rock_sand" & average > 0 & cumsum <= 0.75) %>%
  pivot_longer(ava:avb, names_to = "substrate", values_to = "abundance") %>%
  mutate(abundance = if_else(substrate == "ava", abundance * -1, abundance),
         substrate = if_else(substrate == "ava", "rock", "control")) %>%
  ggplot(aes(interaction(tax, row_id), abundance, fill = substrate)) +
  geom_bar(stat = "identity") +
  ylab("average abundance") +
  coord_flip() + 
  theme(axis.title.y = element_blank()) +
  ggtitle("A. D1 Rock vs. Control Biofilm Simper")

D1_simp_freq_plot <- D1_simper %>%
  filter(id == "rock_sand" & average > 0 & cumsum <= 0.75) %>%
  group_by(tax) %>%
  summarise(frequency = n()) %>%
  ggplot(aes(reorder(tax, frequency), frequency)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  ggtitle("B. Simper Taxa Frequency")

gridExtra::grid.arrange(D1_simper_plot, D1_simp_freq_plot, nrow = 1, widths = c(2, 1))

# rarefaction curves ------------------------------------------------------

# load collated alpha div data
data_path <- "../data/asvs/rarefaction"   # path to the data
files <- dir(data_path, pattern = "*.csv", full.names = T) # get file names

alpha_div = tibble::tibble(File = files) %>%
  tidyr::extract(File, "method", "(?<=/data/asvs/rarefaction/)(.*)(?=[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File, -site, -trip, -note, -date) %>% 
  rename(sample_id = `sample-id`) %>%
  pivot_longer(c(-sample_id, -method), names_to = "iteration", values_to = "value") %>%
  mutate(sample_id = str_replace_all(sample_id, "_", ".")) %>%
  inner_join(metadata %>% filter(coupon_id != "D3T15rep")) %>%
  separate(iteration, into = c("depth", "iteration"), sep = "_") %>%
  filter(!is.na(value)) %>%
  group_by(method, sample_id, depth, site, substrate) %>%
  summarise(value = mean(value, na.rm = T)) %>%
  mutate(depth = as.numeric(str_remove_all(depth, "depth-")))

rarefaction_curve <- alpha_div %>%
  filter(!site == "DuselD") %>%
  ggplot(aes(depth, value, group=sample_id, color=substrate)) +
  geom_line(aes(linetype = site)) + 
  geom_vline(xintercept = 47468, color = "gray", linetype = "dashed") +
  #scale_colour_discrete(guide = 'none') +
  theme(legend.position="bottom") +
  theme_grey() +
  facet_wrap( ~ method, ncol=2, scales = "free") +
  guides(colour = guide_legend(title = "substrate"))

#rarefaction_curve
#plotly::ggplotly(rarefaction_curve)

