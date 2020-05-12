#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate, randomcoloR)

#import otu table, metadata, taxonomy
otu_table <- read_delim("../data/DeMMO_NativeRock_noChimera_otuTable_withTaxa_46822.txt", delim = '\t', comment = "# ")
metadata <- read_csv("../data/metadata.csv") 

#create color dictionary for figure
phylum_color <- read_csv("../orig_data/phylum_color_data.csv")

phylum_color_dict <- phylum_color$hex.color
names(phylum_color_dict) <- phylum_color$full.name

taxonomy <- otu_table %>%
  select(`#OTU ID`, taxonomy) %>%
  mutate(tax = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy), #fix taxonomy for Beta's,
         taxonomy = str_remove_all(tax, "D_0__| D_1__| D_2__| D_3__| D_4__| D_5__| D_6__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species"))

#get rid of duplicate samples for a given experiment
samples_to_remove <- c()

#create taxon abundance table 
abundance_table <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>% #normalize to relative abundance 
  gather(sample_id, abundance, `9.DeMMO3.steri.11June2019`:`17.DeMMO1.tube8.12Dec2019`) 


taxon_abundance <- function(level, name){
    otu_table %>%
    select(`#OTU ID`, taxonomy) %>%
    mutate(taxonomy = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy),
           taxa = str_extract(taxonomy, level),
           taxa = if_else(is.na(taxa), taxonomy, taxa)) %>%
    right_join(abundance_table) %>%
    ungroup() %>%
    group_by(sample_id, taxa) %>%
    summarise(abundance = sum(abundance)) %>%
    left_join(metadata) %>% #add metadata
    group_by(site, `Date Sampled`, substrate, type, taxa) %>% 
    summarise(abundance = sum(abundance))
}

family_level <- "(.*)(?=; D_5__)"
class_level <- "(.*)(?=; D_3__)"
phylum_level <- "(.*)(?=; D_2__)"

less_abundant_taxa <- taxon_abundance(family_level, "family") %>%
  group_by(taxa) %>%
  filter(max(abundance) < 5) %>% #filter for families that represent less than 5% 
  group_by(site, `Date Sampled`, substrate, type) %>%
  summarise(abundance = sum(abundance)) %>%
  mutate(taxa = 'Less abundant taxa')


taxon_abundance_table <- taxon_abundance(family_level, "family") %>%
  group_by(taxa) %>%
  filter(max(abundance) >= 5) %>%
  bind_rows(less_abundant_taxa)

#palette
palette <- distinctColorPalette(k = length(unique(bar_plot$family)), altCol = FALSE, runTsne = FALSE)
names(palette) <- unique(unique(bar_plot$family)) 


#bar plot figure 
bar_plot <- taxon_abundance_table %>%
  ungroup() %>%
  mutate(taxonomy = str_remove_all(taxa, "D_0__| D_1__| D_2__| D_3__| D_4__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family")) %>%
  mutate(family = if_else(is.na(family) | str_detect(family, "uncultured"), order, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), class, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), phylum, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), domain, family),
         family = if_else(family == "GWF2-40-263", order, family),
         type = factor(type, levels = c("rep", "exp")),
         substrate = factor(substrate, levels = rev(c("fluid", "sand", "Yates", "Homestake", "Poorman", "Ellison")))) %>%
         #family = factor(family, levels = family_color$family)) %>% 
  group_by(site, `Date Sampled`, substrate, type, family) %>%
  summarise(abundance = sum(abundance) *100) %>%
  ggplot(aes(fill=family, y=abundance, x=type, label = `Date Sampled`)) +
  geom_bar(stat='identity', position='fill') +
  scale_fill_manual(values=palette) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() + 
  theme(axis.title = element_blank(),
        legend.title = ggplot2::element_blank(), 
        legend.text = ggplot2::element_text(size = 8),
        legend.key.size = unit(0.25, "cm")) +
  facet_grid(cols = vars(site), rows = vars(substrate)) + 
  guides(fill = guide_legend(ncol = 1))


#bar plot figure 
subset_bar_plot <- taxon_abundance_table %>%
  ungroup() %>%
  mutate(taxonomy = str_remove_all(taxa, "D_0__| D_1__| D_2__| D_3__| D_4__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family")) %>%
  mutate(family = if_else(is.na(family) | str_detect(family, "uncultured"), order, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), class, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), phylum, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), domain, family),
         family = if_else(family == "GWF2-40-263", order, family),
         type = factor(type, levels = c("rep", "exp")),
         substrate = factor(substrate, levels = c("fluid", "sand", "Yates", "Homestake", "Poorman", "Ellison"))) %>%
  #family = factor(family, levels = family_color$family)) %>% 
  filter(type == "exp") %>%
  group_by(site, `Date Sampled`, substrate, type, family) %>%
  summarise(abundance = sum(abundance) *100) %>%
  ggplot(aes(fill=family, y=abundance, x=substrate, label = `Date Sampled`)) +
  geom_bar(stat='identity', position='fill') +
  scale_fill_manual(values=palette) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() + 
  theme(axis.title = element_blank(),
        legend.title = ggplot2::element_blank(), 
        legend.text = ggplot2::element_text(size = 8),
        legend.key.size = unit(0.25, "cm")) +
  facet_grid(cols = vars(site)) + 
  guides(fill = guide_legend(ncol = 1))
