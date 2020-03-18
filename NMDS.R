#load dependencies 
pacman::p_load(plotly, vegan, tidyverse, readr, plyr,lubridate)

#import otu table, metadata, taxonomy
otu_table <- read_delim("../data/DeMMO_NativeRock_noChimera_otuTable_withTaxa_46822.txt", delim = '\t', comment = "# ")
metadata <- read_csv("../data/metadata.csv") 

taxonomy <- otu_table %>%
  select(`#OTU ID`, taxonomy) %>%
  mutate(tax = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy), #fix taxonomy for Beta's,
         taxonomy = str_remove_all(tax, "D_0__| D_1__| D_2__| D_3__| D_4__| D_5__| D_6__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species"))


# NMDS on OTU data with vectors representing phyla  -----------------------
otu_norm <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>%
  gather(sample_id, abundance, `9.DeMMO3.steri.11June2019`:`17.DeMMO1.tube8.12Dec2019`) %>% 
  spread(key = `#OTU ID`,value = 'abundance') %>%
  right_join(metadata %>% select(sample_id)) %>%
  column_to_rownames("sample_id")


NMDS_ord <- otu_norm %>%
  metaMDS(k=2)

#extract eigenvectors for each OTU
taxa_vectors <-envfit(NMDS_ord, otu_norm, perm=1000)

#pull out ordination and vector coordinates for plotting
NMDS_coords <- NMDS_ord[["points"]] %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(metadata)

vector_coords <- data.frame(taxa_vectors[["vectors"]][["arrows"]]*sqrt(taxa_vectors[["vectors"]][["r"]])) %>%
  as_tibble(rownames = "#OTU ID") %>%
  bind_cols(as_tibble(taxa_vectors[["vectors"]][["pvals"]])) %>%
  rename(pval = value) %>%
  filter(pval == min(pval)) %>% #filter only for lowest pval, produces 477 OTUs
  left_join(taxonomy) %>%
  mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum),
         phylum = if_else(is.na(phylum), domain, phylum)) 



#NMDS plot with controls 
NMDS_plot <- NMDS_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(shape=site, color=substrate, label = sample_id)) +
  geom_segment(data=vector_coords,inherit.aes = FALSE, aes(x=0,xend=NMDS1,y=0,yend=NMDS2, color=phylum, label = `#OTU ID`), alpha=0.3)+
  #geom_text(aes(label=paste0(site, '.', date),hjust = 1, vjust = 1),  size=2, color="black") +  stat_ellipse(data=NMDS_coords, aes(color=site)) +
  theme(legend.key.size = unit(.5, "cm"))

#visualize interactive plot
plotly::ggplotly(NMDS_plot)


D1_otus <- abundance_table %>%
  left_join(taxonomy) %>%
  left_join(metadata) %>%
  filter(phylum %in% c("Latescibacteria", "Patescibacteria") & site == "D1")
