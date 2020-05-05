#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate, Hmisc)

#load otu table and metadata
metadata <- read_csv("../../data/metadata.csv") 
otu_table <- read_delim("../../data/DeMMO_NativeRock_noChimera_otuTable_withTaxa_46822.txt", delim = '\t', comment = "# ")


# import data -------------------------------------------------------------


# load summary report data 
directories <- list.dirs("/Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data/", full.names = T , recursive =T)
directories <- directories[str_detect(directories, "reports")]
files <- list.files(directories, full.names = T)


elements = tibble::tibble(File = files[str_detect(files, "element_summary")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=/reports/)(.*)(?=_element_summary[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File)  %>%
  left_join(metadata %>% dplyr::select(sample, sample_id))

cells = tibble::tibble(File = files[str_detect(files, "cell_summary")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=/reports/)(.*)(?=_cell_summary[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File) %>%
  left_join(metadata)

models = tibble::tibble(File = files[str_detect(files, "models.csv")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=/reports/)(.*)(?=_cell_distribution_models[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File) %>%
  left_join(metadata)

total_corr = tibble::tibble(File = files[str_detect(files, "total_correlation.csv")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=/reports/)(.*)(?=_total_correlation[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File) %>%
  left_join(metadata)




# cell stats --------------------------------------------------------------
cell_stats <- cells %>%
  ggplot(aes(substrate, value, shape = experiment_type, color = site)) +
  geom_point() +
  facet_wrap(~observation, scales = "free")


ann_vs_density <- cells %>%
  spread(observation, value, fill = 0) %>%
  ggplot(aes(`cell ANN`, `cell density (cells/cm^2)`, shape = substrate, color = site, label = coupon_id)) +
  geom_point()

#plot shows exponential relationship between density and ANN - first order control on "clustering" is density
plotly::ggplotly(ann_vs_density)


cells_vs_biogenic <- cells %>%
  spread(observation, value, fill = 0) %>%
  ggplot(aes(`coverage cell area (%)`, `coverage biogenic area (%)`, shape = substrate, color = site, label = coupon_id)) +
  geom_point()

plotly::ggplotly(cells_vs_biogenic)


# total correlation -------------------------------------------------------

total_corr_plot <- total_corr %>%
  gather(corr_type, corr_val, c(cell_cor,biogenic_cor)) %>%
  ggplot(aes(element, corr_val)) +
  geom_boxplot(aes(group = element)) +
  geom_point(aes(color =site, shape = substrate, label = coupon_id)) +
  coord_flip() +
  facet_wrap(~corr_type)

plotly::ggplotly(total_corr_plot)

# element color palette ---------------------------------------------------
# Set color palette

#palette source: https://sciencenotes.org/molecule-atom-colors-cpk-colors/

element_colors <- c("#FFFFFF", "#D9FFFF", "#CC80FF", "#C2FF00", "#FFB5B5", "#909090", "#3050F8",
                    "#FF0D0D", "#90E050", "#B3E3F5", "#AB5CF2", "#8AFF00", "#BFA6A6", "#F0C8A0",
                    "#FF8000", "#FFFF30", "#1FF01F", "#80D1E3", "#8F40D4", "#3DFF00", "#E6E6E6",
                    "#BFC2C7", "#A6A6AB", "#8A99C7", "#9C7AC7", "#E06633", "#F090A0", "#50D050",
                    "#C88033", "#7D80B0", "#C28F8F", "#668F8F", "#BD80E3", "#FFA100", "#A62929",
                    "#5CB8D1", "#702EB0", "#00FF00", "#94FFFF", "#94E0E0", "#73C2C9", "#54B5B5",
                    "#3B9E9E", "#248F8F", "#0A7D8C", "#006985", "#C0C0C0", "#FFD98F", "#A67573",
                    "#668080", "#9E63B5", "#D47A00", "#940094", "#429EB0", "#57178F", "#00C900",
                    "#70D4FF", "#FFFFC7", "#D9FFC7", "#C7FFC7", "#A3FFC7", "#8FFFC7", "#61FFC7",
                    "#45FFC7", "#30FFC7", "#1FFFC7", "#00FF9C", "#00E675", "#00D452", "#00BF38",
                    "#00AB24", "#4DC2FF", "#4DA6FF", "#2194D6", "#267DAB", "#266696", "#175487",
                    "#D0D0E0", "#FFD123", "#B8B8D0", "#A6544D", "#575961", "#9E4FB5", "#AB5C00",
                    "#754F45", "#428296", "#420066", "#007D00", "#70ABFA", "#00BAFF", "#00A1FF",
                    "#008FFF", "#0080FF", "#006BFF", "#545CF2", "#785CE3", "#8A4FE3", "#A136D4",
                    "#B31FD4", "#B31FBA", "#B30DA6", "#BD0D87", "#C70066", "#CC0059", "#D1004F",
                    "#D90045", "#E00038", "#E6002E", "#EB0026")
names(element_colors) <- c("H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
                           "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
                           "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo",
                           "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",
                           "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                           "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
                           "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
                           "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt")

# otu vs element correlation ----------------------------------------------


abundance_table <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>% #normalize to relative abundance 
  gather(sample_id, abundance, -`#OTU ID`) %>% 
  spread(`#OTU ID`, abundance)  %>% 
  inner_join(elements %>% select(-transect) %>% spread(element, overview)) %>%
  select(-sample_id) %>%
  column_to_rownames("sample")


otu_element_cor <- cor(abundance_table, method = "spearman")
otu_element_cor<-rcorr(as.matrix(abundance_table))

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

results <- flattenCorrMatrix(otu_element_cor$r, otu_element_cor$P)


# element composition -----------------------------------------------------

element_composition <- elements %>%
  left_join(metadata %>% select(sample, substrate)) %>%
  ggplot() +
  geom_bar(aes(substrate, overview, fill = element), stat = "identity", position = "fill") +
  scale_fill_manual(values = element_colors) +
  coord_flip()


# rock element composition NMDS -------------------------------------------



# community NDMS with chemistry  ------------------------------------------


