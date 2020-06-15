#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate, Hmisc, vegan, heatmaply, htmltools, raster)

#load otu table and metadata

metadata <- read_csv("../../../data/metadata.csv") 
otu_table <- read_delim("../../../data/DeMMO_NativeRock_noChimera_otuTable_withTaxa_46822.txt", delim = '\t', comment = "# ")
#xrf <- read_csv("../../../data/XRF.csv")
fluid_geochem <- read_csv("../../../data/site_geochem.csv") %>%
  mutate(sulfide = sulfide/1000) %>%
  gather(substrate, abundance, methane:DOC) %>%
  group_by(Site, substrate) %>%
  summarise(abundance = mean(na.omit(abundance))) %>%
  mutate(abundance = abundance/100) %>%#convert mg/L to %m/v
  pivot_wider(names_from = substrate, values_from = abundance, values_fill = list(abundance = 0)) %>%
  rename(site = "Site")

taxonomy <- otu_table %>%
  select(`#OTU ID`, taxonomy) %>%
  mutate(tax = gsub("Gammaproteobacteria; D_3__Betaproteobacteriales", "Betaproteobacteria; D_3__Betaproteobacteriales", taxonomy), #fix taxonomy for Beta's,
         taxonomy = str_remove_all(tax, "D_0__| D_1__| D_2__| D_3__| D_4__| D_5__| D_6__")) %>%
  separate(taxonomy ,sep=';',c("domain", "phylum", "class", "order", "family", "genus", "species"))

otu_norm <- otu_table %>%
  select(-taxonomy) %>%
  mutate_at(vars(-`#OTU ID`), funs(./sum(.)*100)) %>%
  gather(sample_id, abundance, `9.DeMMO3.steri.11June2019`:`17.DeMMO1.tube8.12Dec2019`) %>% 
  spread(key = `#OTU ID`,value = 'abundance') %>%
  left_join(metadata %>% select(sample_id) %>% distinct()) 


# import data -------------------------------------------------------------


# load summary report data 
directories <- list.dirs("../../../data", full.names = T , recursive =T)
directories <- directories[str_detect(directories, "DeMMO1|DeMMO3")]
files <- list.files(directories, full.names = T)


elements = tibble::tibble(File = files[str_detect(files, "element_summary")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=Dec2019/)(.*)(?=_element_summary[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File)  

cells = tibble::tibble(File = files[str_detect(files, "cell_summary")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=Dec2019/)(.*)(?=_cell_summary[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File) 

models = tibble::tibble(File = files[str_detect(files, "models.csv")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=Dec2019/)(.*)(?=_cell_distribution_models[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File) 

total_corr = tibble::tibble(File = files[str_detect(files, "total_correlation.csv")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=Dec2019/)(.*)(?=_total_correlation[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File) 



# bulk element composition from 70X scans ---------------------------------




# cell stats --------------------------------------------------------------
cell_stats <- cells %>%
  left_join(metadata %>% dplyr::select(coupon_id, substrate, experiment_type, site)) %>%
  ggplot(aes(substrate, value, shape = experiment_type, color = site)) +
  geom_point() +
  facet_wrap(~observation, scales = "free")


ann_vs_density <- cells %>%
  left_join(metadata %>% dplyr::select(coupon_id, substrate, experiment_type, site)) %>%
  spread(observation, value, fill = 0) %>%
  ggplot(aes(`cell ANN`, `cell density (cells/cm^2)`, shape = substrate, color = site, label = coupon_id)) +
  geom_point()

#plot shows exponential relationship between density and ANN - first order control on "clustering" is density
plotly::ggplotly(ann_vs_density)


cells_vs_biogenic <- cells %>%
  left_join(metadata %>% dplyr::select(coupon_id, substrate, experiment_type, site)) %>%
  spread(observation, value, fill = 0) %>%
  ggplot(aes(`coverage cell area (%)`, `coverage biogenic area (%)`, shape = substrate, color = site, label = coupon_id)) +
  geom_point()

plotly::ggplotly(cells_vs_biogenic)

ann_vs_sim <- cells %>%
  left_join(metadata %>% dplyr::select(coupon_id, substrate, experiment_type, site)) %>%
  spread(observation, value, fill = 0) %>%
  mutate(ann_dif = `cell ANN` - `mean random ANN`) %>%
  ggplot(aes(ann_dif, `cell density (cells/cm^2)`, shape = substrate, color = site, label = coupon_id)) +
  geom_point()

ann_data <- read_csv("../../../data/ann_pvals.csv") %>%
  left_join(metadata) %>%
  left_join(cells) %>%
  filter(!is.na(observation)) %>%
  pivot_wider(names_from = observation, values_from = value, values_fill = list(value=0)) %>%
  ggplot(aes(ANN_pval, `cell density (cells/cm^2)`, shape = substrate, color = site, label = coupon_id)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "gray") +
  theme(legend.position = "none") +
  geom_point()

dens_vs_elements <- cells %>%
  filter(observation == "cell density (cells/cm^2)") %>%
  inner_join(elements %>% 
               left_join(metadata %>% select(sample_id, coupon_id, site)) %>% 
               select(-transect) %>% 
               spread(element, overview)) %>%
  filter(site == "D3") %>%
  select(-sample_id, -observation, -site) %>%
  rename(cell_dens = "value") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("coupon_id") %>%
  select_if(function(col) sum(col >0) > (0.3*length(col))) 

res <- cor(dens_vs_elements)
round(res, 2)

#heatmaply(res[1, 2:ncol(res)])
heatmaply(res)
  

# models ------------------------------------------------------------------
selected_elements <- elements %>%
  group_by(coupon_id) %>%
  filter(transect >= 1) %>%
  select(coupon_id, element, transect)

selected_models <- models %>%
  mutate(n_elements = str_count(model, "[+]") + 1,
          model = str_extract(model, "(?<=ppm[(]cell_centroids_ppp ~ )(.*)(?<= [)])"),
         model = str_remove(model, " [)]"),
         model_vars = model) 
         #model_elements = if_else(str_detect(model, "[+]"), as.list(strsplit(model, "[+]")), list(model)[[id]])) %>%
selected_models <- selected_models %>%  
  separate(model_vars, letters[1:max(na.omit(selected_models$n_elements))],"[+]") %>%
  gather(elem_id, element, a:e, na.rm = T) %>% 
  full_join(selected_elements, by = c("coupon_id", "element")) %>%
  group_by(coupon_id, model) %>%
  filter(!any(is.na(transect))) %>%
  select(coupon_id, model, pval, n_elements) %>%
  distinct() %>%
  mutate(model_vars = model) 

max_elements <- max(na.omit(selected_models$n_elements))

selected_models <- selected_models %>%  
  separate(model_vars, letters[1:max_elements],"[+]") %>%
  group_by(coupon_id) %>%
  #filter_if(na.omit(vars(a:e)), any_vars(!. %in% unlist(element))) %>%
  filter(pval == min(na.omit(pval))) %>%
  filter(n_elements == min(n_elements)) %>%
  left_join(cells %>% filter(observation == "cell density (cells/cm^2)") %>% dplyr::select(coupon_id, value)) %>%
  rename(cell_dens_sq_cm = "value") %>%
  left_join(metadata) 

model_plot <- selected_models %>%
  gather(elem_id, element, a:tail(letters[1:max_elements], n=1)) %>%
  filter(!is.na(element)) %>%
  group_by(site, element) %>%
  summarise(freq = n()) %>%
  ggplot(aes(reorder(element, freq), freq, fill=element)) + 
  scale_fill_manual(values = element_colors) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  theme(axis.title.y =  element_blank()) +
  theme(legend.position = "none") +
  facet_wrap(~site) 
  
model_plot <- plotly::ggplotly(model_plot)

model_pval_vs_dens <- selected_models %>%
  ggplot(aes(pval, cell_dens_sq_cm, color = site, shape = substrate)) +
  geom_point() +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "gray") +
  theme(legend.position = "none")

model_pval_vs_dens <- plotly::ggplotly(model_pval_vs_dens)

model_table <- selected_models %>% 
  select(site, substrate, model, pval, n_elements, cell_dens_sq_cm, coupon_id, experiment_type) %>%
  DT::datatable(options = list(lengthMenu = c(5, 10, 20), pageLength = 5))


browsable(
  tagList(list(
    tags$div(
      style = 'width:66%;display:block;float:left;',
      model_plot
    ),
    tags$div(
      style = 'width:34%;display:block;float:left;',
      model_pval_vs_dens
    ),
    tags$div(
      style = 'width:100%;display:block;float:left;',
      model_table
    )
  ))
)

# total correlation -------------------------------------------------------

total_corr_plot <- total_corr %>%
  left_join(metadata %>% select(coupon_id, site, substrate)) %>%
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
                    "#D90045", "#E00038", "#E6002E", "#EB0026", "black")
names(element_colors) <- c("H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si",
                           "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
                           "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo",
                           "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba",
                           "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                           "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
                           "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf",
                           "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "LE")

# otu vs element correlation ----------------------------------------------


element_composition <- elements %>%
  select(-transect) %>%
  left_join(metadata %>% select(coupon_id, substrate)) %>%
  group_by(substrate, element) %>%
  summarise(mean_abundance = mean(na.omit(overview))) %>%
  filter(!is.na(mean_abundance)) 

taxa_selector <- function(taxa_level, cutoff){
  trim = "D_0__Bacteria; |D_0__Archaea; |D_1__|D_2__|D_3__|D_4__|D_5__|D_6__"
    otu_norm  %>% 
    gather(`#OTU ID`, abundance, OTU_1:OTU_999) %>%
    group_by(sample_id, `#OTU ID`) %>%
    summarise(abundance = sum(abundance)) %>%
    left_join(taxonomy) %>%
    ungroup() %>%
    mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum),
           phylum = if_else(is.na(phylum), domain, phylum),
           taxa_level = taxa_level, 
           taxa = "",
           taxa = if_else(taxa_level == "phylum", phylum, taxa),
           taxa = if_else(taxa_level == "class", class, taxa),
           taxa = if_else(taxa_level == "order", paste(class, order, sep = "; "), taxa),
           taxa = if_else(taxa_level == "family", paste(phylum, class, order, family, sep = "; "), taxa),
           taxa = str_remove_all(taxa, "; NA"), 
           taxa = if_else(taxa=="NA", str_remove_all(tax, trim), taxa),
           taxa = if_else(taxa_level=="otu", `#OTU ID`, taxa)) %>% ##need to make this line dynamic
    #filter(!domain == "Unassigned") %>% #remove unassigned taxa 
    group_by(sample_id, taxa) %>%
    summarise(abundance = sum(abundance)) %>% 
    filter(abundance >= cutoff) %>%
    pivot_wider(names_from = taxa, values_from = abundance, values_fill = list(abundance=0)) 
}

n_fam <- taxa_selector("family", 1) %>%
  inner_join(elements %>% select(coupon_id) %>% left_join(metadata %>% select(sample_id, coupon_id, site))) %>%
  ungroup() %>%
  filter(site == "D1") %>%
  select(-sample_id, -coupon_id, -site) %>%
  #select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  select_if(function(col) sum(col >0) > (0.3*length(col))) %>%
  ncol()

fam_abundance_table <- taxa_selector("family", 1) %>%
  inner_join(elements %>% 
               left_join(metadata %>% select(sample_id, coupon_id, site)) %>% 
               select(-transect) %>% 
               spread(element, overview)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("coupon_id") %>%
  filter(site == "D1") %>%
  select(-sample_id, -site) %>%
  #select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  select_if(function(col) sum(col >0) > (0.3*length(col))) #select variables present in threshold percent of samples 

res <- cor(fam_abundance_table)
round(res, 2)

heatmaply(res[1:n_fam, (n_fam+1):ncol(res)], k_row = 3, k_col = 2)


# otu_element_cor <- cor(abundance_table, method = "spearman")
# otu_element_cor<-rcorr(as.matrix(abundance_table))
# 
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }
# 
# results <- flattenCorrMatrix(otu_element_cor$r, otu_element_cor$P)


# element composition -----------------------------------------------------

xrf_plot <- xrf %>%
  mutate(substrate = c("Homestake", "Ellison", "Poorman", "Yates")) %>%
  gather(element, overview, Ti:Ca) %>%
  select(substrate, element, overview) %>%
  group_by(substrate, element) %>%
  summarise(mean_abundance = (overview/1000000)*100) %>%
  mutate(data = "xrf")

transect_compositions <- elements %>%
  left_join(metadata %>% select(coupon_id, substrate)) %>%
  group_by(substrate, element) %>%
  summarise(mean_abundance = mean(na.omit(transect))) %>%
  filter(!is.na(mean_abundance)) %>%
  mutate(data = "sem_transect")
  
element_composition <- elements %>%
  left_join(metadata %>% select(coupon_id, substrate)) %>%
  group_by(substrate, element) %>%
  summarise(mean_abundance = mean(na.omit(overview))) %>%
  filter(!is.na(mean_abundance)) %>%
  mutate(data = "sem_overview") %>%
  bind_rows(transect_compositions) %>%
  bind_rows(xrf_plot) %>%
  mutate(data = factor(data, levels = c("sem_transect", "sem_overview", "xrf"))) %>%
  ggplot() +
  geom_bar(aes(substrate, mean_abundance, fill = element), stat = "identity") +
  scale_fill_manual(values = element_colors) +
  coord_flip() +
  facet_wrap(~data)

plotly::ggplotly(element_composition)


# rock element composition NMDS -------------------------------------------

xrf_data <- xrf %>%
  mutate(id = c("Homestake.xrf", "Ellison.xrf", "Poorman.xrf", "Yates.xrf")) %>%
  gather(element, overview, Ti:Ca) %>%
  group_by(id, element) %>%
  summarise(abundance = (overview/1000000)*100) %>%
  pivot_wider(names_from = element, values_from = abundance, values_fill = list(abundance = 0))



xeds_data <- elements %>%
  gather(type, abundance, overview:transect) %>%
  mutate(id = paste0(coupon_id, ".", type)) %>%
  dplyr::select(-coupon_id, -type) %>%
  mutate(abundance = if_else(is.na(abundance), 0, abundance)) %>%
  pivot_wider(names_from = element, values_from = abundance, values_fill = list(abundance = 0))

all_element_data <- xrf_data %>%
  bind_rows(xeds_data) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("id")

element_NMDS <- all_element_data %>%
  metaMDS(k=2)

element_vectors <-envfit(element_NMDS, all_element_data, perm=1000)


NMDS_coords <- element_NMDS[["points"]] %>%
  as_tibble(rownames = "id") %>%
  separate(id, c("coupon_id", "data"),"[.]") %>%
  left_join(metadata)

vector_coords <- data.frame(element_vectors[["vectors"]][["arrows"]]*sqrt(element_vectors[["vectors"]][["r"]])) %>%
  as_tibble(rownames = "element") %>%
  bind_cols(as_tibble(element_vectors[["vectors"]][["pvals"]])) %>%
  rename(pval = value) %>%
  dplyr::filter(!NMDS1 == 0 & !NMDS2 == 0)


NMDS_plot <- NMDS_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(shape = data, color=substrate, label = coupon_id)) +
  #stat_ellipse(data = NMDS_coords, aes(color=substrate), alpha = 0.5) +
  geom_segment(data=vector_coords,inherit.aes = FALSE, aes(x=0,xend=NMDS1,y=0,yend=NMDS2, label = element), alpha=0.3)+
  geom_text(data=vector_coords, aes(x=NMDS1,y=NMDS2, label = element) ,hjust = 1, vjust = 1,  size=2, color="black") +  
  theme(legend.key.size = unit(.5, "cm"))

#visualize interactive plot
plotly::ggplotly(NMDS_plot)


# community NDMS with chemistry  ------------------------------------------
NMDS_ord <- otu_norm %>%
  column_to_rownames("sample_id") %>%
  metaMDS(k=2)

#extract eigenvectors for each OTU
#taxa_vectors <-envfit(NMDS_ord, otu_norm, perm=1000)

#pull out ordination and vector coordinates for plotting
NMDS_coords <- NMDS_ord[["points"]] %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(metadata %>% filter(!coupon_id == "D3T15rep"))

# vector_coords <- data.frame(taxa_vectors[["vectors"]][["arrows"]]*sqrt(taxa_vectors[["vectors"]][["r"]])) %>%
#   as_tibble(rownames = "#OTU ID") %>%
#   bind_cols(as_tibble(taxa_vectors[["vectors"]][["pvals"]])) %>%
#   rename(pval = value) %>%
#   filter(pval == min(pval)) %>% #filter only for lowest pval, produces 477 OTUs
#   left_join(taxonomy) %>%
#   mutate(phylum = if_else(phylum == "Proteobacteria", class, phylum),
#          phylum = if_else(is.na(phylum), domain, phylum)) 



#NMDS plot with controls 
NMDS_plot <- NMDS_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(shape=substrate, color=site, label = sample_id)) +
  #geom_segment(data=vector_coords,inherit.aes = FALSE, aes(x=0,xend=NMDS1,y=0,yend=NMDS2, color=phylum, label = `#OTU ID`), alpha=0.3)+
  #geom_text(aes(label=paste0(site, '.', date),hjust = 1, vjust = 1),  size=2, color="black") +  stat_ellipse(data=NMDS_coords, aes(color=site)) +
  theme(legend.key.size = unit(.5, "cm"))

#visualize interactive plot
plotly::ggplotly(NMDS_plot)

fam_NMDS <- taxa_selector("family", 0) %>%
  inner_join(elements %>% 
               left_join(metadata %>% select(sample_id, coupon_id, site)) %>% 
               select(-transect) %>% 
               spread(element, overview)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("coupon_id") %>%
  select(-sample_id, -site) %>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  metaMDS(k=2)

fam_coords <- fam_NMDS[["points"]] %>%
  as_tibble(rownames = "coupon_id") %>%
  left_join(metadata )

fam_NMDS_plot <- fam_coords %>%
  ggplot(aes(MDS1, MDS2)) +
  geom_point(size=2, alpha=0.8, aes(shape=substrate, color=site, label = coupon_id)) +
  #geom_segment(data=vector_coords,inherit.aes = FALSE, aes(x=0,xend=NMDS1,y=0,yend=NMDS2, color=phylum, label = `#OTU ID`), alpha=0.3)+
  #geom_text(aes(label=paste0(site, '.', date),hjust = 1, vjust = 1),  size=2, color="black") +  stat_ellipse(data=NMDS_coords, aes(color=site)) +
  theme(legend.key.size = unit(.5, "cm"))

#visualize interactive plot
plotly::ggplotly(fam_NMDS_plot)



#example https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
data(dune.env)
attach(dune.env)
data(dune)
ord <- metaMDS(dune)
plot(ord)
plot(ord, disp="sites", type="n")
ordihull(ord, Management, col=1:4, lwd=3)
ordiellipse(ord, Management, col=1:4, kind = "ehull", lwd=3)
ordiellipse(ord, Management, col=1:4, draw="polygon")
ordispider(ord, Management, col=1:4, label = TRUE)
points(ord, disp="sites", pch=21, col="red", bg="yellow", cex=1.3)
ord.fit <- envfit(ord ~ A1 + Management, data=dune.env, perm=999)
plot(ord.fit)

