#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate, Hmisc, vegan, heatmaply, htmltools)

#load asv data and metadata

asv_abundance <- read_delim("../data/asvs/rarefied_asv_abundance.txt", delim = '\t')

metadata <- read_csv("../data/metadata.csv") 

# import data -------------------------------------------------------------


# load summary report data 
directories <- list.dirs("../data/reports", full.names = T , recursive =T)
files <- list.files(directories, full.names = T)


models = tibble::tibble(File = files[str_detect(files, "significant_models.csv")]) %>%
  tidyr::extract(File, "coupon_id", "(?<=reports/)(.*)(?=/.*_significant_models[.]csv)", remove = FALSE) %>%
  mutate(size = file.info(File)$size) %>%
  filter(size > 100) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File) 


# bulk element composition from 70X scans ---------------------------------
bulk_element_composition <- read_delim("../data/bulkElementComposition.txt", delim="\t")


# models ------------------------------------------------------------------

best_models <- models %>%
  filter(n_elements <= 5) %>%
  #mutate(p.value = round(p.value, 3)) %>%
  arrange(coupon_id, p.value, n_elements, desc(Deviance)) %>% 
  group_by(coupon_id) %>% 
  slice(1) %>%
  left_join(metadata %>% dplyr::select(coupon_id, substrate, site)) 

best_models %>%   
  dplyr::select(substrate, coupon_id, model, p.value, Deviance) %>%
  write_csv("../data/model_summary.csv")


model_plot <- best_models %>%
  separate(model, letters[1:max(na.omit(best_models$n_elements))],"[+]") %>%
  ungroup() %>%
  dplyr::select(site, a,b,c,d,e) %>%
  pivot_longer(a:e, names_to = "id", values_to = "element") %>%
  filter(!is.na(element)) %>%
  group_by(site, element) %>%
  summarise(freq = n()) %>%
  ggplot(aes(reorder(interaction(site, element), freq), freq, fill=element)) +
  scale_fill_manual(values = element_colors) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(axis.title.y =  element_blank()) +
  theme(legend.position = "none") +
  facet_wrap(~site, scales = "free_y") +
  scale_x_discrete(labels = function(i) substring(i, 3))


# asv vs element correlation ----------------------------------------------

heatmapper <- function(site_name){
data <- asv_abundance %>%
  filter(abundance > 1 & site == site_name & !substrate %in% c("fluid", "sand")) %>% 
  mutate(family = if_else(is.na(family) | str_detect(family, "uncultured"), order, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), class, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), phylum, family),
         family = if_else(is.na(family) | str_detect(family, "uncultured"), domain, family),
         family = if_else(family == "GWF2-40-263", order, family)) %>%
         group_by(sample_id, family) %>%
  summarise(abundance = sum(abundance)) %>%
  pivot_wider(names_from = family, values_from = abundance, values_fill = list(abundance=0)) %>%
  left_join(metadata %>% filter(coupon_id!= "D3T15rep") %>% dplyr::select(sample_id, coupon_id)) %>%
  column_to_rownames("coupon_id") %>%
  dplyr::select(-sample_id) %>%
  select_if(function(col) sum(col >0) > (0.5*length(col)))

n_fam <- data %>% ncol()

chemistry <- bulk_element_composition %>% dplyr::select(substrate, element, rel_abundance) %>%
  pivot_wider(names_from = element, values_from = rel_abundance)

fam_abundance_table <- data %>%
  rownames_to_column("coupon_id") %>%
  left_join(metadata %>% filter(!coupon_id == "D3T15rep") %>% dplyr::select(substrate, coupon_id)) %>%
  inner_join(chemistry) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>%
  column_to_rownames("coupon_id") %>%
  dplyr::select(-substrate) %>%
  #select_if(negate(function(col) is.numeric(col) && sum(col) == 0)) %>%
  select_if(function(col) sum(col >0) > (0.5*length(col))) #select variables present in threshold percent of samples 

res <- cor(fam_abundance_table)

#write_csv(data.frame(res[1:n_fam, (n_fam+1):ncol(res)]) %>% rownames_to_column("taxa"), "../data/", site_name, "_taxa_corr.csv")

heatmaply(res[1:n_fam, (n_fam+1):ncol(res)], k_row = 3, k_col = 2, file = paste0("../data/", site_name, "_heatmaply_plot.png"))
}


plotly::subplot(heatmapper("D1"), heatmapper("D3"), nrows = 1, margin = 0.15)
