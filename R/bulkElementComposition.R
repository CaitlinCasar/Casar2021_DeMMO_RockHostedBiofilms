pacman::p_load(tidyverse, raster)

directories <- list.dirs("../data/DeMMO2019_lowMagXEDS", full.names = T , recursive =T)
directories <- directories[2:length(directories)]

#function to read in raster stack into dataframe, convert from rgb to binary
read_files <- function(file){
  element <- str_extract(file, "(?<=DeMMO2019_lowMagXEDS/)(.*)(?=/)")
  file %>% 
    stack() %>% 
    calc(sum) %>%
    as.data.frame(xy=T) %>%
    rename(!!quo_name(element) := 3) %>%
    mutate_at(vars(-x, -y), ~if_else(.==0, 0, 1))
}

#store names of coupons in a list
coupons <- c("D1T1exp", "D1T1rep", "D1T2exp", "D1T2rep",
             "D1T3exp", "D1T3rep", "D1T4exp", "D1T4rep",
             "D1T5exp", "D1T5rep", "D1T6exp", "D1T6rep",
             "D1T7exp", "D1T7rep", "D1T8exp", "D1T8rep")

metadata <- data.frame(coupon = coupons) %>%
  mutate(substrate = c(rep("Poorman", 4), rep("Ellison", 4), rep("Homestake", 4), rep("Yates", 4)))

#add site names to each coupon 
coupon_sites <- c()
for(coupon in coupons){
  coupon_sites <- append(coupon_sites, paste0(coupon, "_site1"))
  coupon_sites <- append(coupon_sites, paste0(coupon, "_site2"))
}

#function to summarize element compositions from coupons
element_calculator <- function(coupon){
  files <- list.files(directories, full.names = T, pattern = coupon, ignore.case = T)
  if(length(files) > 0){
    file_list = lapply(files, read_files)
  data <- reduce(file_list, full_join) 
  n_pixels = nrow(data)
  data %>%
    dplyr::select(-x,-y) %>%
    summarise_each(funs(sum)) %>%
    mutate(coupon = coupon,
           n_pixels = n_pixels) 
  }
}

coupon_list <- lapply(coupon_sites, element_calculator)

chemistry_data <- reduce(coupon_list, bind_rows) %>%
  dplyr::select(-n_pixels) %>%
  pivot_longer(cols = (-coupon), names_to = "element", values_to = "abundance") %>%
  separate(coupon, c("coupon", "scan_id")) %>%
  full_join(metadata) %>%
  group_by(substrate, element) %>%
  filter(!element == "Os") %>%
  summarise(abundance = sum(na.omit(abundance))) %>%
  filter(!is.na(abundance)) %>%
  group_by(substrate) %>%
  mutate(total = sum(abundance),
         rel_abundance = (abundance/total)*100)

#create color palette
element_colors <- c("#383636", "#f9d2d2", "#d76384", "#d19fc9", "#993721", "#f8c045", "#fbebac", "#dfe452", 
                    "#ead934", "#9ccb3c", "#4da146", "#5cbf91", "#5cbee0", "#5d6ab1", "#9782bc", "#b476b2",
                    "#bab9b9", "#848589", "#8eb397", "#4c7070")
names(element_colors) <- c("Al", "C", "Ca", "Eu", "Fe", "Ge", "Hg", "Ir",
                           "S", "Mg", "Mn", "Mo", "Na", "O", "P", "K", 
                           "Si", "Ti", "Y", "Zr")

rock_composition_plot <- chemistry_data %>%
  mutate(element = factor(element, levels = names(element_colors))) %>%
  ggplot(aes(substrate, abundance, fill = element)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = element_colors)

#plotly::ggplotly(plot)

#write_delim(chemistry_data, "../data/bulkElementComposition.txt", delim = "\t")
