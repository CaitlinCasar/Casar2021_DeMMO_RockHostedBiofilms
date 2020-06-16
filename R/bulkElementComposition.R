pacman::p_load(tidyverse, raster)

directories <- list.dirs("/Volumes/Elements/DeMMO/DeMMO_Publications/DeMMO_NativeRock/data/Dec2019/DeMMO2019_lowMagXEDS", full.names = T , recursive =T)
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

#
coupon_list <- lapply(coupon_sites, element_calculator)

data <- reduce(coupon_list, bind_rows) %>%
  dplyr::select(-n_pixels) %>%
  pivot_longer(cols = (-coupon), names_to = "element", values_to = "abundance") %>%
  separate(coupon, c("coupon", "scan_id")) %>%
  left_join(metadata) %>%
  group_by(substrate, element) %>%
  filter(!element == "Os") %>%
  summarise(abundance = sum(na.omit(abundance))) %>%
  filter(!is.na(abundance)) %>%
  group_by(substrate) %>%
  mutate(total = sum(abundance),
         rel_abundance = (abundance/total)*100)

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

plot <- data %>%
  ggplot(aes(substrate, abundance, fill = element)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = element_colors)

plotly::ggplotly(plot)
