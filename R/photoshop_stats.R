pacman::p_load(tidyverse, ggridges, cowplot)


image_sizes <- read_csv("../data/photoshop_data/image_sizes.csv") %>%
  mutate(width = width/22,
         height = height/22,
         transect_area = if_else(str_detect(coupon_id, "cont"), 5*width * height, width * height)) %>%
  select(coupon_id, transect_area)
files <- list.files("../data/photoshop_data", full.names = T, pattern = ".txt")

photoshop_stats = tibble::tibble(File = files) %>%
  tidyr::extract(File, "coupon_id", "(?<=photoshop_data/)(.*)(?=_Dec2019_)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_delim, delim = "\t")) %>%
  tidyr::unnest(Data) %>%
  dplyr::select(-File) %>%
  filter(is.na(Count) | Count == 1) %>%
  select(coupon_id, Area, Perimeter, Circularity, Height, Width) %>%
  filter(!coupon_id == "D3T15exp") %>%
  mutate(type = if_else(str_detect(coupon_id, "fungal"), "fungal", "cell"),
         coupon_id = str_remove(coupon_id, "_fungal"),
         filament = if_else(Circularity <= 0.1, 1, 0),
         cocci = if_else(Circularity > 0.7, 1, 0),
         rod = if_else(Circularity <= 0.7 & Circularity > 0.1, 1, 0)) %>%
  left_join(metadata %>% select(coupon_id, site, substrate)) 

summarised_photoshop_stats = photoshop_stats %>%
  group_by(coupon_id, type) %>%
  summarise(mean_area = mean(Area), 
            cocci = sum(cocci), rod = sum(rod), filament = sum(filament), total_observations = n(),
            total_area = sum(Area)) %>%
  left_join(metadata %>% select(coupon_id, site, substrate)) %>%
  left_join(image_sizes) %>%
  mutate(coupon_id = str_remove(coupon_id, "exp|rep|cont"),
         cocci = cocci/total_observations*100, rod = rod/total_observations*100, filament = filament/total_observations*100,
         coverage_area = total_area/transect_area*100)


cell_size_plot <- photoshop_stats %>%
  select(coupon_id, site, substrate, Area, Circularity) %>%
  mutate(coupon_id = str_remove(coupon_id, "exp|rep|cont")) %>%
  mutate(id = paste0(substrate, "_", coupon_id)) %>%
  pivot_longer(Area:Circularity, names_to = "observation", values_to = "value") %>%
  mutate(observation = if_else(observation == "Area", "Area (sq. micron)", observation)) %>%
  ggplot(aes(value, id, fill = substrate)) + 
  geom_density_ridges() +
  theme_ridges() + 
  xlim(0, 1) +
  theme(legend.position = "none") +
  facet_grid(vars(site),vars(observation), scales = "free") +
  geom_vline(data = data.frame(xint=c(0.1, 0.7),observation="Circularity"), aes(xintercept = xint), linetype = "dotted") +
  theme(axis.title = element_blank(),
        strip.background.y = element_blank(), strip.text.y = element_blank())

cell_stats_table <- summarised_photoshop_stats %>%
  select(type, total_observations, coverage_area, substrate, site) %>%
  mutate(coverage_area = round(coverage_area, 1)) %>%
  pivot_wider(id_cols = c(coupon_id, substrate, site), names_from = "type", values_from = c("total_observations", "coverage_area"), values_fill = list(total_observations = 0, coverage_area = 0)) %>%
  ungroup() %>%
  mutate(coupon_id = str_remove(coupon_id, "exp|rep|cont"),
         id = paste0(substrate, "_", coupon_id)) %>%
  select(id, site, total_observations_cell, coverage_area_cell, coverage_area_fungal) %>%
  pivot_longer(total_observations_cell:coverage_area_fungal, names_to = "observation", values_to = "value") %>%
  mutate(observation= if_else(observation == "coverage_area_cell", "cell coverage %", 
                              if_else(observation == "coverage_area_fungal", "fungal coverage %", 
                                      if_else(observation == "total_observations_cell", "# cells", observation)))) %>%
  ggplot(aes(x=observation,y=id, label = value)) + 
  geom_tile(fill = "white") + geom_text(colour = "black") +
  facet_grid(vars(site),vars(observation), scales = "free") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        strip.background.y = element_blank(), strip.text.y = element_blank(),
        axis.ticks=element_blank())


plot_grid(cell_size_plot, cell_stats_table, nrow=1, align="h", rel_widths = c(3,1.5))

