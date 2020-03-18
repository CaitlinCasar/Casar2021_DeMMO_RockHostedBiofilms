#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate)

#load metadata
metadata <- read_csv("../data/metadata.csv") 

# load collated alpha div data
data_path <- "../data/collated_alpha"   # path to the data
files <- dir(data_path, pattern = "*.txt") # get file names
files <- paste(data_path, '/', files, sep="")

alpha_div = tibble::tibble(File = files) %>%
  tidyr::extract(File, "method", "(?<=/data/collated_alpha/)(.*)(?=[.]txt)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_delim, delim = '\t')) %>%
  tidyr::unnest(Data) %>%
  select(-File, -X1, -iteration) %>%
  gather(sample_id, value, `9.DeMMO3.steri.11June2019`:`17.DeMMO1.tube8.12Dec2019`) %>%
  mutate(value = as.numeric(value)) %>%
  group_by(method, `sequences per sample`, sample_id) %>%
  summarise(value = mean(value)) %>%
  left_join(metadata)

#plot rarefaction curves 
alpha_plot <- alpha_div %>%
  ggplot(aes(`sequences per sample`, value, group=sample_id, color=substrate)) +
  geom_line() + 
  scale_colour_discrete(guide = 'none') +
  theme(legend.position="bottom") +
  theme_grey() +
  facet_wrap( ~ method, ncol=2, scales = "free_y") +
  guides(colour = guide_legend(title = "Site")) 

rarefied_alpha_plot <- alpha_div %>%
  filter(`sequences per sample` == 46500 & method == "observed_otus" & !is.na(site)) %>%
  ggplot(aes(substrate, value, color=substrate, group = substrate)) + 
  geom_line(size = 2) +
  stat_summary(fun.y=mean, geom="point", size=1, color="white") +
  coord_flip() +
  theme(legend.position = "none", 
        axis.title = element_blank()) +
  facet_wrap(~site)
