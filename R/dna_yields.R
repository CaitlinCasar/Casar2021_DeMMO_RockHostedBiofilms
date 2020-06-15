pacman::p_load(tidyverse, readxl, lubridate, plotly)


DNA_yields <- read_excel("../../../data/DNA_concentrations.xlsx", sheet = "samples")

DNA_yields_plot <- DNA_yields %>%
  mutate(Date = paste0(month(mdy(Date), label=T), '.', year(mdy(Date)))) %>%
  filter(!Date == "Apr.2018" & !rock == "braille" & !is.na(site)) %>%
  group_by(site, rock, modification) %>%
  summarise(`concentration (ng/mL)` = max(`concentration (ng/mL)`)) %>%
  ggplot(aes(reorder(rock, `concentration (ng/mL)`), `concentration (ng/mL)`, group = modification, color = modification)) +
  geom_line() +
  coord_flip() +
  theme(axis.title.y = element_blank()) +
  facet_wrap(~site)

ggplotly(DNA_yields_plot)


