pacman::p_load(tidyverse, readxl, lubridate, plotly)


DNA_yields <- read_excel("data/DNA_concentrations.xlsx", sheet = "samples")

DNA_yields_plot <- DNA_yields %>%
  mutate(Date = paste0(month(mdy(Date), label=T), '.', year(mdy(Date)))) %>%
  filter(!Date == "Apr.2018") %>%
  ggplot(aes(rock, `concentration (ng/mL)`, group = modification, color = modification, label = `sample id`)) +
  geom_line() +
  coord_flip() +
  facet_wrap(~site)

ggplotly(DNA_yields_plot)


