#load dependencies 
pacman::p_load(tidyverse, readr, plyr, plotly, lubridate)

#load otu table and metadata
metadata <- read_csv("../../data/metadata.csv") 
otu_table <- read_delim("../../data/DeMMO_NativeRock_noChimera_otuTable_withTaxa_46822.txt", delim = '\t', comment = "# ")


# load summary report data 
directories <- list.dirs("../../data", full.names = T , recursive =T)
directories <- directories[str_detect(directories, "reports")]
files <- list.files(directories, full.names = T)


elements = tibble::tibble(File = files[str_detect(files, "element_summary")]) %>%
  tidyr::extract(File, "sample", "(?<=/reports/)(.*)(?=_element_summary[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File) %>%
  left_join(metadata)

cells = tibble::tibble(File = files[str_detect(files, "cell_summary")]) %>%
  tidyr::extract(File, "sample", "(?<=/reports/)(.*)(?=_cell_summary[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File) %>%
  left_join(metadata)

models = tibble::tibble(File = files[str_detect(files, "models.csv")]) %>%
  tidyr::extract(File, "sample", "(?<=/reports/)(.*)(?=_cell_distribution_models[.]csv)", remove = FALSE) %>%
  mutate(Data = lapply(File, readr::read_csv)) %>%
  tidyr::unnest(Data) %>%
  select(-File) %>%
  left_join(metadata)
