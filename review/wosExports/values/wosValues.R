# paste together wos export files - values

library(readxl)
library(tidyverse)

read_xlsx("savedrecs1.xlsx", col_types = "text") %>%
  bind_rows(read_xlsx("savedrecs2.xlsx", col_types = "text")) %>%
  bind_rows(read_xlsx("savedrecs3.xlsx", col_types = "text")) %>%
  bind_rows(read_xlsx("savedrecs4.xlsx", col_types = "text")) %>%
  bind_rows(read_xlsx("savedrecs5.xlsx", col_types = "text")) %>%
  bind_rows(read_xlsx("savedrecs6.xlsx", col_types = "text")) %>%
  bind_rows(read_xlsx("savedrecs7.xlsx", col_types = "text")) %>%
  bind_rows(read_xlsx("savedrecs8.xlsx", col_types = "text")) %>%
  # keep only journal articles
  filter(PT == "J") %>%
  # calculate citations per year published
  mutate(citesPerYear = parse_number(Z9) / (2021 - parse_number(PY))) %>%
  # arrange by citations per year published
  arrange(desc(citesPerYear)) %>%
  # keep relevant columns
  select(PY, AU, TI, SO, VL, IS, DI, Z9, citesPerYear) %>%
  # write
  write_excel_csv(path = "wosValues.csv")
