library(tidyverse)
library(lubridate)
library(magrittr)
library(RColorBrewer)

database <- read_rds("database.rds")
consortium_uk <- read_rds("consortium_uk.rds")
mutations_uk <- read_rds("mutations_uk.rds")

epi_levels <- min(consortium_uk$epi_week):max(consortium_uk$epi_week) %>% as.character

consortium_uk %<>%
  mutate(across(epi_week, as_factor)) %>%
  mutate(epi_week = fct_expand(epi_week, epi_levels) %>% fct_inseq)

mutations_uk %<>%
  mutate(across(c(epi_week, gene, position), as_factor)) %>%
  mutate(epi_week = fct_expand(epi_week, epi_levels) %>% fct_inseq)

mutations_s_uk <- mutations_uk %>% filter(gene == "S")






dataset_date <- ymd("2021-01-13") #TODO derive from filename

sample_date_28 <- max(consortium_uk$sample_date) - days(27) # calculate 28 day period up to and including latest sample date

#TODO Pre-load key mutations

# Construct a regular expression to match the sublineages of a lineage
sublineage_regex <- function(lineage){
  str_replace_all(lineage, "\\.", "\\\\\\.") %>% 
    paste0("^", ., "\\.")
}

sum_key_mutations_uk <- function(date_from = NULL){
  if(!is_null(date_from)){
    date_from %<>% ymd()
  }

  if(is.Date(date_from)){
    consortium_uk %<>% filter(sample_date >= date_from)
  }
  
  consortium_uk %>%
  group_by(lineage) %>%
  summarise(sequences = n(),
            D614G = sum(d614g == "G"),
            A222V = sum(a222v == "V"),
            N439K = sum(n439k == "K"),
            N501Y = sum(n501y == "Y"),
            Y453F = sum(y453f == "F"),
            `∆69-70` = sum(del_21765_6 == "del"),
            `N439K + ∆69-70` = sum(n439k == "K" & del_21765_6 == "del"),
            `N501Y + ∆69-70` = sum(n501y == "Y" & del_21765_6 == "del"),
            `Y453F + ∆69-70` = sum(y453f == "F" & del_21765_6 == "del")
  )
}

sum_key_mutations_by_lineage_uk <- function(lineages = NULL, date_from = NULL){
  if(is_character(lineages)){
    n_uk_lineages <- sum_key_mutations_uk(date_from)
  
    lapply(lineages, function(x){
      n_uk_lineages %>% 
        filter(lineage == x | str_detect(lineage, sublineage_regex(x))) %>% 
        select(-lineage) %>% 
        summarise_all(funs(sum)) %>% #TODO replace deprecated funs
        mutate(lineage = x, .before = 1)
    }) %>% 
      bind_rows() %>% 
      select(-sequences) %>% 
      gather(key = "variant", value = "n_sequences", D614G:`Y453F + ∆69-70`)
  }
}
