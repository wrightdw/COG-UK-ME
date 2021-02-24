library(tidyverse)
library(lubridate)
library(magrittr)
library(RColorBrewer)

database <- read_rds("2021-02-16/database.rds")
consortium_uk <- read_rds("2021-02-16/consortium_uk.rds")
mutations_uk <- read_rds("2021-02-16/mutations_uk.rds")
mutation_reference_counts <- read_rds("2021-02-16/mutation_reference_counts.rds") # precomputed mutation counts 

mutations_s_uk <- 
  mutations_uk %>% 
  filter(gene == "S") %>% 
  mutate(across(c(variant, position), fct_drop))

dataset_date <- ymd("2021-02-16") #TODO derive from filename
sample_date_28 <- max(consortium_uk$sample_date) - days(27) # calculate 28 day period up to and including latest sample date

total_sequences <- n_distinct(consortium_uk$sequence_name)
total_sequences_28 <- 
  consortium_uk %>% 
  filter(sample_date >= sample_date_28) %$% 
  n_distinct(sequence_name)

lineages_t2 <- c("B.1", "B.1.177", "B.1.141", "B.1.258", "B.1.1", "B.1.1.7", "B.1.1.70", "B.1.351", "B.1.1.298", 
                 "P.2", "P.1", "B.1.222", "A.23.1", "B.1.1.119", "B.1.177.4", "B.1.525")

lineages_t3 <- 
  c("B.1.1.7" = "UK associated variant. Has 17 mutations (14 replacements and 3 deletions) including: T1001I, A1708D, I2230T, SGF 3675-3677 del In the ORF1ab; 69-70 del, Y144 del, N501Y, A570D, P681H, T716I, S982A and D1118H in the Spike; Q27stop, R52I and Y73C in ORF8; D3L and S235F in the N. Noteworthily, N501Y enhances ACE2 binding affinity, and P681H occurs at the furin cleavage site, known for biological significance in membrane fusion.", 
    "B.1.351" = "Variant associated with South Africa. Has eight mutations in the Spike: D80A, D215G, E484K, N501Y, A701V, L18F, R246I and K417N. Three of these in the RBM, K417N, E484K and N501Y. K417N and E484K have been shown to escape some mAbs.", 
    "P.1" = "Variant associated with Brazil. Has 10 mutations in the Spike including L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y,H655Y and T1027I. Noteworthy  E484K, N501Y and K417T have biological significance.") %>% 
  enframe("lineage", "reason")

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
            E484K = sum(e484k == "K"),
            `∆69-70` = sum(del_21765_6 == "del"),
            `N439K + ∆69-70` = sum(n439k == "K" & del_21765_6 == "del"),
            `N501Y + ∆69-70` = sum(n501y == "Y" & del_21765_6 == "del"),
            `Y453F + ∆69-70` = sum(y453f == "F" & del_21765_6 == "del"),
            `N501Y + E484K` = sum(n501y == "Y" & e484k == "K")
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
      gather(key = "variant", value = "n_sequences", sequences:`N501Y + E484K`)
  }
}

n_uk_lineages_all <- inner_join(sum_key_mutations_by_lineage_uk(lineages_t2), 
                                sum_key_mutations_by_lineage_uk(lineages_t2, date_from = sample_date_28) %>% 
                                  rename(n_sequences_28 = n_sequences))