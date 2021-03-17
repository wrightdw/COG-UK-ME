library(tidyverse)
library(lubridate)
library(magrittr)
library(RColorBrewer)

dataset_date <- ymd("2021-03-15") #TODO derive from directory name

database <- str_c(dataset_date, "/database.rds") %>% read_rds
consortium_uk <- str_c(dataset_date, "/consortium_uk.rds") %>% read_rds
mutations_uk <- str_c(dataset_date, "/mutations_uk.rds") %>% read_rds # TODO drop unused columns
mutation_reference_counts <- str_c(dataset_date, "/mutation_reference_counts.rds") %>% read_rds # precomputed mutation counts 

mutations_s_uk <- 
  mutations_uk %>% 
  filter(gene == "S") %>% 
  mutate(across(c(variant, position), fct_drop))

sample_date_28 <- max(consortium_uk$sample_date) - days(27) # calculate 28 day period up to and including latest sample date

total_sequences <- n_distinct(consortium_uk$sequence_name)
total_sequences_28 <- 
  consortium_uk %>% 
  filter(sample_date >= sample_date_28) %$% 
  n_distinct(sequence_name)

lineages_t2 <- c("B.1", "B.1.177", "B.1.141", "B.1.258", "B.1.1", 
                 "B.1.1.7", "B.1.1.70", "B.1.351", "B.1.1.298", "P.2", "P.1",
                 "B.1.222", "A.23.1", "B.1.1.119", "B.1.177.4", "B.1.525")

lineages_t3 <- 
  c("B.1.1.7" = "UK associated variant. Has 17 mutations (14 replacements and 3 deletions) including: T1001I, A1708D, I2230T, SGF 3675-3677 del In the ORF1ab; 69-70 del, Y144 del, N501Y, A570D, P681H, T716I, S982A and D1118H in the Spike; Q27stop, R52I and Y73C in ORF8; D3L and S235F in the N. Noteworthily, N501Y enhances ACE2 binding affinity, and P681H occurs at the furin cleavage site, known for biological significance in membrane fusion.", 
    "B.1.351" = "Variant associated with South Africa. Has eight mutations in the Spike: D80A, D215G, E484K, N501Y, A701V, L18F, R246I and K417N. Three of these in the RBM, K417N, E484K and N501Y. K417N and E484K have been shown to escape some mAbs.", 
    "P.1" = "Variant associated with Brazil. Has 10 mutations in the Spike including L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y,H655Y and T1027I. Noteworthy  E484K, N501Y and K417T have biological significance.",
    "A.23.1" = "International variant with mutations of biological significance F157L, V367F, Q613H and P681R. Q613H is predicted to be functionally equivalent to the D614G mutation that arose early in 2020.",
    "B.1.525" = "International variant with mutations of biological significance E484K, Q677H, F888L and a similar suite of deletions to B.1.1.7.",
    
    "B.1.429" = "Variant associated with California, USA.",
    "B.1.526" = "Variant associated with New York, USA.",
    "A.27" = "Variant associated with France.",
    # "B.1.324.1" = "UK associated variant.",
    "P.2" = "Variant associated with Brazil.",
    "B.1.1.318" = "Variant associated with England.", 
    "B.1.1.28" = "Variant associated with the Philippines."
    ) %>% 
  enframe("lineage", "reason")

lineages_t2 %<>% c(lineages_t3$lineage) %>% unique # combine table 2 and 3 lineages for counting

# Construct a regular expression to match the sublineages of a lineage
sublineage_regex <- function(lineage){
  str_replace_all(lineage, "\\.", "\\\\\\.") %>% 
    paste0("^", ., "\\.")
}

sum_key_mutations_uk <- function(..., date_from = NULL){
  if(!is_null(date_from)){
    date_from %<>% ymd()
  }

  if(is.Date(date_from)){
    consortium_uk %<>% filter(sample_date >= date_from)
  }
  
  consortium_uk %>%
    group_by(...) %>%
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
            `N501Y + E484K` = sum(n501y == "Y" & e484k == "K"),
            .groups = "keep"
  )
}

sum_key_mutations_by_lineage_uk <- function(lineages = NULL, date_from = NULL){
  if(is_character(lineages)){
    n_nations_lineages <- sum_key_mutations_uk(lineage, adm1, date_from = date_from) # grouped by lineage, adm1
    
    n_uk_lineages <- 
      n_nations_lineages %>% 
      group_by(lineage) %>% 
      summarise(across(sequences:`N501Y + E484K`, sum)) %>% 
      mutate(adm1 = "UK", .after = lineage)
    
    n_uk_lineages_uk_nations <- bind_rows(n_nations_lineages, n_uk_lineages) # grouped by lineage, adm1
  
    lapply(lineages, function(x){
      n_uk_lineages_uk_nations %>% 
        filter(lineage == x | str_detect(lineage, sublineage_regex(x))) %>% 
        group_by(adm1) %>% 
        select(-lineage) %>% 
        summarise_all(funs(sum)) %>% #TODO replace deprecated funs
        mutate(lineage = x, .before = 1)
    }) %>% 
      bind_rows() %>% 
      gather(key = "variant", value = "n_sequences", sequences:`N501Y + E484K`)
  }
}

lineage_plus_variant <- function(lineage, variant){
  mutations_s_uk_lv <- 
    mutations_s_uk %>%
      filter(variant == !!variant) %>%
      filter(lineage == !!lineage | str_detect(lineage, sublineage_regex(!!lineage))) 
  
  mutations_s_uk_lv_28 <- 
    mutations_s_uk_lv %>%
    filter(sample_date >= sample_date_28)
  
    left_join(
      mutations_s_uk_lv %>%
        group_by(adm1) %>% 
        summarise(n_sequences = n_distinct(sequence_name)) %>% 
        bind_rows(summarise(., n_sequences = sum(n_sequences)) %>% 
                    mutate(adm1 = "UK")),
      
        mutations_s_uk_lv_28 %>%
          group_by(adm1) %>% 
          summarise(n_sequences_28 = n_distinct(sequence_name)) %>% 
          bind_rows(summarise(., n_sequences_28 = sum(n_sequences_28)) %>% 
                  mutate(adm1 = "UK"))
    ) %>%
      mutate(adm1 = recode(adm1, 
                        `UK-ENG` = "England",
                        `UK-NIR` = "Northern_Ireland",
                        `UK-SCT` = "Scotland",
                        `UK-WLS` = "Wales")) %>% 
      pivot_wider(names_from = adm1, values_from = c(n_sequences, n_sequences_28)) %>%
      mutate(lineage = !!lineage, variant = !!variant, .before = 1)
}

# TODO precompute and include lineage/variant combinations
n_uk_lineages_all <-
  left_join(
    sum_key_mutations_by_lineage_uk(lineages_t2),
    sum_key_mutations_by_lineage_uk(lineages_t2, date_from = sample_date_28) %>%
      rename(n_sequences_28 = n_sequences)
  ) %>% 
  pivot_wider(names_from = adm1, values_from = c(n_sequences, n_sequences_28)) %>% 
  mutate(across(everything(), ~replace_na(.x, 0L)))