library(tidyverse)
library(lubridate)
library(magrittr)
library(RColorBrewer)

dataset_date <- ymd("2021-04-08") #TODO derive from directory name

database <- str_c(dataset_date, "/database.rds") %>% read_rds
consortium_uk <- str_c(dataset_date, "/consortium_uk.rds") %>% read_rds
mutations_uk <- str_c(dataset_date, "/mutations_uk.rds") %>% read_rds # TODO drop unused columns
mutation_reference_counts <- str_c(dataset_date, "/mutation_reference_counts.rds") %>% read_rds # precomputed mutation counts 
antigenic_mutations_lineages <- str_c(dataset_date, "/antigenic_mutations_lineages.rds") %>% read_rds # antigenic mutation counts by lineage

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
  c(
    "B.1.1.7" = "UK. L18F, Δ69-70, Δ144, N501Y, A570D, P681H, T716I, S982A and D1118H.", 
    "B.1.351" = "South Africa. L18F, D80A, D215G, Δ242-244, R246I, K417N, E484K, N501Y and A701V.", 
    "P.1" = " Japan ex Brazil. L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, H655Y and T1027I",
    "A.23.1" = "UK. R102I, F157L, V367F, E484K, Q613H and P681R.",
    "B.1.525" = "UK ex West Africa. Q52R, A67V, Δ69-70, Δ144, E484K, Q677H and F888L.",
    "B.1.1.318" = "UK ex West Africa. T95I, Δ144, E484K, P681H and D796H.", 
    "B.1.526" = "New York, USA. L5F, T95I, D253G, E484K or S477N and A701V",
    "A.27" = "Mayotte. L18F, L452R, N501Y, A653V, H655Y, Q677H, D796Y and G1219V.",
    "B.1.1.28" = "The Philippines. Δ141-143, Δ243-244, E484K, N501Y, P681H, E1092K, H1101Y and V1176F.",
    
    "B.1.429" = "California, USA. D614G, G1251V, L452R, P26S, S13I, S1252C and W152C.",
    # "B.1.324.1" = "UK associated variant. E484K, S494P, N501Y, D614G, P681H and E1111K in the Spike. ",
    "P.2" = "Brazil. E484K and V1176F.",
    "P.3" = "The Philippines. 141-143del, E484K, N501Y, P681H, E1092K, H1101Y, V1176F and in some cases 243-244del.",
    "B.1.617" = "India. G142D, E154K, L452R, E484K, P681R and Q1071H."
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