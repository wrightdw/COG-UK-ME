library(tidyverse)
library(lubridate)
library(magrittr)
library(RColorBrewer)

source("helpers.R")

dataset_date <- ymd("2021-06-08") #TODO derive from directory name

database <- str_c(dataset_date, "/database.rds") %>% read_rds # spike database
consortium_uk <- str_c(dataset_date, "/consortium_uk.rds") %>% read_rds
mutations_uk <- str_c(dataset_date, "/mutations_uk.rds") %>% read_rds # TODO drop unused columns
mutation_reference_counts <- str_c(dataset_date, "/mutation_reference_counts.rds") %>% read_rds # precomputed mutation counts 
database_tcell_predictions <- str_c(dataset_date, "/database_tcell_predictions.rds") %>% read_rds # spike T cell info and predictions
deletions <- str_c(dataset_date, "/deletions.rds") %>% read_rds # deletions (genomic coordinates)
vui_voc <- read_rds("vui_voc.rds") # VUI/VOC defining mutations in spike protein
wt <- read_rds(str_c(dataset_date, "/wt.rds")) # spike protein wild type amino acid counts
lineages_weeks_uk <- read_rds(str_c(dataset_date, "/lineages_weeks_uk.rds")) # VUI/VOC/other counts

mutations_s_uk <- 
  mutations_uk %>% 
  filter(gene == "S") %>% 
  select(-gene) %>% 
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
    "B.1.1.7" = "UK. L18F, Δ69-70, Δ144, N501Y, A570D, P681H, T716I, S982A and D1118H. WHO label: <strong>Alpha</strong>.", 
    "B.1.351" = 
      "South Africa. D80A, D215G, Δ242-244, K417N, E484K, N501Y and A701V. WHO label: <strong>Beta</strong>.", 
    "P.1" = " Japan ex Brazil. L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, H655Y and T1027I. WHO label: <strong>Gamma</strong>.",
    "A.23.1" = "UK. R102I, F157L, V367F, E484K, Q613H and P681R.",
    "B.1.525" = "UK ex West Africa. Q52R, A67V, Δ69-70, Δ144, E484K, Q677H and F888L. WHO label: <strong>Eta</strong>.",
    "B.1.1.318" = "UK ex West Africa. T95I, Δ144, E484K, P681H and D796H.", 
    "B.1.526" = "New York, USA. L5F, T95I, D253G, E484K or S477N and A701V. WHO label: <strong>Iota</strong>.",
    "A.27" = "Mayotte. L18F, L452R, N501Y, A653V, H655Y, Q677H, D796Y and G1219V.",
    "B.1.1.28" = "The Philippines. Δ141-143, Δ243-244, E484K, N501Y, P681H, E1092K, H1101Y and V1176F.",
    "B.1.429" = "California, USA. D614G, G1251V, L452R, P26S, S13I, S1252C and W152C. WHO label: <strong>Epsilon</strong>.",
    "B.1.324.1" = "UK associated variant. E484K, S494P, N501Y, D614G, P681H and E1111K in the Spike. ",
    "P.2" = "Brazil. E484K and V1176F. WHO label: <strong>Zeta</strong>.",
    "P.3" = "The Philippines. Δ141-143, E484K, N501Y, P681H, E1092K, H1101Y, V1176F and in some cases Δ243-244. WHO label: <strong>Theta</strong>.",
    # "B.1.617" = "India. G142D, E154K, L452R, E484Q, P681R and Q1071H.",
    "B.1.617.1" = "India. E154K, L452R, E484Q and P681R. WHO label: <strong>Kappa</strong>.",
    "B.1.617.2" = "India. T19R, Δ156-157, R158G, L452R, T478K, D614G, P681R and D950N. WHO label: <strong>Delta</strong>.",
    "B.1.617.3" = "India. T19R, Δ156-158, L452R, E484Q, D614G, P681R and D950N.",
    "AV.1" = "UK, Greece and Chad. D80G, T95I, G142D, Δ144, N439K, E484K, D614G, P681H, I1130V and D1139H.",
    "C.36.3" = "Egypt. S12F, Δ69-70, W152R, R346S, L452R, Q677H and A899S."
    ) %>% 
  enframe("lineage", "reason")

lineages_t2 %<>% c(lineages_t3$lineage) %>% unique # combine table 2 and 3 lineages for counting

# TODO precompute and include lineage/mutation combinations
n_uk_lineages_all <-
  left_join(
    sum_key_mutations_by_lineage_uk(lineages_t2),
    sum_key_mutations_by_lineage_uk(lineages_t2, date_from = sample_date_28) %>%
      rename(n_sequences_28 = n_sequences)
  ) %>% 
  pivot_wider(names_from = adm1, values_from = c(n_sequences, n_sequences_28)) %>% 
  mutate(across(everything(), ~replace_na(.x, 0L)))

# remove VOCs/VUIs with zero counts
vui_voc %<>% 
  semi_join(n_uk_lineages_all %>% filter(variant == "sequences" & n_sequences_UK > 0), by = "lineage") %>% 
  mutate(lineage = fct_drop(lineage))