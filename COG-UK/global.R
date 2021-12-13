library(tidyverse)
library(lubridate)
library(magrittr)

# Work out choice of dataset based on hour of day or most recent dataset
get_dataset_date <- function(rollover = 7){
  # list dataset dates removing today
  dirs <- 
    list.dirs() %>% 
    ymd(quiet = TRUE) %>% 
    .[. != today()]
  
  # if before 7am, remove yesterday as well
  if(hour(now()) < rollover){
    dirs %<>% .[. != (today() - days(1))]
  }
  
  # get next most recent dataset
  max(dirs, na.rm = TRUE)
}

# Dataset directories are named according to date e.g. "2021-11-24".
# if after midnight, use yesterday
# else before midnight, use 2 days ago
dataset_date <- get_dataset_date(0)

# Alternatively, set date here instead to switch to specific dataset.
# dataset_date <- as.Date("2021-12-11")

database_genome <- str_c(dataset_date, "/database_genome.rds") %>% read_rds # mutation database
consortium_uk <- str_c(dataset_date, "/consortium_uk.rds") %>% read_rds
mutations_uk <- str_c(dataset_date, "/mutations_uk.rds") %>% read_rds # TODO drop unused columns
mutation_reference_counts <- str_c(dataset_date, "/mutation_reference_counts.rds") %>% read_rds # precomputed mutation counts 
database_tcell_predictions <- str_c(dataset_date, "/database_tcell_predictions.rds") %>% read_rds # spike T cell info and predictions
deletions <- str_c(dataset_date, "/deletions.rds") %>% read_rds # deletions (genomic coordinates)
vui_voc <- read_rds("vui_voc.rds") # VUI/VOC defining mutations in spike protein
wt <- read_rds(str_c(dataset_date, "/wt.rds")) # spike protein wild type amino acid counts
lineages_weeks_uk_all <- read_rds(str_c(dataset_date, "/lineages_weeks_uk_all.rds")) # lineage counts by epiweek
lineages_days_uk_all <- read_rds(str_c(dataset_date, "/lineages_days_uk_all.rds")) # lineage counts by sample date
therapeutics <- read_rds(str_c(dataset_date, "/therapeutics.rds")) # antiviral drug resistance mutations
insertions <- str_c(dataset_date, "/insertions.rds") %>% read_rds # deletions (genomic coordinates)
spike_tab <- read_rds(str_c(dataset_date, "/spike_table.rds"))

source("helpers.R")

database <- 
  database_genome %>% 
  filter(gene == "S") %>% 
  select(-gene) %>% 
  mutate(across(where(is.factor), fct_drop))

mutations_s_uk <- 
  mutations_uk %>% 
  filter(gene == "S") %>% 
  select(-gene) %>% 
  mutate(across(c(variant, position), fct_drop)) # drop non-spike mutations from factor levels

sample_date_28 <- max(consortium_uk$sample_date) - days(27) # calculate 28 day period up to and including latest sample date

total_sequences <- n_distinct(consortium_uk$sequence_name)
total_sequences_28 <- 
  consortium_uk %>% 
  filter(sample_date >= sample_date_28) %$% 
  n_distinct(sequence_name)

lineages_t3 <- 
  c(
    "B.1.1.7" = "UK. L18F, Δ69-70, Δ144, N501Y, A570D, P681H, T716I, S982A and D1118H. WHO label: <strong>Alpha</strong>.", 
    "B.1.351" = 
      "South Africa. D80A, D215G, Δ242-244, K417N, E484K, N501Y and A701V. WHO label: <strong>Beta</strong>.", 
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
    
    "B.1.617.1" = "India. E154K, L452R, E484Q and P681R. WHO label: <strong>Kappa</strong>.",
    "B.1.617.2" = "India. T19R, G142D, Δ156-157, R158G, L452R, T478K, D614G, P681R and D950N. WHO label: <strong>Delta</strong>.",
    "B.1.617.3" = "India. T19R, Δ156-158, L452R, E484Q, D614G, P681R and D950N.",
    "AV.1" = "UK, Greece and Chad. D80G, T95I, G142D, Δ144, N439K, E484K, D614G, P681H, I1130V and D1139H.",
    "C.36.3" = "Egypt. S12F, Δ69-70, W152R, R346S, L452R, Q677H and A899S.",
    
    "C.37" = "South America. Δ246-252, G75V, T76I, L452Q, F490S, D614G, and T859N. WHO label: <strong>Lambda</strong>.",
    "B.1.621" = "Colombia. T95I, R346K, E484K, N501Y and P681H. WHO label: <strong>Mu</strong>.",
    
    "AY.4" = "Alias of B.1.617.2.4, UK. WHO label: <strong>Delta</strong>.",
    "AY.4.2" = "Subineage of AY.4 with the addition of Y145H and A222V. Alias of  B.1.617.2.4.2. WHO label: <strong>Delta</strong>.",
    "AY.4.2.1" = "Sublineage of AY.4.2 with the addition of V36F. Alias of  B.1.617.2.4.2.1. WHO label: <strong>Delta</strong>.",
    
    "P.1" = " Japan ex Brazil. L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, H655Y and T1027I. WHO label: <strong>Gamma</strong>.",
    "P.1.8" = "Brazil. S: Gamma + T470N, P681R, C1235F;
NSP3: Gamma + I441V; NSP4: A446V; ORF3a: S216L; ORF8: G8*STOP; N: TRS insertion. WHO label: <strong>Gamma</strong>.",

    "BA.1" = "Southern Africa. Full Spike profile: A67V, Δ69-70, T95I, G142D/Δ143-145, Δ211/L212I, ins214EPE, G339D, S371L, S373P, S375F, K417N, N440K, G446S, S477N, T478K, E484A, Q493R, G496S, Q498R, N501Y, Y505H, T547K, D614G, H655Y, N679K, P681H, N764K, D796Y, N856K, Q954H, N969K, L981F. WHO label: <strong>Omicron</strong>.",
  "BA.2" = "Southern Africa. Full Spike profile: T19I, Δ24-26, A27S, G142D, V213G, G339D, S371F, S373P, S375F, T376A, D405N, R408S, K417N, N440K, S477N, T478K, E484A, Q493R, Q498R, N501Y, Y505H, D614G, H655Y, N679K, P681H, N764K, D796Y, Q954H, N969K. WHO label: <strong>Omicron</strong>."
    ) %>% 
  enframe("lineage", "reason")

lineages_t2 <- c(vui_voc %>% levels, lineages_t3$lineage) %>% unique # lineages for counting

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

# VUI/VOC lineages with AY.x removed
vui_voc_lineages <- 
  vui_voc %>% 
  filter(!str_starts(lineage, fixed("AY."))) %>% 
  mutate(across(c(lineage, lineage_display), fct_drop)) 

vui_voc_lineages <- 
  levels(vui_voc_lineages$lineage) %>% 
  setNames(levels(vui_voc_lineages$lineage_display))

geo_all <- str_c(dataset_date, "/geo_all.rds") %>% read_rds # geographical NUTS1 counts
mapdata <- read_rds("mapdata.rds") # UK map NUTS1 topology as dataframe
# antigenic_mutations_lineages_all <- read_rds(str_c(dataset_date, "/antigenic_mutations_lineages_all.rds"))