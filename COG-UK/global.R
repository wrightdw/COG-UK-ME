library(tidyverse)
library(lubridate)
library(magrittr)
library(waiter)
library(stringi)

# Work out choice of dataset based on hour of day or most recent dataset
get_dataset_date <- function(rollover = 7){
  # list dataset dates removing today
  dirs <- 
    list.dirs() %>% 
    ymd(quiet = TRUE) %>% 
    .[. != today()]
  
  # if before rollover time, remove yesterday as well
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
# dataset_date <- as.Date("2022-06-14")

consortium_uk <- str_c(dataset_date, "/consortium_uk.rds") %>% read_rds
deletions <- str_c(dataset_date, "/deletions.rds") %>% read_rds # deletions (genomic coordinates) # TODO remove dependency
mutation_reference_counts <- str_c(dataset_date, "/mutation_reference_counts.rds") %>% read_rds # precomputed mutation counts
database_deletions <- read_rds(str_c(dataset_date, "/database_deletions.rds"))
mutations_indels_uk_28 <- read_rds(str_c(dataset_date, "/mutations_indels_uk_28.rds"))

# # not required for dashboard app
# deletions_mapping <- read_rds(str_c(dataset_date, "/deletions_mapping.rds"))
# insertions <- str_c(dataset_date, "/insertions.rds") %>% read_rds # insertions (genomic coordinates)
# insertions_mapping <- read_rds(str_c(dataset_date, "/insertions_mapping.rds"))
# mutations_uk <- str_c(dataset_date, "/mutations_uk.rds") %>% read_rds

database_genome <- str_c(dataset_date, "/database_genome.rds") %>% read_rds # mutation database
mutations_s_uk <- str_c(dataset_date, "/mutations_s_escape_uk.rds") %>% read_rds # S gene escape mutations
database_tcell_predictions <- str_c(dataset_date, "/database_tcell_predictions.rds") %>% read_rds # spike T cell info and predictions
vui_voc <- read_rds("vui_voc.rds") # VUI/VOC defining mutations in spike protein
wt <- read_rds(str_c(dataset_date, "/wt.rds")) # spike protein wild type amino acid counts
lineages_weeks_uk_all <- read_rds(str_c(dataset_date, "/lineages_weeks_uk_all.rds")) # lineage counts by epiweek
lineages_days_uk_all <- read_rds(str_c(dataset_date, "/lineages_days_uk_all.rds")) # lineage counts by sample date
therapeutics <- read_rds(str_c(dataset_date, "/therapeutics.rds")) # antiviral drug resistance mutations
spike_tab <- read_rds(str_c(dataset_date, "/spike_table.rds"))
database_insertions <- read_rds(str_c(dataset_date, "/database_insertions.rds"))

source("helpers.R")

dms_antigenic <- read_rds("dms_integrated_data@13.rds")
dms_antigenic %<>% select(-c("real",	"convergent",	"bloom_bind_avg",	"bloom_expr_avg",	"color")) %>%
  filter(!is.na(semantic_score)) %>%
  filter(!is.na(grammaticality)) %>%
  filter(!is.na(evolutionary_index)) %>%
  filter(!is.na(entropy))

ref_wuhan = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
ref_wuhan <- unlist(strsplit(ref_wuhan, split = "", fixed=T)) # Split into individual amino acid characters

database <- 
  database_genome %>% 
  filter(gene == "S") %>% 
  select(-gene) %>% 
  mutate(across(where(is.factor), fct_drop))

database_deletions %<>%
  rename(
    mutation = del_id,
    `numSeqs UK` = UK,
    `numSeqs UK 28 days` = UK_28,
    `numSeqs Eng 28 days` = England_28,
    `numSeqs Scotland 28 days` = Scotland_28,
    `numSeqs Wales 28 days` = Wales_28,
    `numSeqs NI 28 days` = Northern_Ireland_28)

database_insertions %<>%
  rename(
    mutation = insertion_id,
    `numSeqs UK` = UK,
    `numSeqs UK 28 days` = UK_28,
    `numSeqs Eng 28 days` = England_28,
    `numSeqs Scotland 28 days` = Scotland_28,
    `numSeqs Wales 28 days` = Wales_28,
    `numSeqs NI 28 days` = Northern_Ireland_28)

database_genome %<>%
  bind_rows(database_deletions) %>%
  bind_rows(database_insertions)

sample_date_28 <- max(consortium_uk$sample_date) - days(27) # calculate 28 day period up to and including latest sample date

total_sequences <- n_distinct(consortium_uk$sequence_name)
total_sequences_28 <- 
  consortium_uk %>% 
  filter(sample_date >= sample_date_28) %$% 
  n_distinct(sequence_name)

# Lineage names and descriptions for lineage table
# Counted without sublineages TODO count sublineages for everything
# Omicron BA.x dealt with separately to include sublineages
lineages_t3 <- 
  c(
    "B.1.1.7" = "UK. L18F, Δ69-70, Δ144, N501Y, A570D, P681H, T716I, S982A and D1118H. WHO label: <strong>Alpha</strong>.", 
    "B.1.351" = 
      "South Africa. D80A, D215G, Δ242-244, K417N, E484K, N501Y and A701V. WHO label: <strong>Beta</strong>.", 
    "B.1.525" = "UK ex West Africa. Q52R, A67V, Δ69-70, Δ144, E484K, Q677H and F888L. WHO label: <strong>Eta</strong>.",
    "B.1.526" = "New York, USA. L5F, T95I, D253G, E484K or S477N and A701V. WHO label: <strong>Iota</strong>.",
    "B.1.429" = "California, USA. D614G, G1251V, L452R, P26S, S13I, S1252C and W152C. WHO label: <strong>Epsilon</strong>.",
    "B.1.324.1" = "UK associated variant. E484K, S494P, N501Y, D614G, P681H and E1111K in the Spike. ",
    "P.2" = "Brazil. E484K and V1176F. WHO label: <strong>Zeta</strong>.",
    "P.3" = "The Philippines. Δ141-143, E484K, N501Y, P681H, E1092K, H1101Y, V1176F and in some cases Δ243-244. WHO label: <strong>Theta</strong>.",
    
    "C.37" = "South America. Δ246-252, G75V, T76I, L452Q, F490S, D614G, and T859N. WHO label: <strong>Lambda</strong>.",
    "B.1.621" = "Colombia. T95I, R346K, E484K, N501Y and P681H. WHO label: <strong>Mu</strong>.",
    "B.1.617.1" = "India. E154K, L452R, E484Q and P681R. WHO label: <strong>Kappa</strong>.",
    
    # Delta
    "B.1.617.2" = "India. T19R, G142D, Δ156-157, R158G, L452R, T478K, D614G, P681R and D950N. WHO label: <strong>Delta</strong>.",
    "AY.4" = "Alias of B.1.617.2.4, UK. WHO label: <strong>Delta</strong>.",
    "AY.4.2" = "Subineage of AY.4 with the addition of Y145H and A222V. Alias of  B.1.617.2.4.2. WHO label: <strong>Delta</strong>.",
    "AY.4.2.1" = "Sublineage of AY.4.2 with the addition of V36F. Alias of  B.1.617.2.4.2.1. WHO label: <strong>Delta</strong>.",
    
    # Gamma
    "P.1" = " Japan ex Brazil. L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, H655Y and T1027I. WHO label: <strong>Gamma</strong>.",
    "P.1.8" = "Brazil. S: Gamma + T470N, P681R, C1235F;
NSP3: Gamma + I441V; NSP4: A446V; ORF3a: S216L; ORF8: G8*STOP; N: TRS insertion. WHO label: <strong>Gamma</strong>.",

  # Omicron
  "Unassigned" = "Omicron sequences not yet assigned lineages. WHO label: <strong>Omicron</strong>.",
  "BA.2.75" = "Sublineage of BA.2. 33 non-synonymous mutations
in Spike. Reversion in Spike relative to BA.2: R493Q. K147E, W152R, F157L, I210V and G257S in the N-
terminal domain. G339H, G446S, and N460K in the receptor binding domain. WHO label: <strong>Omicron</strong>."
  ) %>% 
  enframe("lineage", "reason")

# Recombinant lineages and descriptions
lineages_recomb <- 
  c("XA" = "Recombinant lineage with parental lineages B.1.1.7 and B.1.177, UK.",
    "XB" = "Recombinant lineage with parental lineages B.1.634 and B.1.631, Central and North America.",
    "XE" = "Recombinant lineage of BA.1 and BA.2, UK.",
    "XN" = "Recombinant lineage of BA.1 and BA.2, UK.",
    "XL" = "Recombinant lineage of BA.1 and BA.2, UK.",
    "XM" = "Recombinant lineage of BA.1.1 and BA.2, Europe.",
    "XQ" = "Recombinant lineage of BA.1.1 and BA.2, UK.",
    "XR" = "Recombinant lineage of BA.1.1 and BA.2, UK.",
    "XP" = "Recombinant lineage of BA.1.1 and BA.2, UK.",
    "XF" = "Recombinant lineage of Delta and BA.1, UK.",
    "XG" = "Recombinant lineage of BA.1 and BA.2, Denmark",
    "XH" = "Recombinant lineage of BA.1 and BA.2, Denmark",
    "XJ" = "Recombinant lineage of BA.1 and BA.2, Finland"
  ) %>% 
  enframe("lineage", "reason")

lineages_t2 <- c(vui_voc %$% levels(lineage), lineages_t3$lineage) %>% unique # lineages for counting

sum_lineages <- function(lineages, use_regex = FALSE){
  left_join(
    sum_key_mutations_by_lineage_uk(lineages, use_regex = use_regex),
    sum_key_mutations_by_lineage_uk(lineages, date_from = sample_date_28, use_regex = use_regex) %>%
      rename(n_sequences_28 = n_sequences)
  ) %>% 
    pivot_wider(names_from = adm1, values_from = c(n_sequences, n_sequences_28)) %>% 
    mutate(across(where(is.numeric), ~replace_na(.x, 0L)))
}

# TODO precompute and include lineage/mutation combinations
## count VOC/VUI
n_uk_lineages_all <- sum_lineages(lineages_t2) # count non-omicron VOC/VUI

# TODO refactor
# count Omicron sublineages by filtering full lineage names to include renamed sublineages
n_uk_lineages_ba_1 <- sum_lineages(
  consortium_uk %>% 
    distinct(lineage, lineage_full) %>% 
    filter(lineage == "BA.1" | str_starts(lineage_full, fixed("B.1.1.529.1."))) %$% 
    lineage # sum BA.1 and sublineages
)

n_uk_lineages_ba_2 <- sum_lineages(
  consortium_uk %>% 
    distinct(lineage, lineage_full) %>% 
    filter(lineage == "BA.2" | str_starts(lineage_full, fixed("B.1.1.529.2."))) %$% 
    lineage # sum BA.2 and sublineages
)

n_uk_lineages_ba_3 <- sum_lineages(
  consortium_uk %>% 
    distinct(lineage, lineage_full) %>% 
    filter(lineage == "BA.3" | str_starts(lineage_full, fixed("B.1.1.529.3."))) %$% 
    lineage # sum BA.1 and sublineages
)
n_uk_lineages_ba_4 <- sum_lineages(
  consortium_uk %>% 
    distinct(lineage, lineage_full) %>% 
    filter(lineage == "BA.4" | str_starts(lineage_full, fixed("B.1.1.529.4."))) %$% 
    lineage # sum BA.1 and sublineages
)

n_uk_lineages_ba_5 <- sum_lineages(
  consortium_uk %>% 
    distinct(lineage, lineage_full) %>% 
    filter(lineage == "BA.5" | str_starts(lineage_full, fixed("B.1.1.529.5."))) %$% 
    lineage # sum BA.1 and sublineages
)

n_uk_recombinants <- sum_lineages(
  consortium_uk %>% 
    distinct(lineage) %>% 
    filter(str_starts(lineage, "X")) %$% 
    lineage
) %>% 
  bind_rows(n_uk_lineages_all %>% slice(0)) # dirty hack to ensure all nations columns are included for missing values

# remove VOCs/VUIs with zero counts
vui_voc %<>% 
  semi_join(n_uk_lineages_all %>% filter(variant == "sequences" & n_sequences_UK > 0), by = "lineage") %>% 
  mutate(across(c(lineage, lineage_display), fct_drop)) 

# VUI/VOC lineages with AY.x removed
vui_voc_lineages <- 
  vui_voc %>% 
  filter(!str_starts(lineage, fixed("AY."))) %>% 
  mutate(across(c(lineage, lineage_display), fct_drop)) 

vui_voc_lineages <- 
  levels(vui_voc_lineages$lineage) %>% 
  setNames(levels(vui_voc_lineages$lineage_display)) %>% 
  append(list("B.1.177/B.1.177.x" = "B.1.177", 
              "BA.4/BA.4.x (Omicron)" = "BA.4",
              "BA.5/BA.5.x (Omicron)" = "BA.5",
              "Unassigned (Omicron)" = "Unassigned")) # special cases, add only to lineage bar chart

geo_all <- str_c(dataset_date, "/geo_all.rds") %>% read_rds # geographical NUTS1 counts
mapdata <- read_rds("mapdata.rds") # UK map NUTS1 topology as dataframe

# Waiter loading data screen
loading_screen <- tagList(
  h3("COG-UK Mutation Explorer", style = "color:#333"),
  spin_plus(),
  h4("Loading data . . .", style = "color:#333")
) 