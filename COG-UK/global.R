library(tidyverse)
library(lubridate)
library(magrittr)
library(RColorBrewer)

dataset_date <- ymd("2021-04-01") #TODO derive from directory name

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
    "A.23.1" = "International variant with mutations of biological significance: F157L, V367F, Q613H and P681R. Q613H is predicted to be functionally equivalent to the D614G mutation that arose early in 2020.",
    "A.27" = "Variant associated with France. Comprises several mutations and in particular, in the Spike, L452R N501Y and Q677H.",
    "B.1.1.28" = "Variant associated with the Philippines.",
    "B.1.1.318" = "Variant associated with England. Has D614G, D796H, E484K, P681H, T95I and Y144del in the Spike.", 
    "B.1.1.7" = "UK associated variant. Has 17 mutations (14 replacements and 3 deletions) including: T1001I, A1708D, I2230T, SGF 3675-3677del in the ORF1ab; 69-70del, Y144del, N501Y, A570D, P681H, T716I, S982A and D1118H in the Spike; Q27stop, R52I and Y73C in ORF8; D3L and S235F in the N. Noteworthily, N501Y enhances ACE2 binding affinity, and P681H occurs at the furin cleavage site, known for biological significance in membrane fusion.", 
    "B.1.351" = "Variant associated with South Africa. Has eight mutations in the Spike: D80A, D215G, E484K, N501Y, A701V, L18F, R246I and K417N. Three of these in the RBM: K417N, E484K and N501Y. K417N and E484K have been shown to escape some mAbs.", 
    "P.1" = "Variant associated with Brazil. Has 10 mutations in the Spike, including L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, H655Y and T1027I. Noteworthy  E484K, N501Y and K417T have biological significance.",
    "B.1.525" = "International variant with mutations of biological significance E484K, Q677H, F888L and a similar suite of deletions to B.1.1.7.",
    "B.1.429" = "Variant associated with California, USA. Comprises D614G, G1251V, L452R, P26S, S13I, S1252C and W152C in the Spike.",
    "B.1.526" = "Variant associated with New York, USA. Has A701V, D253G, D614G, E484K, G1251V, L5F, S982A, S1252C and T95I in the Spike.",
    # "B.1.324.1" = "UK associated variant. Has E484K, S494P, N501Y, D614G, P681H and E1111K in the Spike. ",
    # "P.3" = It has 141-143del, E484K, N501Y, D614G, P681H, E1092K, H1101Y, V1176F and in some cases 243-244del
    "P.2" = "Variant associated with Brazil. It has D614G, E484K and V1176F in the Spike."
    # "B.1.617" = "Originally identified in India it has six amino acid substitution in the Spike protein: G142D, E154K, L452R, E484K, P681R and Q1071H."
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