suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(Cairo)
suppressPackageStartupMessages(library(dendextend))
library(seriation)
library(magick)
library(colorspace)

# Epidemic week / Sunday date conversion
epi_lookup <-
  tibble(
    epi_date = seq(from = ymd("2020-01-26"), to = consortium_uk %$% max(sample_date), by = "week"),
    epi_week = consortium_uk %$% levels(epi_week) %>% as.integer %>% sort
  )

# https://slowkow.com/notes/pheatmap-tutorial/#quantile-breaks
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

antibody_complex_heatmap <- function(mutations_lineages_epi_weeks){
  
  RBD1_class <- c(403, 405, 406, 408, 409, 414, 415, 416, 417, 420, 421, 449, 453, 455, 456, 457, 458, 459, 460, 473, 474, 475, 476, 477, 484, 486, 487, 489, 490, 492, 493, 494, 495, 496, 498, 500, 501, 502, 503, 504, 505)
  RBD2_class <- c(338, 339, 342, 343, 346, 351, 368, 371, 372, 373, 374, 403, 405, 406, 417, 436, 444, 445, 446, 447, 448, 449, 450, 452, 453, 455, 456, 470, 472, 473, 475, 478, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505)
  RBD3_class <- c(333, 334, 335, 337, 339, 340, 341, 342, 343, 344, 345, 346, 354, 356, 357, 358, 359, 360, 361, 438, 439, 440, 441, 442, 443, 446, 499, 500)
  RBD4_class <- c(369, 370, 371, 372, 374, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 390, 430, 431)
  NTD_class <- c(15, 18, 19, 22, 28, 74, 77, 80, 123, 136, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 151, 152, 157, 158, 164, 244, 246, 247, 248, 249, 250, 251, 252, 253, 255, 257, 258)
  
  horz_heat <-
    mutations_lineages_epi_weeks %>% 
    inner_join(database %>% 
                 select(position, mutation, mab, plasma, vaccine_sera, support, domain) %>% 
                 add_row(position = 243, mutation = "del243-244", mab = TRUE, plasma = NA, vaccine_sera = NA, support = "lower", domain = "NTD"), 
               by = c("variant" = "mutation")) %>%   
    mutate(across(where(is.logical), ~na_if(.x, FALSE))) %>%
    arrange(position, variant) %>% 
    mutate(across(domain, as_factor)) %>% 
    rename(confidence = support) %>% 
    column_to_rownames("variant") 
  
  qb <- 
    horz_heat %>% 
    select(-(position:domain)) %>% 
    unlist(use.names = FALSE) %>% 
    quantile_breaks(101) 
  
  col_fun <- colorRamp2(qb, sequential_hcl(qb %>% length, palette = "Greens 3", rev = TRUE))
    
  horz_heat$RBD1 <- ifelse(horz_heat$position %in% RBD1_class, TRUE, NA)
  horz_heat$RBD2 <- ifelse(horz_heat$position %in% RBD2_class, TRUE, NA)
  horz_heat$RBD3 <- ifelse(horz_heat$position %in% RBD3_class, TRUE, NA)
  horz_heat$RBD4 <- ifelse(horz_heat$position %in% RBD4_class, TRUE, NA)
  horz_heat$NTD.1 <- ifelse(horz_heat$position %in% NTD_class, TRUE, NA)

  input <- data.matrix(horz_heat) 
  
  # annotation row
  row_ha = rowAnnotation(
    `Effect mAb` = horz_heat$mab,
    `Effect plasma` = horz_heat$plasma,
    `Effect vaccine` = horz_heat$vaccine_sera,
    Confidence = horz_heat$confidence,
    na_col = 'white',
    col = list(
      Confidence = c(
        "lower" = "lightgoldenrod",
        "medium" = "lightgoldenrod3",
        "high" = "lightgoldenrod4"
      ),
      `Effect mAb` = c("TRUE" = "black"),
      `Effect plasma` = c("TRUE" = "black"),
      `Effect vaccine` = c("TRUE" = "black")
    ),
    annotation_legend_param = list(Confidence = list (at = c(
      "high", "medium", "lower"
    )))
  )
  
  domain_palette <- qualitative_hcl(6, palette = "Set 3")
  
  # domain
  row_ha2 = rowAnnotation(
    Domain = horz_heat$domain,
    `Ab class 1` = horz_heat$RBD1,
    `Ab class 2` = horz_heat$RBD2,
    `Ab class 3` = horz_heat$RBD3,
    `Ab class 4` = horz_heat$RBD4,
    `Ab class 5` = horz_heat$NTD.1,
    col = list(
      Domain = c(
        "SP" = domain_palette[4],
        "NTD" = domain_palette[2],
        "FP" = domain_palette[3],
        "RBD" = domain_palette[1],
        "RBM" = domain_palette[6],
        "S2" = domain_palette[5]
      ),
      `Ab class 1` = c ("TRUE" = "lightgreen"),
      `Ab class 2` = c("TRUE" = "goldenrod1"),
      `Ab class 3` = c("TRUE" = "cornflowerblue"),
      `Ab class 4` = c("TRUE" = "tomato"),
      `Ab class 5` = c("TRUE" = "magenta")
    ),
    show_legend = c(
      Domain = FALSE,
      `Ab class 1` = FALSE,
      `Ab class 2` = FALSE,
      `Ab class 3` = FALSE,
      `Ab class 4` = FALSE,
      `Ab class 5` = FALSE
    )
  )
  
  input <- subset(input, select = -c(position:NTD.1))
  
  heatmap <- Heatmap(
    input,
    name = "Percentage %",
    column_title = "Sample date",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 20),
    
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_by_magick = TRUE,
    
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_order = order((horz_heat$domain)),
    row_split = (horz_heat$domain),
    column_names_rot = 90,
    row_gap = unit(3, "mm"),
    border = TRUE,
    width = ncol(input) * unit(4.5, "mm"),
    height = nrow(input) * unit(4.5, "mm"),
    col = col_fun,
    na_col = 'white',
    column_names_gp = grid::gpar(fontsize = 9),
    row_names_gp = grid::gpar(fontsize = 11),
    right_annotation = row_ha,
    left_annotation = row_ha2
  )
  heatmap
}

antigenic_mutations_lineages <- function(nation = c("UK", "England", "Scotland", "Wales", "Northern_Ireland"), lineage = "B.1.1.7", defining = "N501Y"){
  nation = match.arg(nation)
  
  if(lineage %in% c("AY.4", "AY.4.2", "BA.1")) { # alias lineage name - include sublineages of alias name
    mutations_s_uk %<>% 
      filter( (lineage == !!lineage | str_starts(lineage, fixed(str_c(!!lineage, "."))) ) & !(variant %in% !!defining))
    
    consortium_uk %<>% 
      filter( lineage == !!lineage | str_starts(lineage, fixed(str_c(!!lineage, "."))) )
  } else { # lineage name is not an alias - include sublineages by full lineage name
    mutations_s_uk %<>% 
      filter((lineage == !!lineage | str_starts(lineage_full, fixed(str_c(!!lineage, "."))) )& !(variant %in% !!defining))
    
    consortium_uk %<>% 
      filter( lineage == !!lineage | str_starts(lineage_full, fixed(str_c(!!lineage, "."))) )
    
  } 

  if(nation != "UK"){
    mutations_s_uk %<>% 
      filter(adm1 == nation) 
    
    consortium_uk %<>%
      filter(adm1 == nation) 
  }
  
  del_22289_6_samples <- 
    deletions %>% 
    filter(
      (ref_start == 22289 & length == 6) 
    ) %$% samples 
  
  ### Antigenic deletions by lineage 22289-22294 (6nt)
  del_22289_6 <- 
    consortium_uk %>% 
    filter(sequence_name %in% del_22289_6_samples) %>% 
    dplyr::count(epi_week) %>% 
    mutate(variant = "del243-244", .before = 1)
  
  sequences_by_week_lineages <- 
    consortium_uk %>% 
    dplyr::count(epi_week, name = "n_sequences_lineage")
  
  escape_mutations <-
    database %>%
    filter(!is.na(escape)) %$% 
    mutation 
  
  antigenic_mutations <- 
    mutations_s_uk %>% 
    filter(variant %in% escape_mutations) %>% 
    filter(sequence_name != "England/NEWC-2729532/2021") %>% # remove Delta outlier from January 2021
    dplyr::count(variant, epi_week, sort = TRUE) %>% 
    bind_rows(del_22289_6)
  
  # no non-defining antigenic mutations in lineage
  if(plyr::empty(antigenic_mutations)){
    return(NULL)
  }
  
  antigenic_mutations_lineages_all <- 
    inner_join(antigenic_mutations, sequences_by_week_lineages) %>% 
    mutate(percentage = n / n_sequences_lineage * 100 ) 
  
  antigenic_mutations_lineages <- 
    antigenic_mutations_lineages_all %>% 
    complete(epi_week, nesting(variant), fill = list(n = 0, n_sequences_lineage = 0, percentage = 0)) %>%
    mutate(epi_week = epi_week %>% as.character %>% as.integer)

  # remove epiweeks before first occurrence
  first_occurrence <- 
    antigenic_mutations_lineages %>% 
    filter(n > 0) %$% 
    min(epi_week)
  
  antigenic_mutations_lineages %<>%
    filter(epi_week >= first_occurrence) %>%
    inner_join(epi_lookup) %>% 
    pivot_wider(names_from = epi_date, values_from = percentage, names_sort = TRUE, values_fill = 0, id_cols = variant) 
  
  antigenic_mutations_lineages
}

sum_key_mutations_by_lineage_uk <- function(lineages = NULL, date_from = NULL, use_regex = FALSE){
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
        when(
          use_regex ~filter(., lineage == x | str_detect(lineage, sublineage_regex(x))),
          ~filter(., lineage == x)
        ) %>%
        group_by(adm1) %>% 
        select(-lineage) %>% 
        summarise_all(funs(sum)) %>% #TODO replace deprecated funs
        mutate(lineage = x, .before = 1)
    }) %>% 
      bind_rows() %>% 
      gather(key = "variant", value = "n_sequences", sequences:`N501Y + E484K`)
  }
}

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