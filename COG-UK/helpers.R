library(tidyverse)

library(ComplexHeatmap)
library(circlize)

library(Cairo)
library(dendextend)
library(seriation)

antibody_complex_heatmap <- function(mutations_lineages_epi_weeks, spike_database){
  
  RBD1_class <- c(403, 405, 406, 408, 409, 414, 415, 416, 417, 420, 421, 449, 453, 455, 456, 457, 458, 459, 460, 473, 474, 475, 476, 477, 484, 486, 487, 489, 490, 492, 493, 494, 495, 496, 498, 500, 501, 502, 503, 504, 505)
  RBD2_class <- c(338, 339, 342, 343, 346, 351, 368, 371, 372, 373, 374, 403, 405, 406, 417, 436, 444, 445, 446, 447, 448, 449, 450, 452, 453, 455, 456, 470, 472, 473, 475, 478, 481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505)
  RBD3_class <- c(333, 334, 335, 337, 339, 340, 341, 342, 343, 344, 345, 346, 354, 356, 357, 358, 359, 360, 361, 438, 439, 440, 441, 442, 443, 446, 499, 500)
  RBD4_class <- c(369, 370, 371, 372, 374, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 390, 430, 431)
  NTD_class <- c(15, 18, 19, 22, 28, 74, 77, 80, 123, 136, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 151, 152, 157, 158, 164, 244, 246, 247, 248, 249, 250, 251, 252, 253, 255, 257, 258)
  
  horz_heat <-
    mutations_lineages_epi_weeks %>% 
    filter(lineage == "B.1.1.7" & variant != "N501Y") %>% 
    select(-lineage) %>% 
    inner_join(spike_database %>% 
                 select(position, mutation, mab, plasma, vaccine_sera, support, domain) %>% 
                 add_row(position = 243, mutation = "del243-244", mab = TRUE, plasma = NA, vaccine_sera = NA, support = "lower", domain = "NTD"), 
               by = c("variant" = "mutation")) %>%   
    mutate(across(where(is.logical), ~na_if(.x, FALSE))) %>%
    arrange(position, variant) %>% 
    mutate(across(domain, as_factor)) %>% 
    rename(confidence = support) %>% 
    column_to_rownames("variant")

  horz_heat$RBD1 <- ifelse(horz_heat$position %in% RBD1_class, TRUE, NA)
  horz_heat$RBD2 <- ifelse(horz_heat$position %in% RBD2_class, TRUE, NA)
  horz_heat$RBD3 <- ifelse(horz_heat$position %in% RBD3_class, TRUE, NA)
  horz_heat$RBD4 <- ifelse(horz_heat$position %in% RBD4_class, TRUE, NA)
  horz_heat$NTD.1 <- ifelse(horz_heat$position %in% NTD_class, TRUE, NA)
  
  input <- data.matrix(horz_heat)
  
  # define colour heatmap for frequency
  col_fun = colorRamp2(c( 0, 0.015, 0.5, 2), c("white", "darkolivegreen1","darkolivegreen3","forestgreen"))
  
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
        "FP" = "seashell2",
        "NTD" = "navajowhite",
        "RBD" = "pink",
        "RBM" = "plum1",
        "SP" = "lightblue1"
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
    column_title = "Epiweeks",
    column_title_side = "bottom",
    column_title_gp = gpar(fontsize = 20),
    
    use_raster = TRUE,
    raster_device = "CairoPNG",
    #TODO raster_by_magick
    
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_order = order((horz_heat$domain)),
    row_split = (horz_heat$domain),
    column_names_rot = 0,
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
