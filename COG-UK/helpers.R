library(tidyverse)

library(ComplexHeatmap)
library(circlize)

library(Cairo)
library(dendextend)
library(seriation)

antibody_complex_heatmap <- function(mutations_lineages_epi_weeks, spike_database){
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
    select(-position) %>% 
    rename(confidence = support) %>% 
    column_to_rownames("variant")
  
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
  row_ha2 = rowAnnotation(Domain = (horz_heat$domain),
                          col = list(
                            Domain = c(
                              "FP" = "seashell2",
                              "NTD" = "navajowhite",
                              "RBD" = "pink",
                              "RBM" = "plum1",
                              "SP" = "lightblue1"
                            )
                          ))
  
  heatmap <- Heatmap(
    subset(input, select = -c(mab:domain)),
    name = "Percentage %",
    column_title = "Antigenic Mutations in Lineage B.1.1.7",
    column_title_gp = gpar(fontsize = 18),
    
    use_raster = TRUE,
    raster_device = "CairoPNG",
    #TODO raster_by_magick
    
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    row_order = order((horz_heat$domain)),
    row_split = (horz_heat$domain),
    column_names_rot = 0,
    row_gap = unit(4, "mm"),
    border = TRUE,
    width = ncol(input) * unit(3.6, "mm"),
    height = nrow(input) * unit(3.6, "mm"),
    col = col_fun,
    na_col = 'white',
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 8),
    right_annotation = row_ha,
    left_annotation = row_ha2
  )
  heatmap
}