library(ggplot2)

selected_variants <- c("B.1.1.7","AY.4.2", "B.1.617.2")
vui_voc_lineages <- 
  levels(vui_voc$lineage) %>% 
  append("Other Delta", after = 9) %>% 
  append("Other")

# variants by week

  variants_other_week <- 
    lineages_weeks_uk_all %>% 
    filter(!(lineage %in% selected_variants)) %>% 
    group_by(epi_date) %>% 
    summarise(n_week = sum(n_week)) %>% 
    mutate(lineage = "Other", .before = epi_date) %>% 
    ungroup
  
  
  lineages_weeks_uk <- 
    lineages_weeks_uk_all %>% 
    filter(lineage %in% selected_variants) %>% 
    bind_rows(variants_other_week) %>% 
       mutate(lineage = recode_factor(lineage, # recode WHO Greek display names as factor and order levels to define colour/legend order
                                   "AV.1" = "AV.1", 
                                   "B.1.1.318" = "B.1.1.318",
                                   "B.1.1.7" = "B.1.1.7 (Alpha)",
                                   "B.1.351" = "B.1.351 (Beta)",
                                   "B.1.525" = "B.1.525 (Eta)",
                                   "B.1.617.1" = "B.1.617.1 (Kappa)",
                                   "B.1.617.2" = "B.1.617.2/AY.x (Delta)",
                                   "AY.4" = "AY.4/AY4.x (Delta)",
                                   "AY.4.2" = "AY.4.2 (Delta)",
                                   "B.1.617.3" = "B.1.617.3",
                                   "P.1" = "P.1 (Gamma)",
                                   "P.2" = "P.2 (Zeta)",
                                   "P.3" = "P.3 (Theta)",
                                   "Other Delta" = "Other Delta",
                                   "Other" = "Other"
    )) %>% 
    rename(Variant = lineage, `Start date` = epi_date, Sequences = n_week)
  
  
  vui_voc_plot <- 
    lineages_weeks_uk %>%
    filter(`Start date` >= "2020-03-15" & `Start date` <= "2021-11-09") %>% 
    ggplot(aes(fill = Variant, y = Sequences, x = `Start date`) ) +
    theme_classic() +
    
    scale_fill_discrete_qualitative(palette = "Dynamic", 
                                    nmax = vui_voc_lineages %>% length , # extra colour for Other
                                    
                                    # fix variant/colour combos plus extra colour for Other
                                    order = match(selected_variants, 
                                                  vui_voc_lineages) %>% c(vui_voc_lineages %>% length)
    ) +
    
    scale_x_date(breaks = date_breaks("1 month"),
                 labels = date_format("%b %y")) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Sample date",
         y = "Sequences"
    )
  
  ymax <- 
    ymax <- 
    lineages_weeks_uk %>% 
    group_by(`Start date`) %>% 
    summarise(total_week = sum(Sequences)) %$% 
    max(total_week)
  
  vui_voc_plot <- 
    vui_voc_plot + 
    geom_bar(position="stack", stat="identity") +
    ylim(0, ymax)
  
  
  vui_voc_plot
