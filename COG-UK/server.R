library(shiny)
library(plotly)
library(scales)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(ggseqlogo)
library(RColorBrewer)

library(maps)
library(mapdata)
library(maptools)
library(ggmap)
library(ggplot2)
# library(plyr)
library(viridis)

lineage_plus_variant <- function(lineage, variant, variant2 = NULL, use_regex = FALSE){
  
  mutations_s_uk_lv <- 
    mutations_s_uk %>%
    when(
      use_regex ~filter(., lineage == !!lineage | str_detect(lineage, sublineage_regex(!!lineage))),
      ~filter(., lineage == !!lineage)
    ) 
    
  if(is.null(variant2)){
    mutations_s_uk_lv %<>%  
      filter(variant == !!variant) 
  } else {
    mutations_s_uk_lv1 <-  
      mutations_s_uk_lv %>% 
      filter(variant == !!variant) %>% select(-variant, -position)
    
    mutations_s_uk_lv2 <-  
      mutations_s_uk_lv %>% 
      filter(variant == !!variant2) %>% select(-variant, -position)
    
    mutations_s_uk_lv <- 
      intersect(mutations_s_uk_lv1, mutations_s_uk_lv2)
  }
  
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
    pivot_wider(names_from = adm1, values_from = c(n_sequences, n_sequences_28)) %>%
    mutate(lineage = !!lineage, variant = !!variant, .before = 1) %>% 
    when(
      is.null(variant2) ~mutate(., variant = !!variant),
      ~mutate(., variant = str_c(!!variant, " + ", !!variant2))
    ) 
  
}

## Table functions
# TODO table caching
# 
# Mutations
table_1 <- function(){
  database_genome %>% 
    arrange(desc(`numSeqs UK`)) %>% 
    filter(`numSeqs UK` >= 5) %>% 
    mutate(mutation = mutation %>% fct_drop %>% fct_inorder) %>% 
    select(gene, mutation, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`, earliest) %>% 
    mutate(`Cumulative sequences in UK (%)` = `numSeqs UK` / total_sequences,
           .after = `numSeqs UK`) %>%
    mutate(`Sequences over the last 28 days in UK (%)` = `numSeqs UK 28 days` / total_sequences_28,
           .after = `numSeqs UK 28 days`) %>%
    rename(Gene = gene,
           `Amino acid replacement` = mutation, 
           `Cumulative sequences in UK` = `numSeqs UK`, 
           `Sequences over the last 28 days in UK` = `numSeqs UK 28 days`,
           `Sequences over the last 28 days in England` = `numSeqs Eng 28 days`,
           `Sequences over the last 28 days in Scotland` = `numSeqs Scotland 28 days`,
           `Sequences over the last 28 days in Wales` = `numSeqs Wales 28 days`,
           `Sequences over the last 28 days in Northern Ireland` = `numSeqs NI 28 days`,
           `Date of first detection in UK` = earliest) 
}

# no longer used
table_2 <- function(){
  n_uk_lineages_all %>% 
    filter(    
      (variant == "D614G" & lineage == "B.1" ) | 
                 
     (variant == "A222V" & lineage == "B.1.177" ) | 
     
     (variant == "N439K" & lineage == "B.1.141" ) | 
     (variant == "N439K" & lineage == "B.1.258" ) |
     (variant == "N439K + ∆69-70" & lineage == "B.1.258" ) |
     
     (variant == "∆69-70" & lineage == "B.1.1" ) |
     (variant == "∆69-70" & lineage == "B.1.258" ) |
     
     (variant == "N501Y + ∆69-70" & lineage == "B.1.1.7" ) |
     (variant == "N501Y" & lineage == "B.1.1.70" ) |
     
     (variant == "Y453F" & lineage == "B.1.1" ) |
     (variant == "Y453F" & lineage == "B.1.1.298" ) |
     
     (variant == "E484K" & lineage == "B.1.351" ) |
     (variant == "E484K" & lineage == "P.2" ) |
     (variant == "E484K" & lineage == "P.1" ) |
     (variant == "E484K" & lineage == "A.23.1" ) |
     (variant == "E484K" & lineage == "B.1.1.119" ) |
     (variant == "E484K" & lineage == "B.1.177.4" ) |
     (variant == "E484K" & lineage == "B.1.222" ) |
     (variant == "E484K" & lineage == "B.1.177" ) |
     (variant == "E484K" & lineage == "B.1" ) |
     (variant == "E484K" & lineage == "B.1.525" ) |
     
     (variant == "N501Y + E484K" & lineage == "B.1.351") |
     (variant == "N501Y + E484K" & lineage == "B.1.1.7")) %>% 
    relocate(variant) %>% 
    relocate(n_sequences_UK, .after = lineage) %>% 
    mutate(`UK (%)` = n_sequences_UK / total_sequences,
           .after = n_sequences_UK) %>%
    relocate(n_sequences_28_UK, .after = `UK (%)`) %>% 
    mutate(`UK 28 days (%)` = n_sequences_28_UK / total_sequences_28,
           .after = n_sequences_28_UK) %>%
    relocate(n_sequences_28_England, .after = n_sequences_England) %>% 
    relocate(n_sequences_28_Northern_Ireland, .after = n_sequences_Northern_Ireland) %>% 
    relocate(n_sequences_28_Scotland, .after = n_sequences_Scotland) %>% 
    relocate(n_sequences_28_Wales, .after = n_sequences_Wales) %>% 
    arrange(variant, lineage) %>% 
    rename(Mutation = variant, 
           `Lineage` = lineage, 
           `UK` = n_sequences_UK, 
           `UK 28 days` = n_sequences_28_UK,
           
           England = n_sequences_England,
           `Northern Ireland` = n_sequences_Northern_Ireland,
           Scotland = n_sequences_Scotland,
           Wales = n_sequences_Wales,
           
           `England 28 Days` = n_sequences_28_England,
           `Northern Ireland 28 Days` = n_sequences_28_Northern_Ireland,
           `Scotland 28 Days` = n_sequences_28_Scotland,
           `Wales 28 Days` = n_sequences_28_Wales
    )
}

# Variants
table_3 <- function(){
  bind_rows(
    n_uk_lineages_all %>% 
      filter(lineage %in% lineages_t3$lineage & variant == "sequences") %>% 
      inner_join(lineages_t3) %>% 
      select(-variant) %>% 
      relocate(reason, .after = lineage),
    
    n_uk_lineages_all %>%
      filter(variant == "E484K" & lineage == "A.23.1") %>%
      mutate(lineage = str_c(lineage, " + ", variant), .keep = "unused")  %>%
      mutate(reason = "As A.23.1, with the addition of E484K."),
    
    n_uk_lineages_all %>%
      filter(variant == "E484K" & lineage == "B.1.1.7") %>%
      mutate(lineage = str_c(lineage, " + ", variant), .keep = "unused")  %>%
      mutate(reason = "As B.1.1.7, with the addition of E484K."),
    
    lineage_plus_variant(lineage = "B.1.1.7", variant = "S494P") %>% # for non-key mutations
      mutate(lineage = str_c(lineage, " + ", variant), .keep = "unused") %>%
      mutate(reason = "As B.1.1.7, with the addition of S494P."),
    
    # lineage_plus_variant(lineage = "AY.4", variant = "A222V", variant2 = "Y145H") %>% # for non-key double mutations
    #   mutate(lineage = 'AY.4.2', .keep = "unused") %>%
    #   mutate(reason = "Sublineage of interest carrying a further set of mutations. As AY.4, with the addition of Y145H and A222V. <strong>Note</strong>: Lineage is temporarily undercounted due to sequencing aftefact"),

    n_uk_lineages_all %>%
      filter(variant == "E484K" & lineage == "B.1.324.1") %>%
      mutate(lineage = str_c(lineage, " + ", variant), .keep = "unused")  %>%
      mutate(reason = "As B.1.324.1, with the addition of E484K.")
  ) %>% 
    mutate(across(everything(), ~replace_na(.x, 0L))) %>% 
    filter(n_sequences_UK > 0) %>%
    relocate(n_sequences_UK, .after = reason) %>% 
    mutate(`UK (%)` = n_sequences_UK / total_sequences,
           .after = n_sequences_UK) %>%
    relocate(n_sequences_28_UK, .after = `UK (%)`) %>% 
    mutate(`UK 28 days (%)` = n_sequences_28_UK / total_sequences_28,
           .after = n_sequences_28_UK) %>%
    arrange(desc(n_sequences_28_UK), desc(n_sequences_UK), lineage) %>% 
    relocate(n_sequences_28_England, .after = n_sequences_England) %>% 
    relocate(n_sequences_28_Northern_Ireland, .after = n_sequences_Northern_Ireland) %>% 
    relocate(n_sequences_28_Scotland, .after = n_sequences_Scotland) %>% 
    relocate(n_sequences_28_Wales, .after = n_sequences_Wales) %>% 
    rename(`Variant` = lineage,	
           UK = n_sequences_UK,	 
           `UK 28 days` = n_sequences_28_UK,                    
           `Country first detected & lineage defining mutations in Spike` = reason, 
           
           England = n_sequences_England,
           `Northern Ireland` = n_sequences_Northern_Ireland,
           Scotland = n_sequences_Scotland,
           Wales = n_sequences_Wales,
           
           `England 28 Days` = n_sequences_28_England,
           `Northern Ireland 28 Days` = n_sequences_28_Northern_Ireland,
           `Scotland 28 Days` = n_sequences_28_Scotland,
           `Wales 28 Days` = n_sequences_28_Wales
    ) 
}

# T cell table
table_5 = function(assay_filter = NULL){
  if(!is.null(assay_filter)){
    if(assay_filter == "recognition"){
      database_tcell_predictions %<>% 
        filter(assay %in% c("Reduced T-cell recognition (full)", "Reduced T-cell recognition (partial)"))
    } else { # epitope_studies
      database_tcell_predictions %<>% 
        filter(!assay %in% c("Reduced T-cell recognition (full)", "Reduced T-cell recognition (partial)"))
    }
    
    database_tcell_predictions %<>% mutate(across(c(gene, mutation, Epitope, CD4_CD8, HLA, assay), fct_drop)) # drop unused factor levels
  }
  
  database_tcell_predictions %>% 
    filter(`numSeqs UK` > 0) %>%
    select(gene, mutation, Epitope:Fold, `numSeqs UK`, `numSeqs UK 28 days`, -`End position`) %>%
    arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), desc(Fold)) %>%
    mutate(across(c(Epitope, HLA, assay, CD4_CD8), ~fct_relevel(.x, sort))) %>% 
    rename(`Gene` = gene,
           `Amino acid replacement` = mutation,
           `Cumulative sequences in UK` = `numSeqs UK`,
           `Sequences over 28 days` = `numSeqs UK 28 days`,
           `Assay` = assay,
           `CD4 CD8` = CD4_CD8,
           Start = `Start position`,
           `Mut Percentile Rank Value` = IC50_mutation,
           `WT Percentile Rank Value` = `IC50 WT`,
           `Fold difference` = Fold,
           DOI = doi)
}

table_therapeutics <- function(){
  therapeutics %>% 
    arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`)) %>% 
    select(gene, variant, Protein, resistance, drug, assay, detail, `quantification (fold)`, note, `numSeqs UK`, `numSeqs UK 28 days`, anchor) %>%
    rename_with(str_to_title, c(gene, resistance, drug, assay, detail, `quantification (fold)`, note)) %>% 
    rename(`Amino acid replacement` = variant, 
           `Cumulative sequences in UK` = `numSeqs UK`,
           `Sequences over 28 days` = `numSeqs UK 28 days`,
           `Reference` = anchor
  ) 
}

shinyServer(function(input, output, session) {

    output$table_1 <- renderDataTable({
      table_1() %>% 
        datatable(filter = "top", rownames = FALSE, 
                  options = list(lengthMenu = c(20, 50, 100, 200), pageLength = 20, scrollX = TRUE)) %>% 
        formatPercentage(c("Cumulative sequences in UK (%)", "Sequences over the last 28 days in UK (%)"), digits = 2)
    })
    
    ######## Data download inputs ########
    
    ## Mutations
    # Reactive value to generate downloadable table for selected table 1 mutation metadata
    datasetInput <- reactive({
      mutations_uk %>% 
        filter(gene == input$dataset_gene) %>% 
        filter(variant == input$dataset) %>% 
        filter(sample_date >= sample_date_28) %>% 
        select(sequence_name, sample_date, epi_week, epi_date, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete table 1 data
    table1Input <- reactive({
      table_1() %>% 
        mutate(across(ends_with("(%)"), ~.x * 100)) # convert decimal fraction to percentage
    })
    
    ## Variants
    # Reactive value to generate downloadable table for selected lineage + mutation
    # TODO regex switch
    concernInput <- reactive({
      if(input$concern == "B.1.1.7 + E484K"){
        concern_download <- 
          consortium_uk %>% 
          filter(e484k == "K") %>% 
          filter(lineage == "B.1.1.7") 
        # filter(lineage == "B.1.1.7" | str_detect(lineage, sublineage_regex("B.1.1.7")))
      } else if(input$concern == "A.23.1 + E484K") {
        concern_download <- 
          consortium_uk %>% 
          filter(e484k == "K") %>% 
          filter(lineage == "A.23.1") 
          # filter(lineage == "A.23.1" | str_detect(lineage, sublineage_regex("B.1.1.7"))) 
      } else if(input$concern == "B.1.1.7 + S494P") {
        concern_download <- 
          mutations_s_uk %>% 
          filter(variant == "S494P") %>% 
          filter(lineage == "B.1.1.7")
          # filter(lineage == "B.1.1.7" | str_detect(lineage, sublineage_regex("B.1.1.7")))
      } else {
        concern_download <- 
          consortium_uk %>% 
          filter(lineage == input$concern)
          # filter(lineage == input$concern | str_detect(lineage, sublineage_regex(input$concern)))
      }
      
      concern_download %>% 
        select(sequence_name, sample_date, epi_week, epi_date, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete table 3 data
    table3Input <- reactive({
      table_3() %>% 
        mutate(across(ends_with("(%)"), ~.x * 100)) # convert decimal fraction to percentage
    })
    
    ## Antigenic mutations
    # Reactive value to generate downloadable table for selected mutation
    escapeInput <- reactive({
      mutations_s_uk %>% 
        filter(variant == input$selectEscape) %>% 
        select(sequence_name, sample_date, epi_week, epi_date, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete T cell table data
    table5Input <- reactive({
      table_5() 
    })
    
    ######## Download handlers ########
    
    ## Table 1
    # Downloadable CSV of selected mutation metadata
    output$downloadData <- downloadHandler(
      filename = function() {
        str_c(input$dataset_gene, "_", input$dataset, "_UK_28_days_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(datasetInput(), file)
      },
      contentType = "text/csv"
    )
    
    # Downloadable CSV of complete table 1 data
    output$downloadTable1 <- downloadHandler(
      filename = function() {
        str_c("table_1_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(table1Input(), file)
      },
      contentType = "text/csv"
    )
    
    ## Table 3
    # Downloadable CSV of selected mutation
    output$downloadConcern <- downloadHandler(
      filename = function() {
        str_c(input$concern, "_UK_cumulative_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(concernInput(), file)
      },
      contentType = "text/csv"
    )
    
    # Downloadable CSV of complete table 1 data
    output$downloadTable3 <- downloadHandler(
      filename = function() {
        str_c("table_3_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(table3Input(), file)
      },
      contentType = "text/csv"
    )
    
    ## Antigenic
    # Downloadable CSV of selected mutation metadata
    output$downloadEscape <- downloadHandler(
      filename = function() {
        str_c(input$selectEscape, "_UK_cumulative_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(escapeInput(), file)
      },
      contentType = "text/csv"
    )
    
    # Downloadable CSV of complete T cell data
    output$downloadTable5 <- downloadHandler(
      filename = function() {
        str_c("tcell_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(table5Input(), file)
      },
      contentType = "text/csv"
    )
    
    # Variants
    output$table_3 <- renderDT({
      table_3() %>% 
        datatable(filter = "none", escape = FALSE, rownames = FALSE, 
                  options = list(dom = 't', paging = FALSE, scrollX = TRUE)) %>% 
        formatPercentage(c("UK (%)", "UK 28 days (%)"), digits = 2)
    })
    
    # Antigenic Mutations
    output$table_4 <- renderDT({
      # database %<>%
      #   filter(`numSeqs UK` > 0) # filter zero counts from predicted antibodies not observed in mutations
      
      if("monoclonal" %in% input$escape){
        database %<>% filter(mab == TRUE)
      }
      
      if("convalescent" %in% input$escape){
        database %<>% filter(plasma == TRUE)
      }
      
      if("vaccine" %in% input$escape){
        database %<>% filter(vaccine_sera == TRUE)
      }
      
      database %>%
        filter(!is.na(escape)) %>% 
          mutate(mutation = fct_drop(mutation)) %>%
          select(mutation, domain, escape, `numSeqs UK`, `numSeqs UK 28 days`, support, anchor) %>%  
          arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), mutation) %>% 
          rename(`Amino acid replacement` = mutation, 
                 `Cumulative sequences in UK` = `numSeqs UK`,
                 `Sequences over 28 days` = `numSeqs UK 28 days`,
                 `Escape mutations details` = escape,
                 `References` = anchor,
                 Confidence = support, 
                 Domain = domain) %>% 
        datatable(filter = "top", escape = FALSE, rownames = FALSE,
                  options = list(lengthMenu = c(20, 50, 100, 200), pageLength = 20, scrollX = TRUE)) %>% 
        formatStyle(
          'Confidence',
          target = 'row',
          backgroundColor = styleEqual(c("lower", "medium", "high"), c('LemonChiffon', 'DarkOrange', 'FireBrick')), 
          color = styleEqual(c("lower", "medium", "high"), c('DarkSlateGray', 'White', 'Snow'))) # TODO anchor colour
        })
    
    output$table_5 <- renderDT({
      table_5(input$t_cell_experiment) %>% 
        mutate(Reference = str_c("<a href='", DOI, "'target='_blank'>", Reference,"</a>"), .keep = "unused", .after = `Supporting references`) %>% # hyperlink to citation DOI
        datatable(filter = "top", escape = FALSE, rownames = FALSE,
                    options = list(lengthMenu = c(10, 20, 50, 100, 200), pageLength = 10, scrollX = TRUE))
    })
    
    output$table_therapeutics <- renderDT({
      table_therapeutics() %>% 
      datatable(filter = "none", escape = FALSE, rownames = FALSE, 
                options = list(dom = 't', paging = FALSE, scrollX = TRUE)) 
        
    })
    
    # output$functional <- renderDT({
    #   functional %>% 
    #     arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`)) %>% 
    #     select(gene, variant, Epitope, `CD4/CD8`, `HLA type`, `funtional observation`, `numSeqs UK`, `numSeqs UK 28 days`, anchor) %>%
    #     rename_with(str_to_title, c(gene)) %>%
    #     rename(`Amino acid replacement` = variant,
    #            `Cumulative sequences in UK` = `numSeqs UK`,
    #            `Sequences over 28 days` = `numSeqs UK 28 days`,
    #            `Reference` = anchor
    #     ) %>% 
    #     datatable(filter = "none", escape = FALSE, rownames = FALSE, 
    #             options = list(dom = 't', paging = FALSE, scrollX = TRUE)) 
    # })
    
    # always display wild type on percentage chart
    observeEvent(input$percentage, {
      if(input$percentage){
        disable("ref")
        updatePrettySwitch(
          session = session,
          inputId = "ref",
          value = as.logical(input$percentage)
        )
      } else {
        enable("ref")
      }
    })

    # Filter mutation data and create plot
    mutation_plot <- reactive({
      mutation_reference_counts %<>% 
        filter(gene == input$gene & position == input$position) %>%
        arrange(desc(n)) %>%
        mutate(variant = 
                 variant %>% 
                 fct_drop %>% 
                 fct_inorder %>% 
                 fct_relevel("WT", after = 0)
               ) %>%   # fix variant colour WT first then by frequency
        filter(adm1 == input$nation) 
        
      variants <- mutation_reference_counts %$% levels(variant)
      
      if(!input$ref){
        mutation_reference_counts %<>% filter(variant != "WT") 
      }
      
      colour_order = match(mutation_reference_counts %$% fct_drop(variant) %>% levels, variants)
      
      mutation_reference_plot <- 
        mutation_reference_counts %>% 
        filter(epi_date >= input$mutation_range[1] & epi_date <= input$mutation_range[2]) %>% 
        rename(`Sample date` = epi_date, Sequences = n, Mutation = variant) %>% # display names
        ggplot(aes(fill = Mutation, y = Sequences, x = `Sample date`) ) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(title = str_c(c("Gene", "Position"), c(input$gene, input$position), sep = " : ", collapse = "\n")) +
        scale_fill_discrete_qualitative(palette = "Set 3", 
                                        nmax = mutation_reference_counts %$% levels(variant) %>% length,
                                        order = colour_order,
                                        rev = FALSE) +
        scale_x_date(breaks = date_breaks("2 month"),
                     labels = date_format("%b %y"))
  
      mutation_reference_plot
    }) 
    
    # mutation plot by percentage or count
    mutation_plot_bar <- reactive({
      
      max_date <- mutation_reference_counts %$% max(epi_date)
      gg_bar <- mutation_plot()
      
      # upper range slider is in most recent 2 weeks
      if (input$mutation_range[2] >= max_date - days(7)){
        if(input$mutation_range[2] >= max_date ){
          xmax = max_date + days(3)
        } else{
          xmax = max_date + days(-4)
        }
        
        if(input$percentage){
          ymax <- 1
        } else {
          # TODO cache maximum weekly counts for UK and nations
          if(input$nation == "UK"){
            ymax <-
              consortium_uk %>%
              dplyr::count(epi_date) %$%
              max(n)
          } else {
            ymax <-
              consortium_uk %>%
              filter(adm1 == input$nation) %>% 
              dplyr::count(epi_date) %$%
              max(n)
          }
        }
        
        gg_bar <- 
          gg_bar + 
          annotate("rect",
                   xmax = xmax,
                   xmin = max_date - days(10),
                   ymin = 0,
                   ymax = ymax,
                   alpha = 0.2)
      }
      
      if(input$percentage){
        gg_bar <- 
          gg_bar +
          geom_bar(position="fill", stat="identity") +
          scale_y_continuous(labels = scales::percent_format())
      } else {
        gg_bar <- 
          gg_bar + 
          geom_bar(position="stack", stat="identity")
      }
      
      gg_bar
    }) %>% debounce(500) # allow 500ms to update percentage switch so don't display plot immediately
    
    # Display mutation plot with percentage option
    output$mutation_time <- renderPlotly({
      mutation_plot_bar() %>% ggplotly
    })
    
    values <- reactiveValues()

    observe({
      defining <- 
        vui_voc %>% 
        filter(lineage == input$lineage_antigenic) %$% 
        mutation
      
      antigenic_mutations <- antigenic_mutations_lineages(nation = input$nation_antigenic, lineage = input$lineage_antigenic, defining = defining)
      
      if(is.null(antigenic_mutations)){
        values$antigenic <- NULL
        values$antigenic_title <- str_c(input$lineage_antigenic, " (", input$nation_antigenic %>% str_replace_all("_", " "), ")", ": no antigenic mutations")
      } else {
        values$antigenic_title <- str_c("Antigenic mutations on the top of ", input$lineage_antigenic, " defining mutations (", 
                                        input$nation_antigenic %>% str_replace_all("_", " "), ")")
        values$antigenic <- 
          antibody_complex_heatmap(antigenic_mutations)
      }
    })
    
    output$title_heatmap <- renderText({
      values$antigenic_title
    })
    
    # Display antibody heatmap
    # TODO plot caching
    output$antibody_heatmap <- renderPlot({
      if(is.null(values$antigenic)){
        ""  
      } else {
        draw(values$antigenic) 
      }
    }, height = function(){
      if(is.null(values$antigenic)){
        1
      } else {
        values$antigenic %>% nrow * 13 + 100
      }
    })
    
    observeEvent(input$gene, {
      updateSelectInput(session, "position",
                        choices = mutations_uk %>%
                          filter(gene == input$gene) %>%
                          distinct(position) %>%
                          arrange(position))
    })
    
    observeEvent(input$dataset_gene, {
      updateSelectizeInput(session, 
                           "dataset",
                           choices = database_genome %>% 
                             filter(gene == input$dataset_gene) %>% 
                             distinct(mutation))
    })
    
    
    observeEvent(input$dataset,{
      toggleState("downloadData", input$dataset != "")
    })
    
    output$epitope_sequence <- renderPlot({
      if(!input$epitope_ref){
        wt %<>% mutate(`numSeqs UK` = 0)
      }
      
      if(input$method){
        method = "prob"
      } else {
        method = "bits"
      }
      
      database_logo <- 
        database %>% 
        mutate(WT = str_sub(mutation, 1, 1), AA = str_sub(mutation, -1), .after = mutation) %>% 
        filter(`numSeqs UK` > 0) %>% 
        select(position, AA, `numSeqs UK`) %>% 
        bind_rows(wt) %>%
        pivot_wider(names_from = position, values_from = `numSeqs UK`, values_fill = 0, names_sort = TRUE) %>%
        arrange(AA) %>% 
        column_to_rownames(var = "AA") %>% 
        as.matrix
      
      epitopes_positions <- 
        database_tcell_predictions %>% 
        filter(`Start position` <= input$epitope_position & `End position` >= input$epitope_position) %>% 
        distinct(Epitope, `Start position`, `End position`)
      
      if(nrow(epitopes_positions) > 0){
        epitopes_positions %>% 
          pmap(function(Epitope, `Start position`, `End position`) {
            database_logo[,`Start position`:`End position`]
          }) %>% purrr::set_names(epitopes_positions %$% str_c(Epitope, " ", `Start position`, ":", `End position`)) %>% 
        ggseqlogo(method = method, ncol = 2) # TODO separate reactive for method to update plot only
      } else {
        ""
      }
    })
    
    variant_plot <- reactive({
      lineages_days_uk_all
      selected_variants <- input$variant_vui_voc
      vui_voc_lineages <- 
        levels(vui_voc$lineage) %>% 
        append("Other Delta", after = 9) %>% 
        append("Other")
      
      if( "B.1.617.2" %in% input$variant_vui_voc && input$variant_delta != "B.1.617.2"){
        selected_variants <- replace(selected_variants, selected_variants == "B.1.617.2", input$variant_delta)
        vui_voc_lineages <- replace(vui_voc_lineages, vui_voc_lineages == "B.1.617.2", input$variant_delta)
      }
      
      if(input$variant_day){
        
        if( "B.1.617.2" %in% input$variant_vui_voc){
          # remove unselected Deltas so they aren't counted in Other
          delta_options <- c("B.1.617.2", "AY.4", "AY.4.2")
          delta_options <- delta_options[delta_options != input$variant_delta]
          
          lineages_days_uk_all %<>% 
            filter(!(lineage %in% delta_options)) 
          
          # filter out Delta minus counts
          if(input$variant_delta == "AY.4"){
            lineages_days_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4.2")
            
            # add Other Delta after AY.4 in colour ordering
            ay_index <- match("AY.4", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4", after = ay_index) # don't include in Other
          } else if(input$variant_delta == "AY.4.2"){
            lineages_days_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4")
            
            # add Other Delta after AY.4.2 in colour ordering
            ay_index <- match("AY.4.2", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4.2", after = ay_index) # don't include in Other
          } else {
            lineages_days_uk_all %<>% 
              filter(!str_starts(lineage, fixed("Delta_minus_")))
          }
        } else { # Delta not selected, count only B.1.617.2 in Other
          lineages_days_uk_all %<>% 
            filter(!str_starts(lineage, fixed("AY."))) %>% 
            filter(!str_starts(lineage, fixed("Delta_minus_")))
        } 
        
        variants_other_day <- 
          lineages_days_uk_all %>% 
          filter(!(lineage %in% selected_variants)) %>% 
          group_by(sample_date, adm1) %>% 
          summarise(n_day = sum(n_day)) %>% 
          mutate(lineage = "Other", .before = sample_date) %>% 
          ungroup
        
        lineages_days_uk <- 
          lineages_days_uk_all %>% 
          filter(lineage %in% selected_variants) %>% 
          bind_rows(variants_other_day) %>% 
          mutate(lineage = recode(lineage,
                                  "Delta_minus_AY.4.2" = "Other Delta",
                                  "Delta_minus_AY.4" = "Other Delta")) %>%
          mutate(lineage = recode_factor(lineage, # recode WHO Greek display names as factor and order levels to define colour/legend order
                                  "AV.1" = "AV.1", 
                                  "B.1.1.318" = "B.1.1.318",
                                  "B.1.1.7" = "B.1.1.7 (Alpha)",
                                  "B.1.351" = "B.1.351 (Beta)",
                                  "B.1.525" = "B.1.525 (Eta)",
                                  "B.1.617.1" = "B.1.617.1 (Kappa)",
                                  "B.1.617.2" = "B.1.617.2/AY.x (Delta)",
                                  "Other Delta" = "Other Delta",
                                  "AY.4" = "AY.4/AY4.x (Delta)",
                                  "AY.4.2" = "AY.4.2 (Delta)",
                                  "B.1.617.3" = "B.1.617.3",
                                  "P.1" = "P.1 (Gamma)",
                                  "P.2" = "P.2 (Zeta)",
                                  "P.3" = "P.3 (Theta)",
                                  "Other" = "Other"
                                  )) %>% 
          rename(Variant = lineage, `Sample date` = sample_date, Sequences = n_day, adm1 = adm1)
        
        selected_variants <- replace(selected_variants, selected_variants %in% c("Delta_minus_AY.4.2", "Delta_minus_AY.4"), "Other Delta") 
        selected_variants <- replace(selected_variants, selected_variants == "B.1.617.2", input$variant_delta)
        
        if (input$nations_vui_voc!="UK"){
          lineages_days_uk %<>% 
            filter(adm1 == input$nations_vui_voc)}
       else
         lineages_days_uk %<>% 
          
            
        vui_voc_plot <- 
          lineages_days_uk %>% 
                   filter(`Sample date` >= input$variant_range[1] & `Sample date` <= input$variant_range[2]) %>% 
          ggplot(aes(fill = Variant, y = Sequences, x = `Sample date`) ) +
          theme_classic() +
          
          scale_fill_discrete_qualitative(palette = "Dynamic", 
                                          nmax = vui_voc_lineages %>% length , # extra colour for Other
                                          
                                          # fix variant/colour combos plus extra colour for Other
                                          order = match(selected_variants, 
                                                        vui_voc_lineages) %>% c(vui_voc_lineages %>% length)
                                          ) +
          
          scale_x_date(breaks = date_breaks("1 month"),
                       labels = date_format("%b %y")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          labs(x = "Sample date",
               y = "Sequences"
          ) 
        
        # display annotation from Sunday on 2nd last epiweek for consistency with weeks plot
        max_epi_date <- lineages_weeks_uk_all %$% max(epi_date) 
        
        # Set height of y-axis and annotation box according to highest weekly total sequences 
        if(input$variant_percentage){
          ymax <- 1
        } else {
          ymax <- 
            ymax <- 
            lineages_days_uk %>% 
            group_by(`Sample date`) %>% 
            summarise(total_week = sum(Sequences)) %$% 
            max(total_week)
        }
        
        # display annotation if upper slider is in latest 2 weeks
        if (input$variant_range[2] >= max_epi_date - days(7)){
          vui_voc_plot <- 
            vui_voc_plot +
            annotate("rect", 
                     xmax = input$variant_range[2] + days(1), 
                     xmin = max_epi_date - days(8), 
                     ymin = 0, 
                     ymax = ymax,
                     alpha = 0.2)
        }
        
        if(input$variant_percentage){
          vui_voc_plot <- 
            vui_voc_plot +
            geom_bar(position="fill", stat="identity", width = 1) +
            scale_y_continuous(labels = scales::percent_format())
        } else {
          vui_voc_plot <- 
            vui_voc_plot + 
            geom_bar(position="stack", stat="identity", width = 1) +
            ylim(0, ymax)
        }
        
        vui_voc_plot
      } else {  # variants by week
        
        if( "B.1.617.2" %in% input$variant_vui_voc){
          # remove unselected Deltas so they aren't counted in Other
          delta_options <- c("B.1.617.2", "AY.4", "AY.4.2")
          delta_options <- delta_options[delta_options != input$variant_delta]
          
          lineages_weeks_uk_all %<>%
            filter(!(lineage %in% delta_options)) 
          
          # filter out Delta minus counts
          if(input$variant_delta == "AY.4"){
            lineages_weeks_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4.2")
            
            # add Other Delta after AY.4 in colour ordering
            ay_index <- match("AY.4", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4", after = ay_index) # don't include in Other
          } else if(input$variant_delta == "AY.4.2"){
            lineages_weeks_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4")
            
            # add Other Delta after AY.4.2 in colour ordering
            ay_index <- match("AY.4.2", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4.2", after = ay_index) # don't include in Other
          } else {
            lineages_weeks_uk_all %<>% 
              filter(!str_starts(lineage, fixed("Delta_minus_")))
          }
        } else { # Delta not selected, count only B.1.617.2 in Other
          lineages_weeks_uk_all %<>% 
            filter(!str_starts(lineage, fixed("AY."))) %>% 
            filter(!str_starts(lineage, fixed("Delta_minus_")))
        } 
        
        variants_other_week <- 
          lineages_weeks_uk_all %>% 
          filter(!(lineage %in% selected_variants)) %>% 
          group_by(epi_date, adm1) %>% 
          summarise(n_week = sum(n_week)) %>% 
          mutate(lineage = "Other", .before = epi_date) %>% 
          ungroup
        
        
        lineages_weeks_uk <- 
          lineages_weeks_uk_all %>% 
          filter(lineage %in% selected_variants) %>% 
          bind_rows(variants_other_week) %>% 
          mutate(lineage = recode(lineage,
                                  "Delta_minus_AY.4.2" = "Other Delta",
                                  "Delta_minus_AY.4" = "Other Delta")) %>%
          mutate(lineage = recode_factor(lineage, # recode WHO Greek display names as factor and order levels to define colour/legend order
                                         "AV.1" = "AV.1", 
                                         "B.1.1.318" = "B.1.1.318",
                                         "B.1.1.7" = "B.1.1.7 (Alpha)",
                                         "B.1.351" = "B.1.351 (Beta)",
                                         "B.1.525" = "B.1.525 (Eta)",
                                         "B.1.617.1" = "B.1.617.1 (Kappa)",
                                         "B.1.617.2" = "B.1.617.2/AY.x (Delta)",
                                         "Other Delta" = "Other Delta",
                                         "AY.4" = "AY.4/AY4.x (Delta)",
                                         "AY.4.2" = "AY.4.2 (Delta)",
                                         "B.1.617.3" = "B.1.617.3",
                                         "P.1" = "P.1 (Gamma)",
                                         "P.2" = "P.2 (Zeta)",
                                         "P.3" = "P.3 (Theta)",
                                         "Other" = "Other"
          )) %>% 
          rename(Variant = lineage, `Start date` = epi_date, adm1 = adm1, Sequences = n_week)
        
        selected_variants <- replace(selected_variants, selected_variants %in% c("Delta_minus_AY.4.2", "Delta_minus_AY.4"), "Other Delta") 
        selected_variants <- replace(selected_variants, selected_variants == "B.1.617.2", input$variant_delta)
        
        if (input$nations_vui_voc!="UK"){
          lineages_weeks_uk %<>% 
          filter(adm1 == input$nations_vui_voc)}
        
        vui_voc_plot <- 
          lineages_weeks_uk %>% 
         
          filter(`Start date` >= input$variant_range[1] & `Start date` <= input$variant_range[2]) %>% 
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
          theme(plot.title = element_text(hjust = 0.5)) +
          labs(x = "Sample date",
               y = "Sequences"
          ) 
        
          

        max_date <- lineages_weeks_uk %$% max(`Start date`)
        
        # Set height of annotation box according to highest weekly total sequences
        if(input$variant_percentage){
          ymax <- 1
        } else {
          ymax <- 
            ymax <- 
            lineages_weeks_uk %>% 
            group_by(`Start date`) %>% 
            summarise(total_week = sum(Sequences)) %$% 
            max(total_week)
        }
        
        # display annotation if upper slider is in latest 2 weeks
        if (input$variant_range[2] >= max_date - days(7)){
          if(input$variant_range[2] >= max_date ){
            xmax = max_date + days(3) 
          } else{
            xmax = max_date + days(-4) 
          }
          
          vui_voc_plot <- 
            vui_voc_plot +
            annotate("rect", 
                     xmax = xmax, 
                     xmin = max_date - days(10), 
                     ymin = 0, 
                     ymax = ymax,
                     alpha = 0.2)
        }
        
        # TODO position_nudge
        if(input$variant_percentage){
          vui_voc_plot <- 
            vui_voc_plot +
            geom_bar(position="fill", stat="identity") +
            scale_y_continuous(labels = scales::percent_format())
        } else {
          vui_voc_plot <- 
            vui_voc_plot + 
            geom_bar(position="stack", stat="identity") +
            ylim(0, ymax)
        }
        vui_voc_plot
      } # end else by day
      vui_voc_plot
    }) %>% debounce(500)
    
    output$variant_time <- renderPlotly({
      variant_plot() %>% ggplotly
    })
    
    ########### Ronapreve plot
    # always display wild type on percentage chart
    observeEvent(input$ronapreve_28, {
      if(input$ronapreve_28){
        output$title_ronapreve <- renderText("28 days to latest UK sequence date")
        output$ronapreve_plot <- renderImage({
          list(src = str_c(dataset_date, "/Ronapreve_28.png"),
               alt = "Ronapreve plot 28 days"
          )
        }, deleteFile = FALSE)
      } else {
        output$title_ronapreve <- renderText("All time")
        output$ronapreve_plot <- renderImage({
          list(src = str_c(dataset_date, "/Ronapreve.png"),
               alt = "Ronapreve plot all time"
          )
        }, deleteFile = FALSE)
      }
    })
    
    ########### Map and geographical distribution of variants
    #
    # map_weekInput <- reactive({ # find epiweek from date selected by user
    #   x<-input$variant_date
    #   map_week<- tibble(x) %>% rename("epi_date" = "x") %>%
    #     plyr::join(epi_lookup, by = "epi_date") %>% select(epi_week)
    #   map_week$epi_week
    # })
    
      
    geo_all_1<- geo_all %>% filter(epi_week == "97" & lineage == "B.1.617.2")
    
    geo_all_1<-dplyr::rename(geo_all_1, "value" = "Count")
    geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "objectid")]
    geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "epi_week")]
    geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "lineage")]
    geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "Proportion")]
   
    # mapdata_1<- dplyr::rename(mapdata_1, "NUTS1" = "id")
    #Join mydata with mapdata
    df <- plyr::join(mapdata_1, geo_all_1, by= c("NUTS1"))
    # df<- df[, -which(names(df) == "objectid")]
    
    
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(magrittr))
    
    gg <- ggplot() + geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = value), color = "#FFFFFF", size = 0.25)
    gg <- gg + scale_fill_gradient2(low = "blue", mid = "red", high = "yellow", na.value = "white")
    gg <- gg + coord_fixed(1)
    gg <- gg + theme_minimal()
    gg <- gg +  scale_fill_viridis(trans = "log", breaks=c(1,5,10,20,50,100), name="Number of Sequences", guide = guide_legend( keyheight = unit(3, units = "mm"), keywidth=unit(12, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) )
    gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'right')
    gg <- gg + theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    gg <- gg + theme(axis.title.y=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    # }) 
    
    output$map <- renderPlot({

      gg

    })
    
    
})