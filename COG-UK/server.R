library(shiny)
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(scales)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(ComplexHeatmap)

source("helpers.R")

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

lineage_plus_variant <- function(lineage, variant, use_regex = FALSE){
  mutations_s_uk_lv <- 
    mutations_s_uk %>%
    filter(variant == !!variant) %>%
    when(
      use_regex ~filter(., lineage == !!lineage | str_detect(lineage, sublineage_regex(!!lineage))),
      ~filter(., lineage == !!lineage)
    ) 
  
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
    mutate(adm1 = recode(adm1, 
                         `UK-ENG` = "England",
                         `UK-NIR` = "Northern_Ireland",
                         `UK-SCT` = "Scotland",
                         `UK-WLS` = "Wales")) %>% 
    pivot_wider(names_from = adm1, values_from = c(n_sequences, n_sequences_28)) %>%
    mutate(lineage = !!lineage, variant = !!variant, .before = 1)
}

# TODO caching
table_1 <- function(){
  database %>% 
    arrange(desc(`numSeqs UK`)) %>% 
    filter(`numSeqs UK` >= 5) %>% 
    mutate(mutation = mutation %>% fct_drop %>% fct_inorder) %>% 
    select(mutation, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`, earliest) %>% 
    mutate(`Cumulative sequences in UK (%)` = `numSeqs UK` / total_sequences,
           .after = `numSeqs UK`) %>%
    mutate(`Sequences over the last 28 days in UK (%)` = `numSeqs UK 28 days` / total_sequences_28,
           .after = `numSeqs UK 28 days`) %>%
    rename(`Amino acid replacement` = mutation, 
           `Cumulative sequences in UK` = `numSeqs UK`, 
           `Sequences over the last 28 days in UK` = `numSeqs UK 28 days`,
           `Sequences over the last 28 days in England` = `numSeqs Eng 28 days`,
           `Sequences over the last 28 days in Scotland` = `numSeqs Scotland 28 days`,
           `Sequences over the last 28 days in Wales` = `numSeqs Wales 28 days`,
           `Sequences over the last 28 days in Northern Ireland` = `numSeqs NI 28 days`,
           `Date of first appearance in UK` = earliest) 
}

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
    
    lineage_plus_variant("B.1.1.7", "S494P") %>% # for non-key mutations
      mutate(lineage = str_c(lineage, " + ", variant), .keep = "unused") %>%
      mutate(reason = "As B.1.1.7, with the addition of S494P."),
    
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
    arrange(lineage) %>% 
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

table_5 = function(){
  database_tcell_predictions %>% 
    select(mutation, Epitope:Fold, `numSeqs UK`, `numSeqs UK 28 days`, -`End position`) %>%
    arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), desc(Fold)) %>%
    mutate(across(c(Epitope, HLA, assay, CD4_CD8), ~fct_relevel(.x, sort))) %>% 
    rename(`Amino acid replacement` = mutation,
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

# TODO precompute and include lineage/variant combinations
n_uk_lineages_all <-
  left_join(
    sum_key_mutations_by_lineage_uk(lineages_t2),
    sum_key_mutations_by_lineage_uk(lineages_t2, date_from = sample_date_28) %>%
      rename(n_sequences_28 = n_sequences)
  ) %>% 
  pivot_wider(names_from = adm1, values_from = c(n_sequences, n_sequences_28)) %>% 
  mutate(across(everything(), ~replace_na(.x, 0L)))

shinyServer(function(input, output, session) {

    output$table_1 <- renderDataTable({
      table_1() %>% 
        datatable(filter = "top", rownames = FALSE, 
                  options = list(lengthMenu = c(20, 50, 100, 200), pageLength = 20, scrollX = TRUE)) %>% 
        formatPercentage(c("Cumulative sequences in UK (%)", "Sequences over the last 28 days in UK (%)"), digits = 2)
    })
    
    ######## Data download inputs ########
    ## Table 1
    
    # Reactive value to generate downloadable table for selected table 1 mutation metadata
    datasetInput <- reactive({
      mutations_s_uk %>% 
        filter(variant == input$dataset) %>% 
        filter(sample_date >= sample_date_28) %>% 
        select(sequence_name, sample_date, epi_week, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete table 1 data
    table1Input <- reactive({
      table_1() %>% 
        mutate(across(ends_with("(%)"), ~.x * 100)) # convert decimal fraction to percentage
    })
    
    ## Table 3
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
        select(sequence_name, sample_date, epi_week, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete table 3 data
    table3Input <- reactive({
      table_3() %>% 
        mutate(across(ends_with("(%)"), ~.x * 100)) # convert decimal fraction to percentage
    })
    
    # Reactive value to generate downloadable table for selected mutation
    escapeInput <- reactive({
      mutations_s_uk %>% 
        filter(variant == input$selectEscape) %>% 
        select(sequence_name, sample_date, epi_week, lineage) %>% 
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
        str_c(input$dataset, "_UK_28_days_", dataset_date, ".csv")
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
    
    
    output$table_2 <- renderDT({
      table_2() %>% 
          datatable(filter = "none", rownames = FALSE, 
                    options = list(dom = 't', paging = FALSE, scrollX = TRUE)) %>% 
          formatPercentage(c("UK (%)", "UK 28 days (%)"), digits = 2)
    })
    
    output$table_3 <- renderDT({
      table_3() %>% 
        datatable(filter = "none", rownames = FALSE, 
                  options = list(dom = 't', paging = FALSE, scrollX = TRUE)) %>% 
        formatPercentage(c("UK (%)", "UK 28 days (%)"), digits = 2)
    })
    
    output$table_4 <- renderDT({
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
      table_5() %>% 
        mutate(Reference = str_c("<a href='", DOI, "'target='_blank'>", Reference,"</a>"), .keep = "unused", .after = `Supporting references`) %>% # hyperlink to citation DOI
        datatable(filter = "top", escape = FALSE, rownames = FALSE,
                    options = list(lengthMenu = c(10, 20, 50, 100, 200), pageLength = 10, scrollX = TRUE))
        
    })
    
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
      if(!input$ref){
        mutation_reference_counts %<>% filter(variant != "WT") 
      }
      
      mutation_reference_counts %>% 
        filter(gene == input$gene & position == input$position) %>%
        mutate(variant = variant %>% fct_infreq) %>% # fix variant colour by frequency
        filter(adm1 == input$nation) %>%
        filter(epi_week %in% c(input$epi_week[1]:input$epi_week[2])) %>% # match because epi_week is factor
        mutate(epi_week = fct_drop(epi_week, only = {.} %$% 
                                                    levels(epi_week) %>% 
                                                    as.numeric %>% 
                                                    keep(~ .x < input$epi_week[1] | .x > input$epi_week[2]) %>% 
                                                    as.character)) %>% # drop filtered epi_weeks to exclude from x-axis
        ggplot(aes(fill=variant, y=n, x=epi_week) ) +
        scale_x_discrete(drop=FALSE) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(x = "Epidemic Week", 
             y = "Sequences", 
             fill = "Variant",
             title = str_c(c("Gene", "Position"), c(input$gene, input$position), sep = " : ", collapse = "\n")) +
        scale_fill_manual(values = brewer.pal(name = "Set2", n = 8)) 
    }) 
    
    mutation_plot_bar <- reactive({
      if(input$percentage){
        gg_bar <- 
          mutation_plot() +
          geom_bar(position="fill", stat="identity") +
          scale_y_continuous(labels = scales::percent_format())
      } else {
        gg_bar <- 
          mutation_plot() +
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
    
    observeEvent(input$dataset,{
      toggleState("downloadData", input$dataset != "")
    })
})