library(shiny)
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(scales)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(ComplexHeatmap)
library(circlize)

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
      mutate(reason = "As B.1.1.7, with the addition of E484K, which is located in the RBM and has been shown to escape some mAbs."),
    
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
           `Reason for tracking` = reason, 
           
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
    concernInput <- reactive({
      if(input$concern == "B.1.1.7 + E484K"){
        concern_download <- 
          consortium_uk %>% 
          filter(e484k == "K") %>% 
          filter(lineage == "B.1.1.7" | str_detect(lineage, sublineage_regex("B.1.1.7"))) 
      } else if(input$concern == "A.23.1 + E484K") {
        concern_download <- 
          consortium_uk %>% 
          filter(e484k == "K") %>% 
          filter(lineage == "A.23.1" | str_detect(lineage, sublineage_regex("B.1.1.7"))) 
      } else if(input$concern == "B.1.1.7 + S494P") {
        concern_download <- 
          mutations_s_uk %>% 
          filter(variant == "S494P") %>% 
          filter(lineage == "B.1.1.7" | str_detect(lineage, sublineage_regex("B.1.1.7")))
      } else {
        concern_download <- 
          consortium_uk %>% 
          filter(lineage == input$concern | str_detect(lineage, sublineage_regex(input$concern)))
      }
      
      concern_download %>% 
        select(sequence_name, sample_date, epi_week, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete table 1 data
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
    
    ## Antibody sites
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
      database %>% 
        filter(!is.na(Epitopes) & `numSeqs UK` > 0) %>% 
        mutate(mutation = fct_drop(mutation)) %>%
        select(mutation, Epitopes:Assays, `numSeqs UK`, `numSeqs UK 28 days`) %>% 
        arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), mutation) %>% 
        rename(`Amino acid replacement` = mutation, 
               `Cumulative sequences in UK` = `numSeqs UK`,
               `Sequences over 28 days` = `numSeqs UK 28 days`) %>% 
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

    # Display antibody heatmap
    output$antibody_heatmap <- renderPlot({
      
      horz_heat <-
        antigenic_mutations_lineages %>% 
        filter(lineage == "B.1.1.7" & variant != "N501Y") %>% 
        select(-lineage) %>% 
        inner_join(database %>% 
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
        Effect_mab = horz_heat$mab,
        Effect_plasma = horz_heat$plasma,
        Effect_vaccine = horz_heat$vaccine_sera,
        confidence = horz_heat$confidence,
        na_col = 'white',
        col = list(
          confidence = c(
            "lower" = "lightgoldenrod",
            "medium" = "lightgoldenrod3",
            "high" = "lightgoldenrod4"
          ),
          Effect_mab = c("TRUE" = "black"),
          Effect_plasma = c("TRUE" = "black"),
          Effect_vaccine = c("TRUE" = "black")
        ),
        annotation_legend_param = list(confidence = list (at = c(
          "high", "medium", "lower"
        )))
      )
      
      # domain
      row_ha2 = rowAnnotation(domain = (horz_heat$domain),
                              col = list(
                                domain = c(
                                  "FP" = "seashell2",
                                  "NTD" = "navajowhite",
                                  "RBD" = "pink",
                                  "RBM" = "plum1",
                                  "SP" = "lightblue1"
                                )
                              ))
      
      Heatmap(
        subset(input, select = -c(mab:domain)),
        name = "Percentage %",
        column_title = "Antigenic mutations in lineage B.1.1.7",
        column_title_gp = gpar(fontsize = 18),
                use_raster = TRUE,
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