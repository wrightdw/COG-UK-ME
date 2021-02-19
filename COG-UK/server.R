library(shiny)
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(scales)
library(formattable)
library(shinyWidgets)
library(shinyjs)
library(DT)

shinyServer(function(input, output, session) {

    output$table_1 <- renderDataTable({
      database %>% 
        arrange(desc(`numSeqs UK`)) %>% 
        filter(`numSeqs UK` >= 5) %>% 
        select(mutation, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`, earliest) %>% 
        mutate(`Cumulative sequences in UK (%)` = formattable::percent(`numSeqs UK` / total_sequences, digits = 1L) %>% as.character,
               .after = `numSeqs UK`) %>% 
        mutate(`Sequences over the last 28 days in UK (%)` = formattable::percent(`numSeqs UK 28 days` / total_sequences_28, digits = 1L) %>% as.character,
               .after = `numSeqs UK 28 days`) %>% 
        rename(`Amino acid replacement` = mutation, 
               `Cumulative sequences in UK` = `numSeqs UK`, 
               `Sequences over the last 28 days in UK` = `numSeqs UK 28 days`,
               `Sequences over the last 28 days in England` = `numSeqs Eng 28 days`,
               `Sequences over the last 28 days in Scotland` = `numSeqs Scotland 28 days`,
               `Sequences over the last 28 days in Wales` = `numSeqs Wales 28 days`,
               `Sequences over the last 28 days in Northern Ireland` = `numSeqs NI 28 days`,
               `Date of first appearance in UK` = earliest)
    }, options = list(lengthMenu = c(20, 50, 100, 200), pageLength = 20))
    
    # Reactive value to generate downloadable table for selected mutation
    datasetInput <- reactive({
      mutations_s_uk %>% 
        filter(variant == input$dataset) %>% 
        filter(sample_date >= sample_date_28) %>% 
        select(sequence_name, sample_date, epi_week, lineage, uk_lineage, phylotype) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for selected mutation
    escapeInput <- reactive({
      mutations_s_uk %>% 
        filter(variant == input$escape) %>% 
        select(sequence_name, sample_date, epi_week, lineage, uk_lineage, phylotype) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for selected mutation
    concernInput <- reactive({
      
      if(input$concern == "B.1.1.7 + E484K"){
        consortium_uk %<>% 
          filter(lineage == "B.1.1.7" | str_detect(lineage, sublineage_regex("B.1.1.7"))) %>% 
          filter(e484k == "K") 
      } else {
        consortium_uk %<>% 
          filter(lineage == input$concern | str_detect(lineage, sublineage_regex(input$concern)))
      }
      
      consortium_uk %>% 
        select(sequence_name, sample_date, epi_week, lineage, uk_lineage, phylotype) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Downloadable CSV of selected mutation
    output$downloadData <- downloadHandler(
      filename = function() {
        str_c(input$dataset, "_UK_28_days_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(datasetInput(), file)
      },
      contentType = "text/csv"
    )
    
    # Downloadable CSV of selected mutation
    output$downloadEscape <- downloadHandler(
      filename = function() {
        str_c(input$escape, "_UK_cumulative_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(escapeInput(), file)
      },
      contentType = "text/csv"
    )
    
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
    
    output$table_2 <- renderTable({
        table_2 <- 
          n_uk_lineages_all %>% 
            filter(    (variant == "D614G" & lineage == "B.1" ) | 
                           
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
                       
                       (variant == "E484K" & lineage == "A" ) |
                       (variant == "E484K" & lineage == "B.1.1.119" ) |
                       (variant == "E484K" & lineage == "B.1.177.4" ) |
                       (variant == "E484K" & lineage == "B.1.222" ) |
                       (variant == "E484K" & lineage == "B.1.177" ) |
                       (variant == "E484K" & lineage == "B.1" ) |
                       
                       (variant == "N501Y + E484K" & lineage == "B.1.351") |
                       (variant == "N501Y + E484K" & lineage == "B.1.1.7"))
        
        table_2 %<>% 
            relocate(variant) %>% 
            arrange(variant, lineage) %>% 
            rename(Mutation = variant, 
                   `Lineage(s) in which detected` = lineage, 
                   `Cumulative UK sequences` = n_sequences, 
                   `UK Sequences over 28 days` = n_sequences_28)
        table_2
    })
    
    output$table_3 <- renderTable({
          bind_rows(
            n_uk_lineages_all %>% 
              filter(lineage %in% lineages_t3$lineage & variant == "sequences") %>% 
              inner_join(lineages_t3) %>% 
              select(-variant) %>% 
              relocate(reason, .after = lineage),
            
            n_uk_lineages_all %>% 
              filter(variant == "E484K" & lineage == "B.1.1.7") %>% 
              mutate(lineage = str_c(lineage, " + ", variant), .keep = "unused")  %>% 
              mutate(reason = "As above, with the addition of E484K, which is located in the RBM and has been shown to escape some mAbs.")
          ) %>% 
                arrange(lineage) %>% 
                rename(`Variant/ lineage` = lineage,	
                         `Cumulative sequences in UK` = n_sequences,	 
                         `Sequences over 28 days` = n_sequences_28,                    
                         `Reason for tracking` = reason)
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
                 # `Monoclonal Ab` = mab,
                 # `Convalescent sera` = plasma,
                 # `Vaccine sera` = vaccine_sera,
                 Confidence = support, 
                 Domain = domain) %>% 
        datatable(filter = "top", escape = FALSE) %>% 
        formatStyle(
          'Confidence',
          target = 'row',
          backgroundColor = styleEqual(c("lower", "medium", "high"), c('LemonChiffon', 'DarkOrange', 'FireBrick')), 
          color = styleEqual(c("lower", "medium", "high"), c('DarkSlateGray', 'White', 'Snow'))) # TODO anchor colour
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
        mutate(variant = variant %>% as_factor %>% fct_infreq) %>% # fix variant colour by frequency
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

    observeEvent(input$gene, {
      updateSelectInput(session, "position",
                        choices = mutations_uk %>%
                          filter(gene == input$gene) %>%
                          distinct(position) %>%
                          arrange(position))
    })
})