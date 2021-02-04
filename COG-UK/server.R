library(shiny)
library(tidyverse)
library(formattable)
library(plotly)
library(RColorBrewer)
library(scales)
library(shinyWidgets)
library(shinyjs)

lineages_t2 <- c("B.1", "B.1.177", "B.1.141", "B.1.258", "B.1.1", "B.1.1.7", "B.1.1.70", "B.1.351", "B.1.1.298")

lineages_t3 <- 
    c("B.1.1.7" = "UK associated variant. Has 17 mutations (14 replacements and 3 deletions) including: T1001I, A1708D, I2230T, SGF 3675-3677 del In the ORF1ab; 69-70 del, Y144 del, N501Y, A570D, P681H, T716I, S982A and D1118H in the Spike; Q27stop, R52I and Y73C in ORF8; D3L and S235F in the N. Noteworthily, N501Y enhances ACE2 binding affinity, and P681H occurs at the furin cleavage site, known for biological significance in membrane fusion.", 
      "B.1.351" = "Variant associated with South Africa. Has eight mutations in the Spike: D80A, D215G, E484K, N501Y, A701V, L18F, R246I and K417N. Three of these in the RBM, K417N, E484K and N501Y. K417N and E484K have been shown to escape some mAbs.", 
      "P.1" = "Variant associated with Brazil. Has 10 mutations in the Spike including L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y,H655Y and T1027I. Noteworthy  E484K, N501Y and K417T have biological significance.") %>% 
    enframe("lineage", "reason")

shinyServer(function(input, output, session) {

    output$table_1 <- renderTable(
        database %>% 
            slice_max(`numSeqs UK`, n = 20) %>% 
            select(mutation, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`) %>% 
            mutate(`Sequences over the last 28 days in UK (%)` = percent(`numSeqs UK 28 days` / `numSeqs UK`) %>% as.character, .after = `numSeqs UK`) %>% 
            rename(`Amino acid replacement` = mutation, 
                   `Cumulative sequences in UK` = `numSeqs UK`, 
                   `Sequences over 28 days` = `numSeqs UK 28 days`,
                   `Sequences over the last 28 days in England` = `numSeqs Eng 28 days`,
                   `Sequences over the last 28 days in Scotland` = `numSeqs Scotland 28 days`,
                   `Sequences over the last 28 days in Wales` = `numSeqs Wales 28 days`,
                   `Sequences over the last 28 days in Northern Ireland` = `numSeqs NI 28 days`)
    )
    
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
    
    output$table_2 <- renderTable({
        n_uk_lineages <- sum_key_mutations_by_lineage_uk(lineages_t2)  
        
        n_uk_lineages_28 <- 
            sum_key_mutations_by_lineage_uk(lineages_t2, date_from = sample_date_28)  %>% 
            rename(n_sequences_28 = n_sequences)
        
        table_2 <- 
            inner_join(n_uk_lineages, n_uk_lineages_28) %>% 
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
                       (variant == "Y453F" & lineage == "B.1.1.298" ))
        
        lineages_E484K <- c("B.1.351", "P.2", "P.1", "B.1.1.7")
        table_2 %<>% bind_rows( 
            inner_join(
                lapply(lineages_E484K, function(x){
                    mutations_s_uk %>% 
                        filter(variant == "E484K") %>% 
                        filter(lineage == x | str_detect(lineage, sublineage_regex(x))) %>% 
                        summarise(!!x := n_distinct(sequence_name))
                }) %>% 
                    bind_cols() %>% 
                    gather(key = "lineage", value = "n_sequences"),
                
                lapply(lineages_E484K, function(x){
                    mutations_s_uk %>% 
                        filter(variant == "E484K") %>% 
                        filter(sample_date >= sample_date_28) %>% 
                        filter(lineage == x | str_detect(lineage, sublineage_regex(x))) %>% 
                        summarise(!!x := n_distinct(sequence_name))
                }) %>% 
                    bind_cols() %>% 
                    gather(key = "lineage", value = "n_sequences_28")
            ) %>% mutate(variant = "E484K") 
        )
        
        n_N501Y_E484K <-
            intersect(
                mutations_s_uk %>% 
                    filter(variant == "E484K") %>% 
                    filter(lineage == "B.1.351" | str_detect(lineage, sublineage_regex("B.1.351"))) %>% 
                    select(sequence_name),
                
                consortium_uk %>% 
                    filter(n501y == "Y") %>% 
                    filter(lineage == "B.1.351" | str_detect(lineage, sublineage_regex("B.1.351"))) %>% 
                    select(sequence_name)
            ) %>% nrow
        
        n_N501Y_E484K_28 <- 
         intersect(
            mutations_s_uk %>% 
                filter(variant == "E484K") %>% 
                filter(lineage == "B.1.351" | str_detect(lineage, sublineage_regex("B.1.351"))) %>% 
                filter(sample_date >= sample_date_28) %>% 
                select(sequence_name),
            
            consortium_uk %>% 
                filter(n501y == "Y") %>% 
                filter(lineage == "B.1.351" | str_detect(lineage, sublineage_regex("B.1.351"))) %>% 
                filter(sample_date >= sample_date_28) %>% 
                select(sequence_name)
        ) %>% nrow
        
        table_2 %<>% 
            add_row(lineage = "B.1.351", variant = 'N501Y + E484K', 
                             n_sequences = n_N501Y_E484K, n_sequences_28 = n_N501Y_E484K_28) %>% 
            relocate(variant) %>% 
            arrange(variant, lineage) %>% 
            rename(Mutation = variant, 
                   `Lineage(s) in which it has been detected` = lineage, 
                   `Cumulative sequences in UK` = n_sequences, 
                   `Sequences over 28 days` = n_sequences_28)
        table_2
    })
    
    output$table_3 <- renderTable({
        key_lineages <- 
        lapply(lineages_t3$lineage, function(x){
            sum_key_mutations_uk() %>% 
                filter(lineage == x | str_detect(lineage, sublineage_regex(x)) ) %>% 
                select(-lineage) %>%
                summarise(sequences_sum = sum(sequences)) %>% 
                mutate(lineage = x, .before = 1)
        }) %>% bind_rows()
        
        key_lineages_28 <- 
            lapply(lineages_t3$lineage, function(x){
                sum_key_mutations_uk(date_from = sample_date_28) %>% 
                    filter(lineage == x | str_detect(lineage, sublineage_regex(x)) ) %>% 
                    select(-lineage) %>%
                    summarise(sequences_sum = sum(sequences)) %>% 
                    mutate(lineage = x, .before = 1) %>% 
                    rename(sequences_sum_28 = sequences_sum)
            }) %>% bind_rows() 
        
        inner_join(key_lineages, key_lineages_28) %>% 
            inner_join(lineages_t3) %>% 
            relocate(reason, .after = lineage) %>% 
            rename(`Variant/ lineage` = lineage,	
                   `Cumulative sequences in UK` = sequences_sum,	 
                   `Sequences over 28 days` = sequences_sum_28,                    
                   `Reason for tracking` = reason)
    })
    output$table_4 <- renderDataTable({
        database %>%
            filter(!is.na(escape)) %>% 
            select(mutation, escape, `numSeqs UK`, `numSeqs UK 28 days`, citation, doi) %>%  
            arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), mutation) %>% 
            rename(`Amino acid replacement` = mutation, 
                   `Cumulative sequences in UK` = `numSeqs UK`,
                   `Sequences over 28 days` = `numSeqs UK 28 days`,
                   `Escape mutations details` = escape,
                   `References` = citation,
                   DOI = doi)
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
        filter(adm1 == input$nation) %>%
        filter(gene == input$gene & position == input$position) %>%
        filter(epi_week %in% c(input$epi_week[1]:input$epi_week[2])) %>% # match because epi_week is factor
        mutate(epi_week = fct_drop(epi_week)) %>% # drop filtered epi_weeks to exclude from x-axis
        select(-position, -gene) %>%
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