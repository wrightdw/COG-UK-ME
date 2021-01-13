library(shiny)
library(tidyverse)

sample_date_28 <- "2020-12-02"
lineages <- c("B.1", "B.1.177", "B.1.141", "B.1.258", "B.1.1", "B.1.1.7", "B.1.1.70", "B.1.351", "B.1.1.298")

shinyServer(function(input, output, session) {

    output$table_1 <- renderTable(
        database %>% 
            slice_max(`numSeqs UK`, n = 15) %>% 
            select(replacement, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`) %>% 
            rename(`Amino acid replacement` = replacement, 
                   `Cumulative sequences in UK` = `numSeqs UK`, 
                   `Sequences over 28 days` = `numSeqs UK 28 days`,
                   `Sequences over the last 28 days in England` = `numSeqs Eng 28 days`,
                   `Sequences over the last 28 days in Scotland` = `numSeqs Scotland 28 days`,
                   `Sequences over the last 28 days in Wales` = `numSeqs Wales 28 days`,
                   `Sequences over the last 28 days in Northern Ireland` = `numSeqs NI 28 days`)
    )
    
    output$table_2 <- renderTable({
        n_uk_lineages <- sum_key_mutations_by_lineage_uk(lineages)  
        
        n_uk_lineages_28 <- 
            sum_key_mutations_by_lineage_uk(lineages, date_from = sample_date_28)  %>% 
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
        
        n_N501Y_E484K <-
            intersect(
                mutations_s_uk %>% 
                    filter(variant == "E484K") %>% 
                    filter(lineage == "B.1.351" | str_detect(lineage, "B\\.1\\.1\\.351\\.")) %>% 
                    distinct(sequence_name),
                
                consortium_uk %>% 
                    filter(n501y == "Y") %>% 
                    filter(lineage == "B.1.351" | str_detect(lineage, "B\\.1\\.1\\.351\\.")) %>% 
                    distinct(sequence_name)
            ) %>% nrow
        
        n_N501Y_E484K_28 <- 
         intersect(
            mutations_s_uk %>% 
                filter(variant == "E484K") %>% 
                filter(lineage == "B.1.351" | str_detect(lineage, "B\\.1\\.1\\.351\\.")) %>% 
                filter(sample_date >= sample_date_28) %>% 
                distinct(sequence_name),
            
            consortium_uk %>% 
                filter(n501y == "Y") %>% 
                filter(lineage == "B.1.351" | str_detect(lineage, "B\\.1\\.1\\.351\\.")) %>% 
                filter(sample_date >= sample_date_28) %>% 
                distinct(sequence_name)
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
        sum_key_mutations_uk() %>% 
            filter(lineage %in% c("B.1.1.7", "B.1.351"))

    })
    output$table_4 <- renderTable({})
    
    # output$mutation_time <- renderPlot({
    #     
    #     mutations %>% filter(gene == input$gene 
    #                          & variant == input$variant
    #                          & sample_date >= input$date_range[1] 
    #                          & sample_date <= input$date_range[2]) %>% ggplot(aes(x = sample_date)) + geom_bar() + theme_minimal()
    # })
    # 
    # output$raw_mutations <- DT::renderDataTable(
    #     mutations %>% 
    #     filter(gene == input$gene 
    #            & position == input$position 
    #            & variant == input$variant 
    #            & sample_date >= input$date_range[1] 
    #            & sample_date <= input$date_range[2])
    # )
    
    
    # output$summary <- renderTable({
    #     mutations %>% 
    #         filter(gene == input$gene 
    #                & position == input$position 
    #                & variant == input$variant
    #                & sample_date >= input$date_range[1] 
    #                & sample_date <= input$date_range[2]) %>% count(country)
    # })
    
    # database %>% 
    #     slice_max(`numSeqs UK`, n = 15) %>% 
    #     select(replacement, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`)
    
    
# 
#     observeEvent(input$gene, {
#         updateSelectInput(session, "position",
#                           choices = mutations %>% filter(gene == input$gene) %>% distinct(position) %>% arrange(position))
#     })
#     
#     observeEvent(input$position, {
#         updateSelectInput(session, "variant",
#                           choices = mutations %>% filter(gene == input$gene & position == input$position) %>% distinct(variant) %>% arrange(variant))
#     })
    
})
