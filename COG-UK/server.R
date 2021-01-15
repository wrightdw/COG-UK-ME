library(shiny)
library(tidyverse)

lineages_t2 <- c("B.1", "B.1.177", "B.1.141", "B.1.258", "B.1.1", "B.1.1.7", "B.1.1.70", "B.1.351", "B.1.1.298")
lineages_t3 <- c("B.1.1.7", "B.1.351", "P.1")
escape_t4 <- c(
"G446V",
"L452R",
"E484Q",
"K444R",
"E484K",
"Y508H",
"N440K",
"L455F",
"A831V",
"A475V",
"F490S",
"V483A",
"R346K",
"K378N",
"K444N",
"G446S",
"N450D",
"K150R",
"R346S",
"K150T",
"V445A",
"G446A",
"Y449H",
"E484A",
"E484R",
"Q493R")

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
        key_lineages <- 
        lapply(lineages_t3, function(x){
            sum_key_mutations_uk() %>% 
                filter(lineage == x | str_detect(lineage, sublineage_regex(x)) ) %>% 
                select(-lineage) %>%
                summarise(sequences_sum = sum(sequences)) %>% 
                mutate(lineage = x, .before = 1)
        }) %>% bind_rows()
        
        key_lineages_28 <- 
            lapply(lineages_t3, function(x){
                sum_key_mutations_uk(date_from = sample_date_28) %>% 
                    filter(lineage == x | str_detect(lineage, sublineage_regex(x)) ) %>% 
                    select(-lineage) %>%
                    summarise(sequences_sum = sum(sequences)) %>% 
                    mutate(lineage = x, .before = 1) %>% rename(sequences_sum_28 = sequences_sum)
            }) %>% bind_rows() 
        
        inner_join(key_lineages, key_lineages_28) %>% rename(`Variant/ lineage` = lineage,	
                                                             `Cumulative sequences in UK` = sequences_sum,	 
                                                             `Sequences over 28 days` = sequences_sum_28)
    })
    output$table_4 <- renderTable({
        database %>%
            filter(replacement %in% escape_t4) %>%
            select(replacement, `numSeqs UK`, `numSeqs UK 28 days`) %>%  
            mutate(replacement = factor(replacement, levels = escape_t4)) %>%
            arrange(replacement) %>% 
            rename(`Amino acid replacement` = replacement, 
                   `Cumulative sequences in UK` = `numSeqs UK`,
                   `Sequences over 28 days` = `numSeqs UK 28 days`)
    })
    
    output$mutation_time <- renderPlot({

        mutations_s_uk %>% 
            filter(variant == input$variant
                                  & sample_date >= input$date_range[1]
                                  & sample_date <= input$date_range[2]) %>% 
            ggplot(aes(x = sample_date)) + geom_bar() + theme_minimal()
    })

    observeEvent(input$position, {
        updateSelectInput(session, 
                          "variant",
                          choices = 
                              mutations_s_uk %>% 
                              filter(position == input$position) %>% 
                              distinct(variant) %>% 
                              arrange(variant)
        )
    })
    
})
