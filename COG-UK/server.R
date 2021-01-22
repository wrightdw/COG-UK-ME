library(shiny)
library(tidyverse)
library(formattable)
library(plotly)
library(RColorBrewer)

lineages_t2 <- c("B.1", "B.1.177", "B.1.141", "B.1.258", "B.1.1", "B.1.1.7", "B.1.1.70", "B.1.351", "B.1.1.298")

lineages_t3 <- 
    c("B.1.1.7" = "UK associated variant. Has 17 mutations (14 replacements and 3 deletions) including: T1001I, A1708D, I2230T, SGF 3675-3677 del In the ORF1ab; 69-70 del, Y144 del, N501Y, A570D, P681H, T716I, S982A and D1118H in the Spike; Q27stop, R52I and Y73C in ORF8; D3L and S235F in the N. Noteworthily, N501Y enhances ACE2 binding affinity, and P681H occurs at the furin cleavage site, known for biological significance in membrane fusion.", 
      "B.1.351" = "Variant associated with South Africa. Has eight mutations in the Spike: D80A, D215G, E484K, N501Y, A701V, L18F, R246I and K417N. Three of these in the RBM, K417N, E484K and N501Y. K417N and E484K have been shown to escape some mAbs.", 
      "P.1" = "Variant associated with Brazil. Has 10 mutations in the Spike including L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y,H655Y and T1027I. Noteworthy  E484K, N501Y and K417T have biological significance.") %>% 
    enframe("lineage", "reason")

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

sequences_by_week <- consortium_uk %>% count(epi_week, name = "n_sequences")

positions_by_week <- 
  mutations_uk %>%
  group_by(epi_week, gene, position, .drop = FALSE) %>% # TODO don't count gene/position combinations that don't exist
  summarise(n_variant_sequences = n_distinct(sequence_name))

reference_counts <- inner_join(sequences_by_week, positions_by_week) %>% 
  mutate( across(n_variant_sequences, ~replace_na(.x, 0L)) ) %>% 
  mutate(n = n_sequences - n_variant_sequences, variant = "REF", .keep = "unused") 

mutation_counts <- mutations_uk %>%
  count(epi_week, gene, position, variant)

mutation_reference_counts <- bind_rows(mutation_counts, reference_counts)

shinyServer(function(input, output, session) {

    output$table_1 <- renderTable(
        database %>% 
            slice_max(`numSeqs UK`, n = 15) %>% 
            select(replacement, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`) %>% 
            mutate(`Sequences over the last 28 days in UK (%)` = percent(`numSeqs UK 28 days` / `numSeqs UK`) %>% as.character, .after = `numSeqs UK`) %>% 
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
    output$table_4 <- renderTable({
        database %>%
            filter(replacement %in% escape_t4) %>%
            select(replacement, `numSeqs UK`, `numSeqs UK 28 days`) %>%  
            arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), replacement) %>% 
            rename(`Amino acid replacement` = replacement, 
                   `Cumulative sequences in UK` = `numSeqs UK`,
                   `Sequences over 28 days` = `numSeqs UK 28 days`)
    })
    
    
    
    output$mutation_time <- renderPlotly({
# 
#         mutations_s_uk %>% 
#             filter(variant == input$variant
#                                   & sample_date >= input$date_range[1]
#                                   & sample_date <= input$date_range[2]) %>% 
#             ggplot(aes(x = sample_date)) + geom_bar() + theme_minimal()
      
      
      
      # gg_bar <- mutations_s_uk %>%
      #   # mutate(across(epi_week, as_factor)) %>%
      #   ggplot(aes(x = epi_week)) + geom_bar() + theme_minimal() + scale_x_discrete(drop=FALSE)
      # 
      # gg_bar_plotly <- gg_bar %>% ggplotly
      # gg_bar_plotly
      
      
      
      
      # mutations_s_uk %>%
      #   count(epi_week, position, variant) %>%
      #   filter(position == 614) %>%
      #   select(-position) %>%
      #   pivot_wider(names_from = variant, values_from = n) %>%
      #   # mutate(across(D614N, D614V), ~replace_na(.x, 0L)) %>%
      #   plot_ly(x = ~epi_week, y = ~D614G, type = 'bar', name = 'D614G') %>%
      #   add_trace(y = ~D614N, name = 'D614N') %>%
      #   add_trace(y = ~D614V, name = 'D614V') %>%
      #   layout(yaxis = list(title = 'Sequences'), barmode = 'stack')
      # 
      # 
      # 
      # 
      
      gg_bar <- 
        mutation_reference_counts %>%
        filter(gene == input$gene & position == input$position
                # & sample_date >= input$date_range[1]
                # & sample_date <= input$date_range[2]
               ) %>%
        select(-position, -gene) %>%
        ggplot(aes(fill=variant, y=n, x=epi_week)) +
        geom_bar(position="stack", stat="identity") +
        scale_x_discrete(drop=FALSE) +
        theme_classic() +
        # xlab("SNP") +
        # ylab("Participants") +
        # ggtitle(title) +
        scale_fill_manual(values = brewer.pal(name = "Dark2", n = 7))
      
      gg_bar_plotly <- gg_bar %>% ggplotly
      gg_bar_plotly
      
    })

    observeEvent(input$gene, {
      updateSelectInput(session, "position",
                        choices = mutations_uk %>%
                          filter(gene == input$gene) %>%
                          distinct(position) %>%
                          arrange(position))
    })

})
