library(shiny)
library(tidyverse)

shinyServer(function(input, output, session) {

    output$table_1 <- renderTable(
        database %>% 
            slice_max(`numSeqs UK`, n = 15) %>% 
            select(replacement, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`)
        
    )
    
    output$table_2 <- renderTable({
        n_uk_lineages <- 
            consortium_uk %>% 
            # filter(sample_date >= params$sample_date_28) %>% 
            group_by(lineage) %>%
            summarise(sequences = n(), 
                      D614G = sum(d614g == "G"), 
                      A222V = sum(a222v == "V"), 
                      N439K = sum(n439k == "K"), 
                      N501Y = sum(n501y == "Y"), 
                      Y453F = sum(y453f == "F"),
                      DEL_69_70 = sum(del_21765_6 == "del"), 
                      N439K_DEL_69_70 = sum(n439k == "K" & del_21765_6 == "del"), 
                      N501Y_DEL_69_70 = sum(n501y == "Y" & del_21765_6 == "del"),
                      Y453F_DEL_69_70 = sum(y453f == "F" & del_21765_6 == "del")
            )
        
        n_uk_lineages %>% 
            filter(lineage == "B.1" |str_detect(lineage, "^B\\.1\\.")) %>% 
            select(-lineage) %>% 
            summarise_all(funs(sum))
        
    })
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
