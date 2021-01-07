library(shiny)
library(tidyverse)

shinyServer(function(input, output, session) {

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
    
    output$raw_database <- renderTable(
        database %>% 
            slice_max(`numSeqs UK`, n = 15) %>% 
            select(replacement, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`)
        
    )
    
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
