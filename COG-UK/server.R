#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

shinyServer(function(input, output, session) {

    output$mutation_time <- renderPlot({
        
        mutations %>% filter(gene == input$gene 
                             & variant == input$variant
                             & sample_date >= input$date_range[1] 
                             & sample_date <= input$date_range[2]) %>% ggplot(aes(x = sample_date)) + geom_bar() + theme_minimal()
    })
    
    output$raw_table <- DT::renderDataTable(
        mutations %>% 
        filter(gene == input$gene 
               & position == input$position 
               & variant == input$variant 
               & sample_date >= input$date_range[1] 
               & sample_date <= input$date_range[2])
    )
    
    output$summary <- renderTable({
        mutations %>% 
            filter(gene == input$gene 
                   & position == input$position 
                   & variant == input$variant
                   & sample_date >= input$date_range[1] 
                   & sample_date <= input$date_range[2]) %>% count(country)
    })

    observeEvent(input$gene, {
        updateSelectInput(session, "position",
                          choices = mutations %>% filter(gene == input$gene) %>% distinct(position) %>% arrange(position))
    })
    
    observeEvent(input$position, {
        updateSelectInput(session, "variant",
                          choices = mutations %>% filter(gene == input$gene & position == input$position) %>% distinct(variant) %>% arrange(variant))
    })
    
})
