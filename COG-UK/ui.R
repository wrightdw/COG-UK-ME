library(shiny)
library(DT)

library(shinydashboard)

dashboardPage(
    dashboardHeader(title = "COG-UK"),
    dashboardSidebar(
        # selectInput("gene", "Gene:", mutations %>% distinct(gene) %>% arrange(gene), selected = "S"),
        # selectInput("position", "Position:", mutations %>% distinct(position) %>% arrange(position)),
        # selectInput("variant", "Variant:", mutations %>% distinct(variant) %>% arrange(variant)),
        # dateRangeInput('date_range',
        #                label = 'Date range:',
        #                start = min(mutations$sample_date), end = max(mutations$sample_date)
        # ),
        
        sidebarMenu(
            menuItem("Report", tabName = "report", icon = icon("table"))
            # menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
            # menuItem("Raw Data", tabName = "raw", icon = icon("database"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "report",
                    fluidRow(
                        h2(paste0("COG-UK report on SARS-CoV-2 Spike mutations of interest")),
                        h4("2nd January 2021"),
                        h3("Table 1: Spike mutations (top 15) present in the UK at high frequency"),
                        tableOutput("table_1"),
                        
                        h3("Table 2: Spike gene mutations of potential importance"),
                        tableOutput("table_2")
                        
                    )
            )
        ) # end tabItems
            
            # tabItem(tabName = "dashboard",
            #         fluidRow(
            #             box(plotOutput("mutation_time"), width = 12),
            #             box(h4("Counts of variant by country"),
            #                 tableOutput("summary"), width = 6)
            #         ), 
            # ),
            
            # tabItem(tabName = "raw",
            #         h2("Raw Data"), 
            #         h3("Mutations"),
            #         DT::dataTableOutput("raw_mutations"), 
            #         h3("Spike Database"),
            #         DT::dataTableOutput("raw_database")
            # )
     #   )
    ) # end dashboardBody
)

