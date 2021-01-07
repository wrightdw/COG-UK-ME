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
                        tableOutput("raw_database")
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

