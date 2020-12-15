library(shiny)
library(DT)


library(shinydashboard)

dashboardPage(
    dashboardHeader(title = "COG-UK"),
    dashboardSidebar(
        selectInput("gene", "Gene:", mutations %>% distinct(gene) %>% arrange(gene), selected = "S"),
        selectInput("position", "Position:", mutations %>% distinct(position) %>% arrange(position)),
        selectInput("variant", "Variant:", mutations %>% distinct(variant) %>% arrange(variant)),
        dateRangeInput('date_range',
                       label = 'Date range:',
                       start = min(mutations$sample_date), end = max(mutations$sample_date)
        ),
        
        sidebarMenu(
            menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
            menuItem("Raw Data", tabName = "raw", icon = icon("table"))
        )
    ),
    dashboardBody(
        tabItems(
            # First tab content
            tabItem(tabName = "dashboard",
                    fluidRow(
                        box(plotOutput("mutation_time"), width = 12)
                    ), 
            ),
            
            # Second tab content
            tabItem(tabName = "raw",
                    h2("Raw Data"), 
                    DT::dataTableOutput("raw_table")
            )
        )
    )
)

