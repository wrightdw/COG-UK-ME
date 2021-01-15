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
                        h1(paste0("COG-UK report on SARS-CoV-2 Spike mutations of interest")),
                        h3("2nd January 2021"),
                        
                        h2("Data source and processing"),
                        p("The analysis described in this report is based on # UK-derived genomes UK-derived genomes sequenced by COG-UK (complete data in MRC-CLIMB database to (add date)."),
                        p("A report of the geographic distribution and prevalence of SARS-CoV-2 lineages in general, and global variants of interest, can be found here (link to https://cov-lineages.org/global_report.html). Amino acid replacement, insertion and deletion counts for all SARS-CoV-2 genes in the global GISAID database can be found here (http://cov-glue.cvr.gla.ac.uk/)."),
                        
                        tabBox(
                            title = "Report tables", width = 12,
                            id = "tabs_report",
                            tabPanel("Table 1", 
                                     h3("1. High frequency spike gene mutations"),
                                     p("Individual amino acid replacements detected at the fifteen highest frequencies in UK
                                     genomes are shown in Table 1. Table S1 (link to table S1) shows all amino acid
                                     replacements detected in complete SARS-CoV-2 genomes in the COG-UK dataset when
                                     present in 5 or more sequences. Insertions or deletions are not shown. "),
                                     h4("Table 1. Spike mutations (top 15) present in the UK at high frequency"),
                                     p(em("NB number of genomes is not equal to number of COVID-19 cases as data have not been deduplicated.")),
                                     tableOutput("table_1")
                            ),
                            tabPanel("Table 2", 
                                     
                                     
                                     h3("2. Spike gene mutations of potential importance"),
                                     
                                     p("Single spike gene mutations of potential or clinical and public health importance based on
                                     current evidence are listed in Table 2."),
                                     
                                     p("CAVEATS:"),
                                     tags$ul(
                                         tags$li("The table aims to provide information on individual mutations, but this is rapidly
                                         becoming an over-simplification because mutations are increasingly arising in a
                                         range of combinations."),
                                         
                                         tags$li("The term ‘escape’ is used in the table as shorthand to mean weaker neutralisation of
                                         the virus by convalescent plasma from people who have been infected with SARS-
                                             COV-2, and/or some monoclonal antibodies (mAbs) that may be given to patients
                                         with COVID-19.")
                                     ),
                                     
                                     h4("Table 2. Spike S gene mutations, lineage associations and reason for interest, UK
                                     lineages"),
                                     p(em("NB Numbers are lower than in Table 1 because Table 2 only considers specific lineages.")),
                                     tableOutput("table_2")
                            ), 
                            tabPanel("Table 3",                         
                                     h3("3. Variants of concern being monitored by UK PHAs"),
                                     tableOutput("table_3")
                            ), 
                            tabPanel("Table 4",                         
                                     h3("4. Mutations with known antigenic role"),
                                     tableOutput("table_4")
                            )
                        ),
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

