library(shiny)
library(DT)
library(shinydashboard)
library(formattable)
library(plotly)

dashboardPage(
    dashboardHeader(title = tags$a(href='http://cogconsortium.uk', target = "_blank",
                                   tags$img(src='image2.png', height = "50px"))),
    dashboardSidebar(
        
        sidebarMenu(id="sidebar_menu",
            menuItem("Report", tabName = "report", icon = icon("table")), 
            menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard"))
        ),
        
        conditionalPanel(condition =  "input.sidebar_menu == 'dashboard'",
                         selectInput("gene", "Gene:", mutations_uk %>% distinct(gene) %>% arrange(gene), selected = "S"),
                         selectInput("position", "Position:", mutations_uk %>% distinct(position) %>% arrange(position), selected = "614")
                         # selectInput("variant", "Variant:", mutations_s_uk %>% distinct(variant) %>% arrange(variant), selected = "D614G"),
                         # dateRangeInput('date_range',
                         #                label = 'Date range:',
                         #                start = min(mutations_s_uk$sample_date), end = max(mutations_s_uk$sample_date)
                         # )
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "report",
                fluidRow(
                    h1("COG-UK report on SARS-CoV-2 Spike mutations of interest"),
                    h3(dataset_date %>% format("%A %d %B %Y")),
                    
                    h2("Data source and processing"),
                    p("The analysis described in this report is based on ", n_distinct(consortium_uk$sequence_name) %>% comma(format = "d")
," UK-derived genomes sequenced by COG-UK (complete data in MRC-CLIMB database to ", max(consortium_uk$sample_date) %>% format("%d/%m/%Y"), "."),
                    p("A report of the geographic distribution and prevalence of SARS-CoV-2 lineages in general, and global variants of interest, can be found ", a(href = "https://cov-lineages.org/global_report.html", target = "_blank", "here"), ". Amino acid replacement, insertion and deletion counts for all SARS-CoV-2 genes in the global GISAID database can be found ", a(href = "http://cov-glue.cvr.gla.ac.uk/", target = "_blank", "here"), "."),
                    
                    tabBox(
                        title = "Analysis", width = 12,
                        id = "tabs_report",
                        tabPanel("Table 1", 
                                 h3("1. High frequency spike gene mutations"),
                                 p("Individual amino acid replacements detected at the fifteen highest frequencies in UK
                                 genomes are shown in Table 1. Table S1 (link to table S1) shows all amino acid
                                 replacements detected in complete SARS-CoV-2 genomes in the COG-UK dataset when
                                 present in 5 or more sequences. Insertions or deletions are not shown. "),
                                 
                                 h4("Table 1. Spike mutations (top 15) present in the UK at high frequency"),
                                 p(em("NB Number of genomes is not equal to number of COVID-19 cases as data have not been deduplicated.")),
                                 tableOutput("table_1")
                        ),
                        tabPanel("Table 2", 
                                 h3("2. Spike gene mutations of potential importance"),
                                 p("Single spike gene mutations of potential or clinical and public health importance based on
                                 current evidence are listed in Table 2."),
                                 
                                 h5("CAVEATS:"),
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
                                 h3("3. Global variants of interest"),
                                 h4("Table 3. Variants of concern being monitored by UK PHAs"),
                                 tableOutput("table_3")
                        ), 
                        tabPanel("Table 4",                         
                                 h3("4. Spike gene mutations of potential importance"),
                                 p("Table 4 lists those mutations in the spike gene identified in the UK dataset that have been
                                 associated with weaker neutralisation of the virus by convalescent plasma from people who
                                 have been infected with SARS-COV-2, and/or some mAbs that may be given to patients with
                                 COVID-19 (referred to below as ‘escape’)."),
                                 p(strong("There is no evidence at the time of writing for this impacting on the efficacy of current
                                 vaccines or the immune response to natural SARS-CoV-2 infection.")),
                                 
                                 h4("Table 4. Reported ‘escape’ mutations in the spike gene present in the UK"),
                                 tableOutput("table_4")
                        )
                    ), # end tabBox
                    
                    h2("Appendix 1"),
                    h3("Background"),
                    p("Mutations arise naturally in the SARS-CoV-2 genome as the virus replicates and circulates in the
                    human population. As a result of this on-going process, many thousands of mutations have already
                    arisen in the SARS-CoV-2 genome since the virus emerged in late 2019. As mutations continue to
                    arise, novel combinations of mutations are increasingly observed. The vast majority of mutations have
                    no apparent effect on the virus. Only a very small minority are likely to be important and change the
                    virus in any appreciable way. This could include a change in the ability to infect/transmit between
                    people; a change in disease severity; or a change in the way the virus interacts with the immune
                    system (including the response generated by a vaccine). We pay most attention to mutations in the
                    gene that encodes the Spike protein, which is associated with viral entry into cells and it is relevant in
                    the context of immunity and vaccine efficacy."),
                    h3("Definitions"),
                    tags$ul(
                        tags$li(em("Mutation"), " is used to describe a change of a nucleotide in the virus RNA genome, a subset of which
                        results in a change in amino acid (sometimes referred to as a substitution or replacement), or a
                        mutation can refer to a deletion or insertion event in the virus genome. By convention an amino acid
                        change is written N501Y to denote the wildtype (N, asparagine) and replacement amino acid (Y, 
                        tyrosine) at site 501 in the amino acid sequence."),
                        
                        tags$li(em("Viral variant"), " refers to a genetically distinct virus with different mutations to other viruses. Variant can
                        also refer to the founding virus of a cluster/lineage and used to refer collectively to the resulting
                        variants that form the lineage."),
                        
                        tags$li(em("Lineages"), " are assigned combining genetic and, in the case of SARS-CoV-2 due to weak phylogenetic
                        signals, also with epidemiological data. COG-UK uses the nomenclature system introduced by
                        Rambaut et al. (2020), see https://cov-lineages.org."),
                        
                        tags$li(em("VUI"), " is used by Public Health England to indicate Variant Under Investigation."),
                        
                        tags$li(em("VOC"), " is used by Public Health England to indicate Variant of Concern. Lineage B.1.1.7 was initially
                        named by Public Health England as VUI 202012/01 (Variant Under Investigation, year 2020, month
                        12, variant 01) and subsequently redesignated as VOC-202012/01 (Variant of Concern 202012/01).")
                    )
                ) # end fluidRow
            ), # end tabItem
            
            tabItem(tabName = "dashboard",
                    fluidRow(
                        box(plotlyOutput("mutation_time"), width = 12)
                    )
            )
        ) # end tabItems
    ) # end dashboardBody
)

