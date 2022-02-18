library(shiny)
library(DT)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyjs)
library(formattable)
library(plotly)

dashboardPage(
    title = "COG-UK/Mutation Explorer",
    skin = "black",
    
    header = dashboardHeader(
        title = tags$a(
        href = '.',
        tags$img(src = 'COG_ME_LOGO2-ORIZCOLOR.png', width = "200px", title = "COG-UK/Mutation Explorer")
    ),
    
    tags$li(class = "dropdown",
            div(
                dashboardLabel(
                    "FOR RESEARCH PURPOSES ONLY",
                    status = "primary",
                    style = "square"
                ),
                style = "padding:15px"
            ))),
    
    sidebar = dashboardSidebar(
        disable = FALSE, 
        minified = FALSE, 
        collapsed = FALSE,
        
        sidebarMenu(
            id = "sidebar_menu",
            menuItem("VOCs/VUIs in the UK", tabName = "vui_voc", selected = TRUE, icon = icon("viruses")),
            menuItem("VOC Spike Structures", tabName = "strucure_voc", icon = icon("virus")),
            menuItem("Antigenic Mutations", tabName = "immunology", icon = icon("shield-virus")),
            menuItem("VOCs/VUIs + Antigenicity", tabName = "figure_1", icon = icon("fire-alt")),
            menuItem("T Cell Epitope Mutations", tabName = "t_cell", icon = icon("disease")),
            # menuItem("Mutation Counts", tabName = "report", icon = icon("virus")),
            menuItem("Mutations by Week", icon = icon("eye"), tabName = "dashboard"),
            menuItem("Spike Profiles", icon = icon("chart-line"), tabName = "spike_profiles"),
            menuItem("Drug Resistance", icon = icon("prescription-bottle-alt"), tabName = "therapeutics"),
            menuItem("Ronapreve",  tabName = "ronapreve", icon = icon("pills")),
            menuItem("Geographical Distribution", tabName = "map", icon = icon("map")),
            menuItem("Omicron and mAb", tabName = "omicron", icon = icon("disease")),
            menuItem("About", tabName = "about", icon = icon("info-circle"))
        ),
        
        conditionalPanel(
            condition =  "input.sidebar_menu == 'vui_voc'",
            hr(),
            prettySwitch("variant_percentage", "Percentage", FALSE, status = "info", fill = TRUE),
            prettySwitch("variant_day", "By day", FALSE, status = "info", fill = TRUE),
            prettySwitch("other_switch", "Exclude Other", FALSE, status = "info", fill = TRUE),
            prettyCheckboxGroup("variant_vui_voc", "Variant:",
                                vui_voc_lineages,
                                selected = c("B.1.1.7", "B.1.617.2", "BA.1", "BA.2"),
                                shape = "curve",
                                status = "info",
                                fill = TRUE),
            conditionalPanel(
                condition = "input.variant_vui_voc.includes('B.1.617.2')",
                prettyRadioButtons(inputId = "variant_delta",  label = "Delta sublineage:",
                                   choices = c("B.1.617.2", "AY.4", "AY.4.2"),
                                   shape = "round",
                                   status = "info",
                                   selected = "AY.4.2",
                                   fill = TRUE)
            ),
            
            prettyRadioButtons(
              inputId = "nations_vui_voc",
              label = "UK nation:", 
              choices = c("UK", "England", "Northern Ireland" = "Northern_Ireland", "Scotland", "Wales"),
              inline = FALSE, 
              status = "info",
              fill = TRUE,
              selected = "UK"
            ),
            
        ),
        
        conditionalPanel(
            condition =  "input.sidebar_menu == 'dashboard'",
            hr(),
            selectInput("gene", "Gene:", mutation_reference_counts %>% distinct(gene) %>% arrange(gene), 
                        selected = "S"),
            selectInput("position", "Position:", 
                        mutation_reference_counts %>% 
                            filter(gene == "S") %>% 
                            distinct(position) %>% 
                            arrange(position), 
                        selected = "614", selectize = TRUE),
            prettySwitch("percentage", "Percentage", FALSE, status = "info", fill = TRUE),
            prettySwitch("ref", "Wild type / other", TRUE,  status = "info", fill = TRUE),
            
            prettyRadioButtons(
                inputId = "nation",
                label = "UK nation:", 
                choices = c("UK", "England", "Northern Ireland" = "Northern_Ireland", "Scotland", "Wales"),
                inline = FALSE, 
                status = "info",
                fill = TRUE,
                selected = "UK"
            ),
            
            chooseSliderSkin("Modern", "#5bc0de"), # Bootstrap info colour
        ),
        
        conditionalPanel(
            condition = "input.sidebar_menu == 'immunology'", 
            hr(),
            prettyCheckboxGroup("escape", "Escape:",
                                c("Monoclonal Ab" = "monoclonal",
                                  "Convalescent sera" = "convalescent",
                                  "Vaccine sera" = "vaccine"),
                                shape = "curve",
                                status = "info")
        ),

        conditionalPanel(
            condition =  "input.sidebar_menu == 'ronapreve'",
            hr(),
            prettySwitch("ronapreve_28", "Latest 28 days", TRUE, status = "info", fill = TRUE)
        ),
        
        conditionalPanel(
            condition =  "input.sidebar_menu == 't_cell'",
            hr(),
            prettyRadioButtons("t_cell_experiment", "Type of experiment:",
                                c("Reduced T cell recognition" = "recognition",
                                  "Epitope studies" = "epitope_studies"),
                                shape = "round",
                                status = "info",
                                selected = "recognition")
        ),
    conditionalPanel(
        condition =  "input.sidebar_menu == 'map'",
        hr(),
        prettySwitch("percentage_map", "Percentage", TRUE , status = "info", fill = TRUE)
    )),
    
    body = dashboardBody(
        useShinyjs(), # set up the dashboard to use shinyjs 
        
        tags$head(
            tags$link(rel = "shortcut icon", href = "favicon.png"),
            tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
            tags$link(rel = "stylesheet", type = "text/css", href = "lightbox2-2.11.3/dist/css/lightbox.min.css"),
            tags$script(src = "lightbox2-2.11.3/dist/js/lightbox.min.js"),
            includeHTML("google-analytics.html")
        ),
        
        fluidRow(
            infoBox("Latest UK Sequence", max(consortium_uk$sample_date) %>% format("%A %d %B %Y"), icon = icon("calendar")),
            infoBox("Data Compiled on", dataset_date %>% format("%A %d %B %Y"), icon = icon("calendar-day"), color = "light-blue"),
            infoBox("UK Genomes", total_sequences %>% comma(format = "d"), icon = icon("virus"), color = "teal")
        ),
        
        tabItems(
            tabItem(tabName = "strucure_voc",
                    fluidRow(
                        box(title = "Spike protein visualisations", 
                            closable = FALSE, 
                            width = 12, 
                            status = "orange", 
                            collapsible = FALSE, 
                            icon = icon("microscope"), 
                            
                            p("Spike protein structures showing locations of amino acid residues that are mutated in each
                              variant of concern (VOC). The spike protein protrudes from the surface of the SARS-CoV-2
                              virus, is responsible for initating binding to and entry into host cells, and is also the
                              primary target for antibodies that recognise the virus."),
                            p("Each spike consists of three identical protein chains (shown in teal, blue and gold). Here,
                              spike is shown in its 'open' conformation in which the receptor-binding domain of the teal
                              chain is 'up' exposing the binding site that recognises the human ACE2 receptor."),
                            p("On each chain, the locations of amino acid substitutions, deletions (Δ), and insertions
                              (ins) which distinguish each VOC from the original genotype (Wuhan-Hu-1), are highlighted as
                              opaque-surface spheres coloured in red (where they are labelled), blue and gold. The
                              substitution D614G which is shared by common descent by all lineage B.1 descendants is
                              italicised."),
                            p("Visualisations are made using the ectodomain of a complete spike model (Woo et al., 2020)
                              which is in turn based upon a partial cryo-EM structure (RCSB Protein Data Bank (PDB) ID:
                              6VSB (Wrapp et al., 2020)).")
                    )), # End of box and fluid row
                    
                    fluidRow(
                          box(title = "Spike protein mutations (Omicron: BA.1)", 
                                   closable = FALSE, 
                                   width = 6, 
                                   status = "orange", 
                                   collapsible = TRUE, 
                                   icon = icon("microscope"),
                                   height = 480,
                                   tags$a(
                                       href = "mutants_BA.1_ME_web.png",
                                       `data-lightbox` = "structure",
                                       `data-title` = "Spike protein mutations (Omicron: BA.1)",
                                       `data-alt` = "Omicron BA.1 spike structure with mutations",
                                       img(src = "mutants_BA.1_ME_web.png",
                                           class = "center-block img-responsive", 
                                           alt = "Omicron BA.1 spike structure with mutations")
                                   )),
                          
                          box(title = "Spike protein mutations (Delta: B.1.617.2)", closable = FALSE, width = 6, 
                              status = "orange", collapsible = TRUE, icon = icon("microscope"),
                              tags$a(
                                  href = "mutants_B.1.617.2_ME_web.png",
                                  `data-lightbox` = "structure",
                                  `data-title` = "Spike protein mutations (Delta: B.1.617.2)",
                                  `data-alt` = "Delta B.1.617.2 spike structure with mutations",
                                  img(src = "mutants_B.1.617.2_ME_web.png", 
                                      class = "center-block img-responsive", 
                                      alt = "Delta B.1.617.2 spike structure with mutations")
                              ))
                        ),
                    
                    fluidRow(
                        box(title = "Spike protein mutations (Alpha: B.1.1.7)", 
                            closable = FALSE, 
                            width = 6, 
                            status = "orange", 
                            collapsible = TRUE, 
                            icon = icon("microscope"),
                            tags$a(
                                href = "mutants_B.1.1.7_ME_web.png",
                                `data-lightbox` = "structure",
                                `data-title` = "Spike protein mutations (Alpha: B.1.1.7)",
                                `data-alt` = "Alpha B.1.1.7 spike structure with mutations",
                                img(src = "mutants_B.1.1.7_ME_web.png", 
                                    class = "center-block img-responsive", 
                                    alt = "Alpha B.1.1.7 spike structure with mutations")
                            )),
                        
                        box(title = "Spike protein mutations (Gamma: P.1)", 
                            closable = FALSE, 
                            width = 6, 
                            status = "orange", 
                            collapsible = TRUE, 
                            icon = icon("microscope"),
                            tags$a(
                                href = "mutants_P.1_ME_web.png",
                                `data-lightbox` = "structure",
                                `data-title` = "Spike protein mutations (Gamma: P.1)",
                                `data-alt` = "Gamma P.1 spike structure with mutations",
                                img(src = "mutants_P.1_ME_web.png", 
                                    class = "center-block img-responsive", 
                                    alt = "Gamma P.1 spike structure with mutations")
                            ))
                    ), 
                    
                    fluidRow(
                        box(title = "Spike protein mutations (Beta: B.1.351)", 
                            closable = FALSE, 
                            width = 6, 
                            status = "orange", 
                            collapsible = TRUE, 
                            icon = icon("microscope"), 
                            
                            tags$a(
                                href = "mutants_B.1.351_ME_web.png",
                                `data-lightbox` = "structure",
                                `data-title` = "Spike protein mutations (Beta: B.1.351)",
                                `data-alt` = "Beta B.1.351 spike structure with mutations",
                                img(src = "mutants_B.1.351_ME_web.png", 
                                    class = "center-block img-responsive", 
                                    alt = "Beta B.1.351 spike structure with mutations")
                            ))
                    ), 
                    
                    ), # end of Spike structure tab
            
            tabItem(tabName = "about",
                    fluidRow(
                        box(title = "About COG-UK/Mutation Explorer", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("info-circle"),  
                            
                            h3("Overview"),
                            p("The COG-UK/Mutation Explorer (COG-UK/ME) provides information and structural context on mutations and associated variants in the genes encoding SARS-COV-2 proteins that have been identified from sequence data generated by the COVID-19 Genomics (COG-UK) Consortium.",
                              "We focus on SARS-CoV-2 spike gene mutations of potential or known importance based on epidemiological, clinical and/or experimental observations."),
                            p("The Mutation Explorer comprises of:"),
                            tags$ol(
                                tags$li("the designated global variants of concern and their structural contexts"),
                                tags$li("high frequency individual amino acid replacements, a subset of which may be important"),
                                tags$li("heatmap of antigenic mutations accumumulating on top of lineage-defining mutations of VOC/VUI"),
                                tags$li("frequency plots for mutations at specific residue in SARS-CoV-2 ORFs (Mutation Visualiser)"),
                                tags$li("mutations of potential antigenic significance as indicated by experimental studies: shown to lead to weaker neutralisation of the virus by convalescent plasma from people who have been infected with SARS-CoV-2 and/or demonstrated escape from some monoclonal antibodies (mAbs) that may be given to patients with COVID-19 (Antigenic Information: Antibody Sites)"),
                                tags$li("mutations in T cell epitopes as indicated by experimental studies (Antigenic Information: T Cell Epitopes).")
                            ),
                            
                            h3("Data source and processing"),
                            p(.noWS = c("after-begin", "before-end"), "The analysis described in this report is based on ", strong(total_sequences %>% comma(format = "d"), .noWS = "outside"), " UK-derived genomes after dedeuplication, sequenced by COG-UK: complete data in the MRC-CLIMB database to ", strong(dataset_date %>% format("%d/%m/%Y"), .noWS = "outside"), ", with the latest sequence from ", strong(max(consortium_uk$sample_date) %>% format("%d/%m/%Y"), .noWS = "outside"),  "."),
                            p("A report of the geographic distribution and prevalence of SARS-CoV-2 lineages in general, and global variants of interest, can be found ", a(href = "https://cov-lineages.org/global_report.html", target = "_blank", "here", .noWS = "outside"), ". Amino acid replacement, insertion and deletion counts for all SARS-CoV-2 genes in the global GISAID database can be found ", a(href = "http://cov-glue.cvr.gla.ac.uk/", target = "_blank", "here", .noWS = "outside"), ".", .noWS = c("after-begin", "before-end")),
                            
                            h3("Limitations"),
                            tags$ol(
                                tags$li("This report is for information only. The clinical and public health importance of any single mutation, or combination of mutations cannot be determined from sequence data alone."),
                                tags$li("Putative evidence for the importance of any single mutation, or combination of mutations can be derived from computational biology and further evaluated by laboratory experiments. Genomic and laboratory evidence then need to be combined with clinical datasets that are designed to allow detection of increased transmissibility, change in disease severity, drug resistance or altered vaccine efficacy. For this reason, surveillance and risk assessment of mutations and variants is a multi-agency process involving UK Public Health Agencies who have access to detailed information on patients and populations, and other groups including NERVTAG (New and Emerging Respiratory Virus Threats Advisory Group)."),
                                tags$li("COG-UK generates around 10,000 genomes a week, which will rise to 20,000 per week by March 2021. When COVID-19 infection rates are high, not all viruses from infected people will be sequenced and some mutations at low frequency will not be detected, but COG-UK aims to take representative samples from across the UK.")
                            ),
                            
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
                            
                            h3("Glossary"),
                            tags$ul(
                                tags$li(em("Mutation"), "is used to describe a change of a nucleotide in the virus RNA genome, a subset of which
                                        results in a change in amino acid (sometimes referred to as a substitution or replacement), or a
                                        mutation can refer to a deletion or insertion event in the virus genome. By convention an amino acid
                                        change is written N501Y to denote the wildtype (N, asparagine) and replacement amino acid (Y, 
                                        tyrosine) at site 501 in the amino acid sequence."),
                                
                                tags$li(em("Viral variant"), "refers to a genetically distinct virus with different mutations to other viruses. Variant can
                                        also refer to the founding virus of a cluster/lineage and used to refer collectively to the resulting
                                        variants that form the lineage."),
                                
                                tags$li(em("Lineages"), "are assigned combining genetic and, in the case of SARS-CoV-2 due to weak phylogenetic
                                        signals, also with epidemiological data. COG-UK uses the nomenclature system introduced by
                                        Rambaut et al. (2020), see https://cov-lineages.org."),
                                
                                tags$li(em("VUI"), "is used by Public Health England to indicate Variant Under Investigation."),
                                
                                tags$li(em("VOC"), "is used by Public Health England to indicate Variant of Concern.")
                            ),
                            
                            h3("Disclaimer"),
                            p("Dashboard reports are not advice.
                          They capture research findings which are always necessarily provisional.
                          They are for research use only.
                          Commercial use/resale is not permitted."),
                          
                          h3("Credits"),
                          p("COG-UK/ME is developed within and funded by the COVID-19 Genomics UK Consortium by Derek W. Wright, Joseph Hughes, William Harvey, MacGregor Cox, Rachel Colquhoun, Ben Jackson, Andrew Rambaut, Thomas Peacock, David L. Robertson, Alessandro M. Carabelli.",
                            "COG-UK/ME is based on the CLIMB framework, and maintained by the ", a(href = "https://www.gla.ac.uk/researchinstitutes/iii/cvr/", target = "_blank", .noWS = "outside", "MRC-University of Glasgow Centre for Virus Research"), ".",
                            "Follow", a(href ="https://twitter.com/CovidGenomicsUK", target = "_blank", "COG-UK"), "to be notified of updates.", .noWS = c("after-begin", "before-end")),
                          
                          h3("Contact Us"),
                          p("To request features or report issues, contact us on ", a(href = "https://github.com/wrightdw/COG-UK-ME/issues", target = "_blank", .noWS = "outside", "GitHub"), ".")
                        )
                    )),
            
            tabItem("vui_voc", 
                    fluidRow(
                        box(title = "Variants of concern (VOC) and under investigation (VUI) and any other variant by weeks and days", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("chart-line"),
                            
                            plotlyOutput("variant_time", height = "70vh"),
                            
                            chooseSliderSkin("Modern", "#5bc0de"), # Bootstrap info colour
                            sliderInput(
                                inputId = "variant_range", 
                                label = "Date range:",
                                min = lineages_days_uk_all %$% min(sample_date), 
                                max = lineages_days_uk_all %$% max(sample_date), 
                                value = c(lineages_weeks_uk_all %>% filter(lineage %in% vui_voc_lineages) %$% min(epi_date),
                                          lineages_days_uk_all %$% max(sample_date)),
                                step = 1,
                                ticks = FALSE,
                                animate = TRUE,
                                timeFormat = "%d %b %y"
                            ),
                            
                            p("Variant sequence counts are grouped either by week starting on Sunday or by day.
                              The most recent sequence data (approx. the last two weeks) have low sample numbers,
                              so are highlighted with a grey box for the last two weeks of the weekly chart 
                              or from the second-to-last Sunday onwards for the daily chart.")
                        ))
                    ,

                    fluidRow(
                        box(title = "Variants of concern (VOC) and under investigation (VUI) detected in the UK data", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("table"),
                            p("DISCLAIMER: COG-UK uses curated sequences for determining the counts of a given lineage. Other sources of information may be reporting cases with partial sequence information or other forms of PCR testing."),
                            dataTableOutput("table_3")
                        )
                    ),
                    
                    fluidRow(
                        column(width = 6, 
                               box(title = "Download metadata", closable = FALSE, width = NULL, 
                                   status = "info", collapsible = FALSE, icon = icon("file-download"),
                                   p("Download a CSV file, for each variant, containing COG-UK sequence name, sample date, epidemiological week, epidemiological week start date and global lineage. Cumulative UK sequences are filtered by the selected lineage of concern."), 
                                   selectizeInput("concern", "Choose lineage:",
                                                  choices = lineages_t3$lineage %>% sort),
                                   downloadButton("downloadConcern", "Download", class = "btn-info")),
                               
                               box(title = "Download table", closable = FALSE, width = NULL, 
                                   status = "info", collapsible = FALSE, icon = icon("file-download"), 
                                   p("Download a CSV file comprising complete table data."),
                                   downloadButton("downloadTable3", "Download", class = "btn-info"))
                        ))
            ),
            
            tabItem(tabName = "report",
                    
            ), # end tabItem
            
            tabItem(tabName = "figure_1", 
                    fluidRow(
                        box(width = 12, closable = FALSE,  status = "warning", collapsible = FALSE, icon = icon("fire-alt"),
                            title = "Antigenic amino acid replacements in variants of concern (VOC) and variants under investigation (VUI) in addition to their defining mutations",
                            fluidRow(
                                column(
                                    width = 2,
                                    prettyRadioButtons(
                                        inputId = "nation_antigenic",
                                        label = "UK nation:",
                                        choices = c("UK", "England", "Northern Ireland" = "Northern_Ireland", "Scotland", "Wales"),
                                        inline = FALSE,
                                        status = "info",
                                        fill = TRUE,
                                        selected = "UK"
                                    ),
                                    
                                    pickerInput(
                                        inputId = "lineage_antigenic",
                                        label = "Lineage:",
                                        choices = (function(){
                                          picks <- vui_voc %$% levels(lineage)
                                          names(picks) <- vui_voc %$% levels(lineage_display)
                                          picks
                                        })(),
                                        selected = "BA.1"
                                    ),
                                    
                                    sliderInput(
                                        inputId = "percentage_range",
                                        label = "Percentage range:",
                                        min = 0,
                                        max = 100,
                                        value = c(0, 100))
                                ), 
                                column(width = 10,                                      
                                       h4(textOutput("title_heatmap", inline = TRUE), class = "text-center"),
                                       plotOutput("antibody_heatmap", height = "auto")))
                        )
                    )
            ), # end tabItem
            
            tabItem(tabName = "dashboard",
                    fluidRow(box(
                        plotlyOutput("mutation_time", height = "80vh"), 
                        
                        sliderInput(
                          inputId = "mutation_range", 
                          label = "Date range:",
                          min = mutation_reference_counts %$% min(epi_date), 
                          max = mutation_reference_counts %$% max(epi_date), 
                          value = c(mutation_reference_counts %$% min(epi_date),
                                    mutation_reference_counts %$% max(epi_date)),
                          step = 7,
                          ticks = FALSE,
                          animate = TRUE,
                          timeFormat = "%d %b %y"
                        ),
                        
                        p("Mutation counts are grouped by week, starting on Sunday.
                              The most recent sequence data (approx. the last two weeks) have low sample numbers
                              so are highlighted with a grey box."),
                        
                        width = 12,
                        status = "info",
                        collapsible = FALSE,
                        closable = FALSE, 
                        title = "Amino acid replacement counts and percentages by week in the UK data",
                        icon = icon("chart-bar")
                    )),
                    fluidRow(
                        box(title = "Amino acid replacements detected in the UK data: counts, percentages, nations and date of first detection", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("table"),   
                            
                            p("Individual amino acid replacements detected in UK genomes are shown (sequences ≥ 5). Neither insertions nor deletions, nor synonymous mutations are included."),
                            p(em("NB Number of genomes is not equal to number of COVID-19 cases as data have not been deduplicated.")),
                            dataTableOutput("table_1")
                        ) # end box
                    ),
                    
                    fluidRow(
                        box(title = "Download metadata", closable = FALSE, width = 6, height = 500,
                            status = "info", collapsible = FALSE, icon = icon("file-download"),
                            fluidRow(
                                column(
                                    width = 8,
                                    p("Download a CSV file, for each amino acid replacement/del/ins, comprising COG-UK sequence name, sample date, epidemiological week, epidemiological week start date and global lineage. UK sequences are filtered by a 28 day period up to and including the most recent UK sequence date.")
                                ),
                                
                                (function() { 
                                    database_genome %<>% filter(`numSeqs UK` >= 5 & `numSeqs UK 28 days` > 0) 
                                    
                                    column(
                                        width = 4, 
                                        selectizeInput(
                                            inputId = "dataset_gene", 
                                            label = "Gene:",
                                            choices = database_genome$gene %>% levels, 
                                            selected = "S",
                                            options = list(dropdownParent = 'body')# prevent dropdown opening behind footer
                                        ),
                                        
                                        selectizeInput(
                                            inputId = "dataset", 
                                            label = "Amino acid replacement (type to search):",
                                            choices = NULL, 
                                            size = 10,
                                            options = list(dropdownParent = 'body') # prevent dropdown opening behind footer
                                        ),
                                        downloadButton("downloadData", "Download", class = "btn-info"))
                                })()
                            ) 
                        ), 
                        
                        box(title = "Download table", closable = FALSE, width = 6, height = 500,
                            status = "info", collapsible = FALSE, icon = icon("file-download"),
                            fluidRow(
                                column(
                                    width = 8,
                                    p("Download a CSV file comprising complete table data.")
                                ),
                                
                                column(
                                    width = 4, 
                                    downloadButton("downloadTable1", "Download", class = "btn-info"))
                            )
                        )
                    )
                    ),
            
            tabItem(tabName = "immunology",
                    value = "antibody",
                    box(title = "Spike amino acid replacements reported to confer antigenic change relevant to antibodies, detected in the UK data", 
                        closable = FALSE, width = 12,
                        status = "info", collapsible = FALSE, icon = icon("table"), 
                        p('The table lists those mutations in the spike gene identified in the UK dataset that have been
                                         associated with weaker neutralisation of the virus by convalescent plasma from people who
                                         have been infected with SARS-CoV-2, and/or monoclonal antibodies (mAbs) that recognise the SARS-CoV-2 spike protein (referred to below as "escape").'),
                        p(strong("There is no evidence at the time of writing for this impacting on the efficacy of current
                                         vaccines or the immune response to natural SARS-CoV-2 infection.")),
                        
                        DTOutput("table_4"),
                        
                        h4("Table key"),
                        h5("Confidence"),
                        tags$ul(
                            tags$li(
                                span(style = "background-color:firebrick; color:snow", "High", .noWS = "after"),
                                ": Antigenic role of mutation is supported by multiple studies including at least one that reports an effect observed with (post-infection serum) convalescent plasma.",
                                .noWS = c("after-begin", "before-end")
                            ),
                            tags$li(
                                span(style = "background-color:darkorange; color:white", "Medium", .noWS = "after"),
                                ": Antigenic role of mutation is supported by multiple studies.",
                                .noWS = c("after-begin", "before-end")
                            ),
                            tags$li(
                                span(style = "background-color:lemonchiffon; color:darkslategrey", "Lower", .noWS = "after"),
                                ": Mutation is supported by a single study.",
                                .noWS = c("after-begin", "before-end")
                            )
                        ),
                        
                        h5("Spike protein domain definitions"),
                        tags$ul(
                            tags$li("SP, signal protein (residues 1-13)"),
                            tags$li("NTD, N-terminal domain (14-303)"),
                            tags$li(
                                "RBD, receptor-binding domain (331-527) which includes the RBM, receptor-binding motif (437-508)"
                            ),
                            tags$li("FP, fusion peptide (815-834)"),
                            tags$li(
                                "Residues outside of these specific domains are labelled by subunit, S1 (residues 1-685) or S2 (residues 686-1173)"
                            )
                        )
                    ), # end box
                    
                    box(title = "Download data", closable = FALSE, width = 12, height = 500,
                        status = "info", collapsible = FALSE, icon = icon("file-download"),
                        fluidRow(
                            column(
                                width = 6,
                                p("Download a CSV file containing COG-UK sequence name, sample date, epidemiological week, epidemiological week start date and global lineage.
                                        Cumulative UK sequences are filtered by the selected amino acid replacement.")
                            ),
                            
                            column(
                                width = 6, 
                                selectizeInput("selectEscape", "Choose amino acid replacement:",
                                               choices = database %>% 
                                                   filter(!is.na(escape)) %>% 
                                                   filter(`numSeqs UK` > 0) %>% 
                                                   arrange(desc(`numSeqs UK`), mutation) %$% 
                                                   mutation, 
                                               options = list(dropdownParent = 'body')), # prevent dropdown opening behind footer
                                
                                downloadButton("downloadEscape", "Download", class = "btn-info")
                            )
                        ) # end fluidRow
                    ), # end box
            ), # end tabItem
            
            tabItem(tabName = "t_cell",
                    fluidRow(
                        box(title = "Spike amino acid replacements in T cell epitopes, detected in the UK data", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("table"),   
                            
                            p("T-cell epitope data have been compiled by Dhruv Shah, Sharon Hsu and Thushan de Silva, University of Sheffield."),
                            p("Data are filtered depending on the experiments that have been used either showing 'Reduced T cell recognition' or through 'Epitope studies' (options on the left hand side)."),
                            p("Predicted binding percentile rank values have been calculated by Morten Nielsen, The Technical University of Denmark."),
                            DTOutput("table_5"),
                            
                            h4("Table Key"),
                            tags$ul(
                                tags$li("WT Percentile Rank Value and Mut Percentile Rank Value: predicted IC50 nM for the corresponding reported restricting allele. 
                                    Predictions were performed using the NetMHCpan BA 4.1 algorithm, hosted by the IEDB."),
                                tags$li("Fold difference indicates Increase/decrease in affinity defined by a two-fold difference in predicted IC50 nM."),
                                tags$li("Binding is reported as a percentile rank value (as described ",a("here", href = "http://www.cbs.dtu.dk/services/NetMHCpan/", target = "_blank", .noWS = "outside"),"), the lower the value the stronger the binding.", 
                                        tags$ul(tags$li("For HLA-I, values less then 2 are binders and values less than 0.5 strong binders."), 
                                                tags$li("For HLA-II, values less then 5 are binders and values less than 1 strong binders.")),
                                        .noWS = c("after-begin"))
                            )
                        )
                    ),
                    
                    fluidRow(
                        box(title = "T cell epitope sequence viewer", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("disease"),
                            
                            sliderInput("epitope_position", "Position:",
                                        min = 1, max = wt %$% max(position),
                                        value = 484, step = 1),
                            
                            switchInput(
                                inputId = "epitope_ref", 
                                label = "Ref AA", 
                                offLabel = "Exclude",
                                onLabel = "Include",
                                value = FALSE,  
                                onStatus = "success", 
                                offStatus = "danger",
                                inline = TRUE,
                                size = "large"
                            ),
                            
                            switchInput(
                                inputId = "method",
                                offLabel = "Bits",
                                onLabel = "Probability",
                                value = FALSE,
                                label = "Method",
                                inline = TRUE,
                                onStatus = "info",
                                offStatus = "warning",
                                size = "large"
                            ),
                            
                            plotOutput("epitope_sequence"),
                            
                            p(.noWS = c("after-begin", "before-end"), 
                              "Move the slider to see sequence logs showing amino acid replacements in any epitope that overlaps on a specific position in the spike protein sequence.
                              Each letter represents an amino acid replacement present in a specific epitope. 
                              The number below the sequence logo shows the position relative to the start position of the epitope.
                              The height of a letter gives a measure of frequency of a mutation, whereas colour indicates amino acid chemistry.
                              Frequencies are normalised within each epitope on a scale of entropy (0 to 4.3 bits).
                              The wild-type epitope sequence and the start and end positions of the epitope are displayed above each sequence logo."),
                            p(.noWS = c("after-begin", "before-end"), 
                              "We thank Wagih, Omar.", 
                              a(.noWS = c("after-begin", "before-end"), href = "https://doi.org/10.1093/bioinformatics/btx469", target = "_blank", "ggseqlogo: a versatile R package for drawing sequence logos."), em("Bioinformatics (2017)"), ".")
                        )
                    ),
                    
                    fluidRow(
                        column(width = 10, offset = 1,
                               box(title = "Download table", closable = FALSE, width = 12,
                                   status = "info", collapsible = FALSE, icon = icon("file-download"),
                                   fluidRow(
                                       column(
                                           width = 8,
                                           p("Download a CSV file comprising complete table data.")
                                       ),
                                       
                                       column(
                                           width = 4,
                                           downloadButton("downloadTable5", "Download", class = "btn-info"))
                                   )
                               )
                        )
                    )
            ), # end tabItem t_cell
            
            
            # Spike Profiles tab
            tabItem(tabName = "spike_profiles",
                    fluidRow(
                        box(title = "Spike profile expansion and contraction", width = 12,
                            collapsible = FALSE, icon = icon("chart-line"),
                            p("Each spike profile is a set of amino acid substitutions listed relative to the original genotype (Wuhan-Hu-1).
                              Spike profiles sampled within 7 days of the latest UK sequence are plotted below, with count of sequences in the
                              latest 28-day period on the x-axis and a statistic estimating recent significant change in profile frequency on
                              the y-axis."),
                            p("Profiles of the Delta or Omicron variants of concern (VOCs) are described as amino acid substitutions relative
                              to the core VOC profiles listed below."),
                            p("Hovering the cursor over a point reveals the substitutions defining a spike profile, count in the latest 28-day
                              period and associated pango lineages. The plot and table below can be re-drawn including data for each of the
                              UK nations by selecting a dataset below."),
                            
                            prettyRadioButtons(
                                inputId = "spike_geo",
                                label = h5("Dataset:"),
                                choices = list("United Kingdom" = "United Kingdom" ,
                                               "England" = "England",
                                               "Northern Ireland" = "Northern Ireland",
                                               "Scotland" = "Scotland",
                                               "Wales" = "Wales"),
                                selected = "United Kingdom"
                            ),
                            
                            plotlyOutput("spikePlot_count28", height = 550),
                            br(),
                            p(strong("VOC core spike profiles:"), "'+' indicates additional substitutions and '-' marks the absence of a
                              substitution present in the core profiles below"),
                            p("Delta: T19R, G142D, Δ156-157/R158G, L452R, T478K, D614G, P681R, D950N"),
                            p("Omicron (BA.1): A67V, Δ69-70, T95I, G142D/Δ143-145, Δ211/L212I, ins214EPE, G339D, S371L, S373P, S375F, K417N,
                              N440K, G446S, S477N, T478K, E484A, Q493R, G496S, Q498R, N501Y, Y505H, T547K, D614G, H655Y, N679K, P681H, N764K,
                              D796Y, N856K, Q954H, N969K, L981F"),
                            p("Omicron (BA.2): T19I, L24S/Δ25-27, G142D, V213G, G339D, S371F, S373P, S375F, T376A, D405N, R408S, K417N, N440K,
                              S477N, T478K, E484A, Q493R, Q498R, N501Y, Y505H, D614G, H655Y, N679K, P681H, N764K, D796Y, Q954H, N969K"),
                            p(em("Substitutions in VOC core profiles that are sometimes not identified VOC sequences (due to amplicon dropout
                              during sequencing) are not listed as absent. For Delta profiles, G142D, is not listed as absent. For Omicron,
                              due to widespread undercalling of core substitutions, the absence of core substitutions is currently not
                              shown")),
                            
                            p(strong("Expansion/contraction:"), "For each profile,", em("i"), ", the absolute value for this statistic is
                              calculated using the observed frequency", em("O", tags$sub("i")), "of profile", em("i"), "in each of the most
                              recent 2-week periods", em("j"), "according to:",
                              tags$div(HTML(paste("Sum{(O", tags$sub("i,j"), " - E", tags$sub("i"), ")", tags$sup("2"), "/ E", tags$sub("i"), "}", sep = ""))),
                              "where", em("E", tags$sub("i")), "is the frequency of profile",  em("i"), "over the full 8-week period. This value
                              is shown as negative or positive according to direction of change and represents both the rate of change in
                              frequency and the overall frequency of a profile.")
                            ) # end of box
                    ), # end of fluid row
                    
                    fluidRow(
                        box(title = "Spike profiles detected in the UK during the last week", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("table"),
                            dataTableOutput("spikeTable"),
                            p(strong("Profile"), "lists the amino acid substitutions (currently deletions and insertions are not included) in the spike protein relative to the original genotype (Wuhan-Hu-1). Note: Incomplete spike profiles may be called where the underlying sequence data is incomplete. An example of this is the substitution G142D which is present in Delta sequences but often not called due to an amplicon dropout."),
                            p(strong("Amino acid substitutions"), "is the count of spike amino acid substitutions relative to the orginal genotype (Wuhan-Hu-1)."),
                            p(strong("Frequency change vs. prev 28 days (%)"), "28 day periods are calculated relative to the date of the most recent UK sequence. A blank cell indicates that a spike profile was not sequenced in the 28-day period preceding the most recent 28-day period. It shows the percentage change in the frequency of a profile among all sequenced genomes in the most recent 28-day period.")
                            )
                        )
            ), # end tabitem spike_profiles
            
            
            tabItem(tabName = "therapeutics",
                    fluidRow(
                        box(title = "Amino acid mutations reported to confer resistance to antiviral therapies, detected in the UK data", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("prescription-bottle-alt"),
                            DTOutput("table_therapeutics"),
                            br(),
                            p("The table lists those mutations in the SARS-CoV-2 genome identified in the UK dataset that have been associated with resistance of the virus to antiviral treatments. There is variation in the detail of the viral assays between the different studies displayed here.")
                        )
                    )
            ), # end tabItem therapeutics
            
            tabItem(tabName = "ronapreve",
                    fluidRow(
                        box(title = "Combinations of spike amino acid substitutions that may confer resistance to antibodies in the therapeutical antibody cocktail, Ronapreve.", closable = FALSE, width = 12,
                            status = "info", collapsible = FALSE, icon = icon("pills"),
                            h4(textOutput("title_ronapreve", inline = TRUE), class = "text-center"),
                            imageOutput("ronapreve_plot", 
                                        width = "100%", 
                                        height = "100%"),
                            br(),
                            p(.noWS = c("after-begin", "before-end"), "Plot showing the frequency of mutations affecting Ronapreve constituent monoclonal antibodies (mAbs) and their combinations (shown as lines) in cumulative UK SARS-CoV-2 genome sequence data. Spike amino acid substitutions known to affect either ", em("casirivimab", .noWS = "outside"), " or ", em("imdevimab", .noWS = "outside"), " mAbs were considered. The upper histogram shows the number of sequences per combination whereas the bottom left histogram shows the number of sequences with each specific substitution. Rows are coloured according to the mAb to which the greatest fold-decrease in binding was recorded (blue = ", em("casirivimab", .noWS = "outside"), ", orange = ", em("imdevimab", .noWS = "outside"), "), with a lighter shade indicating a fold-decrease of less than 100 and darker shade indicating 100 or greater."),
                            p(.noWS = c("after-begin", "before-end"), "The plot is generated using data from ", tags$a("here.", href="https://www.fda.gov/drugs/drug-safety-and-availability/fda-authorizes-revisions-fact-sheets-address-sars-cov-2-variants-monoclonal-antibody-products-under", target = "_blank", .noWS = "outside"))
                        )
                    )
            ), # end tabItem ronapreve
            
            tabItem(tabName = "map",
                    fluidRow( box(title = "Geographical distribution", closable = FALSE, width = 12,
                                       status = "info", collapsible = FALSE, icon = icon("map"),
                                  fluidRow(
                                      column(width = 2,
                                             prettyRadioButtons(
                                                 inputId = "variant_map",
                                                 label = "Lineage:",
                                                 choices = (function(){
                                                     picks <- 
                                                         vui_voc %$% 
                                                         levels(lineage)
                                                     
                                                     names(picks) <- vui_voc %$% levels(lineage_display)
                                                     
                                                     picks %>% c("Delta non-AY.4" = "Delta_minus_AY.4",
                                                                 "Delta non-AY.4.2" = "Delta_minus_AY.4.2")
                                                 })(),
                                                 selected = "BA.1",
                                             ),
                                             
                                             # selectizeInput(
                                             #     inputId = "antigenic_mutation",
                                             #     label = "Antigenic mutation:",
                                             #     choices = NULL, 
                                             #     size = 10,
                                             #     selected = ""
                                             # )
                                      ),
                                      
                                      column( width = 10,
                                              plotlyOutput("map", height = "70vh")
                                      )), # end of fluidRow
                                  
                                           sliderInput(
                                               inputId = "variant_date",
                                               label = "Date:",
                                               min = geo_all %$% min(epi_date),
                                               max = geo_all %$% max(epi_date) - 7,
                                               value = c(
                                                   geo_all %$% max(epi_date)) - 7,
                                               step = 7,
                                               ticks = FALSE,
                                               animate = TRUE,
                                               timeFormat = "%d %b %y"
                                           ),
                                  p("Map showing the geographical distribution of variants, as either number of sequences or percentage relative to a specific region. N.B. Sequences without geographical information have been excluded from this analysis, so overall counts may be slightly lower than reported in VOCs/VUIs in the UK.")
                                           
                                  ) # end of box
                            )

),
tabItem(tabName = "omicron",
         fluidRow(
             box(title = "Possible effect of Omicron against mAbs", closable = FALSE, width = 12,
                 status = "info", collapsible = FALSE, icon = icon("map")
                 ,fluidRow(
                     p("The table shows the fold reduction in neutralisation by monoclonal antibodies for circulating variants and single mutation spike profiles. Comparison between single mutation and full variant data indicates the role of each mutation in the evasion of the full variant. Some suggestion is made of the likely evasion profile of the Omicron (B.1.1.529) variant based on the mutations it contains. However, more definitive conclusions await clinical data and neutralisation assays involving the full virus."),
                     htmlOutput("omicron_table"))
             )
         )
) 
        ) # end map
    ), # end dashboardBody ##8a7967
    
    footer = dashboardFooter(
        left = fluidRow(
            column(2, tags$a(img(src = "MRC.png", alt = "Medical Research Council logo", class = "img-responsive center-block footer-logo"), href="https://mrc.ukri.org/", target = "_blank")), # TODO remove left and right padding
            column(2, tags$a(img(src = "UOG.png", alt = "University of Glasgow logo", class = "img-responsive center-block footer-logo"), href="https://www.gla.ac.uk/", target = "_blank")),
            column(2, tags$a(img(src = "CVR.png", alt = "Centre for Virus Research logo", class = "img-responsive center-block footer-logo"), href="https://www.gla.ac.uk/researchinstitutes/iii/cvr/", target = "_blank")),
            column(2, tags$a(img(src = "cog-uk.svg", alt = "COG-UK Consortium logo", class = "img-responsive center-block footer-logo"), href="https://www.cogconsortium.uk", target = "_blank")),
            column(2, tags$a(img(src = "climb_big_data.png", alt = "CLIMB logo", class = "img-responsive center-block footer-logo"), href="https://www.climb.ac.uk/", target = "_blank")),
            column(2, tags$a(img(src = "University of Cambridge-colour.gif", alt = "University of Cambridge logo", class = "img-responsive center-block footer-logo"), href="https://www.cam.ac.uk"))
        )
    )
 
)
