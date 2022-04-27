library(shiny)
library(plotly)
library(scales)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(ggseqlogo)
library(RColorBrewer)
library(viridis)

## Table functions
# TODO table caching
# 
# Mutations
table_1 <- function(){
  
database_genome %>% 
    select(gene, mutation, `numSeqs UK`, `numSeqs UK 28 days`, `numSeqs Eng 28 days`, `numSeqs Scotland 28 days`, `numSeqs Wales 28 days`, `numSeqs NI 28 days`, earliest) %>%
    # bind_rows(database_deletions) %>% 
    arrange(desc(`numSeqs UK`)) %>% 
    filter(`numSeqs UK` >= 5) %>% 
    mutate(mutation = mutation %>% fct_drop %>% fct_inorder) %>% 
    mutate(`Cumulative sequences in UK (%)` = `numSeqs UK` / total_sequences,
           .after = `numSeqs UK`) %>%
    mutate(`Sequences over the last 28 days in UK (%)` = `numSeqs UK 28 days` / total_sequences_28,
           .after = `numSeqs UK 28 days`) %>%
    rename(Gene = gene,
           `Amino acid replacement/ indel` = mutation, 
           `Cumulative sequences in UK` = `numSeqs UK`, 
           `Sequences over the last 28 days in UK` = `numSeqs UK 28 days`,
           `Sequences over the last 28 days in England` = `numSeqs Eng 28 days`,
           `Sequences over the last 28 days in Scotland` = `numSeqs Scotland 28 days`,
           `Sequences over the last 28 days in Wales` = `numSeqs Wales 28 days`,
           `Sequences over the last 28 days in Northern Ireland` = `numSeqs NI 28 days`,
           `Date of first detection in UK` = earliest) #%T>% print
}

# Variants
# TODO precompute
table_3 <- function(){
  bind_rows(
    n_uk_lineages_all %>% 
      filter(lineage %in% lineages_t3$lineage & variant == "sequences") %>% 
      inner_join(lineages_t3) %>% # join descriptions
      select(-variant) %>% 
      relocate(reason, .after = lineage),
    
    n_uk_lineages_all %>%
      filter(variant == "E484K" & lineage == "B.1.324.1") %>%
      mutate(lineage = str_c(lineage, " + ", variant), .keep = "unused")  %>%
      mutate(reason = "As B.1.324.1, with the addition of E484K.")
  ) %>% 
    lineages_table()
}

table_recombinants <- function(){
  n_uk_recombinants %>% 
    filter(variant == "sequences") %>% 
    left_join(lineages_recomb) %>% # descriptions 
    select(-variant) %>% 
    relocate(reason, .after = lineage) %>% 
    lineages_table
}

# Format lineages table for display
lineages_table <- function(summed_lineages){
  summed_lineages %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0L))) %>% 
    filter(n_sequences_UK > 0) %>%
    relocate(n_sequences_UK, .after = reason) %>% 
    mutate(`UK (%)` = n_sequences_UK / total_sequences,
           .after = n_sequences_UK) %>%
    relocate(n_sequences_28_UK, .after = `UK (%)`) %>% 
    mutate(`UK 28 days (%)` = n_sequences_28_UK / total_sequences_28,
           .after = n_sequences_28_UK) %>%
    arrange(desc(n_sequences_28_UK), desc(n_sequences_UK), lineage) %>% 
    relocate(n_sequences_28_England, .after = n_sequences_England) %>% 
    relocate(n_sequences_28_Northern_Ireland, .after = n_sequences_Northern_Ireland) %>% 
    relocate(n_sequences_28_Scotland, .after = n_sequences_Scotland) %>% 
    relocate(n_sequences_28_Wales, .after = n_sequences_Wales) %>% 
    rename(`Variant` = lineage,	
           UK = n_sequences_UK,	 
           `UK 28 days` = n_sequences_28_UK,                    
           Description = reason, 
           
           England = n_sequences_England,
           `Northern Ireland` = n_sequences_Northern_Ireland,
           Scotland = n_sequences_Scotland,
           Wales = n_sequences_Wales,
           
           `England 28 Days` = n_sequences_28_England,
           `Northern Ireland 28 Days` = n_sequences_28_Northern_Ireland,
           `Scotland 28 Days` = n_sequences_28_Scotland,
           `Wales 28 Days` = n_sequences_28_Wales
    ) 
  
}

# T cell table
table_5 = function(assay_filter = NULL){
  if(!is.null(assay_filter)){
    if(assay_filter == "recognition"){
      database_tcell_predictions %<>% 
        filter(assay %in% c("Reduced T-cell recognition (full)", "Reduced T-cell recognition (partial)"))
    } else { # epitope_studies
      database_tcell_predictions %<>% 
        filter(!assay %in% c("Reduced T-cell recognition (full)", "Reduced T-cell recognition (partial)"))
    }
    
    database_tcell_predictions %<>% mutate(across(c(gene, mutation, Epitope, CD4_CD8, HLA, assay), fct_drop)) # drop unused factor levels
  }
  
  database_tcell_predictions %>% 
    filter(`numSeqs UK` > 0) %>%
    select(gene, mutation, Epitope:Fold, `numSeqs UK`, `numSeqs UK 28 days`, -`End position`) %>%
    arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), desc(Fold)) %>%
    mutate(across(c(Epitope, HLA, assay, CD4_CD8), ~fct_relevel(.x, sort))) %>% 
    rename(`Gene` = gene,
           `Amino acid replacement` = mutation,
           `Cumulative sequences in UK` = `numSeqs UK`,
           `Sequences over 28 days` = `numSeqs UK 28 days`,
           `Assay` = assay,
           `CD4 CD8` = CD4_CD8,
           Start = `Start position`,
           `Mut Percentile Rank Value` = IC50_mutation,
           `WT Percentile Rank Value` = `IC50 WT`,
           `Fold difference` = Fold,
           DOI = doi)
}

table_therapeutics <- function(){
  therapeutics %>% 
    arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`)) %>% 
    select(gene, variant, Protein, mutation_protein, resistance, drug, assay, detail, quantification, note, `numSeqs UK`, `numSeqs UK 28 days`, anchor) %>%
    rename_with(str_to_title, c(gene, resistance, drug, assay, detail, quantification, note)) %>% 
    rename(`Amino acid replacement (polyprotein)` = variant, 
           `Amino acid replacement (protein)` = mutation_protein, 
           `Cumulative sequences in UK` = `numSeqs UK`,
           `Sequences over 28 days` = `numSeqs UK 28 days`,
           `Reference` = anchor
  ) 
}

### spike profiles - non reactive - start
names(spike_tab) <- c('Profile', 'N_change', 'Lineage', 'Count', 'Count_28',
                      'Growth', 'Expansion', 'VOC', 'Geography')
spike_tab$Profile <- gsub(';', ', ', spike_tab$Profile)
spike_tab$Profile <- gsub('Y505H, T547K', 'Y505H,\nT547K', spike_tab$Profile)
spike_tab$Lineage <- gsub('AY.46.5, ', 'AY.46.5,\n', spike_tab$Lineage)
spike_tab$VOC <- gsub('Omicron ', 'Omicron\n', spike_tab$VOC)

##### colour scale stuff
# identify where 0 is in range
col_max <- max(c(0.02, plyr::round_any(max(abs(spike_tab$Expansion)), 0.005, ceiling)))
col_min <- -col_max
zero_percentile <- (0 - col_min) / (col_max - col_min)

col_max_growth <- plyr::round_any(max(spike_tab$Growth), 10, ceiling)
col_min_growth <- plyr::round_any(min(spike_tab$Growth), 10, floor)
zero_percentile_growth <- (0 - col_min_growth) / (col_max_growth - col_min_growth)

#### Re-formatting of spike table for data table
spike_table <- spike_tab
spike_table$Growth <- round(spike_table$Growth, 2)
spike_table$Expansion <- signif(spike_table$Expansion, 4)

names(spike_table)[2] <- 'Amino acid substitutions'
names(spike_table)[3] <- 'Lineage(s)'
names(spike_table)[4] <- 'Sequences 56 days'
names(spike_table)[5] <- 'Sequences 28 days'
names(spike_table)[6] <- 'Average growth rate (%)'
names(spike_table)[7] <- 'Expansion/ contraction'
### spike profiles - non reactive - end

shinyServer(function(input, output, session) {
    waiter_hide() # hide loading screen
  
    output$table_1 <- renderDataTable({
      table_1() %>% 
        datatable(filter = "top", rownames = FALSE, 
                  options = list(lengthMenu = c(20, 50, 100, 200), pageLength = 20, scrollX = TRUE)) %>% 
        formatPercentage(c("Cumulative sequences in UK (%)", "Sequences over the last 28 days in UK (%)"), digits = 2)
    })
    
    ######## Data download inputs ########
    
    ## Mutations
    # Reactive value to generate downloadable table for selected mutation or indel metadata
    # TODO database query
    datasetInput <- reactive({
      mutations_indels_uk_28 %>% 
        filter(gene == input$dataset_gene) %>% 
        filter(variant == input$dataset) %>% 
        select(sequence_name, sample_date, epi_week, epi_date, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete table 1 data
    table1Input <- reactive({
      table_1() %>% 
        mutate(across(ends_with("(%)"), ~.x * 100)) # convert decimal fraction to percentage
    })
    
    ## Variants
    # Reactive value to generate downloadable table for selected lineage + mutation
    # TODO regex switch
    concernInput <- reactive({
      if(input$concern == "B.1.1.7 + E484K"){
        concern_download <- 
          consortium_uk %>% 
          filter(e484k == "K") %>% 
          filter(lineage == "B.1.1.7") 
        # filter(lineage == "B.1.1.7" | str_detect(lineage, sublineage_regex("B.1.1.7")))
      } else if(input$concern == "A.23.1 + E484K") {
        concern_download <- 
          consortium_uk %>% 
          filter(e484k == "K") %>% 
          filter(lineage == "A.23.1") 
          # filter(lineage == "A.23.1" | str_detect(lineage, sublineage_regex("B.1.1.7"))) 
      } else {
        concern_download <- 
          consortium_uk %>% 
          filter(lineage == input$concern)
          # filter(lineage == input$concern | str_detect(lineage, sublineage_regex(input$concern)))
      }
      
      concern_download %>% 
        select(sequence_name, sample_date, epi_week, epi_date, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete table 3 data
    table3Input <- reactive({
      table_3() %>% 
        mutate(across(ends_with("(%)"), ~.x * 100)) # convert decimal fraction to percentage
    })
    
    ## Antigenic mutations
    # Reactive value to generate downloadable table for selected mutation
    # TODO query from database
    escapeInput <- reactive({
      mutations_s_uk %>% 
        filter(variant == input$selectEscape) %>% 
        select(sequence_name, sample_date, epi_week, epi_date, lineage) %>% 
        arrange(desc(sample_date), lineage)
    })
    
    # Reactive value to generate downloadable table for complete T cell table data
    table5Input <- reactive({
      table_5() 
    })
    
    ######## Download handlers ########
    
    ## Table 1
    # Downloadable CSV of selected mutation metadata
    output$downloadData <- downloadHandler(
      filename = function() {
        str_c(input$dataset_gene, "_", input$dataset, "_UK_28_days_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(datasetInput(), file)
      },
      contentType = "text/csv"
    )
    
    # Downloadable CSV of complete table 1 data
    output$downloadTable1 <- downloadHandler(
      filename = function() {
        str_c("table_1_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(table1Input(), file)
      },
      contentType = "text/csv"
    )
    
    ## Table 3
    # Downloadable CSV of selected mutation
    output$downloadConcern <- downloadHandler(
      filename = function() {
        str_c(input$concern, "_UK_cumulative_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(concernInput(), file)
      },
      contentType = "text/csv"
    )
    
    # Downloadable CSV of complete table 1 data
    output$downloadTable3 <- downloadHandler(
      filename = function() {
        str_c("table_3_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(table3Input(), file)
      },
      contentType = "text/csv"
    )
    
    ## Antigenic
    # Downloadable CSV of selected mutation metadata
    output$downloadEscape <- downloadHandler(
      filename = function() {
        str_c(input$selectEscape, "_UK_cumulative_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(escapeInput(), file)
      },
      contentType = "text/csv"
    )
    
    # Downloadable CSV of complete T cell data
    output$downloadTable5 <- downloadHandler(
      filename = function() {
        str_c("tcell_", dataset_date, ".csv")
      },
      content = function(file) {
        write_csv(table5Input(), file)
      },
      contentType = "text/csv"
    )
    
    # Render VOC/VUI table
    output$table_3 <- renderDT({
      table_3() %>% 
        datatable(filter = "none", escape = FALSE, rownames = FALSE, 
                  options = list(dom = 't', paging = FALSE, scrollX = TRUE)) %>% 
        formatPercentage(c("UK (%)", "UK 28 days (%)"), digits = 2)
    })
    
    # Render recombinants table
    output$table_recomb <- renderDT({
      table_recombinants() %>% 
        datatable(filter = "none", escape = FALSE, rownames = FALSE, 
                  options = list(dom = 't', paging = FALSE, scrollX = TRUE)) %>% 
        formatPercentage(c("UK (%)", "UK 28 days (%)"), digits = 2)
    })
    
    # Antigenic Mutations
    output$table_4 <- renderDT({
      if("monoclonal" %in% input$escape){
        database %<>% filter(mab == TRUE)
      }
      
      if("convalescent" %in% input$escape){
        database %<>% filter(plasma == TRUE)
      }
      
      if("vaccine" %in% input$escape){
        database %<>% filter(vaccine_sera == TRUE)
      }
      
      database %>%
        filter(!is.na(escape)) %>% 
          mutate(mutation = fct_drop(mutation)) %>%
          select(mutation, domain, escape, `numSeqs UK`, `numSeqs UK 28 days`, support, anchor) %>%  
          arrange(desc(`numSeqs UK`), desc(`numSeqs UK 28 days`), mutation) %>% 
          rename(`Amino acid replacement` = mutation, 
                 `Cumulative sequences in UK` = `numSeqs UK`,
                 `Sequences over 28 days` = `numSeqs UK 28 days`,
                 `Escape mutations details` = escape,
                 `References` = anchor,
                 Confidence = support, 
                 Domain = domain) %>% 
        datatable(filter = "top", escape = FALSE, rownames = FALSE,
                  options = list(lengthMenu = c(20, 50, 100, 200), pageLength = 20, scrollX = TRUE)) %>% 
        formatStyle(
          'Confidence',
          target = 'row',
          backgroundColor = styleEqual(c("lower", "medium", "high"), c('LemonChiffon', 'DarkOrange', 'FireBrick')), 
          color = styleEqual(c("lower", "medium", "high"), c('#737272', '#3e3d3d', "#d0cfcf"))) # contrasting colours generated by https://leonardocolor.io/
        })
    
    output$table_5 <- renderDT({
      table_5(input$t_cell_experiment) %>% 
        mutate(Reference = str_c("<a href='", DOI, "'target='_blank'>", Reference,"</a>"), .keep = "unused", .after = `Supporting references`) %>% # hyperlink to citation DOI
        datatable(filter = "top", escape = FALSE, rownames = FALSE,
                    options = list(lengthMenu = c(10, 20, 50, 100, 200), pageLength = 10, scrollX = TRUE))
    })
    
    output$table_therapeutics <- renderDT({
      table_therapeutics() %>% 
      datatable(filter = "none", escape = FALSE, rownames = FALSE, 
                options = list(dom = 't', paging = FALSE, scrollX = TRUE)) 
        
    })
    
    # always display wild type on percentage chart
    observeEvent(input$percentage, {
      if(input$percentage){
        disable("ref")
        updatePrettySwitch(
          session = session,
          inputId = "ref",
          value = as.logical(input$percentage)
        )
      } else {
        enable("ref")
      }
    })

    # Filter mutation and deletion data and create plot
    mutation_plot <- reactive({
      # deletions or mutations only - lump the other category together
      
      mutation_reference_counts %<>%  
        filter(gene == input$gene & position == input$position) 
       
      if(input$mutation_type == "replacement"){ # amino acid replacement
        # replacements with n < 5 over all time in UK
        reps_low <- 
          mutation_reference_counts %>% 
          select(-epi_week, -epi_date, -adm1) %>% 
          filter(variant != "WT" & !str_starts(variant, "del")) %>% 
          group_by(across(-n)) %>% 
          summarise(n = sum(n), .groups = "drop") %>% 
          filter(n < 5) %$% 
          variant %>% 
          as.character %T>% print
        
        # substitute combined deletions by date
        mutation_reference_counts %<>%  
          mutate(variant = fct_collapse(variant, 
                                        "Deletions" = unique(grep("^del", # collapse deletions only
                                                                  variant, value = TRUE)))) %>% 
          mutate(variant = fct_other(variant, drop = reps_low, other_level = "Replacements (n<5)")) %>% # lump low-level replacements
          group_by(across(-n)) %>% 
          summarise(n = sum(n), .groups = "drop")         
      } else { # deletion
        # get deletions with n < 5 in UK over all time
        dels_low <- 
          mutation_reference_counts %>% 
          select(-epi_week, -epi_date, -adm1) %>% 
          filter(str_starts(variant, "del")) %>% 
          group_by(across(-n)) %>% 
          summarise(n = sum(n), .groups = "drop") %>% 
          filter(n < 5) %$% 
          variant %>% 
          as.character
          
        mutation_reference_counts %<>%  
          mutate(variant = fct_collapse(variant, 
                                        "Replacements" = unique(grep("(^del)|(^WT$)", # collapse non-deletions and non-WT
                                                                     variant, value = TRUE, invert = TRUE)))) %>% 
          mutate(variant = fct_other(variant, drop = dels_low, other_level = "Deletions (n<5)")) %>% # lump low-level deletions
          group_by(across(-n)) %>% 
          summarise(n = sum(n), .groups = "drop")        
      }
      
      # calculate variants ordered by frequency over all time in UK
      variants_by_frequency <- 
        mutation_reference_counts %>% 
        select(-epi_week, -epi_date, -adm1) %>% 
        group_by(across(-n)) %>% 
        summarise(n = sum(n), .groups = "drop") %>%
        arrange(desc(n)) %$% 
        variant %>% 
        as.character
                  
      mutation_reference_counts %<>%   
        mutate(variant = 
                 variant %>% 
                 fct_drop %>%
                 fct_relevel(variants_by_frequency) %>% # order colours by frequency over all time
                 fct_relevel("WT", after = 0) # then move WT to always be first colour
               ) %>%   
        filter(adm1 == input$nation) # keep colours the same when switching between nations
        
      variants <- mutation_reference_counts %$% levels(variant)
      
      if(!input$ref){
        mutation_reference_counts %<>% filter(variant != "WT") 
      }
      
      colour_order = match(mutation_reference_counts %$% fct_drop(variant) %>% levels, variants)
      
      mutation_reference_plot <- 
        mutation_reference_counts %>% 
        filter(epi_date >= input$mutation_range[1] & epi_date <= input$mutation_range[2]) %>% 
        rename(`Sample date` = epi_date, Sequences = n, Mutation = variant) %>% # display names
        ggplot(aes(fill = Mutation, y = Sequences, x = `Sample date`) ) +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(title = str_c(c("Gene", "Position"), c(input$gene, input$position), sep = " : ", collapse = "\n")) +
        scale_fill_discrete_qualitative(palette = "Set 3", 
                                        nmax = mutation_reference_counts %$% levels(variant) %>% length,
                                        order = colour_order,
                                        rev = FALSE) +
        scale_x_date(breaks = date_breaks("2 month"),
                     labels = date_format("%b %y"))
  
      mutation_reference_plot
    }) 
    
    # Mutation plot by percentage or count
    mutation_plot_bar <- reactive({
      max_date <- mutation_reference_counts %$% max(epi_date)
      gg_bar <- mutation_plot()
      
      # upper range slider is in most recent 2 weeks
      if (input$mutation_range[2] >= max_date - days(7)){
        if(input$mutation_range[2] >= max_date ){
          xmax = max_date + days(3)
        } else{
          xmax = max_date + days(-4)
        }
        
        if(input$percentage){
          ymax <- 1
        } else {
          # TODO cache maximum weekly counts for UK and nations
          if(input$nation == "UK"){
            ymax <-
              consortium_uk %>%
              dplyr::count(epi_date) %$%
              max(n)
          } else {
            ymax <-
              consortium_uk %>%
              filter(adm1 == input$nation) %>% 
              dplyr::count(epi_date) %$%
              max(n)
          }
        }
        
        gg_bar <- 
          gg_bar + 
          annotate("rect",
                   xmax = xmax,
                   xmin = max_date - days(10),
                   ymin = 0,
                   ymax = ymax,
                   alpha = 0.2)
      }
      
      if(input$percentage){
        gg_bar <- 
          gg_bar +
          geom_bar(position="fill", stat="identity") +
          scale_y_continuous(labels = scales::percent_format())
      } else {
        gg_bar <- 
          gg_bar + 
          geom_bar(position="stack", stat="identity")
      }
      
      gg_bar
    }) %>% debounce(500) # allow 500ms to update percentage switch so don't display plot immediately
    
    # Display mutation plot with percentage option
    output$mutation_time <- renderPlotly({
      mutation_plot_bar() %>% ggplotly
    })
    
    values <- reactiveValues()

    observe({
      defining <- 
        vui_voc %>% 
        filter(lineage == input$lineage_antigenic) %$% 
        mutation
      
      antigenic_mutations <- antigenic_mutations_lineages(nation = input$nation_antigenic, lineage = input$lineage_antigenic, defining = defining, input$percentage_range)

      if(is.null(antigenic_mutations)){
        values$antigenic <- NULL
        values$antigenic_title <- str_c(input$lineage_antigenic, " (", input$nation_antigenic %>% str_replace_all("_", " "), ")", ": no antigenic mutations")
      } else {
        values$antigenic_title <- str_c("Antigenic mutations on the top of ", input$lineage_antigenic, " defining mutations (", 
                                        input$nation_antigenic %>% str_replace_all("_", " "), ")")
        values$antigenic <- 
          antibody_complex_heatmap(antigenic_mutations, input$percentage_range)
      }
    })
    
    output$title_heatmap <- renderText({
      values$antigenic_title
    })
    
    # Display antibody heatmap
    # TODO plot caching
    output$antibody_heatmap <- renderPlot({
      if(is.null(values$antigenic)){
        ""  
      } else {
        draw(values$antigenic) 
      }
    }, height = function(){
      if(is.null(values$antigenic)){
        1
      } else {
        values$antigenic %>% nrow * 13 + 100
      }
    })
    
    # position input for mutations/indels chart
    observeEvent(input$gene, {
      updateSelectInput(session, "position",
                        choices = mutation_reference_counts %>%
                          filter(gene == input$gene) %>%
                          distinct(position) %>%
                          arrange(position))
    })
    
    # variant select for mutations/indels download
    observeEvent(input$dataset_gene, {
      updateSelectizeInput(session, 
                           "dataset",
                           choices = database_genome %>% 
                             filter(gene == input$dataset_gene) %>% 
                             filter(`numSeqs UK 28 days` > 0) %>% 
                             distinct(mutation))
    })
    
    
    observeEvent(input$dataset,{
      toggleState("downloadData", input$dataset != "")
    })
    
    output$epitope_sequence <- renderPlot({
      if(!input$epitope_ref){
        wt %<>% mutate(`numSeqs UK` = 0)
      }
      
      if(input$method){
        method = "prob"
      } else {
        method = "bits"
      }
      
      database_logo <- 
        database %>% 
        mutate(WT = str_sub(mutation, 1, 1), AA = str_sub(mutation, -1), .after = mutation) %>% 
        filter(`numSeqs UK` > 0) %>% 
        select(position, AA, `numSeqs UK`) %>% 
        bind_rows(wt) %>%
        pivot_wider(names_from = position, values_from = `numSeqs UK`, values_fill = 0, names_sort = TRUE) %>%
        arrange(AA) %>% 
        column_to_rownames(var = "AA") %>% 
        as.matrix
      
      epitopes_positions <- 
        database_tcell_predictions %>% 
        filter(`Start position` <= input$epitope_position & `End position` >= input$epitope_position) %>% 
        distinct(Epitope, `Start position`, `End position`)
      
      if(nrow(epitopes_positions) > 0){
        epitopes_positions %>% 
          pmap(function(Epitope, `Start position`, `End position`) {
            database_logo[,`Start position`:`End position`]
          }) %>% purrr::set_names(epitopes_positions %$% str_c(Epitope, " ", `Start position`, ":", `End position`)) %>% 
        ggseqlogo(method = method, ncol = 2) # TODO separate reactive for method to update plot only
      } else {
        ""
      }
    })
    
    variant_plot <- reactive({
      # only generate plot if variant selected and exclude other switch is off
      if(is.null(input$variant_vui_voc) && input$other_switch){
        ggplot() + theme_void()
      } else {
      selected_variants <- input$variant_vui_voc
      vui_voc_lineages <- 
        vui_voc_lineages %>% 
        append("Other Delta", after = which(. == "B.1.617.2")[1]) %>% 
        append("Other")
      
      if( "B.1.617.2" %in% input$variant_vui_voc && input$variant_delta != "B.1.617.2"){
        selected_variants <- replace(selected_variants, selected_variants == "B.1.617.2", input$variant_delta)
        vui_voc_lineages <- replace(vui_voc_lineages, vui_voc_lineages == "B.1.617.2", input$variant_delta)
      }
      
      if(input$variant_day){
                  lineages_days_uk_all %<>% 
            filter(adm1 == input$nations_vui_voc)
        
        if( "B.1.617.2" %in% input$variant_vui_voc){
          # remove unselected Deltas so they aren't counted in Other
          delta_options <- c("B.1.617.2", "AY.4", "AY.4.2")
          delta_options <- delta_options[delta_options != input$variant_delta]
          
          lineages_days_uk_all %<>% 
            filter(!(lineage %in% delta_options)) 
          
          # filter out Delta minus counts
          if(input$variant_delta == "AY.4"){
            lineages_days_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4.2")
            
            # add Other Delta after AY.4 in colour ordering
            ay_index <- match("AY.4", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4", after = ay_index) # don't include in Other
          } else if(input$variant_delta == "AY.4.2"){
            lineages_days_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4")
            
            # add Other Delta after AY.4.2 in colour ordering
            ay_index <- match("AY.4.2", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4.2", after = ay_index) # don't include in Other
          } else {
            lineages_days_uk_all %<>% 
              filter(!str_starts(lineage, fixed("Delta_minus_")))
          }
        } else { # Delta not selected, count only B.1.617.2 in Other
          lineages_days_uk_all %<>% 
            filter(!str_starts(lineage, fixed("AY."))) %>% 
            filter(!str_starts(lineage, fixed("Delta_minus_")))
        } 
        
        variants_other_day <- 
          lineages_days_uk_all %>% 
          filter(!(lineage %in% selected_variants)) %>% 
          group_by(sample_date) %>%  
          summarise(n_day = sum(n_day)) %>% 
          mutate(lineage = "Other", .before = sample_date) %>% 
          ungroup
        
        # TODO legend names from vui_voc$lineage_display
        lineages_days_uk <- 
          lineages_days_uk_all %>% 
          filter(lineage %in% selected_variants) %>% 
          bind_rows(variants_other_day) %>% 
          mutate(lineage = recode(lineage,
                                  "Delta_minus_AY.4.2" = "Other Delta",
                                  "Delta_minus_AY.4" = "Other Delta")) %>%
          mutate(lineage = recode_factor(lineage, # recode WHO Greek display names as factor and order levels to define colour/legend order
                                  "BA.1" = "BA.1/BA.1.x (Omicron)",
                                  "BA.2" = "BA.2/BA.2.x (Omicron)",
                                  "B.1.1.7" = "B.1.1.7 (Alpha)",
                                  "B.1.351" = "B.1.351 (Beta)",
                                  "B.1.617.2" = "B.1.617.2/AY.x (Delta)",
                                  "Other Delta" = "Other Delta",
                                  "AY.4" = "AY.4/AY.4.x (Delta)",
                                  "AY.4.2" = "AY.4.2/AY.4.2.x (Delta)",
                                  "P.1" = "P.1 (Gamma)",
                                  "Other" = "Other"
                                  )) %>% 
          rename(Variant = lineage, `Sample date` = sample_date, Sequences = n_day)
        
        selected_variants <- replace(selected_variants, selected_variants %in% c("Delta_minus_AY.4.2", "Delta_minus_AY.4"), "Other Delta") 
        selected_variants <- replace(selected_variants, selected_variants == "B.1.617.2", input$variant_delta)
        
        if(input$other_switch){
          lineages_days_uk %<>% 
            filter(`Sample date` >= input$variant_range[1] & `Sample date` <= input$variant_range[2] & Variant !="Other") 
        } else {
          lineages_days_uk %<>% 
            filter(`Sample date` >= input$variant_range[1] & `Sample date` <= input$variant_range[2]) 
        }
        
        vui_voc_plot <- 
          lineages_days_uk %>% 
          filter(`Sample date` >= input$variant_range[1] & `Sample date` <= input$variant_range[2]) %>% 
          ggplot(aes(fill = Variant, y = Sequences, x = `Sample date`) ) +
          theme_classic() +
          
          scale_fill_discrete_qualitative(palette = "Dynamic", 
                                          nmax = vui_voc_lineages %>% length , # extra colour for Other
                                          
                                          # fix variant/colour combos plus extra colour for Other
                                          order = match(selected_variants, 
                                                        vui_voc_lineages) %>% c(vui_voc_lineages %>% length)
          ) +
          
          scale_x_date(breaks = date_breaks("2 month"),
                       labels = date_format("%b %y")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          labs(x = "Sample date",
               y = "Sequences"
          ) 

        # display annotation from Sunday on 2nd last epiweek for consistency with weeks plot
        max_epi_date <- lineages_weeks_uk_all %$% max(epi_date) 
        
        # Set height of y-axis and annotation box according to highest weekly total sequences 
        if(input$variant_percentage){
          
          ymax <- 1
        } else {
          ymax <- 
            ymax <- 
            lineages_days_uk %>% 
            group_by(`Sample date`) %>% 
            summarise(total_week = sum(Sequences)) %$% 
            max(total_week)
        }
        
        # display annotation if upper slider is in latest 2 weeks
        if (input$variant_range[2] >= max_epi_date - days(7)){
          vui_voc_plot <- 
            vui_voc_plot +
            annotate("rect", 
                     xmax = input$variant_range[2] + days(1), 
                     xmin = max_epi_date - days(8), 
                     ymin = 0, 
                     ymax = ymax,
                     alpha = 0.2)
        }
        
        if(input$variant_percentage){
          vui_voc_plot <- 
            vui_voc_plot +
            geom_bar(position="fill", stat="identity", width = 1) +
            scale_y_continuous(labels = scales::percent_format())
        } else {
          vui_voc_plot <- 
            vui_voc_plot + 
            geom_bar(position="stack", stat="identity", width = 1) +
            ylim(0, ymax)
        }
        
        vui_voc_plot
      } else {  # variants by week
        
          lineages_weeks_uk_all %<>% 
            filter(adm1 == input$nations_vui_voc)
        
        if( "B.1.617.2" %in% input$variant_vui_voc){
          # remove unselected Deltas so they aren't counted in Other
          delta_options <- c("B.1.617.2", "AY.4", "AY.4.2")
          delta_options <- delta_options[delta_options != input$variant_delta]
          
          lineages_weeks_uk_all %<>%
            filter(!(lineage %in% delta_options)) 
          
          # filter out Delta minus counts
          if(input$variant_delta == "AY.4"){
            lineages_weeks_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4.2")
            
            # add Other Delta after AY.4 in colour ordering
            ay_index <- match("AY.4", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4", after = ay_index) # don't include in Other
          } else if(input$variant_delta == "AY.4.2"){
            lineages_weeks_uk_all %<>% 
              filter(lineage != "Delta_minus_AY.4")
            
            # add Other Delta after AY.4.2 in colour ordering
            ay_index <- match("AY.4.2", selected_variants) 
            selected_variants %<>% append("Delta_minus_AY.4.2", after = ay_index) # don't include in Other
          } else {
            lineages_weeks_uk_all %<>% 
              filter(!str_starts(lineage, fixed("Delta_minus_")))
          }
        } else { # Delta not selected, count only B.1.617.2 in Other
          lineages_weeks_uk_all %<>% 
            filter(!str_starts(lineage, fixed("AY."))) %>% 
            filter(!str_starts(lineage, fixed("Delta_minus_")))
        } 
        
        variants_other_week <- 
          lineages_weeks_uk_all %>% 
          filter(!(lineage %in% selected_variants)) %>% 
          group_by(epi_date) %>% 
          summarise(n_week = sum(n_week)) %>% 
          mutate(lineage = "Other", .before = epi_date) %>% 
          ungroup
        
        # TODO legend names from vui_voc$lineage_display
        lineages_weeks_uk <- 
          lineages_weeks_uk_all %>% 
          filter(lineage %in% selected_variants) %>% 
          bind_rows(variants_other_week) %>% 
          mutate(lineage = recode(lineage,
                                  "Delta_minus_AY.4.2" = "Other Delta",
                                  "Delta_minus_AY.4" = "Other Delta")) %>%
          mutate(lineage = recode_factor(lineage, # recode WHO Greek display names as factor and order levels to define colour/legend order
                                         "BA.1" = "BA.1/BA.1.x (Omicron)",
                                         "BA.2" = "BA.2/BA.2.x (Omicron)",
                                         "B.1.1.7" = "B.1.1.7 (Alpha)",
                                         "B.1.351" = "B.1.351 (Beta)",
                                         "B.1.617.2" = "B.1.617.2/AY.x (Delta)",
                                         "Other Delta" = "Other Delta",
                                         "AY.4" = "AY.4/AY.4.x (Delta)",
                                         "AY.4.2" = "AY.4.2/AY.4.2.x (Delta)",
                                         "P.1" = "P.1 (Gamma)",
                                         "Other" = "Other"
          )) %>% 
          rename(Variant = lineage, `Start date` = epi_date, Sequences = n_week)
        
        selected_variants <- replace(selected_variants, selected_variants %in% c("Delta_minus_AY.4.2", "Delta_minus_AY.4"), "Other Delta") 
        selected_variants <- replace(selected_variants, selected_variants == "B.1.617.2", input$variant_delta)
        
        if(input$other_switch){
          lineages_weeks_uk %<>% 
            filter(`Start date` >= input$variant_range[1] & `Start date` <= input$variant_range[2] & Variant !="Other")
        } else {
          lineages_weeks_uk %<>% 
            filter(`Start date` >= input$variant_range[1] & `Start date` <= input$variant_range[2])
        }
        
        vui_voc_plot <- 
          lineages_weeks_uk %>% 
          filter(`Start date` >= input$variant_range[1] & `Start date` <= input$variant_range[2]) %>% 
          ggplot(aes(fill = Variant, y = Sequences, x = `Start date`) ) +
          theme_classic() +
          
          scale_fill_discrete_qualitative(palette = "Dynamic", 
                                          nmax = vui_voc_lineages %>% length , # extra colour for Other
                                          
                                          # fix variant/colour combos plus extra colour for Other
                                          order = match(selected_variants, 
                                                        vui_voc_lineages) %>% c(vui_voc_lineages %>% length) 
          ) +
          
          scale_x_date(breaks = date_breaks("2 month"),
                       labels = date_format("%b %y")) +
          theme(plot.title = element_text(hjust = 0.5)) +
          labs(x = "Sample date",
               y = "Sequences"
          ) 
          
          max_date <- lineages_weeks_uk %$% max(`Start date`)
        
        # Set height of annotation box according to highest weekly total sequences
        if(input$variant_percentage){
          ymax <- 1
        } else {
          ymax <- 
            ymax <- 
            lineages_weeks_uk %>% 
            group_by(`Start date`) %>% 
            summarise(total_week = sum(Sequences)) %$% 
            max(total_week)
        }
        
        # display annotation if upper slider is in latest 2 weeks
        if (input$variant_range[2] >= max_date - days(7)){
          if(input$variant_range[2] >= max_date ){
            xmax = max_date + days(3) 
          } else{
            xmax = max_date + days(-4) 
          }
          
          vui_voc_plot <- 
            vui_voc_plot +
            annotate("rect", 
                     xmax = xmax, 
                     xmin = max_date - days(10), 
                     ymin = 0, 
                     ymax = ymax,
                     alpha = 0.2)
        }
        
        # TODO position_nudge
        if(input$variant_percentage){
          vui_voc_plot <- 
            vui_voc_plot +
            geom_bar(position="fill", stat="identity") +
            scale_y_continuous(labels = scales::percent_format())
        } else {
          vui_voc_plot <- 
            vui_voc_plot + 
            geom_bar(position="stack", stat="identity") +
            ylim(0, ymax)
        }
        vui_voc_plot
      } # end else by day
      vui_voc_plot
      
      } # end if other switch
    }) %>% debounce(500)
    
    output$variant_time <- renderPlotly({
      variant_plot() %>% ggplotly
    })
    
    # if no variants selected, disable and turn off Exclude Other
    # else variants are selected, enable
    observeEvent(input$variant_vui_voc, {
      state <- !is.null(input$variant_vui_voc)
      toggleState("other_switch", condition = state)
      if(!state){
        updatePrettySwitch(
          session = session,
          inputId = "other_switch",
          value = state
        )
      }
    }, ignoreNULL = FALSE)
    
    # Observer for variant Percentage switch to enable/disable and change state of Exclude Other switch
    observeEvent(input$variant_percentage, {
      
      # if percentage on, OR if percentage off and no variants selected, disable Exlude Other
      # else if percentage off and variants selected, enable Exclude Other
      # toggleState("other_switch", condition = !input$variant_percentage)
      if(input$variant_percentage || (!input$variant_percentage && is.null(input$variant_vui_voc)) ){
        disable("other_switch")
      } else {
        enable("other_switch")
      }
      
      # percentage off, also turn off Exclude Other
      if(input$variant_percentage){
        updatePrettySwitch(
          session = session,
          inputId = "other_switch",
          value = FALSE
        )
      }
    })
    
    ########### Ronapreve plot
    # always display wild type on percentage chart
    observeEvent(input$ronapreve_28, {
      css_class <- "center-block img-responsive"
      if(input$ronapreve_28){
        x<-paste("28 days to latest UK sequence date (# sequences:", total_sequences_28,")")
        output$title_ronapreve <- renderText(x
          )
        output$ronapreve_plot <- renderImage({
          list(src = str_c(dataset_date, "/Ronapreve_28.png"),
               alt = "Ronapreve plot 28 days",
               class = css_class
          )
        }, deleteFile = FALSE)
      } else {
        output$title_ronapreve <- renderText("All time")
        output$ronapreve_plot <- renderImage({
          list(src = str_c(dataset_date, "/Ronapreve.png"),
               alt = "Ronapreve plot all time",
               class = css_class
          )
        }, deleteFile = FALSE)
      }
    })
    
    ########### Map and geographical distribution of variants
    
    map_weekInput <- reactive({
      # if((input$antigenic_mutation) == ""){ # if input is empty show only lineages
      
        geo_all_1<- geo_all %>% filter(lineage == input$variant_map)
        
        max_count<- max(geo_all_1$Count)
        max_proportion<- max(geo_all_1$Proportion)
        geo_all_1<-geo_all_1 %>% filter(epi_date == input$variant_date)
        
        if(input$percentage_map == TRUE){ 
          geo_all_1<-dplyr::rename(geo_all_1, "value" = "Proportion")
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "Count")]
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "objectid")]
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "epi_week")]
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "lineage")]
          
          #Join mydata with mapdata
          df <- plyr::join(mapdata, geo_all_1, by= c("NUTS1"))
          c<-"Percentage"
          max_val<-max_proportion
        } else {
          geo_all_1<-dplyr::rename(geo_all_1, "value" = "Count")
          c <-"Number of sequences"
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "objectid")]
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "epi_week")]
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "lineage")]
          geo_all_1<- geo_all_1[, -which(names(geo_all_1) == "Proportion")]
          
          #Join mydata with mapdata
          df <- plyr::join(mapdata, geo_all_1, by= c("NUTS1"))
          max_val <- max_count
        }

      # generate plot
      gg <- ggplot() + geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = value), color = "#FFFFFF", size = 0.25)
      gg <- gg + scale_fill_gradient2(low = "blue", mid = "red", high = "yellow", na.value = "white")
      gg <- gg + coord_fixed(1) # TODO conic conformal map projection?
      gg <- gg + theme_minimal()
      gg <- gg + scale_fill_viridis(name= c, direction = -1, limits = c(0,max_val)) 
      gg <- gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'right')
      gg <- gg + theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
      gg <- gg + theme(axis.title.y=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
      gg
    }) 
    
    output$map <- renderPlotly({
      y <- map_weekInput() 
      ggplotly(y + theme(legend.position = 'bottom')) %>% 
        layout(legend = list(orientation = "h", x = -0.5, y =-1))
    })
  
  output$omicron_table <- renderUI({
      tags$embed(seamless="seamless", 
                  src= "Omicron_5.htm",
                  width=1600, 
                  height=1200)
    })
  
  ### 'Spike profiles' tab - outputs - start
  # Subset spike profiles for plot based on user selected geography (input$spike_geo)
  spike_tab_subset <- reactive({
    spike_tab_sub <- spike_tab %>%
      filter(Geography == input$spike_geo)
    spike_tab_sub <- spike_tab_sub[order(abs(spike_tab_sub$Expansion), decreasing = F),]
    as.data.frame(spike_tab_sub)
    
  })
  
  
  output$spikePlot_count28 <- renderPlotly({
    
    p <- ggplot(data = spike_tab_subset(),
                aes(text = paste('</br>Profile: ', Profile,
                                 '</br>Count 28 days = ', Count_28,
                                 '</br>Lineage(s) = ', Lineage))) +
      geom_point(aes_string(x = "Count_28", y = input$spike_y, col = input$spike_col)) +
      # scale_y_continuous(breaks = seq(-20, 20, ifelse(col_max - col_min > 2, 0.2, 0.1))) +
      scale_x_continuous(trans = 'log', breaks = c(1, 10, 100, 1000, 10000, 100000),
                         labels = c("1", "10", "100", "1,000", "10k", "100k")) +
      xlab('Count of sequences in latest 28-day period') +
      theme_minimal() +
      theme(text = element_text(size = 13))
    
    ## Depending on input$spike_y, need different scale and label
    if (input$spike_y == "Expansion") {
      p <- p +
        scale_y_continuous(breaks = seq(-20, 20, ifelse(col_max - col_min > 2, 0.2, 0.1)),
                           name = "Expansion/contraction")
    } else if (input$spike_y == "Growth") {
      p <- p +
        scale_y_continuous(name = "Average growth rate (%)")
    } else if (input$spike_y == "N_change") {
      p <- p +
        scale_y_continuous(breaks = seq(0, 100, 2),
                           name = "Amino acid substitutions")
    }
    
    ## Depending on input$spike_x, need different colour scale
    if (input$spike_col == "Expansion") {
      p <- p +
        scale_color_gradientn(limits = c(col_min, col_max),
                              values = c(0, zero_percentile - 0.02, zero_percentile, zero_percentile + 0.02, 1),
                              colours = c('dodgerblue4', 'lightblue1', 'grey90', 'coral', 'firebrick4'),
                              name = 'Expansion/\nContraction')
      
    } else if (input$spike_col == "Growth") {
      p <- p + 
        scale_color_gradientn(limits = c(col_min_growth, col_max_growth),
                              values = c(0, zero_percentile_growth - 2, zero_percentile_growth, zero_percentile_growth + 0.02, 1),
                              colours = c('dodgerblue4', 'lightblue1', 'grey90', 'coral', 'firebrick4'),
                              name = 'Growth\nrate')
    } else if (input$spike_col == "VOC") {
      p <- p +
        scale_color_brewer(palette = 'Set2', name = 'Variant')
    }
    
    ggplotly(p,
             tooltip = c("text"), textposition = 'topright')
  })
  
  
  # Subset spike profiles for table based on user selected geography (input$spike_geo)
  spike_table_subset <- reactive({
    spike_table_sub <- spike_table %>%
      filter(Geography == input$spike_geo)
    spike_table_sub <- spike_table_sub %>% select(-c("Geography", "VOC"))
    spike_table_sub <- spike_table_sub[order(spike_table_sub$`Expansion/ contraction`, decreasing = T),]
    as.data.frame(spike_table_sub)
    
  })
  
  output$spikeTable <- renderDataTable({
    datatable(spike_table_subset(),
              options = list(lengthMenu = c(20, 50, 100, 200)), rownames = FALSE) 
  })
  ### 'Spike profiles' tab - outputs - end
})
