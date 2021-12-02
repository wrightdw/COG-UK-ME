### Ronopreve upset plot
### Upset plot showing combinations of mutations known to influence mabs in cocktail

library(UpSetR)
library(RColorBrewer)
library(dplyr)
library(magrittr)

generate_upset <- function(filter_date = NULL){
  ## load info on amino acid substitutions
  
  ## Read in mab data which has mutations assoc. fold differences
  mab_data <- read.csv('COG-UK/data_auxiliary/mab_mutations_data.csv')
  mab_data <- mab_data[which(mab_data$Cas > 0 | mab_data$Imd > 0),]
  mab_data <- mab_data[is.na(mab_data$Cas) == F,]
  mab_data <- mab_data[is.na(mab_data$Imd) == F,]
  
  # use mutations df to get to viruses df
  mutations_uk <- 
    mutations_uk %>% 
    filter(gene == 'S')
  
  if(!is.null(filter_date)){
    mutations_uk %<>% 
      filter(sample_date >= filter_date)
  }
  
  mutations_uk <- 
    mutations_uk %>%
    mutate(epi_week = as.numeric(as.character(epi_week)))
  
  mutations_uk <- mutations_uk[,c('sequence_name', 'sample_date', 'epi_week',
                                  'epi_date', 'lineage', 'adm1', 'gene',
                                  'variant')]
  
  viruses <- dplyr::summarise(group_by(mutations_uk, sequence_name, sample_date, epi_week,
                                       epi_date, lineage, adm1),
                              mutations = paste0(variant, collapse = ';'))
  
  # create mab which evaluates presence/absence of each mutation in 'mab_data'
  # in each spike profile in 'viruses'
  mab <- data.frame(sequence_name = viruses$sequence_name)
  for(i in 1:nrow(mab_data)) {
    curr_mutation <- mab_data$set[i]
    new_col <- as.numeric(grepl(curr_mutation, viruses$mutations))
    mab <- cbind(mab, new_col)
    names(mab)[ncol(mab)] <- curr_mutation
  }
  
  # create copy of 'mab_data' that will be read into plot as metadata
  metadata <- mab_data
  names(metadata)[1] <- 'sets'
  
  # create a variable that summarises phenotype
  # phenotype dictates colour of row assoc. mutation
  # shows mab with greater effect and indication of effect size
  metadata$pheno <- NA
  display_threshold <- 100
  for (i in 1:nrow(metadata)) {
    if (metadata$Cas[i] > metadata$Imd[i]) {
      metadata$pheno[i] <- 'Cas_1'
      
      if (metadata$Cas[i] >= display_threshold) {
        metadata$pheno[i] <- 'Cas_2'
      }
      
    } else if (metadata$Cas[i] < metadata$Imd[i]) {
      metadata$pheno[i] <- 'Imd_1'
      
      if (metadata$Imd[i] >= display_threshold) {
        metadata$pheno[i] <- 'Imd_2'
      }
      
    }
  }
  
  ## Use brewer pal to get some colours - upset() matrix has max of 4
  cas_pal <- brewer.pal(6, name = "Blues")[c(3, 4)]
  imd_pal <- brewer.pal(6, name = "Oranges")[c(3, 4)]
  
  ### Make upset plot
  # mb.ratio sets ration in heights of upper histogram & upset plot but problem with
  # overlap of x-axis of historgram when used
  ronapreve_upset <- upset(mab, keep.order = TRUE, nsets = 13, nintersects = 20, #mb.ratio = c(0.35, 0.65), group.by = 'sets',
                           text.scale = 1.4, 
                           point.size = 2.1, 
                           line.size = 1,
                           mainbar.y.label = 'Combination count',
                           sets.x.label = 'Mutation count',
                           
                           set.metadata = list(data = metadata,
                                               plots = list(
                                                 list(type = "matrix_rows",
                                                      column = "pheno",
                                                      colors = c(Cas_1 = cas_pal[1],
                                                                 Cas_2 = cas_pal[2],
                                                                 Imd_1 = imd_pal[1],
                                                                 Imd_2 = imd_pal[2]),
                                                      alpha = 1))))
  
  ronapreve_upset
}

