## Function: take mutations_uk format, filter to spike and get one row per virus
virus_profiles <- function(mutation_df = NA, filter_date = NULL) {
  mutation_df <- 
    mutation_df %>% 
    filter(gene == 'S')
  
  if(!is.null(filter_date)){
    mutation_df %<>% 
      filter(sample_date >= filter_date)
  }
  
  mutation_df <- 
    mutation_df %>%
    mutate(epi_week = as.numeric(as.character(epi_week)))
  
  mutation_df <- mutation_df[,c('sequence_name', 'sample_date', 'epi_week',
                                'epi_date', 'lineage', 'adm1', 'gene',
                                'variant')]
  
  viruses <- dplyr::summarise(group_by(mutation_df, sequence_name, sample_date, epi_week,
                                       epi_date, lineage, adm1),
                              mutations = paste0(variant, collapse = ';'))
  viruses
}

## Function: take virus_df, a dataframe with list of mutations separated by ';' and sample dates
spike_profiles <- function(viruses = NA, spike_period = 7,
                           do_pretty_profiles = T) {
  
  if (do_pretty_profiles == T) {
    viruses <- viruses[viruses$sample_date > (max(viruses$sample_date) - 56), ]
    profile_lineage_combos <- dplyr::distinct(viruses[, c('lineage', 'mutations')])
    profile_lineage_combos <- get_pretty_profiles(x = profile_lineage_combos)
    # Sort to place 'lineage == 'None' later in order so they'll be removed by duplicated()
    # when they also show up linked to a 
    profile_lineage_combos$sort_by <- ifelse(profile_lineage_combos$lineage == 'None', 2, 1)
    profile_lineage_combos <- profile_lineage_combos[order(profile_lineage_combos$sort_by),]
    profile_lineage_combos <- profile_lineage_combos[!duplicated(profile_lineage_combos$mutations),]
    profile_lineage_combos <- dplyr::distinct(profile_lineage_combos[, c('mutations',
                                                                         'pretty_profile',
                                                                         'N_change')])
    viruses <- dplyr::left_join(viruses, profile_lineage_combos, by = 'mutations')
    viruses$mutations <- viruses$pretty_profile
  }
  
  spike_recent <- unique(viruses$mutations[viruses$sample_date > max(viruses$sample_date) - spike_period])
  
  # Get total N in each latest and previous 28 day periods, fortnights, and 56 days
  N_28 <- sum(viruses$sample_date > max(viruses$sample_date) - 28)
  N_28_prev <- sum(viruses$sample_date > max(viruses$sample_date) - 56 &
                     viruses$sample_date <= max(viruses$sample_date) - 28)
  N_f1 <- sum(viruses$sample_date > max(viruses$sample_date) - 14)
  N_f2 <- sum(viruses$sample_date > max(viruses$sample_date) - 28 &
                viruses$sample_date <= max(viruses$sample_date) - 14)
  N_f3 <- sum(viruses$sample_date > max(viruses$sample_date) - 42 &
                viruses$sample_date <= max(viruses$sample_date) - 28)
  N_f4 <- sum(viruses$sample_date > max(viruses$sample_date) - 56 &
                viruses$sample_date <= max(viruses$sample_date) - 42)
  N_56 <- sum(viruses$sample_date > max(viruses$sample_date) - 56)
  
  
  # subsets of viruses by time period
  viruses_28 <- viruses[viruses$sample_date > (max(viruses$sample_date) - 28), ]
  viruses_28p <- viruses[viruses$sample_date > max(viruses$sample_date) - 56 &
                           viruses$sample_date <= max(viruses$sample_date) - 28, ]
  viruses_f1 <- viruses[viruses$sample_date > max(viruses$sample_date) - 14, ]
  viruses_f2 <- viruses[viruses$sample_date > max(viruses$sample_date) - 28 &
                          viruses$sample_date <= max(viruses$sample_date) - 14, ]
  viruses_f3 <- viruses[viruses$sample_date > max(viruses$sample_date) - 42 &
                          viruses$sample_date <= max(viruses$sample_date) - 28, ]
  viruses_f4 <- viruses[viruses$sample_date > max(viruses$sample_date) - 56 &
                          viruses$sample_date <= max(viruses$sample_date) - 42, ]
  viruses_56 <- viruses[viruses$sample_date > (max(viruses$sample_date) - 56), ]
  
  spike_tab <- data.frame(profile = spike_recent,
                          lineage = NA,
                          N = NA,
                          N_28 = NA,
                          N_28p = NA,
                          N_f1 = NA,
                          N_f2 = NA,
                          N_f3 = NA,
                          N_f4 = NA,
                          N_56 = NA,
                          Freq_28 = NA,
                          Freq_28p = NA,
                          Freq_f1 = NA,
                          Freq_f2 = NA,
                          Freq_f3 = NA,
                          Freq_f4 = NA,
                          Freq_56 = NA,
                          Growth = NA,
                          N_change = NA#,
                          #N_change_antigenic = NA
  )
  
  for (i in 1:nrow(spike_tab)) {
    profile <- spike_tab$profile[i]
    spike_tab$lineage[i] <- paste(unique(viruses_28$lineage[viruses_28$mutations == profile]), collapse = ', ')
    spike_tab$N[i] <- sum(viruses$mutations == profile)
    spike_tab$N_28[i] <- sum(viruses_28$mutations == profile)
    spike_tab$N_28p[i] <- sum(viruses_28p$mutations == profile)
    spike_tab$N_f1[i] <- sum(viruses_f1$mutations == profile)
    spike_tab$N_f2[i] <- sum(viruses_f2$mutations == profile)
    spike_tab$N_f3[i] <- sum(viruses_f3$mutations == profile)
    spike_tab$N_f4[i] <- sum(viruses_f4$mutations == profile)
    spike_tab$N_56[i] <- sum(viruses_56$mutations == profile)
    if (do_pretty_profiles == F) {
      spike_tab$N_change[i] <- stringr::str_count(profile, ';') + 1
    } else {
      spike_tab$N_change[i] <- viruses$N_change[viruses$mutations == profile][1]
    }
    #spike_tab$N_change_antigenic[i] <- stringr::str_count(viruses$muts_ab[viruses$mutations == profile][1], ';') + 1
  }
  
  # Frequencies from counts
  spike_tab$Freq_28 <- spike_tab$N_28 / N_28
  spike_tab$Freq_28p <- spike_tab$N_28p / N_28_prev
  spike_tab$Freq_f1 <- spike_tab$N_f1 / N_f1
  spike_tab$Freq_f2 <- spike_tab$N_f2 / N_f2
  spike_tab$Freq_f3 <- spike_tab$N_f3 / N_f3
  spike_tab$Freq_f4 <- spike_tab$N_f4 / N_f4
  spike_tab$Freq_56 <- spike_tab$N_56 / N_56
  
  spike_tab$Growth <- (spike_tab$Freq_28 / spike_tab$Freq_28p * 100) - 100
  spike_tab$chi_2 <- ((spike_tab$Freq_f1 - spike_tab$Freq_56)^2 / spike_tab$Freq_56) + 
    ((spike_tab$Freq_f2 - spike_tab$Freq_56)^2 / spike_tab$Freq_56) +
    ((spike_tab$Freq_f3 - spike_tab$Freq_56)^2 / spike_tab$Freq_56) +
    ((spike_tab$Freq_f4 - spike_tab$Freq_56)^2 / spike_tab$Freq_56)
  
  # Create simple column of increase or decrease
  spike_tab$Increase <- NA
  spike_tab$Increase <- ifelse(spike_tab$Freq_28 > spike_tab$Freq_28p, 1, spike_tab$Increase)
  spike_tab$Increase <- ifelse(spike_tab$Freq_28p > spike_tab$Freq_28, 0, spike_tab$Increase)
  
  # use column to put take absolute value of change and put sign on it
  spike_tab$expansion <- ifelse(spike_tab$Increase == 0, spike_tab$chi_2 - spike_tab$chi_2*2, spike_tab$chi_2)
  
  spike_tab <- spike_tab[,c('profile', 'N_change', 'lineage', 
                            'N', 'N_28', 'Growth', 'expansion')]
  names(spike_tab) <- c('Profile',
                        'Number of amino acid changes',
                        'lineage(s)',
                        'Count',
                        'Count in latest 28 days',
                        'Change in Frequency vs. previous 28 days (%)',
                        'Expansion')
  spike_tab
}


spike_profiles_nations <- function(viruses = NA) {
  spike_tab_uk <- spike_profiles(viruses)
  spike_tab_eng <- spike_profiles(viruses[viruses$adm1 == 'England',])
  spike_tab_nir <- spike_profiles(viruses[viruses$adm1 == 'Northern_Ireland',])
  spike_tab_sco <- spike_profiles(viruses[viruses$adm1 == 'Scotland',])
  spike_tab_wal <- spike_profiles(viruses[viruses$adm1 == 'Wales',])
  
  spike_tab_uk$Geography <- 'United Kingdom'
  spike_tab_eng$Geography <- 'England'
  spike_tab_nir$Geography <- 'Northern Ireland'
  spike_tab_sco$Geography <- 'Scotland'
  spike_tab_wal$Geography <- 'Wales'
  
  spike_profiles_nations <- rbind(spike_tab_uk, spike_tab_eng,
                                  spike_tab_nir, spike_tab_sco,
                                  spike_tab_wal)
  spike_profiles_nations
  
}


get_pretty_profiles <- function(x = NA) {
  
  # Delta
  profile_delta <- c('T19R', 'G142D', 'L452R', 'T478K', 'D614G', 'P681R', 'D950N')
  N_change_delta <- 8
  # Omicron BA.1/B.1.1.529
  profile_omicron <- c('A67V', 'T95I', 'G142D', 'G339D', 'S371L', 'S373P', 'S375F',
                       'K417N', 'N440K', 'G446S', 'S477N', 'T478K', 'E484A', 'Q493R', 'G496S',
                       'Q498R', 'N501Y', 'Y505H', 'T547K', 'D614G', 'H655Y', 'N679K',
                       'P681H', 'N764K', 'D796Y', 'N856K', 'Q954H', 'N969K', 'L981F')
  N_change_omicron <- 30
  # Omicron BA.2
  profile_omicron_ba2 <- c('T19I', 'G142D', 'V213G', 'G339D', 'S371F', 'S373P', 'S375F',
                           'T376A', 'D405N', 'R408S', 'K417N', 'N440K', 'S477N', 'T478K',
                           'E484A', 'Q493R', 'Q498R', 'N501Y', 'Y505H', 'D614G', 'H655Y',
                           'N679K', 'P681H', 'N764K', 'D796Y', 'Q954H', 'N969K')
  N_change_omicron_ba2 <- 28
  
  
  x$pretty_profile <- NA
  x$N_change <- NA
  for (i in 1:nrow(x)) {
    if (grepl('AY.', x$lineage[i]) |
        grepl('B.1.617.2', x$lineage[i])) {
      
      profile <- str_split(x$mutations[i], pattern = ';')[[1]]
      
      voc_extra <- setdiff(profile, profile_delta)
      voc_absent <- setdiff(profile_delta, profile)
      # Correction to account for most common amplicon dropout
      voc_absent <- setdiff(voc_absent, 'G142D')
      
      # use collected info to generate 'pretty_profile'
      pretty_profile <- 'Delta'
      pretty_profile <- ifelse(length(voc_extra) > 0,
                               paste0(pretty_profile, ' +', paste(voc_extra, collapse = ',')),
                               pretty_profile)
      pretty_profile <- ifelse(length(voc_absent) > 0,
                               paste0(pretty_profile, ' -', paste(voc_absent, collapse = ',')),
                               pretty_profile)
      x$pretty_profile[i] <- pretty_profile
      x$N_change[i] <- N_change_delta + length(voc_extra) - length(voc_absent)
      
    } else if (grepl('BA.1', x$lineage[i]) |
               grepl('B.1.1.529', x$lineage[i])){
      
      profile <- str_split(x$mutations[i], pattern = ';')[[1]]
      voc_extra <- setdiff(profile, profile_omicron)
      voc_absent <- setdiff(profile_omicron, profile)
      # Correction to account for most common amplicon dropout errors
      voc_absent <- setdiff(voc_absent, 'G142D')
      voc_absent <- setdiff(voc_absent, 'K417N')
      voc_absent <- setdiff(voc_absent, 'N440K')
      voc_absent <- setdiff(voc_absent, 'G446S')
      voc_absent <- setdiff(voc_absent, 'N764K')
      voc_absent <- c()
      
      # use collected info to generate 'pretty_profile'
      pretty_profile <- 'Omicron (BA.1)'
      pretty_profile <- ifelse(length(voc_extra) > 0,
                               paste0(pretty_profile, ' +', paste(voc_extra, collapse = ',')),
                               pretty_profile)
      pretty_profile <- ifelse(length(voc_absent) > 0,
                               paste0(pretty_profile, ' -', paste(voc_absent, collapse = ',')),
                               pretty_profile)
      x$pretty_profile[i] <- pretty_profile
      x$N_change[i] <- N_change_omicron + length(voc_extra) - length(voc_absent)
      
    } else if (grepl('BA.2', x$lineage[i])){
      
      profile <- str_split(x$mutations[i], pattern = ';')[[1]]
      voc_extra <- setdiff(profile, profile_omicron_ba2)
      voc_absent <- setdiff(profile_omicron_ba2, profile)
      
      # Correction to account for most common amplicon dropout errors
      voc_absent <- c()
      
      # use collected info to generate 'pretty_profile'
      pretty_profile <- 'Omicron (BA.2)'
      pretty_profile <- ifelse(length(voc_extra) > 0,
                               paste0(pretty_profile, ' +', paste(voc_extra, collapse = ',')),
                               pretty_profile)
      pretty_profile <- ifelse(length(voc_absent) > 0,
                               paste0(pretty_profile, ' -', paste(voc_absent, collapse = ',')),
                               pretty_profile)
      x$pretty_profile[i] <- pretty_profile
      x$N_change[i] <- N_change_omicron_ba2 + length(voc_extra) - length(voc_absent)
      
    } else {
      x$pretty_profile[i] <- x$mutations[i]
      x$N_change[i] <- stringr::str_count(x$mutations[i], ';') + 1
    }
    
  }
  x
}

