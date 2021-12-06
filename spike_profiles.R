## Function: take virus_df, a dataframe with list of mutations separated by ';' and sample dates
spike_profiles <- function(viruses = NA, spike_period = 7) {
  
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
                          N_change = NA,
                          N_change_antigenic = NA)
  
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
    spike_tab$N_change[i] <- str_count(profile, ';') + 1
    spike_tab$N_change_antigenic[i] <- str_count(viruses$muts_ab[viruses$mutations == profile][1], ';') + 1
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
  
  # use column to put sign on 
  spike_tab$expansion <- ifelse(spike_tab$Increase == 0, spike_tab$chi_2 - spike_tab$chi_2*2, spike_tab$chi_2)
  
  spike_tab <- spike_tab[,c('profile', 'N_change', 'N_change_antigenic', 'lineage',
                            'N', 'N_28', 'Growth', 'expansion')]
  spike_tab
}