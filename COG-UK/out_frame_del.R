## Out-of-frame deletion
# Generate descriton of codons most likely affected by deletion when out-of-frame
# ref_nt can be provided by function like 
# only desinged to work with deletions of length multiple of 3
out_frame_del <- function(ref_nt = NA, del_S = NA, del_length = NA) {
  
  del_E <- del_S + del_length + 1
  
  # copy the refstrain and then insert the deletion in the climb coordinates
  del_nt <- ref_nt
  del_string <- paste0(rep('-', length = del_length), collapse = "")
  substring(del_nt, del_S, del_E) <- del_string
  
  # Identify the codons affected by the deletion
  codon_n_first <- ceiling(del_S / 3)
  codon_n_last <- ceiling(del_E / 3)
  
  ## Extract the relevant codons of the wt and del mutant nt sequence
  # length of these extracts should should be del_length + 3
  # multiplying by 3 gives last nucleotide in codon so -2 to get start
  ref_nt_ext <- stringr::str_sub(ref_nt,
                                 start = (codon_n_first * 3) - 2,
                                 end = (codon_n_last * 3))
  codon_first <- stringr::str_sub(ref_nt_ext, 1, 3)
  codon_last <- stringr::str_sub(ref_nt_ext, -3, -1)
  
  ref_aa_ext <- bioseq::seq_translate(bioseq::dna(ref_nt_ext))
  ref_aa_ext <- as.character(ref_aa_ext)
  
  # for del sequence, remove gaps before translating
  del_nt_ext <- stringr::str_sub(del_nt,
                                 start = (codon_n_first * 3) - 2,
                                 end = (codon_n_last * 3))
  del_codon <- gsub('-', '', del_nt_ext)
  del_aa <- bioseq::seq_translate(bioseq::dna(del_codon))
  del_aa <- as.character(del_aa)
  
  ## Decision
  decision = F
  # if the del aa matches one but not both of the terminal aa, make decision
  if (del_aa == stringr::str_sub(ref_aa_ext, 1, 1) &
      del_aa != stringr::str_sub(ref_aa_ext, -1, -1)) {
    decision <- 1
  } else if (del_aa != stringr::str_sub(ref_aa_ext, 1, 1) &
             del_aa == stringr::str_sub(ref_aa_ext, -1, -1)) {
    decision <- 2
  } 
  
  # if tied, move to hamming distance between codons - full codon first
  if (decision == F) {
    hamming_first <- sum(strsplit(codon_first,  "")[[1]] != strsplit(del_codon,  "")[[1]])
    hamming_last <- sum(strsplit(codon_last,  "")[[1]] != strsplit(del_codon,  "")[[1]])
    
    # favour smaller distance
    if (hamming_first < hamming_last) {
      decision <- 1
    } else if (hamming_first > hamming_last) {
      decision <- 2
    }
  }
  
  # if tied move to hamming distance between codons - [[1]][2]
  if (decision == F) {
    hamming_first <- sum(strsplit(codon_first,  "")[[1]][2] != strsplit(del_codon,  "")[[1]][2])
    hamming_last <- sum(strsplit(codon_last,  "")[[1]][2] != strsplit(del_codon,  "")[[1]][2])
    
    if (hamming_first < hamming_last) {
      decision <- 1
    } else if (hamming_first > hamming_last) {
      decision <- 2
    }
  }
  
  # finally priority is given to first pos rather than last [[1]][1]
  if (decision == F) {
    hamming_first <- sum(strsplit(codon_first,  "")[[1]][1] != strsplit(del_codon,  "")[[1]][1])
    hamming_last <- sum(strsplit(codon_last,  "")[[1]][1] != strsplit(del_codon,  "")[[1]][1])
    
    if (hamming_first < hamming_last) {
      decision <- 1
    } else if (hamming_first > hamming_last) {
      decision <- 2
    }
  }
  
  ambiguous <- F
  if (decision == F) {
    decision <- 1
    ambiguous <- T
    cat('Warning - ambiguous deletion\n')
  }
  
  # generate output based on decision
  if (del_length == 3 & decision == 1) {
    output <- paste0('del', codon_n_last)
  } else if (del_length > 3 & decision == 1) {
    output <- paste0('del', codon_n_first + 1, '-', codon_n_last)
  } else if (del_length == 3 & decision == 2) {
    output <- paste0('del', codon_n_first)
  } else if (del_length > 3 & decision == 2) {
    output <- paste0('del', codon_n_first, '-', codon_n_last -1)
  }
  
  if (ambiguous == T) {
    cat('Returning', output, '\n\n')
  }
  
  output
  
}

