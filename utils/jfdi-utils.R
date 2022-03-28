#' get.sensitivity.dataframe
#'
#' Takes the JFDI dataframe and generates the total, TP, FN, TPR and FNR columns.
#'
#' @param dat JFDI raw dataframe
#' @returns sensitivity dataframe
get.sensitivity.dataframe <- function(dat){
  dat %>%
    dplyr::filter(!is.na(Observed) & expected_af != "WT") %>%
    dplyr::count(participant, sequencing_platform, library_method,
          contrived_source, dna_source, expected_af, sample_replicate,
          Observed) %>%
    dplyr::group_by(participant, sequencing_platform, library_method,
             contrived_source, dna_source, expected_af, sample_replicate) %>%
    dplyr::mutate(total=sum(n)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from=Observed, values_from = n, values_fill = 0) %>%
    dplyr::rename(TP=Detected, FN=Missed) %>%
    dplyr::mutate(TPR=TP/total, FNR=FN/total)
}

#' get.specificity.dataframe
#'
#' Takes the JFDI dataframe and generates the total, TN, FP, TNR, and FPR columns.
#'
#' @param dat JFDI raw dataframe
#' @returns specificity dataframe
get.specificity.dataframe <- function(dat){
  dat %>%
    dplyr::filter(!is.na(Observed) & expected_af == "WT") %>%
    dplyr::count(participant, sequencing_platform, library_method,
          contrived_source, dna_source, expected_af, sample_replicate,
          Observed) %>%
    dplyr::group_by(participant, sequencing_platform, library_method,
             contrived_source, dna_source, expected_af, sample_replicate) %>%
    dplyr::mutate(total=sum(n)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from=Observed, values_from = n, values_fill = 0) %>%
    dplyr::rename(FP=Detected, TN=Missed) %>%
    dplyr::mutate(FPR=FP/total, TNR=TN/total)
}

#' get.within.conocrdance.df
#'
#' Takes the JFDI dataframe and estimates concordance rate within a participant.
#'
#' @param dat JFDI raw dataframe
#' @returns formatted concordance dataframe
get.within.concordance.df <- function(dat){
  # First, rework the input dataframe
  dat.sel.rep <- dat %>%
    dplyr::select(contrived_source, dna_source, participant,
                  expected_af, variant, sample_replicate, Observed)

  # This is the ugly nested lapply method I used, so bear with me.
  # First is across contrived source
  rep.result.both <- lapply(c("Horizon", "SeraCare"), function(csource){
    # Subset to this contrived source and get the unique participants
    parts <- (dat.sel.rep %>%
                dplyr::filter(contrived_source==csource) %>%
                dplyr::distinct(participant))$participant
    # Now lapply across each participant
    parts.dsource.af.comb.res <- lapply(parts, function(bp){
      # Subset to this contrived source and participant to get the set of dna source
      src <- (dat.sel.rep %>%
                dplyr::filter(contrived_source==csource & participant==bp) %>%
                dplyr::distinct(dna_source))$dna_source

      # Now lapply across each dna_source
      dsource.af.comb.res <- lapply(src, function(dsource){
        # Subset to this contrived source, participant, and dna_source. get expected AF levels
        eaf <- (dat.sel.rep %>%
                  dplyr::filter(contrived_source==csource & participant==bp & dna_source == dsource) %>%
                  dplyr::distinct(expected_af))$expected_af

        # Another gross lapply across each expected AF
        af.comb.res <- lapply(eaf, function(af){
          # Subset to the current working dataframe, that was pivoted to wide format
          curr <- dat.sel.rep %>%
            dplyr::filter(contrived_source==csource & participant==bp & dna_source == dsource & expected_af==af) %>%
            dplyr::select(variant, sample_replicate, Observed) %>%
            tidyr::pivot_wider(names_from=sample_replicate, values_from=Observed)

          # Convert to matrix (i skip first column cause that is variant which is the identifier)
          curr.mat <- as.matrix(curr[,c(2:ncol(curr))])
          row.names(curr.mat) <- curr$variant
          # I drop NA's here
          na.omit(curr.mat)
          # Use combinat package to get all combinations across the replicates
          combs <- combinat::combn(colnames(curr.mat), 2)

          # Now, we actually get to our final lapply across all replicate combinations
          if(is.vector(combs)){
            sel.a <- curr.mat[,1] == "Detected"
            # mask of detected records in second rep in combination
            sel.b <- curr.mat[,2] == "Detected"
            # Use mask of detected records in first rep in comb to find where concordant
            lres <- curr.mat[,1][sel.a] == curr.mat[,2][sel.a]
            # nc: Number of concordant
            nc <- sum(lres, na.rm=TRUE)
            # na: Number of detected in first rep of combination
            na <- sum(sel.a, na.rm=TRUE)
            # nb: Number of detected in second rep of combination
            nb <- sum(sel.b, na.rm=TRUE)
            # rc: Concordant rate
            rc <- nc / mean(c(na, nb), na.rm=TRUE)
            # Return results as a dataframe
            return(data.frame(contrived_source=csource, participant=bp,
                              dna_source=dsource, expected_af=af,
                              nc=nc, na=na, nb=nb, rc=rc))
          } else if(!is.null(dim(combs))) {
            comb.res <- lapply(c(1:dim(combs)[2]), function(comb){
              # mask of detected records in first rep in combination
              sel.a <- curr.mat[,combs[1,comb]] == "Detected"
              # mask of detected records in second rep in combination
              sel.b <- curr.mat[,combs[2,comb]] == "Detected"
              # Use mask of detected records in first rep in comb to find where concordant
              lres <- curr.mat[,combs[1,comb]][sel.a] == curr.mat[,combs[2,comb]][sel.a]
              # nc: Number of concordant
              nc <- sum(lres, na.rm=TRUE)
              # na: Number of detected in first rep of combination
              na <- sum(sel.a, na.rm=TRUE)
              # nb: Number of detected in second rep of combination
              nb <- sum(sel.b, na.rm=TRUE)
              # rc: Concordant rate
              rc <- nc / mean(c(na, nb), na.rm=TRUE)
              # Return results as a dataframe
              return(data.frame(contrived_source=csource, participant=bp,
                                dna_source=dsource, expected_af=af,
                                nc=nc, na=na, nb=nb, rc=rc))
            })
            # All of these do.call are just collapsing the lists of data frames into a single data frame
            return(do.call(rbind, comb.res))
          }
        })
        return(do.call(rbind, af.comb.res))
      })
      return(do.call(rbind, dsource.af.comb.res))
    })
    return(do.call(rbind, parts.dsource.af.comb.res))
  })
  return(do.call(rbind, rep.result.both))
}


#' get.between.concordance.df
#'
#' Takes the JFDI dataframe and estimates concordance rate between participants
#'
#' @param dat JFDI raw dataframe
#' @returns concordance dataframe
get.between.concordance.df <- function(dat){
  # First, rework the input dataframe
  dat.sel.rep <- dat %>%
    dplyr::select(contrived_source, dna_source, participant,
                  expected_af, variant, sample_replicate, Observed)

  results.list <- list()
  # First is across contrived source
  for(csource in c("Horizon", "SeraCare")){
    eafs <- unique(subset(dat.sel.rep, contrived_source==csource)$expected_af)
    # Then across AF
    for(eaf in eafs){
      dsources <- unique(subset(dat.sel.rep, contrived_source==csource & expected_af==eaf)$dna_source)
      # Then across dna sources
      for(dsource in dsources){
        # Subset df
        curr.df <- dat.sel.rep %>%
          filter(contrived_source==csource & expected_af==eaf & dna_source==dsource)
        parts <- unique(curr.df$participant)
        for(pi in parts){
          p.a <- curr.df %>%
            dplyr::filter(participant == pi) %>%
            dplyr::select(variant, sample_replicate, Observed) %>%
            tidyr::pivot_wider(names_from = sample_replicate, values_from=Observed)
          p.a.mat <- as.matrix(p.a[,2:ncol(p.a)])
          row.names(p.a.mat) <- p.a$variant
          for(pj in parts){
            if(pi != pj){
              # do stuff
              p.b <- curr.df %>%
                dplyr::filter(participant == pj) %>%
                dplyr::select(variant, sample_replicate, Observed) %>%
                tidyr::pivot_wider(names_from = sample_replicate, values_from=Observed)
              p.b.mat <- as.matrix(p.b[,2:ncol(p.b)])
              row.names(p.b.mat) <- p.b$variant

              # Loop over each replicate
              for(i in 1:dim(p.a.mat)[2]){
                a.vec <- p.a.mat[which(!is.na(p.a.mat[,i])),i]
                for(j in 1:dim(p.b.mat)[2]){
                  b.vec <- p.b.mat[which(!is.na(p.b.mat[,j])),j]
                  # Get common covered variants
                  ab.comm <- intersect(names(a.vec), names(b.vec))
                  a.vec <- a.vec[ab.comm]
                  b.vec <- b.vec[ab.comm]

                  # Select Detected
                  a.sel <- which(a.vec == "Detected")
                  na <- length(a.sel)
                  b.sel <- which(b.vec == "Detected")

                  nb <- length(b.sel)
                  nc <- sum(a.vec[a.sel] == b.vec[a.sel], na.rm=TRUE)
                  rc <- nc / mean(c(na, nb), na.rm=TRUE)

                  results.list <- append(results.list, list(data.frame(contrived_source=csource, dna_source=dsource,
                                                                       expected_af=eaf, participant.a=pi, participant.b=pj,
                                                                       na=na, nb=nb, nc=nc, rc=rc)))
                }
              }
            }
          }
        }
      }
    }
  }
  return(do.call(rbind, results.list))
}
