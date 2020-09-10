#' Model PTM and/or protein data and make adjustments if needed
#'
#' Takes summarized PTM data from proteinSummarization and models with
#' groupComparisonTMT. Can also take protein level data in the same format
#' and model with groupComparisonTMT. Including protein data allows
#' for adjusting PTM Fold Change by the change in protein abundance
#' without modification.
#'
#' @param data PTM dataset returned by the proteinSummarization function
#' @param protein Protein dataset returned by the proteinSummarization function
#' @param contrast.matrix Comparison between conditions of interests.
#'                        1) default is 'pairwise', which compare all possible pairs between two conditions.
#'                        2) Otherwise, users can specify the comparisons of interest. Based on the levels of conditions,
#'                        specify 1 or -1 to the conditions of interests and 0 otherwise.
#'                        The levels of conditions are sorted alphabetically.
#' @param moderated TRUE will moderate t statistic; FALSE (default) uses ordinary t statistic.
#' @param adj.method Adjusted method for multiple comparison. "BH" is default.
#'
#' @return A list \code{models} of all modeled and adjusted datasets
#'
#' @export
groupComparisonTMTPTM <- function(data, protein = NULL, contrast.matrix = "pairwise",
                                  moderated = FALSE, adj.method = "BH", calc.corr = FALSE) {

  ## save process output in each step
  allfiles <- list.files()
  filenaming <- "msstatstmtptm"

  if (length(grep(filenaming,allfiles)) == 0) {

    finalfile <- "msstatstmtptm.log"
    processout <- NULL

  } else {

    num <- 0
    finalfile <- "msstatstmtptm.log"

    while (is.element(finalfile, allfiles)) {
      num <- num + 1
      lastfilename <- finalfile ## in order to rea
      finalfile <- paste0(paste(filenaming, num, sep="-"), ".log")
    }

    finalfile <- lastfilename
    processout <- as.matrix(read.table(finalfile, header=TRUE, sep="\t"))
  }

  processout <- rbind(processout,
                      as.matrix(c(" ", " ", "MSstatsTMTPTM - groupComparisonTMTPTM function", " "), ncol=1))

  adj.protein = FALSE

  ## Check for missing variables in PTM
  if (is.null(data))
    stop("PTM estimates are missing!")
  required.columns <- c('Run', 'Protein', 'Abundance', 'Channel',
                'BioReplicate', 'Condition', 'TechRepMixture', 'Mixture')
  if (!all(required.columns %in% names(data))) {
    stop("Please include in the PTM list all the following elements: ",
         paste0(sQuote(required.columns), collapse = ", "))
  }

  ## Determine if PTM should be adjusted for protein level
  if (!is.null(protein)) {
    adj.protein = TRUE

    if (!all(required.columns %in% names(protein))) {
      stop("Please include in the Protein list all the following elements: ",
           paste0(sQuote(required.columns), collapse = ", "))
    }
  }

  setDT(data)
  setDT(protein)

  ## MSstatsTMT Modeling
  ptm_model <- MSstatsTMT::groupComparisonTMT(data, contrast.matrix, moderated, adj.method)

  models <- list('PTM.Model' = ptm_model)

  if (adj.protein) {

    ## MSstatsTMT Modeling
    protein_model <- MSstatsTMT::groupComparisonTMT(protein, contrast.matrix, moderated, adj.method)

    ## Parse site from protein name
    regex_protein <- '([^-]+)(?:_[^-]+){1}$'
    regex_site <- '_(?!.*_)([^-]+)'
    ptm_model_site_sep <- ptm_model
    ptm_model_site_sep[, c('Protein', 'Site') := list(stringr::str_match(as.matrix(ptm_model_site_sep[,'Protein']), regex_protein)[,2],
                                                      stringr::str_match(as.matrix(ptm_model_site_sep[,'Protein']), regex_site)[,2])]

    ## adjustProteinLevel function can only compare one label at a time
    comparisons <- unique(ptm_model_site_sep[, Label])

    adjusted_models_list <- lapply(comparisons, apply_ptm_adjustment)
    adjusted_models <- bind_rows(adjusted_models_list)
    setDT(adjusted_models)

    adjusted_models[, Protein:=do.call(paste, c(.SD, sep = "_")), .SDcols=c(1,2)]
    adjusted_models <- adjusted_models[, setdiff(names(adjusted_models), c("Site")), with = FALSE]

    models <- list('PTM.Model' = ptm_model, 'Protein.Model' = protein_model, 'Adjusted.Model' = adjusted_models)

  }

  return(models)
}


## TODO: Grab new code from MSstatsPTM for this function
adjustProteinLevel <- function(diffSite, diffProtein) {
  diffRef <- diffProtein[, c("Protein", "Label", "log2FC", "SE", "DF")]
  names(diffRef)[names(diffRef) == "log2FC"] <- "log2FC_ref"
  names(diffRef)[names(diffRef) == "SE"] <- "SE_ref"
  names(diffRef)[names(diffRef) == "DF"] <- "DF_ref"
  joined <- inner_join(diffSite, diffRef)

  missing_ctrl <- joined[joined$log2FC == Inf, ]  # PTM missing in control
  missing_case <- joined[joined$log2FC == -Inf, ] # PTM missing in case
  res_mctrl <- tibble(Protein = missing_ctrl$Protein,
                      Site = missing_ctrl$Site, Label = missing_ctrl$Label,
                      log2FC = Inf, SE = NA, Tvalue = NA, DF = NA, pvalue = NA)
  res_mcase <- tibble(Protein = missing_case$Protein,
                      Site = missing_case$Site, Label = missing_case$Label,
                      log2FC = -Inf, SE = NA, Tvalue = NA, DF = NA, pvalue = NA)

  idx_full <- abs(joined$log2FC) != Inf & abs(joined$log2FC_ref) != Inf
  full <- joined[idx_full, ]
  log2fc <- full$log2FC - full$log2FC_ref
  s2 <- full$SE ^ 2
  s2_ref <- full$SE_ref ^ 2
  stderr <- sqrt(s2 + s2_ref)
  numer <- (s2 + s2_ref) ^ 2
  denom <- (s2 ^ 2 / full$DF + s2_ref ^ 2 / full$DF_ref)
  df <- numer / denom
  tval <- log2fc / stderr
  pval <- 2 * stats::pt(abs(tval), df, lower.tail = FALSE)
  bh.pval <- p.adjust(pval, method = 'BH')

  res_full <- tibble(Protein = full$Protein,
                     Site = full$Site,
                     Label = full$Label,
                     log2FC = log2fc,
                     SE = stderr,
                     Tvalue = tval,
                     DF = df,
                     pvalue = pval)

  bind_rows(res_full, res_mctrl, res_mcase)
}
