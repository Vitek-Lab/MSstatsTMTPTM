#' Model PTM and/or protein data and make adjustments if needed
#'
#' Takes summarized PTM data from proteinSummarization and models with
#' groupComparisonTMT. Can also take protein level data in the same format
#' and model with groupComparisonTMT. Including protein data allows
#' for adjusting PTM Fold Change by the change in protein abundance
#' without modification.
#'
#' @export
#' @importFrom dplyr filter bind_rows distinct inner_join select tibble %>%
#' @importFrom stats p.adjust xtabs
#' @importFrom utils read.table write.table
#' @importFrom MSstatsTMT groupComparisonTMT
#' @importFrom stringr str_match
#' @param data.ptm Name of the output of proteinSummarization function with PTM data. It should have columns named
#'                   `Protein`, `TechRepMixture`,  `Mixture`, `Run`, `Channel`, `Condition`, `BioReplicate`, `Abundance`.
#' @param data.protein Protein dataset returned by the proteinSummarization function
#' @param contrast.matrix Comparison between conditions of interests.
#'                        1) default is 'pairwise', which compare all possible pairs between two conditions.
#'                        2) Otherwise, users can specify the comparisons of interest. Based on the levels of conditions,
#'                        specify 1 or -1 to the conditions of interests and 0 otherwise.
#'                        The levels of conditions are sorted alphabetically.
#' @param moderated TRUE will moderate t statistic; FALSE (default) uses ordinary t statistic.
#' @param adj.method Adjusted method for multiple comparison. "BH" is default.
#'
#' @return A list \code{models} of all modeled and adjusted datasets
#' @examples
#' # Load summarized datasets from MSstatsTMT proteinSummarization function
#' data(quant.msstats.ptm)
#' data(quant.msstats.protein)
#'
#' # Load specific contrast matrix
#' data(example.comparisons)
#'
#" # test for specified condition comparisons only
#' test.pairwise <- groupComparisonTMTPTM(data.ptm=quant.msstats.ptm,
#'                                       data.protein=quant.msstats.protein,
#'                                       contrast.matrix = example.comparisons)
#'
groupComparisonTMTPTM <- function(data.ptm, data.protein = NULL, contrast.matrix = "pairwise",
                                  moderated = FALSE, adj.method = "BH") {

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

  Protein = Label = Site = NULL

  adj.protein = FALSE

  ## Check for missing variables in PTM
  if (is.null(data.ptm))
    stop("PTM estimates are missing!")
  required.columns <- c('Run', 'Protein', 'Abundance', 'Channel',
                'BioReplicate', 'Condition', 'TechRepMixture', 'Mixture')
  if (!all(required.columns %in% names(data.ptm))) {
    stop("Please include in the PTM list all the following elements: ",
         paste0(sQuote(required.columns), collapse = ", "))
  }

  ## Determine if PTM should be adjusted for protein level
  if (!is.null(data.protein)) {
    adj.protein = TRUE

    if (!all(required.columns %in% names(data.protein))) {
      stop("Please include in the Protein list all the following elements: ",
           paste0(sQuote(required.columns), collapse = ", "))
    }
  }

  ## MSstatsTMT Modeling
  ptm_model <- MSstatsTMT::groupComparisonTMT(data.ptm, contrast.matrix, moderated, adj.method)

  models <- list('PTM.Model' = ptm_model)

  if (adj.protein) {

    ## MSstatsTMT Modeling
    protein_model <- MSstatsTMT::groupComparisonTMT(data.protein, contrast.matrix, moderated, adj.method)

    ## Parse site from protein name
    regex_protein <- '([^-]+)(?:_[^-]+){1}$'
    regex_site <- '_(?!.*_)([^-]+)'
    ptm_model_site_sep <- ptm_model %>% mutate(Site = str_match(
      Protein, regex_site)[,2], Protein = str_match(Protein, regex_protein)[,2])

    ## adjustProteinLevel function can only compare one label at a time
    comparisons <- (ptm_model_site_sep %>% distinct(Label))[[1]]
    adjusted_models <- data.frame()
    for (i in 1:length(comparisons)) {
      temp_adjusted_model <- apply_ptm_adjustment(comparisons[[i]], ptm_model_site_sep, protein_model)
      adjusted_models <- rbind(adjusted_models, temp_adjusted_model)
    }

    adjusted_models$Protein <- paste(adjusted_models$Protein, adjusted_models$Site, sep = '_')
    adjusted_models <- adjusted_models %>% select(-Site)

    models <- list('PTM.Model' = ptm_model, 'Protein.Model' = protein_model, 'Adjusted.Model' = adjusted_models)

  }

  return(models)
}

#' @keywords internal
apply_ptm_adjustment <- function(label, ptm_model, protein_model){

  Label = NULL

  temp_ptm_model <- ptm_model %>% filter(Label == label)
  temp_protein_model <- protein_model %>% filter(Label == label)

  ## Function from MSstatsPTM Compare
  temp_adjusted_model <- adjustProteinLevel(temp_ptm_model, temp_protein_model)
  temp_adjusted_model$adj.pvalue <- p.adjust(temp_adjusted_model$pvalue, method = 'BH')
  temp_adjusted_model
}

## TODO: Replace this with MSstatsPTM function once made external
#' @keywords internal
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
