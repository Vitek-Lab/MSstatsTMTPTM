#' MSstatsTMTPTM: A package for detecting differencially abundant post translational modifications (PTM) in shotgun mass spectrometry-bsed proteomic experiments with tandem mass tag (TMT) labeling.
#'
#' A set of tools for detecting differentially abundant PTMs and proteins in shotgun mass spectrometry-based proteomic experiments with tandem mass tag (TMT) labeling.
#'
#' @section functions :
#' \itemize{
#'   \item \code{\link{dataProcessPlotsTMTPTM}} : Data visualization of PTM and global protein levels. Can plot either Profile plots to identify the potential sources of variation for each protein, or quality control plots to evaluate the systematic bias between MS runs.
#'   \item \code{\link{groupComparisonTMTPTM}} : Tests for significant changes in PTM abundance adjusted for global protein abundance across conditions based on a family of linear mixed-effects models in TMT experiment.
#' }
#'
#' @docType package
#' @name MSstatsTMTPTM
NULL