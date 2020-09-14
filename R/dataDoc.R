
#' Example of input PTM dataset for TMT experiments.
#'
#' It can be the output of PDtoMSstatsTMTFormat or other MSstatsTMT converter functions.
#' It includes peak intensities for a variety of PTMs.
#' The variables are as follows:
#'
#' \itemize{
#'   \item ProteinName : Name of protein with modification site mapped in with an underscore. ie "Protein_4_Y474"
#'   \item PeptideSequence
#'   \item Charge
#'   \item PSM
#'   \item Mixture : Mixture of samples labeled with different TMT reagents, which can be analyzed in
#'   a single mass spectrometry experiment. If the channal doesn't have sample, please add `Empty' under Condition.
#'   \item TechRepMixture : Technical replicate of one mixture. One mixture may have multiple technical replicates.
#'   For example, if `TechRepMixture' = 1, 2 are the two technical replicates of one mixture, then they should match
#'   with same `Mixture' value.
#'   \item Run : MS run ID.
#'   \item Channel : Labeling information (126, ... 131).
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject. If the channal doesn't have sample, please add `Empty' under BioReplicate.
#'   \item Intensity
#' }
#'
#' @format A data frame with 24704 rows and 11 variables.
#' @examples
#' head(raw.ptm)
#'
"raw.ptm"

#' Example of input Protein dataset for TMT experiments.
#'
#' It can be the output of PDtoMSstatsTMTFormat or other MSstatsTMT converter functions.
#' It includes peak intensities for a variety of PTMs.
#' This is the companion file to the raw.ptm dataset, includes unmodified protein data.
#' The variables are as follows:
#'
#' \itemize{
#'   \item ProteinName : Name of protein
#'   \item PeptideSequence
#'   \item Charge
#'   \item PSM
#'   \item Mixture : Mixture of samples labeled with different TMT reagents, which can be analyzed in
#'   a single mass spectrometry experiment. If the channal doesn't have sample, please add `Empty' under Condition.
#'   \item TechRepMixture : Technical replicate of one mixture. One mixture may have multiple technical replicates.
#'   For example, if `TechRepMixture' = 1, 2 are the two technical replicates of one mixture, then they should match
#'   with same `Mixture' value.
#'   \item Run : MS run ID.
#'   \item Channel : Labeling information (126, ... 131).
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject. If the channal doesn't have sample, please add `Empty' under BioReplicate.
#'   \item Intensity
#' }
#'
#' @format A data frame with 620476 rows and 11 variables.
#' @examples
#' head(raw.protein)
#'
"raw.protein"

#' Example of output from proteinSummarizaiton function for PTM data
#'
#' It is made from \code{\link{raw.ptm}}.
#' It is the output of proteinSummarization function from MSstatsTMT.
#' It should include the required columns as below.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID
#'   \item Protein : Protein ID with modification site mapped in. Ex. Protein_1002_S836
#'   \item Abundance: Protein-level summarized abundance
#'   \item Channel : Labeling information (126, ... 131)
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject.
#'   \item TechRepMixture : Unique ID for technical replicate of one TMT mixture.
#'   \item Mixture : Unique ID for TMT mixture.
#' }
#'
#' @format A data frame with 19205 rows and 8 variables.
#' @examples
#' head(quant.msstats.ptm)
#'
"quant.msstats.ptm"

#' Example of output from proteinSummarizaiton function for Protein data
#'
#' It is made from \code{\link{raw.protein}}.
#' It is the output of proteinSummarization function from MSstatsTMT.
#' It should include the required columns as below.
#' The variables are as follows:
#'
#' \itemize{
#'   \item Run : MS run ID
#'   \item Protein : Protein ID
#'   \item Abundance: Protein-level summarized abundance
#'   \item Channel : Labeling information (126, ... 131)
#'   \item Condition : Condition (ex. Healthy, Cancer, Time0)
#'   \item BioReplicate : Unique ID for biological subject.
#'   \item TechRepMixture : Unique ID for technical replicate of one TMT mixture.
#'   \item Mixture : Unique ID for TMT mixture.
#' }
#'
#' @format A data frame with 93258 rows and 8 variables.
#' @examples
#' head(quant.msstats.protein)
#'
"quant.msstats.protein"

#' Example contrast matrix for input into the groupComparisonTMTPTM function
#' 
#' Manually specified comparisons of interest for contrast.matrix
#' arguement of groupComparisonTMTPTM.
#' 
#' \itemize{
#'   \item Condition_1, ... Condition_6 : Column names are conditions in dataset
#'   \item 1-4, ... 5-6 : Row names are comparisons of interest
#' }
#' 
#' @format A data frame with 9 rows and 6 variables.
#' @examples
#' head(example.comparisons)
#'
"example.comparisons"