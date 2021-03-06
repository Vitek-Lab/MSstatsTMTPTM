% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{quant.msstats.protein}
\alias{quant.msstats.protein}
\title{Example of output from proteinSummarizaiton function for Protein data}
\format{
A data frame with 93258 rows and 8 variables.
}
\usage{
quant.msstats.protein
}
\description{
It is made from \code{\link{raw.protein}}.
It is the output of proteinSummarization function from MSstatsTMT.
It should include the required columns as below.
The variables are as follows:
}
\details{
\itemize{
\item Run : MS run ID
\item Protein : Protein ID
\item Abundance: Protein-level summarized abundance
\item Channel : Labeling information (126, ... 131)
\item Condition : Condition (ex. Healthy, Cancer, Time0)
\item BioReplicate : Unique ID for biological subject.
\item TechRepMixture : Unique ID for technical replicate of one TMT
mixture.
\item Mixture : Unique ID for TMT mixture.
}
}
\examples{
head(quant.msstats.protein)

}
\keyword{datasets}
