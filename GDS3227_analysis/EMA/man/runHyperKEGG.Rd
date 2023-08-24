\name{runHyperKEGG}
\alias{runHyperKEGG}

\title{Run KEGG pathway analysis based on hypergeometric test from a probeset list}
\description{
  Run KEGG pathway analysis based on hypergeometric test from a probeset list
}
\usage{
runHyperKEGG(list, pack.annot, categorySize = 1, name = "hyperKEGG",
htmlreport = TRUE, txtreport = TRUE, tabResult = FALSE, pvalue = 0.05)
}

\arguments{
  \item{list}{vector of character with probeset names}
  \item{pack.annot}{character string, annotation package to use}
  \item{categorySize}{integer, minimum size for category, by default = 1}
  \item{name}{character, name for output files, by default "hyperKEGG"}
  \item{htmlreport}{logical, if TRUE, a html report is created, by default TRUE}
  \item{txtreport}{logical, if TRUE, a txt report is created, by default TRUE}
  \item{tabResult}{logical, if TRUE, a list with the results is created, by default FALSE}
  \item{pvalue}{numeric, a cutoff for the hypergeometric test pvalue, by default 0.05}
}

\value{
  Txt and html report

  Data.frame with KEGG Id, pvalue, Odd Ratio, Expected count, Size and KEGG Term

}

\details{
  \code{runHyperKEGG} needs packages \code{GOstats} and \code{KEGG.db} from Bioconductor.
}  

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link[GOstats]{hyperGTest}},\code{\link{runHyperGO}}}
\examples{
\dontrun{
require(hgu133plus2.db)
data(marty)

## Probe list
probeList <- rownames(marty)[1:50]

## Hypergeometric test for KEGG pathway
res <- runHyperKEGG(probeList, htmlreport = FALSE, txtreport = FALSE,
    tabResult = TRUE, pack.annot = "hgu133plus2.db")
}
}

\keyword{htest}

