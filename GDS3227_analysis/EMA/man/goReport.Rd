\name{goReport}
\alias{goReport}

\title{Text report from the 'hyperGTest' function}
\description{
  Create a text report with the significant probesets for each significant GO category from the result of a call to 'hyperGTest'
}
\usage{
goReport(hgOver, fileout = "report.txt", type = c("CC", "MF", "BP"),
pack.annot, pvalue = 0.05, categorySize=1)
}

\arguments{
  \item{hgOver}{result from a call to 'hyperGTest'}
  \item{fileout}{name of txt file, by default "report.txt"}
  \item{type}{type of GO categories, by default = c("CC", "MF", "BP")}
  \item{pack.annot}{the annotation package to use}
  \item{pvalue}{cutoff for significant pvalue, by default 0.05}
  \item{categorySize}{integer, minimum size for category, by default = 1}
}

\details{
  \code{goReport} needs packages \code{GOstats} and \code{GO.db} from Bioconductor.
}

\value{
  Txt report with the results of a call to 'hyperGTest'
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}


\seealso{
  \code{\link[GOstats]{hyperGTest}}, \code{\link{runHyperGO}} }

\keyword{print}
\keyword{internal}

