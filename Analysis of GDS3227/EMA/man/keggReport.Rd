\name{keggReport}
\alias{keggReport}

\title{Text report from the result of the 'hyperGTest' function for KEGG
  pathway analysis}
\description{
 Create a txt report with the significant probesets for each significant KEGG pathway from the result of a call to 'hyperGTest'
}
\usage{
keggReport(hgOver,fileout = "report.txt", pack.annot, pvalue = 0.05)
}

\arguments{
  \item{hgOver}{result from a call to 'hyperGTest' }
  \item{fileout}{name of txt file, by default "report.txt"}
  \item{pack.annot}{the annotation package to use}
  \item{pvalue}{cutoff for significant pvalue, by default 0.05}
}

\details{
  \code{keggReport} needs packages \code{GOstats} and \code{KEGG.db} from Bioconductor.
  
}

\value{
 Txt report with the results of a call to 'hyperGTest'
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link[GOstats]{hyperGTest}},\code{\link{runHyperKEGG}}}
\keyword{print}
\keyword{internal}

