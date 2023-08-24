\name{htmlresult}
\alias{htmlresult}

\title{Html report from the result of the 'hyperGTest' function}
\description{
  Create a html report from the result of a call to 'hyperGTest'
}
\usage{
htmlresult(hgOver, filename, app = FALSE, categorySize = 1, pvalue=0.05)
}

\arguments{
  \item{hgOver}{result from a call to 'hyperGTest'}
  \item{filename}{file name with path for html result}
  \item{app}{if TRUE, overwrite the pre-existing html report file, by default FALSE}
  \item{categorySize}{an integer, minimum size for category, by default  =1}
  \item{pvalue}{the pvalue threshold for significance, default =0.05}
}

\details{
  \code{htmlresult} needs packages \code{GOstats} and \code{GO.db} from Bioconductor.
}

\value{
  An html report with the result from a call to 'hyperGTest'

  Id, pvalue, Odd Ratio, Expected Count, Count, Size and Term are returned 
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link[GOstats]{hyperGTest}},\code{\link{runHyperGO}}}

\keyword{print}
\keyword{internal}
