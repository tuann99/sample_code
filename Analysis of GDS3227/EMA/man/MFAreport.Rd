\name{MFAreport}
\alias{MFAreport}
\title{Function to create a txt and pdf report with the main statistics
  and graphics of the MFA.}
\description{
  This function build an analysis report with the main statistics
  usefull to interpret the results of the MFA and the main graphics.
}
\usage{
MFAreport(resMFA, file.txt, file.pdf)
}
\arguments{
  \item{resMFA}{Output of \code{\link{MFA}} or \code{\link{runMFA}}}
  \item{file.txt}{Character string giving the name of the txt file}
  \item{file.pdf}{Character string giving the name of the pdf file}
}
\details{
  The grahics produced are:
  \itemize{
  \item individuals representation in global analysis, axes 1 to 4
  \item individuals representation in global analysis with partial
    individuals, axes 1 to 4
  \item individuals representation in each separated analysis, axes 1
    to 4
  \item groups representation, axes 1 to 4
}
}
\value{
  A txt file with the main statistics and a pdf file with the main graphics.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{MFA}}, \code{\link{plot.MFA}}, \code{\link{runMFA}}}

\keyword{print}
\keyword{internal}

