\name{FDR.BH}
\Rdversion{1.1}
\alias{FDR.BH}

\title{
  FDR.BH
}
\description{
False Discovery Rate correction
}
\usage{
FDR.BH(pval, TST = FALSE, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{A vector of raw pvalues}
  \item{TST}{If TRUE, the Two-stages procedure is used to correct the data. If FALSE, the standard BH95 procedure is used.}
   \item{q}{The error rate to use for the Two-stages procedure (FDR-TST).}
}


\value{
  A vector with the corrected pvalues.
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
\seealso{
\code{\link[multtest]{mt.rawp2adjp}}}
\keyword{internal}
