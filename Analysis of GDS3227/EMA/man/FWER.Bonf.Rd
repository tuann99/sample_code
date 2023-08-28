\name{FWER.Bonf}
\alias{FWER.Bonf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple testing correction using FWER}
\description{
  Family Wise Error rate - Bonferoni correction.
}
\usage{
FWER.Bonf(pval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{A vector of raw pvalues}
}

\value{
A vector with the correced pvalues.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link[multtest]{mt.rawp2adjp}}}
\keyword{htest}
\keyword{internal}
