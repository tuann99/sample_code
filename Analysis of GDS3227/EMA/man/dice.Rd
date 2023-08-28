\name{dice}
\alias{dice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute Dice distance on a data matrix}
\description{
   Computes Dice distances between the columns of a 0-1-matrix.
}
\usage{
dice(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{0-1-matrix.}
}
\details{
  Missing values are allowed, and are excluded from all computations
  involving the cols within which they occur.
}
\value{
  A distance matrix
 }
\references{Dice, L.R. (1945), Measures of the amount of ecologic
  association between species, J. Ecology, 26, 297-302.}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\keyword{cluster}
\keyword{internal}
