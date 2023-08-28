\name{intersectg}
\alias{intersectg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized version of intersect for n objects}
\description{
  This function returns the intersection for the n objects in argument
}
\usage{
intersectg(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{The objects to compare}
}

\value{
  A vector of intersection between the n objects
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
\examples{
a<-1:10
b<-2:7
c<-5:8
intersectg(a,b,c)
}
\keyword{utilities}

