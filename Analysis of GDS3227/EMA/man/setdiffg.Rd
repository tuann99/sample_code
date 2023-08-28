\name{setdiffg}
\alias{setdiffg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized version of setdiff for n objects}
\description{
  This function returns the difference between the first argument and
  the (n-1) others
}
\usage{
setdiffg(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{The objects to compare}
}

\value{
A vector of differences between the first object and the (n-1) others
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
\examples{
a<-1:10
b<-2:5
c<-7:8
setdiffg(a,b,c)
}
\keyword{utilities}

