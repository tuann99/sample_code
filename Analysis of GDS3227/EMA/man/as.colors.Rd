\name{as.colors}
\alias{as.colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert labels to colors}
\description{
  This function returns a object a colors values according to a palette function
}
\usage{
as.colors(x, col.na="#E6E6E6",palette="rainbow", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A object (vector or matrix) to convert in colors.}
  \item{col.na}{Colors to use for missing values.}
  \item{palette}{The palette function to use}
  \item{...}{Additional argument for the palette function. The 'n'
    parameters must exist in the palette function}
}

\value{
A vector or matrix of color label.
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
\examples{
     lab1<-c(rep("a",5), rep("b",5))
     lab2<-c(rep("c",5), rep("a",5))
     as.colors(lab1)
     as.colors(rbind(lab1, lab2))
     as.colors(lab2, palette="heat.colors", alpha=0.8)
}
\keyword{utilities}

