\name{distrib.plot}
\alias{distrib.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distribution plots of genes expression level}
\description{
   Plot the distribution of the expression level of each gene of interest.
}
\usage{
distrib.plot(data, labels = NULL, plot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Expression matrix, genes on rows and samples on columns}
  \item{labels}{A character string or numeric vector of label}
  \item{plot}{If true, plots are displayed}
  \item{...}{Arguments to be passed to methods, such as graphical parameters (see 'par').}

}
\details{
  For each gene (row of the matrix), the distribution of the expression level
  for all the samples is ploted. The colors are chosen according to the
  label information.
}
\value{
  A list of objects of class 'histogram':
  \item{mids }{The n cell midpoints}
  \item{counts }{n integers; for each cell, the number of 'x[]' inside.}
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{ \code{\link{hist}}}
\examples{

data(marty)

distrib.plot(marty[1:5,], labels=marty.type.cl)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
