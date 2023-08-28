\name{sample.plot}
\alias{sample.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{barplot of genes expression level}
\description{
  Plot a barplot of the expression level of each gene of interest.
}
\usage{
sample.plot(data, labels=NULL, plot=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Expression matrix, genes on rows and samples on columns}
  \item{labels}{A character string or numeric vector of label}
  \item{plot}{Display output barplot if TRUE}
  \item{...}{Arguments to be passed to methods, such as graphical parameters (see 'par').}

}
\details{
  For each gene (row of the matrix), a barplot of the expression level
  for all the samples is ploted. The colors are chosen according to the
  label information.
}


\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\examples{

data(marty)
sample.plot(marty[1:3,], labels=marty.type.cl)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}

