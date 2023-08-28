\name{expFilter}
\alias{expFilter}
\title{Filter expression data}
\description{
  This function takes an expression matrix and filtered out non detected
  genes. 
}
\usage{
expFilter(data, threshold = 3.5, p=0.01, graph = TRUE)
}

\arguments{
  \item{data}{expression matrix, genes on rows and samples on columns.}
  \item{threshold}{minimal value of expression to be reached}
  \item{p}{keep probes with at least p*ncol(data) samples higher than threshold}
  \item{graph}{boolean indicating if an histogram of the data should be
    plotted.}
}
\details{
  The non variant genes are defined gy the \code{threshold} value. A
  gene is kept if at least p*ncol(data) of its values is higher than
  \code{threshold}.

  The graph represents the distribution of all the genes in
  \code{Data}. A line shows the threshold value used.
}
\value{
  An expression matrix.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\examples{

data(marty)
expFilter(marty, threshold = 3.5, graph = TRUE)

}
\keyword{manip}
