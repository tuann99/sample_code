\name{genes.selection}
\alias{genes.selection}

\title{Genes selection}
\description{
Selects the most variant genes of a expression dataset. For each gene,
the difference between the quantile with probability (1 - probs) and the
quantile with probability probs is computed. Default for probs value is
set to 0.25 such that the difference corresponds to the Inter Quartile
Range (IQR).
}
\usage{
genes.selection(data, thres.diff, thres.num, probs=0.25)
}

\arguments{
  \item{data}{A expression matrix, genes on rows and samples on columns}
  \item{thres.diff}{Difference threshold - Genes are selected based on the difference threshold}
  \item{thres.num}{Genes number threshold - Selects the N genes with the
    highest difference value}
  \item{probs}{probability value used to compute both quantiles}  
}
\details{
The difference is computed for each genes. If the thres.diff option is
chosen, the most variant genes are selected according to the difference
threshold. If the thres.num option is chosen, the genes are ordered
according to their difference value level, and the N first genes are selected.
}
\value{
 The name of the selected genes.
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{ \code{\link{IQR}}, \code{\link{quantile}}}
\examples{
data(marty)
data.f<-expFilter(marty, graph=FALSE)

##Select the 50 most variant genes amongst the 1000 first genes
sel<-genes.selection(data.f[1:1000,], thres.num=50)
}
\keyword{manip}

