\name{eval.stability.clustering}
\alias{eval.stability.clustering}
\title{Compares several clustering methods by means of its stability.}
\description{This function compares several clustering methods (link functions and distances) by means of its stability.}
\usage{
eval.stability.clustering(X,nb=c(2:4),f=0.8,nsub=10,s0=0.98,
list_DIS=c("euclidean","pearson"),
list_ALGO=c("average","complete","ward"), pdfname = NULL,
verbose = TRUE)
}
\arguments{
  \item{X}{a data frame with p rows and n columns; if clustering on genes - samples by row and genes by column; if clustering on samples genes by row and samples by column }
  \item{nb}{number of classes for partition; it must start at 2 and be sequential(by default 2,3 and 4)}
  \item{f}{part of the data set which is randomly picked for each subsample in the resampling procedure (by default 0.8)}
  \item{nsub}{half of the number of times the perturbation procedure is applied in the resampling procedure (by default 100)}
  \item{list_DIS}{the list of distances to test}
  \item{list_ALGO}{the list of linkage method to test}
  \item{s0}{similarity threshold, must lie between 0 and 1 (by default 0.98)  }
  \item{pdfname}{pdf file name for saving graphic, by default = NULL}
  \item{verbose}{print results if verbose = TRUE, by default = TRUE}
}

\details{Resampling is done by randomly picking without replacement f of
  the data set; similarity threshold is the value which is
  pertinent to decide that two partitions are similar; see references}

\value{
  \item{stab.methods}{a list containing methods declared stable for each partition}

Returns a graphic containing the frequencies of methods declared stable
}

\references{
  http://bioinfo-out.curie.fr/projects/cgh-clustering/index.html
}


\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}
\seealso{\code{\link{clustering}},\code{\link{clust.dist}}}

\examples{
data(marty)
## Test on a smaller dataset
\dontrun{
example.data<-marty[1:100,]
stab<-eval.stability.clustering(example.data) 
}
}
\keyword{cluster}

