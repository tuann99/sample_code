\name{clust.dist}
\alias{clust.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Computes distances on a data matrix }
\description{
  Computes distances between the columns of an expression matrix using the specified distance measure.
}
\usage{
clust.dist(mat, meth.dis="euclidean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{Expression matrix, genes on rows and samples on columns}
  \item{meth.dis}{Character string specifying the distance measure to be used. This must be one of 'euclidean', 'manhattan', 'pearson', 'pearsonabs', 'spearman', 'spearmanabs', 'jaccard', 'dice'}
}

\details{
  \itemize{
    Available distance measures are (written for two vectors x and y):
    \item 'euclidean' : Usual square distance between the two vectors.
    \item 'manhattan' :  Absolute distance between the two vectors
    \item 'pearson' : Pearson correlation distance. (1 - r)/2
    \item 'pearsonabs' : Absolute Pearson correlation distance. 1 - abs(r)
    \item 'spearman' : Spearman rank correlation distance. (1 - r)/2
    \item 'spearmanabs' :  Absolute Spearlan rnak correlation distance. 1 - abs(r)
    \item 'jaccard' : Jaccard distance on 0-1 matrix
    \item 'dice': Dice distance on 0-1 matrix
  }
}
\value{
 A distance matrix.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{ \code{\link{jaccard}},\code{\link{dice}},\code{\link{daisy}}}
\keyword{cluster}
\keyword{internal}