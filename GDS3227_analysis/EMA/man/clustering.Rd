\name{clustering}
\alias{clustering}
\title{Agglomerative hierarchical clustering}
\description{
  Computes agglomerative hierarchical clustering of the dataset.}
\usage{
  clustering(data, metric="euclidean", method="ward", nb)
}
\arguments{
  \item{data}{Expression matrix, genes on rows and samples on columns}
  \item{metric}{Character string specifying the metric to be used for calculating dissimilarities between the columns of the matrix. This must be one of 'euclidean', 'manhattan', 'pearson', 'pearsonabs', 'spearman', 'spearmanabs', 'jaccard', 'dice'}
  \item{method}{Character string defining the clustering method. This
    must be one of 'average', 'single', 'complete', 'ward'}
  \item{nb}{The number of classes for kmeans and PAM clustering (kcentroids)}
}
\details{
  Available metrics are (written for two vectors x and y):
  \describe{
    \item{euclidean:}{Usual square distance between the two vectors.}
    \item{manhattan:}{Absolute distance between the two vectors}
    \item{pearson:}{Pearson correlation distance. (1 - r)/2}
    \item{pearsonabs:}{Absolute Pearson correlation distance. 1 - abs(r)}
    \item{spearman:}{Spearman rank correlation distance. (1 - r)/2}
    \item{spearmanabs:}{Absolute Spearlan rnak correlation distance. 1 - abs(r)}
    \item{jaccard:}{Jaccard distance on 0-1 matrix}
    \item{dice:}{Dice distance on 0-1 matrix}
  }
  Available agglomerative methods are :
  \describe{
    \item{average:}{The distance between two clusters is the average of the dissimilarities between the points in one cluster and the points in the other cluster.}
    \item{single:}{we use the smallest dissimilarity between a point in the first cluster and a point in the second cluster (nearest neighbor method).}
    \item{complete:}{we use the largest dissimilarity between a point in the first cluster and a point in the second cluster}
    \item{ward:}{Ward's agglomerative method}
    \item{weighted:}{The weighted distance from the agnes package}
    \item{diana:}{computes a divise clustering}
    \item{kcentroids:}{Perform either kmeans clustering if the distance is
      euclidean or PAM clustering. The number of classes nb has to be done.}

  }
}
\value{
   An object of class 'agnes' representing the clustering. See 'agnes.object' for details.
}
\references{ Kaufman, L. and Rousseeuw, P.J. (1990). Finding Groups in Data: An Introduction to Cluster Analysis. Wiley, New York.}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{agnes}},\code{\link{clust.dist}}}
\examples{
data(marty)
c<-clustering(marty, metric="pearson", method="ward")
clustering.plot(c, title="Hierarchical Clustering\nPearson-Ward")
}
\keyword{cluster}
