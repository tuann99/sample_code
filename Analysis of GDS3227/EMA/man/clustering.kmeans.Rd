\name{clustering.kmeans}
\alias{clustering.kmeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kmeans and hierarchical clustering}
\description{
  Kmeans clustering to summarize the genes information and hierarchical clustering on the kmeans' groups
}
\usage{
clustering.kmeans(data, N = 100, iter.max = 20,
title = "Kmeans - Hierarchical Clustering",
dist.s = "pearson", dist.g = "pearsonabs", method = "ward")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Expression matrix, genes on rows and samples on columns}
  \item{N}{The number of a priori clusters for the kmeans}
  \item{iter.max}{The maximum number of iterations allowed for the kmeans clustering}
  \item{title}{The plot title}
  \item{dist.s}{The distance used for the sample clustering}
  \item{dist.g}{The distance used for the genes clustering}
  \item{method}{The linkage used for both clusterings}
  
}
\details{
  The goal of this analysis is to first summarizes the genes information using the kmeans clustering. Then, a two-ways clustering is performed using the center of each kmean groups, and all the samples.
}
\value{
  A list with the kmeans object and the two hierarchical clusterings.
  \item{c.km}{An object of class 'kmeans'.}
  \item{c.sample}{An object of class 'agnes'. The hierarchical clustering on samples}
  \item{c.kcenters}{An object of class 'agnes'. The hierarchical clustering on the kmeans centers}
  
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{kmeans}},\code{\link{agnes}}}
\examples{
data(marty)
##Example on 100 genes for 5 classes
clustering.kmeans(marty[1:100,], N=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}

