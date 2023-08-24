\name{plotVariable}
\alias{plotVariable}

\title{Variable representation for Principal Component Analysis}
\description{
  Variable representation for Principal Component Analysis (PCA)
}
\usage{
plotVariable(acp, axes = c(1, 2), new.plot = FALSE, lab, lim.cos2.var =
0, palette="rainbow", ...)
}

\arguments{
  \item{acp}{result from PCA or do.pca function}
  \item{axes}{axes for variable representation, by default 1 and 2}
  \item{new.plot}{if TRUE, a new graphical device is created, by default = FALSE}
  \item{lab}{variable label}
  \item{palette}{character, name of color palette, by default = "rainbow"}
  \item{lim.cos2.var}{keep variables with cos2 >= lim.cos2.var}
  \item{ ... }{Arguments to be passed to methods, such as graphical parameters (see 'par').}
}

\value{
  Variable representation on axes axes[1] and axes[2]
  
  If PCA is normed, the correlation circle is plotted colored by lab
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{runPCA}},\code{\link{PCA}}}
\examples{
\dontrun{
data(marty)

## PCA on sample on 100 genes
## In practice see genes.selection
##mvgenes<-genes.selection(marty, thres.num=100)

pca <- runPCA(t(marty[1:100,]), verbose = FALSE, plotSample = FALSE,
    plotInertia = FALSE)
\dontrun{
## Variable plot of PCA object
\dontrun{
plotVariable(pca)
}
}
}
}

\keyword{plot}

