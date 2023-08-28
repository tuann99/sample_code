\name{plotSample}
\alias{plotSample}

\title{ Sample representation for Principal Component Analysis}
\description{
   Sample representation for Principal Component Analysis (PCA)
}
\usage{
plotSample(acp, axes = c(1, 2), new.plot = FALSE, lab = "quality",
palette="rainbow", lim.cos2.sample = 0, text = TRUE,
lab.title = NULL, ellipse=FALSE, ...)
}

\arguments{
  \item{acp}{result from PCA or do.pca function }
  \item{axes}{axes for sample representation, by default 1 and 2}
  \item{new.plot}{if TRUE, a new graphical device is created, by default = FALSE}
  \item{lab}{character. Sample label, by default = quality (points are labelled by
  quality index). If lab=NULL, no label is displayed.}
  \item{lim.cos2.sample}{keep samples with cos2 >= lim.cos2.sample, by default = 0}
  \item{palette}{characters. Name of a palette, By default, "rainbow" palette}
  \item{text}{add sample name or not, by default = TRUE}
  \item{lab.title}{title for the legend, by default = NULL}
  \item{ellipse}{if TRUE and lab provided, draw 95$\%$ confidence
  ellipse around barycentre of each group}
  \item{ ... }{Arguments to be passed to methods, such as graphical parameters (see 'par').}
}

\value{
  Sample representation on axes axes[1] and axes[2] colored by quality index (= cos2 of samples) or colored by lab
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}


\seealso{\code{\link{runPCA}},\code{\link{PCA}}}
\examples{
data(marty)

## PCA on sample - example set
example.subset <- marty[1:100,]
pca <- runPCA(t(example.subset), verbose = FALSE, plotInertia = FALSE, plotSample = FALSE)

## Sample plot of PCA object colored by tumour type
perso.colors <- colorRampPalette(c("red", "green"))
\dontrun{
plotSample(pca, lab = marty.type.cl, palette="perso.colors", ellipse=TRUE)
}
}
\keyword{plot}

