\name{runPCA}
\alias{runPCA}
\title{Perform an Principal Component Analysis}
\description{
  This function performs a Principal Component Analysis (PCA) and represents the samples or the variables of the analysis.
}
\usage{
runPCA(X, ncp=5, scale=TRUE, ind.sup=NULL, quanti.sup=NULL,quali.sup=NULL,
sample.qual=TRUE, variable.qual=FALSE, sample.cont=TRUE,variable.cont=FALSE,
plotSample=TRUE, plotVariable=FALSE, plotInertia = TRUE, plotBiplot=FALSE,
lab.sample="quality", lab.var=NULL,palette="rainbow",
lim.cos2.sample=0, lim.cos2.var=0, pdf=FALSE, pdfname= NULL, verbose=FALSE, ...) 
}
\arguments{
  \item{X}{a data frame with n rows (samples) and p columns (variables)}
  \item{ncp}{number of dimensions kept in the results (by default 5)}
  \item{scale}{a boolean, if TRUE (value set by default) then data are scaled to unit variance}
  \item{ind.sup}{a vector indicating the indexes of the supplementary individuals}
  \item{quanti.sup}{a vector indicating the indexes of the quantitative supplementary variables}
  \item{quali.sup}{a vector indicating the indexes of the qualitative supplementary variables}
  \item{sample.qual}{a boolean, if TRUE quality sample is displayed, by default = TRUE}
  \item{variable.qual}{a boolean, if TRUE quality variable is displayed, by default = FALSE}
  \item{sample.cont}{a boolean, if TRUE sample contribution is displayed, by default = TRUE}
  \item{variable.cont}{a boolean, if TRUE variable contribution is displayed, by default = FALSE}
  \item{plotSample}{a boolean, if TRUE samples are displayed, by default = TRUE}
  \item{plotVariable}{a boolean, if TRUE variables are displayed, by default = FALSE}
  \item{plotInertia}{a boolean, if TRUE inertia percentage of components is displayed, by default = TRUE}
  \item{plotBiplot}{a boolean, if TRUE biplot is displayed, by default = FALSE}
  \item{lab.sample}{a vector, sample representation is colored by label.sample, by default = NULL}
  \item{lab.var}{a vector, variable representation is colored by label.var, by default = "quality"}
  \item{palette}{character, name of palette color, by default = "rainbow"}
  \item{lim.cos2.sample}{a numeric, for graphics, keep samples with cos2 >= lim.cos2.sample, by default = 0}
  \item{lim.cos2.var}{a numeric, for graphics, keep variables with cos2 >= lim.cos2.var, by default = 0}
  \item{pdf}{a boolean, if TRUE save all the graphics in a pdf file, by default = FALSE}
  \item{pdfname}{pdf file name for saving graphics}
  \item{verbose}{print results if verbose = TRUE, by default = FALSE}
  \item{ ... }{Arguments to be passed to methods, such as graphical parameters (see 'par').}
}

\value{
  \item{eig}{a matrix containing all the eigenvalues, the percentage of variance and the cumulative percentage of variance}
  \item{var}{a list of matrices containing all the results for the active variables (coordinates, correlation between variables and axes, square cosine, contributions}
  \item{ind}{a list of matrices containing all the results for the active individuals (coordinates, square cosine, contributions)}

Returns the individuals factor map for axes 1 and 2, 1 and 3, 2 and 3
Returns the inertia percentage of components
By default, print sample coordinates, sample quality and sample contribution

}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link{plotSample}}, \code{\link{plotVariable}}, \code{\link{plotInertia}}}
\examples{
data("marty")

## PCA on sample with inertia plot and sample plot colored by tumour type

example.subset <- marty[1:100,]
\dontrun{
pca <- runPCA(t(example.subset), verbose = FALSE, lab.sample = marty.type.cl)
}
}
\keyword{multivariate}

