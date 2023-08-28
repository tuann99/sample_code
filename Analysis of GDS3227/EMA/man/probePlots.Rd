\name{probePlots}
\alias{probePlots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot the expression profiles of the probes corresponding to given probesets }
\description{
  Performs expression plots of probes for giver probesets.
  This function is interesting to show which probes are responsible of
  intensity signal.
}
\usage{
probePlots(abatch, path, pbsList, labAxisProbes=TRUE, labAxisArrays=TRUE,
legendArrays=TRUE, legendProbes=TRUE, cex.axis=0.9, cex.legend=0.8, pdfName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{abatch}{An affybatch object.}
  \item{path}{A character. If no affyBatch object, the path where are
    .CEL files  }
  \item{pbsList}{A character. A vector of charcter listing the
    interesting probesets. }
  \item{labAxisProbes}{A logical. TRUE if probes names have to be
    plotted in the x axis, FALSE otherwise.}
  \item{labAxisArrays}{A logical. TRUE if arrays names have to be
    plotted in the x axis, FALSE otherwise.}
  \item{legendArrays}{A logical. TRUE if legend corresponding to the arrays names have to be
    plotted, FALSE otherwise.}
  \item{legendProbes}{A logical. TRUE if legend corresponding to the probes names have to be
    plotted, FALSE otherwise.}
  \item{cex.axis}{A numeric. The magnification to be used for x axis relative
    to the current setting of 'cex'}
  \item{cex.legend}{A numeric. The magnification to be used for legend relative
    to the current setting of 'cex'}
  \item{pdfName}{A character. Name for a pdf file if needed.}

  
}
\value{
  Three plots for each probesets are generated.
  The first one is an inter chips plot. For each arrays in the
  AffyBatch, the perfect match intensity of probesets' probes are plotted
  The second one is an inter probes plot on the perfect match. For each
  probes in the probeset, the perfect match intensity in each arrays are plotted.
  The third one is an inter probes plot on the mis match. For each
  probes in the probeset, the mis match intensity in each arrays are plotted.
}
\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\examples{
\dontrun{
if (require(affydata)) {
## load example
data(Dilution)
probePlots(Dilution, pbsList=geneNames(Dilution)[1:2])
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
