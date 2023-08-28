\name{normAffy}
\alias{normAffy}
\title{Normalisation of Affymetrix expression arrays}
\description{
  This function converts .cel files into a matrix of normalised data
  using GCRMA, RMA or MAS5 normalisation.
}
\usage{
normAffy(filenames, celfile.path, method = c("GCRMA","RMA","MAS5"),
cdfname = NULL, rmaffx = TRUE, fast=TRUE)
}
\arguments{
  \item{filenames}{The .CEL files to normalize separated by comma, if empty the celfile.path.}
  \item{celfile.path}{Path to the directory containing the cel files}
  \item{method}{Normalisation method to use}
  \item{cdfname}{Used to specify the name of an alternative cdf package. If set to NULL, the usual cdf package based on Affymetrix' mappings will be used. Note that the name should not include the 'cdf' on the end, and that the corresponding probe package is also required to be installed. If either package is missing an error will result.}
  \item{rmaffx}{boolean, indicating wether the probe control should be
    removed}
  \item{fast}{logical, option of GCRMA normalisation. see details}
}
\details{
  This function allows to use several normalisation methods. It is a
  wrapper for functions \code{justRMA},\code{justGCRMA} and \code{mas5} from packages \code{affy}
  and \code{gcrma}.\\
  The \code{fast} parameter of the \code{justRMA} can lead to different
  results. If \code{fast} is true, the Lim et al. correction is
  applied. All expression values smaller than a threshold are put to the
  same values, leading to a sharp pick on the left of the
  distribution. Set \code{fast} to false to obtain a standard GCRMA distribution.

  If the filenames parameter is missing, then all the files of the
  celfile.path are used.
}
\value{
  A matrix containing the normalised data, genes on rows and samples on columns.
}
\examples{
\dontrun{
## GCRMA normalisation
normData <- normAffy(celfile.path="PathToCelFiles", method="GCRMA")
}

}
\keyword{manip}
