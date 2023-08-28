\name{marty}
\Rdversion{1.1}
\alias{marty}
\docType{data}
\title{
marty data
}
\description{
  marty data
}
\usage{data(marty)}
\format{
  A matrix with 54613 rows and 23 columns. Each row represents the
  expression level of a probeset for the 23 samples.
  The original .CEL files (hgu133plus2) were normalized using GCRMA and the AFFX probesets were discarded from the dataset.
}
\details{
The data are available on the ArrayExpress website (E-GEOD-13787) or in
the NCBI-GEO website (GEO - GSE13787).

}
\source{
  http://www.ncbi.nlm.nih.gov/projects/geo/query/acc.cgi?acc=GSE13787
}
\references{
Frequent PTEN genomic alterations and activated phosphatidylinositol
3-kinase pathway in basal-like breast cancer cells. Marty B, Maire V,
Gravier E, Rigaill G, Vincent-Salomon A, Kappler M, Lebigot I, Djelti F,
Tourdes A, Gestraud P, Hupe P, Barillot E, Cruzalegui F, Tucker GC,
Stern MH, Thiery JP, Hickman JA, Dubois T. Breast Cancer Res. 2008;10(6):R101. Epub 2008 Dec 3.
}
\examples{
data(marty)
dim(marty)
class(marty)
}
\keyword{datasets}
