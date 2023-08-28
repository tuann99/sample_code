\name{marty.type.cl}
\Rdversion{1.1}
\alias{marty.type.cl}
\docType{data}
\title{
marty class data for Basal vs HER2 cancer type
}
\description{
  marty class data for Basal vs HER2 cancer type
}
\usage{data(marty)}
\format{
 A vector of 0/1 according to the sample type. The HER2+ samples are
 labeled as 0 and the Basal-like samples as 1.
}
\details{
 The clinical data are available on the ArrayExpress website (E-GEOD-13787) or in
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
## maybe str(marty) ; plot(marty) ...
}
\keyword{datasets}
