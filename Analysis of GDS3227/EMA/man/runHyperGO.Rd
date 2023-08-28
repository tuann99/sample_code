\name{runHyperGO}
\alias{runHyperGO}

\title{Run Gene Ontology analysis based on hypergeometric test from a probeset list}
\description{
  Run Gene Ontology analysis based on hypergeometric test from a probeset list 
}
\usage{
runHyperGO(list, pack.annot, categorySize = 1, verbose = TRUE,
name = "hyperGO", htmlreport = TRUE, txtreport = TRUE,
tabResult = FALSE, pvalue = 0.05)
}

\arguments{
  \item{list}{vector of character with probeset names}
  \item{pack.annot}{annotation package to use}
  \item{categorySize}{integer, minimum size for category, by default = 1}
  \item{verbose}{logical, if TRUE, results are displayed, by default TRUE}
  \item{name}{character, name for output files, by default "hyperGO"}
  \item{htmlreport}{logical, if TRUE, a html report is created, by default TRUE}
  \item{txtreport}{logical, if TRUE, a txt report is created, by default TRUE}
  \item{tabResult}{logical, if TRUE, a list with the results is created, by default FALSE}
  \item{pvalue}{numeric, a cutoff for the hypergeometric test pvalue, by default 0.05}
}

\details{
  The choice of the universe could have a significant impact on the
  results. It is well discussed in the vignette of the GOstats package.
  Here, we decided to apply a non-specific filtering procedure different from the one proposed by Falcon and Gentleman.
  Since not all genes will be expressed under all conditions in our data, we can ask the question of defining the universe only with the expressed genes or with all the genes of the array.
  Actually, we are not able to distinguish the genes which are biologically non expressed, from the ones of low quality.
  That's why we think that the non-expressed probesets could be biologically relevant, as well as the ones with a little variation accross samples, and we decided to first defined the universe with all the genes of the array. 
  Then, we just remove probe sets that have no Entrez Gene identifier in our annotation data or no GO annotation. Finally, the Hypergeometric test is performed on the unique EntrezId of the gene list, and the unique EntrezId of the universe.
  The pvalues in output are not corrected from multiple testing.
  Note that because of the existing dependence structure (between genes, and GO terms) it is difficult to do any multiple testing correction. 
  Moreover the most insteresting genesets are not necessarily the ones with the smallest pvalues. Nodes that are interesting are typically those with a reasonable number of genes (10 or more)
  and small pvalues. 

  \code{runHyperGO} needs packages \code{GOstats} and \code{GO.db} from Bioconductor.

}

\value{
  The R objects or the Txt and html reports
  \item{BP}{Data.frame with results for Biological Process with GO Id, pvalue, Odd Ratio, Expected count, Size and GO Term}
  \item{MF}{Idem for Molecular Function}
  \item{CC}{Idem for Cellular Component}
}

\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link[GOstats]{hyperGTest}},\code{\link{runHyperKEGG}}}
\examples{
\dontrun{
require(hgu133plus2.db)
data(marty)

## Probe list
probeList <- rownames(marty)[1:50]

## Hypergeometric test for GO pathway
res <- runHyperGO(probeList, htmlreport = FALSE, txtreport = FALSE,
    tabResult = TRUE, pack.annot = "hgu133plus2.db")
}
}

\keyword{htest}

