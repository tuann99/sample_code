\name{bioMartAnnot}
\alias{bioMartAnnot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Annotation of probesets using biomaRt}
\description{
  Annot input IDs using the biomaRt package.
}
\usage{
bioMartAnnot(data, inputTypeId, outputTypeId =c("entrezgene","hgnc_symbol",
"ensembl_gene_id", "description", "chromosome_name", "start_position",
"end_position", "band", "strand") , dataset= c("hsapiens_gene_ensembl"),
database = "ensembl",  sort.by = NULL, outfile = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{A data.frame with the input names as rownames or a vector of
probesets}
\item{inputTypeId}{Input Ids type. see details}
\item{outputTypeId}{Output Ids type. if NULL the description of some
genes names are returned. see details}
\item{dataset}{The dataset used for the annotation}
\item{database}{The database used for the annotation}
\item{sort.by}{optional. Allow to sort the output data.frame}
\item{outfile}{optional. Write an html file with the results}
}
\details{
This function is based on the biomaRt package.
First, you have to define the database and the dataset you want to use
for the annotation (default dataset= c("hsapiens_gene_ensembl"), database = "ensembl").
The input IDs type have to be defined as in the biomaRt package. The
listFilters() functions of the biomaRt package lists the available input
type. If the inputs are probe names, use the microarray name :
affy_hg_u95a,affy_hg_u95av2,affy_hg_u133a_2,affy_hg_u133a,
affy_hg_u133b,affy_hg_u133_plus_2, illumina_humanwg_6_v2, ...
Use the biomaRt function 'listAttributes()', to select your output. By
default, the function returns  "description", "chromosome_name",
"start_position", "end_position", "band", "strand", and "ensembl_gene_id".
}


\author{Nicolas Servant, Eleonore Gravier, Pierre Gestraud, Cecile Laurent, Caroline Paccard, Anne Biton, Jonas Mandel, Bernard Asselain, Emmanuel Barillot, Philippe Hupe}

\seealso{\code{\link[biomaRt]{useMart}}}
\examples{
\dontrun{
## load data
data(marty)

##Annotations example
bioMartAnnot(rownames(marty)[1:50], inputTypeId ="affy_hg_u133_plus_2",
   outputTypeId =c("entrezgene","hgnc_symbol"),
   dataset=c("hsapiens_gene_ensembl"),database = "ensembl")
}
}
\keyword{datagen}
