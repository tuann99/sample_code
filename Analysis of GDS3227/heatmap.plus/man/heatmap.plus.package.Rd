\name{heatmap.plus.package}
\alias{heatmap.plus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Heatmap with more sensible behavior }
\description{
  Call exactly as you would the base heatmap() function.  Only two optional parameters' values should differ, RowSideColors and ColSideColors.  They should be matrices instead of vectors, where colnames() on each matrix is a vector of the annotation track names.  One annotation track will be drawn per column.  See help('heatmap') for documentation of arguments not documented here.
}
\usage{
heatmap.plus(x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, distfun = dist, hclustfun = hclust, reorderfun = function(d, w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, verbose = getOption("verbose"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric matrix of the values to be plotted. }
  \item{Rowv}{determines if and how the \emph{row} dendrogram should be
    computed and reordered.  Either a \code{\link{dendrogram}} or a
    vector of values used to reorder the row dendrogram or
    \code{\link{NA}} to suppress any row dendrogram (and reordering) or
    by default, \code{\link{NULL}}, see \emph{Details} below.}
  \item{Colv}{determines if and how the \emph{column} dendrogram should be
    reordered.  Has the same options as the \code{Rowv} argument above and
    \emph{additionally} when \code{x} is a square matrix, \code{Colv =
      "Rowv"} means that columns should be treated identically to the rows.}
  \item{distfun}{function used to compute the distance (dissimilarity)
    between both rows and columns.  Defaults to \code{\link{dist}}.}
  \item{hclustfun}{function used to compute the hierarchical clustering
    when \code{Rowv} or \code{Colv} are not dendrograms.  Defaults to
    \code{\link{hclust}}.}
  \item{reorderfun}{function(d,w) of dendrogram and weights for
    reordering the row and column dendrograms.  The default uses
    \code{\link{reorder.dendrogram}}.}
  \item{add.expr}{expression that will be evaluated after the call to
    \code{image}.  Can be used to add components to the plot.}
  \item{symm}{logical indicating if \code{x} should be treated
    \bold{symm}etrically; can only be true when \code{x} is a square matrix.}
  \item{revC}{logical indicating if the column order should be
    \code{\link{rev}}ersed for plotting, such that e.g., for the
    symmetric case, the symmetry axis is as usual.}
  \item{scale}{character indicating if the values should be centered and
    scaled in either the row direction or the column direction, or
    none.  The default is \code{"row"} if \code{symm} false, and
    \code{"none"} otherwise.}
  \item{na.rm}{logical indicating whether \code{NA}'s should be removed.}
  \item{margins}{numeric vector of length 2 containing the margins
    (see \code{\link{par}(mar= *)}) for column and row names, respectively.}
  \item{ColSideColors}{(optional) character matrix with number of rows
    matching number of rows in x.  Each column is plotted as a row similar
    to heatmap()'s ColSideColors. colnames() will be used for labels if
    present. }
  \item{RowSideColors}{(optional)character matrix with number of rows
    matching number of columns in x.  Each column is plotted as a row
    similar to heatmap()'s RowSideColors. colnames() will be used for labels
    if present. }
  \item{cexRow, cexCol}{positive numbers, used as \code{cex.axis} in
    for the row or column axis labeling.  The defaults currently only
    use number of rows or columns, respectively.}
  \item{labRow, labCol}{character vectors with row and column labels to
    use; these default to \code{rownames(x)} or \code{colnames(x)},
    respectively.}
  \item{main, xlab, ylab}{main, x- and y-axis titles; defaults to none.}
  \item{keep.dendro}{logical indicating if the dendrogram(s) should be
    kept as part of the result (when \code{Rowv} and/or \code{Colv} are
    not NA).}
  \item{verbose}{logical indicating if information should be printed.}
  \item{\dots}{additional arguments passed on to \code{\link{image}},
    e.g., \code{col} specifying the colors.} %


}
\details{ }
\value{ }
\references{ }
\author{ Allen Day }
\note{ }
\seealso{ \code{\link{heatmap}} }
\examples{
  z = matrix(rnorm(30),nrow=5,ncol=6);
  rlab = matrix(as.character(c(1:5,2:6,3:7,4:8)),nrow=5,ncol=4);
  clab = matrix(as.character(c(1:6,6:1)),nrow=6,ncol=2);
  colnames(rlab) = LETTERS[1:dim(rlab)[2]];
  colnames(clab) = 1:dim(clab)[2];
  heatmap.plus(z,ColSideColors=clab,RowSideColors=rlab);
}
\keyword{ data }% at least one, from doc/KEYWORDS
