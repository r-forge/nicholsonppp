\name{fixation.endpoints}
\Rdversion{1.1}
\alias{fixation.endpoints}
\title{ Plot allele frequency by selection type}
\description{Plot gene frequencies for all loci and populations for a given generation, stratified by selection type and coefficient.}
\usage{
fixation.endpoints(lf, main = "Loci fixation depends on selection type and population color", pop.colors = c(pop.colors.default, "black"), par.settings = list(superpose.symbol = list(col = pop.colors, pch = 20)), hilite.locus = NULL, selection.colors = selection.colors.default, sub = deduce.param.label(attr(lf, "parameters")), ...)
}
\arguments{
  \item{lf}{Subset of simulated gene frequency data frame, with just 1 generation.}
  \item{main}{Plot title.}
  \item{pop.colors}{Population color scheme (blue, neutral, red, ancestral).}
  \item{par.settings}{Plot settings as described in trellis.par.get.}
  \item{hilite.locus}{Locus to highlight with a vertical line, or NULL to highlight nothing.}
  \item{selection.colors}{List with element "col" which contains a vector of colors to label the selection types (balancing, none, positive)}
  \item{sub}{Subtitle for the plot.}
  \item{\dots}{
}
}
\details{
}
\value{The lattice plot.}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{}