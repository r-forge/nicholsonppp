\name{anc.est.plot}
\Rdversion{1.1}
\alias{anc.est.plot}
\title{ Ancestral estimate plot}
\description{Plot naive estimates of ancestral allele frequency versus actual values from the simulation, to see if they agree.}
\usage{
anc.est.plot(est.df, hilite.locus = NULL, sub = deduce.param.label(attr(est.df, "parameters")), selection.colors = selection.colors.default, selection.symbols = selection.symbols.default, f = ancestral.est ~ ancestral, xlab = "Simulated blue allele frequency", ylab = "Estimated blue allele frequency", main = "Allele frequency estimates vary with selection type")
}
\arguments{
  \item{est.df}{Data frame describing loci and ancestral estimates. Need columns type ancestral ancestral.est.}
  \item{hilite.locus}{Locus number to highlight on the plot with a circle. NULL means do not highlight.}
  \item{sub}{Subtitle for the plot.}
  \item{selection.colors}{Colors for the different selection types (balancing, neutral, positive)}
  \item{selection.symbols}{Symbols for the different selection types (balancing, neutral, positive)}
  \item{f}{Plot formula for xyplot.}
  \item{xlab}{
}
  \item{ylab}{
}
  \item{main}{
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
\examples{library(nicholsonppp)
sim <- sim.drift.selection()
df <- sim2df(sim)
est <- anc.est.naive(df[df$generation==sim$p$gen,])
anc.est.plot(est)}