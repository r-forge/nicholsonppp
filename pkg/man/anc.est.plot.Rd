\name{anc.est.plot}
\Rdversion{1.1}
\alias{anc.est.plot}
\title{anc.est.plot}
\description{Plot naive estimates of ancestral allele frequency versus actual values from the simulation, to see if they agree.}
\usage{
anc.est.plot(fr, hilite.locus = NULL, sub = deduce.param.label(attr(fr, "parameters")), selection.colors = selection.colors.default, selection.symbols = selection.symbols.default)
}
\arguments{
  \item{fr}{Subset of simulation data frame containing the generation of interest to estimate and plot.}
  \item{hilite.locus}{Locus number to highlight on the plot with a circle. NULL means do not highlight.}
  \item{sub}{Subtitle for the plot.}
  \item{selection.colors}{Colors for the different selection types (balancing, neutral, positive)}
  \item{selection.symbols}{Symbols for the different selection types (balancing, neutral, positive)}
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
\examples{sim <- sim.drift.selection()
df <- sim2df(sim)
anc.est.plot(subset(df,generation==generations))

}