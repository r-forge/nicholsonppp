\name{sim.summary.plot}
\Rdversion{1.1}
\alias{sim.summary.plot}
\title{ Simulation summary plot}
\description{Draw 3 simulation summary plots on the same screen (loci.over.time, anc.est.plot, fixation.endpoints).}
\usage{
sim.summary.plot(fr, g = NULL, hilite.locus = NULL)
}
\arguments{
  \item{fr}{Data frame from a simulation.}
  \item{g}{Generation to plot in anc.est.plot and fixation.endpoints, and generation to emphasize in loci.over.time. NULL means the last generation in the simulation.}
  \item{hilite.locus}{Locus to highlight in the plots. NULL means choose an interesting locus under positive selection.}
}
\details{
}
\value{
}
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
sim.summary.plot(df,ceiling(sim$p$gen/3*2))}