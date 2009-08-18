\name{evolution.animation}
\Rdversion{1.1}
\alias{evolution.animation}
\title{evolution.animation}
\description{Create an animation that summarizes a simulation, using sim.summary.plot for every generation of the simulation.}
\usage{
evolution.animation(subdir, df, ...)
}
\arguments{
  \item{subdir}{Subdirectory for plot files, to be created.}
  \item{df}{Result of sim2df.}
  \item{\dots}{
}
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
\examples{sim <- sim.drift.selection()
df <- sim2df(sim)
evolution.animation(tempfile(),
                    "Allele frequency and ancestral estimate evolution",df)}