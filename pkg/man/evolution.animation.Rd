\name{evolution.animation}
\Rdversion{1.1}
\alias{evolution.animation}
\title{evolution animation}
\description{Create an animation that summarizes a simulation, using sim.summary.plot for every generation of the simulation.}
\usage{
evolution.animation(df, outdir = tempfile(), tit = "Allele frequency and ancestral estimate evolution")
}
\arguments{
  \item{df}{Result of sim2df.}
  \item{outdir}{Subdirectory for plot files, to be created.}
  \item{tit}{Title for the animation.}
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
evolution.animation(df[df$generation\%\%10==1,])}