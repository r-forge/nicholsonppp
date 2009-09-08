\name{nppp.over.time}
\Rdversion{1.1}
\alias{nppp.over.time}
\title{nppp over time}
\description{Fit a series of models over time.}
\usage{
nppp.over.time(sim, generation = seq(0, sim$p$gen, l = 5)[-1])
}
\arguments{
  \item{sim}{Simulation to fit models on.}
  \item{generation}{Vector of generation numbers to fit on. By default, fit 4 equally spaced models.}
}
\details{
}
\value{List of model fit result lists, with attributes from mlply.}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{}