\name{sim.drift.selection}
\Rdversion{1.1}
\alias{sim.drift.selection}
\title{ Simulate drift and selection}
\description{Simulate allele frequency evolution in several populations according to a simple drift and selection model.}
\usage{
sim.drift.selection(populations = 12, popsize = 1000, generations = 50, loci.per.s.value = 10, beta1 = 0.7, beta2 = 0.7, p.neutral = 0.8, s = c(1, 5)/100, adapt.pop = c(blue = 0.4, red = 0.4, neutral = 0.2))
}
\arguments{
  \item{populations}{Number of populations.}
  \item{popsize}{Size of each simulated population.}
  \item{generations}{Generations of evolution to simulate.}
  \item{loci.per.s.value}{Number of loci simulated for each s value.}
  \item{beta1}{Parameter for beta distribution of initial allele frequencies.}
  \item{beta2}{Parameter for beta distribution of initial allele frequencies.}
  \item{p.neutral}{Proportion of neutral alleles to simulate.}
  \item{s}{Vector of selection strength values to simulate.}
  \item{adapt.pop}{Probabilities of color adaptation in populations.}
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
\examples{}