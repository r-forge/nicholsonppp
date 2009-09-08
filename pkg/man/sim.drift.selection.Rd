\name{sim.drift.selection}
\Rdversion{1.1}
\alias{sim.drift.selection}
\title{ Simulate drift and selection}
\description{Simulate allele frequency evolution in several populations according to a simple drift and selection model.}
\usage{
sim.drift.selection(populations = 12, popsize = 1000, generations = 50, loci.per.s.value = 10, beta1 = 0.7, beta2 = 0.7, p.neutral = 0.8, s = c(1, 5)/100, adapt.pop = c(blue = 0.4, red = 0.4, neutral = 0.2), array.fun = array)
}
\arguments{
  \item{populations}{Number of populations.}
  \item{popsize}{Size of simulated populations. Vectors smaller than number of populations will be repeated.}
  \item{generations}{Generations of evolution to simulate.}
  \item{loci.per.s.value}{Number of loci simulated for each s value.}
  \item{beta1}{Parameter for beta distribution of initial allele frequencies.}
  \item{beta2}{Parameter for beta distribution of initial allele frequencies.}
  \item{p.neutral}{Proportion of neutral alleles to simulate.}
  \item{s}{Vector of selection strength values to simulate.}
  \item{adapt.pop}{Probabilities of color adaptation in populations.}
  \item{array.fun}{Function to use to construct the evolution array. By default we use a normal R array, but this can be any function that accepts data as the first argument, and a dim argument, and gives a result that behaves like an array. In particular you can use the ff function from the ff package if the simulation is too large for memory.}
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