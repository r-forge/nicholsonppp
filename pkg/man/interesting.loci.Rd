\name{interesting.loci}
\Rdversion{1.1}
\alias{interesting.loci}
\title{interesting loci}
\description{Find a small subset of loci for each selection type, from a variety of ancestral frequencies. This will be used to show how allele frequency evolution depends on selection type and population color.}
\usage{
interesting.loci(fr)
}
\arguments{
  \item{fr}{Data frame of all the simulated allele frequencies.}
}
\details{
}
\value{Data frame, subset of input data.}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{library(nicholsonppp)
sim <- sim.drift.selection()
s.i <- interesting.loci(sim$s)
print(s.i)
stopifnot(nrow(s.i)==6)}