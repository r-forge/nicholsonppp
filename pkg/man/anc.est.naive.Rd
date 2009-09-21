\name{anc.est.naive}
\Rdversion{1.1}
\alias{anc.est.naive}
\title{ Naive ancestral allele frequency estimates}
\description{Estimate ancestral allele frequencies by simply taking the mean of all the current observations. This is a simple statistic we can use to reality check the model estimates, and which are fast enough to calculate for each generation in the simulation, facilitating animations.}
\usage{
anc.est.naive(fr)
}
\arguments{
  \item{fr}{Subset of simulation data frame containing the generation of interest from which we estimate the ancestral frequency.}
}
\details{
}
\value{Data frame with a row for each locus, summarizing naive ancestral allele frequencies in the ancestral.est column.}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{}