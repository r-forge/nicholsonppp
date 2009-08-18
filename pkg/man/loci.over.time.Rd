\name{loci.over.time}
\Rdversion{1.1}
\alias{loci.over.time}
\title{loci.over.time}
\description{Plot allele frequency evolution over time for each locus, grouping by population color. This shows the difference between selection types and population colors.}
\usage{
loci.over.time(fr, pop.colors = pop.colors.default, generation = NULL, m = "Allele frequency evolution varies with selection type and population color")
}
\arguments{
  \item{fr}{Data frame of simulated allele frequencies.}
  \item{pop.colors}{Colors to use to distinguish populations (blue, neutral, red)}
  \item{generation}{Generation to emphasize, or NULL for no emphasis.}
  \item{m}{Main title.}
}
\details{
}
\value{The ggplot2 plot.}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{sim <- sim.drift.selection()
df <- sim2df(sim)
df.i <- interesting.loci(df)
stopifnot(nrow(df.i)/nrow(df)*attr(df,"parameters")$n.locus==6)
loci.over.time(df.i)}