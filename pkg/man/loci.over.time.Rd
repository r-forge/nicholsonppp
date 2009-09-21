\name{loci.over.time}
\Rdversion{1.1}
\alias{loci.over.time}
\title{ Allele frequency time series plots}
\description{Plot allele frequency evolution over time for each locus, grouping by population color. This shows the difference between selection types and population colors.}
\usage{
loci.over.time(fr, pop.colors = pop.colors.default, generation = NULL, m = "Allele frequency evolution varies with selection type and population color", popsize.enc = NULL)
}
\arguments{
  \item{fr}{Data frame of simulated allele frequencies.}
  \item{pop.colors}{Colors to use to distinguish populations (blue, neutral, red)}
  \item{generation}{Generation to emphasize, or NULL for no emphasis.}
  \item{m}{Main title.}
  \item{popsize.enc}{Aesthetic to encode population size. lty or lwd work well. NULL means do not show population size.}
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
\examples{library(nicholsonppp)
sim <- sim.drift.selection(s=0.5,generations=100)
df <- sim2df(sim)
df.i <- interesting.loci(df)
stopifnot(nrow(df.i)/nrow(df)*sim$p$n.locus==6)
loci.over.time(df.i)}