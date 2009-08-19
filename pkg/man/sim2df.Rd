\name{sim2df}
\Rdversion{1.1}
\alias{sim2df}
\title{ Convert simulation result to data frame}
\description{Convert simulation data in list form to a huge data frame useful for plotting.}
\usage{
sim2df(L)
}
\arguments{
  \item{L}{Result of simulate.drift.selection.}
}
\details{
}
\value{Data frame containing all simulated data and parameter values, 1 line per (locus,population,generation) tuple. Has attribute "parameters" which is the 1-row data frame description of simulation parameters.}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{}