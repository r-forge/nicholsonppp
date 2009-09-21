\name{nicholsonppp}
\Rdversion{1.1}
\alias{nicholsonppp}
\title{ Enhanced Nicholson model}
\description{Fit the Nicholson model and calculate PPP-values, using MCMC implented using fast fortran code. This function is basically just a convenient wrapper around .Fortran() which calls the compiled fortran code, then post-processes the result.  The model is discussed in Nicholson, et. al., Assessing population differentiation and isolation from single-nucleotide polymorphism data. J. R. Statist. Soc. B (2002) 64, Part 4, pp. 695-715.}
\usage{
nicholsonppp(YY, NN = 1000, beta_pi = 1, seed = 4501, nvaleurs = 1000, thin = 50, burn_in = 5000, npilot = 50, pilot_length = 1000, delta_a_init = 0.4, delta_p_init = 1.25, delta_c_init = 0.06, rate_adjust = 1.25, acc_inf = 0.2, acc_sup = 0.45, out_option = 1)
}
\arguments{
  \item{YY}{Matrix of allele counts, one for each (locus,population) pair. Used for nmrk, npop, Y_OBS.}
  \item{NN}{Number of total markers.}
  \item{beta_pi}{Beta parameter for the model's prior.}
  \item{seed}{Random seed for mersenne twister.}
  \item{nvaleurs}{Number of observations to sample in MCMC.}
  \item{thin}{Thinning rate.}
  \item{burn_in}{Burn in period length.}
  \item{npilot}{Max number of pilot runs.}
  \item{pilot_length}{Pilot run length.}
  \item{delta_a_init}{Random walk delta for alpha.}
  \item{delta_p_init}{Random walk delta for pi.}
  \item{delta_c_init}{Random walk delta for c.}
  \item{rate_adjust}{Pilot run adjustment rate.}
  \item{acc_inf}{Target acceptance rate for rejection sampling, min value.}
  \item{acc_sup}{Target acceptance rate for rejection sampling, max value.}
  \item{out_option}{How much detail to print? (deprecated)}
}
\details{
}
\value{List of estimated parameters of the Nicholson model.}
\references{
}
\author{Toby Dylan Hocking <toby.hocking@etu.upmc.fr>}
\note{
}


\seealso{
}
\examples{library(nicholsonppp)

## A really simple test data set:
observed <- matrix(rbeta(100,0.7,0.7)*100,nrow=20,ncol=5)
fit <- nicholsonppp(observed)
df <- cbind(melt(fit$a),melt(observed)$value)
names(df) <- c("locus","population","estimated","observed")
xyplot(estimated~observed,df)

## A data set that comes from an evolution simulation:
sim <- sim.drift.selection(populations=4,loci=1)
ev.obs <- sim$sim[,,sim$p$gen]
ev.fit <- nicholsonppp(ev.obs)
## Compare alpha estimates to simulated values:
head(round(ev.obs*100,2))
head(round(ev.fit$a*100,2))
ev.df <- sim2df(sim)
ev.df <- ev.df[ev.df$generation==sim$p$gen,]
ev.df$estimated <- as.vector(ev.fit$a)
plot(xyplot(simulated~estimated,ev.df))
## Plot ancestral estimates given by model fit:
s <- anc.est.nicholson(sim,ev.fit)
plot(anc.est.plot(s))}