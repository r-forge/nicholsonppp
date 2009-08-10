library(R2WinBUGS)
data(schools)
J <- nrow(schools)
inits <- function(){
  list(theta=rnorm(J, 0, 100),
       mu.theta=rnorm(1, 0, 100),
       sigma.theta=runif(1, 0, 100))
}
ss<-bugs(list(J=J,
              y=schools$estimate,
              sigma.y=schools$sd),
         inits,
         c("theta", "mu.theta", "sigma.theta"),
         "/home/thocking/projects/bugs/schools.bugs")
