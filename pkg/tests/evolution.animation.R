library(nicholsonppp)
sim <- simulate.drift.selection()
df <- sim2df(sim)
evolution.animation(tempfile(),
                    "Allele frequency and ancestral estimate evolution",df)
