# load library
library(BayesianTools)

# load data
d = read.csv('Denv dissemination time course and day 7 data 8-1-23.csv',header=T)

alpha_di_sf = 3.89
beta_di_sf = 0.39
alpha_di_df = 1.76
beta_di_df = 0.24

# 1 - alpha_0
# 2 - alpha_df
# 3 - alpha_wAlb
# 4 - beta_0
# 5 - beta_df
# 6 - beta_wAlb

par.ind = paste(d$Treatment,d$Feeding)
par.ind.name = unique(par.ind)
# "wAlb SF" "wAlb DF" "WT SF"   "WT DF" 
par.ind = as.numeric(sapply(par.ind,function(ii)which(unique(par.ind)==ii)))

ll = function(par){
  sum(dbinom(
    x = d$Positive,
    size = d$Total,
    prob = pgamma(
      q = d$Day.post.infection,
      shape = par[par.ind],
      rate = par[4+par.ind]),
    log=T))
}

# density = function(par){
#   sum(dgamma(
#     x = par[1:4],
#     shape = rep(alpha_di_sf*0.1,4),
#     rate = rep(0.1,4),
#     log=T)) +
#   sum(dlnorm(
#     x = par[4+(1:4)],
#     meanlog = rep(log(beta_di_sf),4),
#     sdlog = rep(1,4),
#     log=T))
# }
# 
# sampler = function(par){
#   c(
#     rgamma(
#       n=4,
#       shape = rep(alpha_di_sf*0.1,4),
#       rate = rep(0.1,4)),
#     rlnorm(
#       n=4,
#       meanlog = rep(log(beta_di_sf),4),
#       sdlog = rep(1,4)))
# }

# load posterior samples from Armstrong et al. analysis for use in priors
posterior.armstrong = cbind(
  unlist(read.csv('armstrongprior_shapessfthin.csv')[,2]), # wAlb SF shape
  unlist(read.csv('armstrongprior_shapesdfthin.csv')[,2]), # wAlb DF shape
  unlist(read.csv('armstrongprior_shapessfthin.csv')[,2]), # WT SF shape
  unlist(read.csv('armstrongprior_shapesdfthin.csv')[,2]), # WT DF shape
  unlist(read.csv('armstrongprior_ratessfthin.csv')[,2]), # wAlb SF rate
  unlist(read.csv('armstrongprior_ratesdfthin.csv')[,2]), # wAlb DF rate
  unlist(read.csv('armstrongprior_ratessfthin.csv')[,2]), # WT SF rate
  unlist(read.csv('armstrongprior_ratesdfthin.csv')[,2])) # WT DF rate
scramble.order = sample(1:nrow(posterior.armstrong),nrow(posterior.armstrong),replace=F)
posterior.armstrong[,1:2] = posterior.armstrong[scramble.order,1:2]
posterior.armstrong[,5:6] = posterior.armstrong[scramble.order,5:6]

# prior on multiplier for shape (which unilaterally affects mean) due to wAlb
# based on https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0003894
# mean(c(rnorm(1e6,6.0,0.58)/rnorm(1e6,4.9,0.21),rnorm(1e6,7.8,0.77)/rnorm(1e6,5.9,0.29)))
prior.walb.mean = 1.276
# sd(c(rnorm(1e6,6.0,0.58)/rnorm(1e6,4.9,0.21),rnorm(1e6,7.8,0.77)/rnorm(1e6,5.9,0.29)))
prior.walb.sd = 0.1469
posterior.armstrong[,1] = posterior.armstrong[,1] *
  rnorm(nrow(posterior.armstrong),prior.walb.mean,prior.walb.sd)
posterior.armstrong[,2] = posterior.armstrong[,2] *
  rnorm(nrow(posterior.armstrong),prior.walb.mean,prior.walb.sd)

# create a prior density based on the Armstrong et al. posterior results
prior.armstrong = createPriorDensity(
  posterior.armstrong, method = 'multivariate', eps = 1e-10,
  lower = rep(1e-4,8), upper = rep(100,8), best = NULL)
# prior = createPrior(
#   density = density, sampler = sampler,
#   lower = rep(1e-4,8), upper = rep(100,8), best = NULL)

bayesianSetup = createBayesianSetup(
  likelihood = ll,
  prior = prior.armstrong)

iter = 1e5
settings = list(iterations = iter, message = FALSE, nrChains = 3,
                adapt = T, DRlevels = 1, gibbsProbabilities = NULL,
                temperingFunction = NULL, optimize = T)
out.full = runMCMC(
  bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)

summary(out.full)

correlationPlot(out.full)

params.map = MAP(out.full)$parametersMAP

ML.full = marginalLikelihood(out.full)$ln.ML
# -70.93477

samples.full = getSample(out.full,start=5e4)

correlationPlot(samples.full)

tracePlot(samples.full)

save(samples.full,file='mcmc_samples.RData')
