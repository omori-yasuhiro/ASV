#---
# Generate simulation data
#--
# set.seed(111)
nobs = 2000;
mu = 0; phi = 0.97; sigma_eta = 0.3; rho = -0.3;

h = 0; Y = c();
for(i in 1:nobs){
  eps = rnorm(1, 0, 1)
  eta = rho*sigma_eta*eps + sigma_eta*sqrt(1-rho^2)*rnorm(1, 0, 1)
  y   = eps * exp(0.5*h)
  h   = mu + phi * (h-mu) + eta
  Y   = append(Y, y)
}
iter = 1500; burn = 500;

# asv_mcmc: function to implement MCMC for SV model with leverage
# asv_mcmc(return vecor, # of iterations, # of burn-in period)
# return: 
out  = asv_mcmc(Y, nsim, burn)
vmu  = out[[1]]; vphi = out[[2]]; vsigma_eta =  out[[3]]; vrho = out[[4]]; 
mh   = out[[5]];

ReportMCMC(cbind(vmu,vphi,vsigma_eta, vrho), vname=c(expression(mu), expression(phi),expression(sigma[eta]), expression(rho)), soutfilename=c("ASV"))
ReportMCMC(cbind(mh[,500], mh[,1000],mh[,1500]), vname=c(expression(h[500]), expression(h[1000]),expression(h[1500])), soutfilename=c("ASV"))

mh_ci = matrix(0, nrow = nobs, ncol = 3)
for(i in 1:nobs){
  mh_ci[i,1:3] = t(quantile(mh[1:iter,i],c(0.025,0.5, 0.975)))
}

library(ggplot2)
mydate = seq(1,nobs)
ggplot()+geom_line(aes(x=mydate, y=mh_ci[,1], color="95%L"),size=1)+geom_line(aes(x=mydate, y=mh_ci[,2], color="Median"),size=1)+geom_line(aes(x=mydate, y=mh_ci[,3], color="95%U"),size=1)+labs(x="Date", y="Log Volatility")

# Particle filter: Log likelihood given theta = (mu, phi, sigma_eta)
mu = 0; phi = 0.97; sigma_eta = 0.3; rho = -0.3;
# Number of Particles
npart = 1000
# repeat m times to compute the standard error of the estimate
m  = 5; loglik = rep(0, m)

for(i in 1:m){
  loglik[i] = asv_pf(mu, phi, sigma_eta, rho, Y, npart)
}
# Estimate and its standard error
loglik
c(mean(loglik), sqrt(var(loglik)/m))


