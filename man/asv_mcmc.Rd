\name{asv_mcmc}
\alias{asv_mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC estimation for stochastic volatility models with leverage
%%  ~~function to do ... ~~
}
\description{
This function estimates model parameters and latent log volatilities for 
stochastic volatility models with leverage:  

y(t) = eps(t)*exp(h(t)/2), h(t+1) = mu + phi*(h(t)-mu) + eta(t)

eps(t)~i.i.d. N(0,1), eta(t)~i.i.d. N(0,sigma_eta^2)

where we assume the correlation betweem eps(t) and eta(t) equals to rho.
Prior distributions are 

mu~N(mu_0,sigma_0^2), (phi+1)/2~Beta(a_0,b_0), (rho+1)/2~Beta(a_1,b_1), sigma_eta^2~IG(n_0/2,S_0/2)

where N, Beta, IG denote normal, beta, inverse gaussian distributions respectively. Note that the probability density function of x ~ IG(a,b) is proportional to (1/x)^(a+1)*exp(-b/x). 

The highly efficient Markov chain Monte Carlo algorithm is based on the mixture 
sampler by Omori, Chib, Shephard and Nakajima (2007), but it further corrects 
the approximation error within the sampling algorithm.
See Takahashi, Watanabe and Omori (2022) for more details.
}
\usage{
asv_mcmc(return_vector, nSim = NULL, nBurn = NULL, vhyper = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{return_vector}{
%%     ~~Describe \code{return_vector} here~~
      T x 1 vector (y(1),...,y(T))' of returns where T is a sample size.
}
  \item{nSim}{
%%     ~~Describe \code{nSim} here~~
Number of iterations for the MCMC estimation. Default value is 5000.
}
  \item{nBurn}{
%%     ~~Describe \code{nBurn} here~~
Number of iterations for the burn-in period. Default value is the maximum integer less than or equal to 2*sqrt(nSim)+1.
}
  \item{vhyper}{
8 x 1 vector of hypermarameters. (mu_0,sigma_0^2,a_0,b_0,a_1,b_1,n_0,S_0). Default values are (0,1000, 1,1,1,1,0.01,0.01).
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A list of a list with components:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{mu_s}{nSim x 1 vector of MCMC samples of mu}
  \item{phi_s}{nSim x 1 vector of MCMC samples of phi}
  \item{sigma_eta_s}{nSim x 1 vector of MCMC samples of sigma_eta}
  \item{rho_s}{nSim x 1 vector of MCMC samples of rho}
  \item{mh_s}{nSim x T matrix of latent log volatilities (h(1),...,h(T)). For example, the first column is a vector of MCMC samples for h(1).}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Omori, Y., Chib, S., Shephard, N., and J. Nakajima (2007), "Stochastic volatility model with leverage: fast and efficient likelihood inference," Journal of Econometrics, 140-2, 425-449.

Takahashi, M., Omori, Y. and T. Watanabe (2022+), Realized volatility and stochastic volatility models. JSS Research Series in Statistics, in press.
}
\author{
Yasuhiro Omori
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
See also \code{\link{ReportMCMC}}, \code{\link{asv_pf}}}
}
\examples{
nobs = 2000;
mu = 0; phi = 0.97; sigma_eta = 0.3; rho = -0.3;
h  = 0;   Y = c();
for(i in 1:nobs){
  eps = rnorm(1, 0, 1)
  eta = rho*sigma_eta*eps + sigma_eta*sqrt(1-rho^2)*rnorm(1, 0, 1)
  y   = eps * exp(0.5*h)
  h   = mu + phi * (h-mu) + eta
  Y   = append(Y, y)
}
nsim = 2000; nburn = 1000;
vhyper = c(0.0,1000,1.0,1.0,1.0,1.0,0.01,0.01)
out  = sv_mcmc(Y, nsim, nburn, vhyper)
vmu = out[[1]]; vphi = out[[2]]; vsigma_eta = out[[3]]; vrho = out[[4]];
mh  = out[[5]];
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
