\name{asv_apf}
\alias{asv_apf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Auxiliary particle filter to compute the log likelihood for stochastic volatility 
models with leverage
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
asv_apf(mu, phi, sigma_eta, rho, Y, I)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{parameter value such as the posterior mean of mu
}
  \item{phi}{parameter value such as the posterior mean of phi
}
  \item{sigma_eta}{parameter value such as the posterior mean of sigma_eta
}
  \item{rho}{parameter value such as the posterior mean of rho
}
  \item{Y}{T x 1 vector (y(1),...,y(T))' of returns where T is a sample size.
}
  \item{I}{Number of particles to approximate the filtering density.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Logarithm of likelihood of Y given parameters (mu, phi, sigma_eta, rho)
using the auxiliary particle filter by Pitt and Shephard (1999).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Pitt, M. K., and N. Shephard (1999), "Filtering via simulation: Auxiliary particle filters." Journal of the American statistical association 94, 590-599.

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
}
\examples{
mu = 0; phi = 0.97; sigma_eta = 0.3; rho = -0.3;
npart = 5000 
# Repeat m times to compute the standard error of the estimate
m  = 10; loglik = rep(0, m)
for(i in 1:m){
  loglik[i] = asv_apf(mu, phi, sigma_eta, rho, Y, npart)
}
# Estimate and its standard error
loglik
c(mean(loglik), sqrt(var(loglik)/m))
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
