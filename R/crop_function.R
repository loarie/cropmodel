#'
#' This function takes responses (y), and 2 environmental variables (x1 & x2) along with prior parameters for beta (b0, Vbcoef) and sigma2 (s1,s2) and produces estimates for beta and sigma using a Gibb's sampler and quadradic Bayesian regression
#' @param b0, Vbcoef, s1, s2, y, x1, x2
#' @keywords crop Bayesian model
#' @export
#' @examples
#' crop_function()

crop_function <- function(b0, Vbcoef, s1, s2, y1, x1, x2){
  library(bayesm)
  library(coda)

  n = length( y1 )

  if(n==1){ #simulate another data point from the model
    y1 = c(y1, sum(b0*c(1, x1, x2, x1^2, x2^2, 115)))
    x1 = c(x1,x1)
    x2 = c(x2,x2)
    n = 2
  }
  y1 = y1+rnorm(n,0,0.0001)
  
  X1 = cbind( rep( 1, times = n ), x1, x2, x1^2, x2^2, 115 )
  p = dim( X1 )[2]

  Data1 = list(y=y1,X=X1)

  # I have manually set the prior information to the function default:
  betabar1 = b0 # prior mean
  A1 = Vbcoef*diag(p) # prior precision matrix (precision = 1/variance)

  nu1 = s1
  ssq1 = s2

  Prior1 = list(betabar=betabar1,A=A1,nu=nu1,ssq=ssq1)

  # This portion specifies MCMC parameters
  R1 = 100000 # This is how many iterations to run the Gibbs sampler 
  keep1 = 1 # This is the thinning parameter. 1 means no thinning.

  MCMC1 = list(R=R1,keep=keep1)

  ### Run the Gibbs sampler! ###
  simOut = runiregGibbs(Data1,Prior1,MCMC1)

  b_mean = apply( simOut$betadraw, 2, mean )
  sigma2_mean = mean(simOut$sigmasqdraw)
  print(paste("{beta:[",paste(b_mean,collapse=","),"], sigma2:",sigma2_mean,"}",sep="")[1])
  #print(sum(b0*c(1, x1[1], x2[1], x1[1]^2, x2[1]^2, 115)))
  #print(sum(b_mean*c(1, x1[1], x2[1], x1[1]^2, x2[1]^2, 115)))
}
