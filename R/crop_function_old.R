#' Crop Function_old 
#'
#' This function takes responses (y), and 2 environmental variables (x1 & x2) along with prior parameters for beta (b0, Vbcoef) and sigma2 (s1,s2) and produces estimates for beta and sigma using a Gibb's sampler and quadradic Bayesian regression
#' @param b0, Vbcoef, s1, s2, y, x1, x2
#' @keywords crop Bayesian model
#' @export
#' @examples
#' crop_function_old ()

crop_function_old <- function(b0, Vbcoef, s1, s2, y, x1, x2){
  library( mvtnorm )
  library( MCMCpack )
  
  n = length( y )
  X = cbind( rep( 1, times = n ), x1, x2, x1^2, x2^2, 2015 )
  p = dim( X )[2]

  #beta prior
  Vb = Vbcoef * diag( p )
  b_val = rmvnorm( 1, b0, Vb )[1,]

  #sigma2 prior
  sigma2_val = rinvgamma( 1, s1, s2 )

  ngibbs = 20000
  b_vals = matrix( NA, ngibbs, p )
  sigma2_vals = rep( NA, times = ngibbs )
  for( g in 1 : ngibbs ){
    #sample sigma2
    u1 = s1 + n / 2
    u2 = s2 + 0.5 * crossprod( y - X %*% b_val )
    sigma2_val = rinvgamma( 1, u1, u2 )
    sigma2_vals[g] = sigma2_val 
    
    #sample b
    v = 1 / ( sigma2_val ) * crossprod( X, y ) + solve( Vb ) %*% b_val
    V_inv = 1/(sigma2_val) * crossprod(X) + solve( Vb )
    b_val = rmvnorm( 1, solve( V_inv ) %*% v, solve( V_inv ) )[1,]
    b_vals[g,] = b_val
  }

  b_mean = apply( b_vals, 2, mean )
  sigma2_mean = mean( sigma2_vals )
  #print(paste("{\"beta\":[",paste(b_mean,collapse=","),"], \"sigma2\":",sigma2_mean,"}",sep="")[1])
  print(paste("{beta:[",paste(b_mean,collapse=","),"], sigma2:",sigma2_mean,"}",sep="")[1])
    
}