#Script that loads a variety of sigmoid functions for creating virtual species curves
library(stats)

#Modified sigmoid function
dsigmoid <- function(x, curve, halfway){
  return(1 / (1 + (exp((-1 * curve) * (x - halfway)))))
}

#Re-scaled beta distribution for producing curves
#See https://stats.stackexchange.com/questions/638026/is-there-a-mathematical-function-with-a-shape-like-an-inverted-sigmoid-pdf-but-x?noredirect=1#comment1192438_638026 for the stackoverflow post that led to these functions


#Produces a decreasing curve
inv_d_sigmoid <- function(x, max_val, curve){
  
  resp_prob <- (1 - qbeta((x / max_val), curve, curve))
  
  resp_prob[is.nan(resp_prob)] <- 0
  
  return(resp_prob)
}

#Produces an increasing curve
inv_i_sigmoid <- function(x, max_val, curve){
  
  resp_prob <- (qbeta((x / max_val), curve, curve))
  
  resp_prob[is.nan(resp_prob)] <- 1
  
  return(resp_prob)
}
