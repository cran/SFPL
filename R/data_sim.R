data_sim <- function(m, M, n, p, K, delta, eta){
  
  comp_prob <- function(beta_true, x, ranking){
    xb = as.matrix(x)%*%beta_true
    lambda = exp(xb)
    lsum = sum(exp(xb[ranking]))
    prob <- 1
    for(j in 1:length(ranking)){
      prob = prob * (lambda[as.numeric(ranking[j])] / lsum)
      lsum = lsum - lambda[as.numeric(ranking[j])] 
    }
    return(prob)
  }
  
  probrank <- function(n, m, M, partial, beta_true, x){
    if(partial == TRUE){
      ranking <- gtools::permutations(M, m, 1:M)
    } else{
      ranking <- gtools::permutations(M, M, 1:M)
    }
    
    probs <- numeric(nrow(ranking))
    for(perm in 1:nrow(ranking)){
      probs[perm] <- comp_prob(beta_true, x, ranking[perm,])
    }
    
    s_idx <- sample(1:nrow(ranking), n, replace = TRUE, prob = probs)
    
    samples <- ranking[s_idx,]
    
    return(samples)
  }
  
  if(m == M){
    partial <- FALSE
  } else{
    partial <- TRUE
  }
  
  beta_true <- matrix(0, nrow = p, ncol = K)
  beta_true[,1] <- runif(p, -1, 1)
  n_sparse <- floor(eta*p)
  sparse_idx <- sample(1:p, n_sparse)
  beta_true[sparse_idx,1] <- 0
  for(k in 2:K){
    n_dif <- floor(delta*p) 
    dif_idx <- sample(1:p, n_dif)
    beta_true[,k] <- beta_true[,1]
    beta_true[dif_idx,k] <- runif(n_dif, -1, 1)
  }
  x <- rnorm(p*M)
  x <- scale(x)
  x <- as.numeric(x)
  x <- matrix(x, nrow = M, byrow = F)
  
  rank_true <- matrix(0, nrow = K, ncol = M)
  for(k in 1:K){
    rank_true[k,] <- order(exp(x%*%beta_true[,k]), decreasing = T) 
  }
  
  rankinglist <- vector("list", K)
  for(k in 1:K){
    if(partial == TRUE){
      rankinglist[[k]] <- probrank(n, m, M, partial = TRUE, beta_true[,k], x)
    } else{
      rankinglist[[k]] <- probrank(n, m, M, partial = FALSE, beta_true[,k], x)
    }
  }
  y <- rankinglist
  return(list(y = y, x = x, beta = beta_true))  
}