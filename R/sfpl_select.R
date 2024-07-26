
sfpl_select <- function(beta_est, x, y, ls_vec, lf_vec){

  loglik <- function(B, x, y){
    loglik = 0
    K <- length(y)
    n_k <- numeric(K)
    for(k in 1:K){
      n_k[k] = nrow(y[[k]])   # number of rankers for group k
      # isolate beta vectors for group k from large multigroup B vector
      if(k == 1){
        beta <- B[1:ncol(x)]
      } else{
        a <- (k-1)*ncol(x)+1
        b <- (k)*ncol(x)
        beta <- B[a:b]
      }
      for(i in 1:n_k[k]){
        ranking = y[[k]][i, ]
        xb = as.matrix(x)%*%beta
        lambda = exp(xb)
        lsum = sum(exp(xb[as.numeric(ranking)]))
        # plackettluce log likelihood, where objects are evaluated in order of rank
        for(j in 1:length(ranking)){
          loglik = loglik + log(lsum) - xb[as.numeric(ranking[j])]
          lsum = lsum - lambda[as.numeric(ranking[j])] # remove already ranked item from denominator in lik
        }
      }
    }
    return(loglik)
  }

  lambdas <- expand.grid(ls_vec, lf_vec)


  lik <- numeric(nrow(lambdas))
  for(l_val in 1:nrow(lambdas)){
    B <- matrix(beta_est[[l_val]], ncol = 1)
    lik[l_val] <- loglik(B, x, y)
  }

  aic <- bic <- numeric(nrow(lambdas))
  for(l_val in 1:nrow(lambdas)){

    nonzero <- length(which(round(beta_est[[l_val]], 3) != 0))
    n_k <- numeric(length(y))
    for(k in 1:length(y)){
      n_k[k] <- nrow(y[[k]])
    }

    aic[l_val] <- 2*nonzero + 2*lik[l_val]
    bic[l_val] <- nonzero*log(sum(n_k)) + 2*lik[l_val]
  }
  aic_sel <- round(beta_est[[which.min(aic)]],3)
  bic_sel <- round(beta_est[[which.min(bic)]],3)

  return(list(model_aic = aic_sel, model_bic = bic_sel))
}


