
sfpl_approx <- function(x, y, ls_vec, lf_vec, epsilon, verbose){

  pl_penloglik <- function(B, x, y, ls, lf, epsilon){
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
    loglik <- (1/sum(n_k))*loglik

    for(k in 1:K){
      if(k == 1){
        beta <- B[1:ncol(x)]
      } else{
        a <- (k-1)*ncol(x)+1
        b <- (k)*ncol(x)
        beta <- B[a:b]
      }
    loglik = loglik + ls*sum(abs(beta) - epsilon*log(1 + (abs(beta)/epsilon)))
    }

    for(kp in 1:(K-1)){
      for(k in (kp+1):K){
        if(k == 1){
          beta_kp <- B[1:ncol(x)]
        } else{
          a <- (kp-1)*ncol(x)+1
          b <- (kp)*ncol(x)
          beta_kp <- B[a:b]
        }
        a <- (k-1)*ncol(x)+1
        b <- (k)*ncol(x)
        beta_k <- B[a:b]
        loglik = loglik + lf*sum(abs(beta_k - beta_kp) - epsilon*log(1 + (abs(beta_k - beta_kp)/epsilon)))
      }
    }
    return(loglik)
  }

  lambdas <- expand.grid(ls_vec, lf_vec)
  K <- length(y)
  p <- ncol(x)

  B <- matrix(rep(0, K*p), ncol = 1)
  beta_est <- vector("list", nrow(lambdas))
  for(l_val in 1:nrow(lambdas)){
    if(verbose){
      cat("Parameter estimation in process:", (l_val/nrow(lambdas))*100, "%", "\n")
    }
    ls <- lambdas[l_val,1]
    lf <- lambdas[l_val,2]

    if(is.null(tryCatch({
      result <- optim(par = B, pl_penloglik, method = "BFGS", x = x, y = y, ls = ls, lf = lf, epsilon = epsilon)
    }, error=function(e){}))){
      result <- nlminb(start = B, pl_penloglik, x = x, y = y, ls = ls, lf = lf, epsilon = epsilon) # use instead of BFGS when numerical issues arise
      mthd <- "nlminb"
    } else{
      result <- optim(par = B, pl_penloglik, method = "BFGS", x = x, y = y, ls = ls, lf = lf, epsilon = epsilon)
      mthd <- "BFGS"
    }
    beta_est[[l_val]] <- matrix(result$par, nrow = ncol(x), byrow = F)
  }
  return(beta_est)
}
