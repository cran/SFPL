sfpl <- function(x, y, ls_vec, lf_vec, epsilon, verbose){

  # some auxilliary functions
  make_vf <- function(B, p, K, epsilon){
    vf <- matrix(0, p*K, p*K)
    valseq <- seq(0, p*K, p)
    for(i in 1:(p*K)){  #rows
      for(j in 1:(p*K)){ # columns
        if(i > p){
          k_val <- min(which(valseq >= i)) - 1
          q_val <- i - (k_val-1)*p
        } else{
          k_val <- 1
          q_val <- i
        }

        if(j > p){
          kp_val <- min(which(valseq >= j)) - 1
          qp_val <- j - (kp_val-1)*p
        } else{
          kp_val <- 1
          qp_val <- j
        }

        if(q_val == qp_val && k_val == kp_val){
          if(k_val < K){
            b_qk <- B[i]
            b_other <- c(1:K)[k_val:(K-1)]
            frac_b = numeric(length(b_other))
            for(l in 1:length(b_other)){
              b_other[l] <- B[(b_other[l]*p) + q_val]
              frac_b[l] <- 1/(abs(b_qk - b_other[l]) + epsilon)
            }
            vf[i,j] <- sum(frac_b)
          }
        }

        if(q_val == qp_val && k_val < kp_val){
          b_qk <- B[i]
          b_qpkp <- B[j]
          vf[i,j] <- -1/(abs(b_qk - b_qpkp) + epsilon)
        }
        #cat(i, j, vf[i,j], "\n")
      }
    }
    return(vf)
  }

  make_vs <- function(B, p, K, epsilon){
    vs <- matrix(0, p*K, p*K)
    valseq <- seq(0, p*K, p)
    for(i in 1:(p*K)){  #rows
      for(j in 1:(p*K)){ # columns
        if(i > p){
          k_val <- min(which(valseq >= i)) - 1
          q_val <- i - (k_val-1)*p
        } else{
          k_val <- 1
          q_val <- i
        }

        if(j > p){
          kp_val <- min(which(valseq >= j)) - 1
          qp_val <- j - (kp_val-1)*p
        } else{
          kp_val <- 1
          qp_val <- j
        }

        if(q_val == qp_val && k_val == kp_val){
          vs[i,j] <- 1/(abs(B[i]) + epsilon)
        }
      }
    }
    return(vs)
  }

  pl_loglik <- function(B, x_data, y){
    # B is a long K*p vector, because the optimization function requires this
    loglik = 0
    K <- length(y)
    n_k <- numeric(K)
    for(k in 1:K){
      n_k[k] = nrow(y[[k]])   # number of rankers for group k
      # isolate beta vectors for group k from large multigroup B vector
      if(k == 1){
        beta <- B[1:ncol(x_data)]
      } else{
        a <- (k-1)*ncol(x_data)+1
        b <- (k)*ncol(x_data)
        beta <- B[a:b]
      }
      for(i in 1:n_k[k]){
        ranking = y[[k]][i, ]
        xb = as.matrix(x_data)%*%beta
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
    return(loglik)
  }

  eval_func <- function(B, x, y, ls, lf){
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
      loglik = loglik + ls*sum(abs(beta))
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
        loglik = loglik + lf*sum(abs(beta_k - beta_kp))
      }
    }
    return(loglik)
  }

  sfpl_inner <- function(x, y, ls, lf, epsilon, B_init){

    K <- length(y)
    p <- ncol(x)

    beta_est <- B_init

    diff <- 1000
    xi <- 0.001
    while(diff > xi){
      vs_mat <- make_vs(beta_est, p, K, epsilon)
      vf_mat <- make_vf(beta_est, p, K, epsilon)

      beta_est_new <- beta_est - solve(pracma::hessian(pl_loglik, beta_est, x_data = x, y = y) + ls*vs_mat + lf*vf_mat)%*%(pracma::grad(pl_loglik, beta_est, x_data = x, y = y) + (ls*vs_mat + lf*vf_mat)%*%beta_est)

      diff <- abs((eval_func(beta_est_new, x, y, ls, lf) - eval_func(beta_est, x, y, ls, lf))/eval_func(beta_est, x, y, ls, lf))
      beta_est <- beta_est_new
    }
    beta_est <- matrix(beta_est, nrow = p)
    return(beta_est)
  }

  lambdas <- expand.grid(ls_vec, lf_vec)

  K <- length(y)
  p <- ncol(x)

  # compute initial estimate B^[0] based on the normal likelihood
  B_init <- matrix(rep(0, K*p), ncol = 1)
  if(is.null(tryCatch({
    result <- optim(par = B_init, pl_loglik, method = "BFGS", x_data = x, y = y)
  }, error=function(e){}))){
    result <- nlminb(start = B_init, pl_loglik, x_data = x, y = y) # use instead of BFGS when numerical issues arise
    mthd <- "nlminb"
  } else{
    result <- optim(par = B_init, pl_loglik, method = "BFGS", x_data = x, y = y)
    mthd <- "BFGS"
  }

  B_init <- result$par

  beta_est <- vector("list", nrow(lambdas))
  for(l_val in 1:nrow(lambdas)){

    if(verbose){
      cat("Parameter estimation in process:", (l_val/nrow(lambdas))*100, "%", "\n")
    }
    ls <- lambdas[l_val,1]
    lf <- lambdas[l_val,2]

    result <- sfpl_inner(x, y, ls, lf, epsilon, B_init)

    beta_est[[l_val]] <- result
  }
  return(beta_est)
}
