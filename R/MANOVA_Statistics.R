## function for calculating test statistics for MANOVA
MANOVA.Stat<- function(data, n, hypo_matrix, iter, alpha, resampling, n.groups, p, CPU){
  
  N <- sum(n)
  H <- hypo_matrix
  x <- data
  
  
  #---------------- useful matrices ---------------------#
  A <-  t(rep(1 / n[1], n[1])) %x% diag(p)
  for (ii in 2:length(n)){
    B <- t(rep(1 / n[ii], n[ii])) %x% diag(p)
    A <- magic::adiag(A, B)
  } 
  # -----------------------------------------------------#
  means <- A %*% x
  
  V <- list(NA)
  n.temp <- cumsum(c(0, n))
  for (i in 1:n.groups){
    y <- matrix(x[(n.temp[i]*p+1):(n.temp[i+1]*p)], ncol = p, byrow = TRUE)
    V[[i]] <- 1 / n[i] * cov(y)
  }
  
  sigma_hat <- V[[1]]
  for (i in 2:n.groups){
    sigma_hat <- magic::adiag(sigma_hat, V[[i]])
  }
  Sn <- N * sigma_hat
  
  # WTS
  T <- t(H) %*% MASS::ginv(H %*% Sn %*% t(H)) %*% H
  WTS <- N * t(means) %*% T %*% means
  df_WTS <- Matrix::rankMatrix(H)[[1]]
  
  # ATS
  C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
  D <- diag(C) * diag(ncol(C))
  spur <- sum(diag(C %*% Sn))
  # Lambda <- diag(1 / (n - 1))
  ATS <- N / spur * t(means) %*% C %*% means
  df_ATS <- spur ^ 2 / sum(diag(C %*% Sn %*% C %*% Sn))
  # df_ATS2 <- spur ^ 2 / sum(diag(D %*% D %*% Sn %*% Sn %*% Lambda))
  df_ATS2 <- Inf
  
  #--------------------------------- parametric bootstrap ---------------------------#
  PBS <- function(i, ...){
    # calculate mvrnorm for each group
    XP <- list()
    meansP <- list()
    for (i in 1:n.groups){
      XP[[i]] <- MASS::mvrnorm(n[i], mu = rep(0, p), Sigma = n[i]*V[[i]])
      meansP[[i]] <- colMeans(XP[[i]])
    }
    meansP <- unlist(meansP)
    
    VP <- list()
    for(i in 1:n.groups){
      VP[[i]] <- 1 / n[i] * cov(XP[[i]])
    }
    
    sigma_hatP <- VP[[1]]
    for (i in 2:n.groups){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hatP
    
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    return(WTPS)
  }
  
  #---------------------------------- Wild bootstrap ---------------------------------#
  WBS <- function(i, ...){
    
    VP <- list(NA)
    xperm <- list(NA)
    for (i in 1:n.groups){
      y <- matrix(x[(n.temp[i]*p+1):(n.temp[i+1]*p)], ncol = p, byrow = TRUE)
      means2 <- rep(colMeans(y), n[i])
      epsi <- 2*rbinom(n[i], 1, 1/2)-1
      xperm[[i]] <- rep(epsi, p)*(y -means2)
      VP[[i]] <- 1 / n[i] * cov(xperm[[i]])
    }
    
    sigma_hatP <- VP[[1]]
    for (i in 2:n.groups){
      sigma_hatP <- magic::adiag(sigma_hatP, VP[[i]])
    }
    SnP <- N * sigma_hat
    
    xperm <- unlist(xperm)
    meansP <- A %*% xperm
    
    # WTS
    TP <- t(H) %*% MASS::ginv(H %*% SnP %*% t(H)) %*% H
    WTPS <- diag(N * t(meansP) %*% TP %*% meansP)
    # ATS
    C <- t(H) %*% MASS::ginv(H %*% t(H)) %*% H
    spur <- sum(diag(C %*% SnP))
    ATS_res <- N / spur * t(meansP) %*% C %*% meansP
    return(list(WTPS, ATS_res))
  }
  
  time1 <- Sys.time()  
  cl <- makeCluster(CPU)
  
  if(resampling == "paramBS"){
    WTPS <- parSapply(cl, 1:iter, FUN = PBS)
    ecdf_WTPS <- ecdf(WTPS)
    p_valueWTPS <- 1-ecdf_WTPS(WTS)    
    p_valueATS_res <- NA       
  } else if(resampling == "WildBS"){
    WTPS <- parSapply(cl, 1:iter, FUN = WBS)
    ecdf_WTPS <- ecdf(unlist(WTPS[1, ]))
    p_valueWTPS <- 1-ecdf_WTPS(WTS)    
    ecdf_ATS_res <- ecdf(unlist(WTPS[2, ]))
    p_valueATS_res <- 1-ecdf_ATS_res(ATS)    
  }
  
  time <- Sys.time() - time1
  parallel::stopCluster(cl)
  
  #------------------------ p-values -------------------------------#
  p_valueWTS <- 1 - pchisq(abs(WTS), df = df_WTS)
  p_valueATS <- 1 - pf(abs(ATS), df1 = df_ATS, df2 = df_ATS2)
  
  #---------------------- CIs -------------------------------------#
  CI_lower <- means - sqrt(diag(Sn) / n) * qt(1 - alpha / 2, df = n)
  CI_upper <- means + sqrt(diag(Sn) / n) * qt(1 - alpha / 2, df = n)
  
  #-------------------- Output ----------------------------------#
  WTS_out <- c(WTS, df_WTS, p_valueWTS)
  ATS_out <- c(ATS, df_ATS, df_ATS2, p_valueATS)
  WTPS_out <- c(p_valueWTPS, p_valueATS_res)
  CI <- cbind(CI_lower, CI_upper)
  result <- list(WTS = WTS_out, WTPS = WTPS_out, ATS = ATS_out,
                 Cov = Sn, Mean = means, CI = CI, time = time)
  return(result)
}