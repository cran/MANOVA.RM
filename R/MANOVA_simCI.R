#' Simultaneous confidence intervals for contrasts in multivariate factorial designs.
#' 
#' @param object A \code{MANOVA} object.
#' @param contrast The contrast matrix of interest, can either be "pairwise" or "user-defined". 
#' @param contmat If contrast = "user-defined", the contrast matrix must be specified here. Note that
#' its rows must sum to zero.
#' @param type If contrast is "pairwise", the type of the pairwise comparison must be specified here. 
#' Calculation is based on the contrMat function in package multcomp, see the corresponding help page 
#' for details on the types of contrasts available.
#' @param base an interger specifying which group is considered the baseline group 
#' for Dunnett contrasts, see \code{\link[multcomp]{contrMat}}.
#' @param ... Not used yet.
#'   
#' @details The simCI() function computes confidence intervals for the chosen contrasts of the multivariate mean vector 
#' based on the sum statistic. Details on the derivation of these confidence intervals can be found in Friedrich and Pauly (2018).
#'    
#' @return Simultaneous confidence intervals for the chosen contrasts. 
#' 
#' @references Friedrich, S., and Pauly, M. (2018). MATS: Inference for potentially singular and
#'   heteroscedastic MANOVA. Journal of Multivariate Analysis, 165, 166-179.
#'   
#' @seealso \code{\link[multcomp]{contrMat}}
#' 
#' @importFrom multcomp contrMat
#'
#' @export
simCI <- function(object, contrast = c("pairwise", "user-defined"), contmat = NULL, type = NULL,
                  base = 1, ...){
  
  if(object$nested){
    stop("The pairwise comparisons cannot be used in nested designs!")
  }
  
  meanvec <- as.vector(t(object$Means))
  n <- object$Descriptive$n
  N <- sum(n)
  p <- object$p
  D <- diag(object$Covariance)*diag(p*length(n))
  factors <- object$factors
  nf <- object$nf
  BSmean <- object$BSMeans
  BSD <- object$BSVar
  alpha <- object$input$alpha
  fl <- object$fl
  lev <- subset(object$Descriptive, select = 1:n)
  lev <- lev[-ncol(lev)]
  
  if(contrast == "user-defined"){
    if(is.null(contmat)){
      stop("Please specify a contrast matrix.")
    }
    if(sum(rowSums(contmat) > 1e-15) != 0){
      stop("The rows of the contrast matrix must sum to zero!")
    }
    if (ncol(contmat) != p*prod(fl)){
      stop(paste("The contrast matrix must have", prod(fl)*p, "columns."))
    }
  }
  
  if(contrast == "pairwise"){
    if(is.null(type)){
      stop("Please specify the type of pairwise comparison, see the multcomp-package for
           details.")
    }
   # one-way
    if(nf == 1){
      names(n) <- lev[, 1]
    } else {
      names(n) <- sort(do.call(paste, c(lev, sep = " ")))
    }
      M <- contrMat(n, type = type, base)
      contmat <- M %x% t(rep(1, p))
}
  
  # calculation of critical value based on sum statistic
  sumstat <- function(mean, var, ...){
    
    HDH <- diag(1/(contmat%*%var%*%t(contmat)))*diag(nrow(contmat))
    S_N_star <- N*t(contmat %*% mean)%*% HDH %*% contmat %*% mean 
    
    return(S_N_star)
  }
  
  S_star <- mapply(sumstat, BSmean, BSD)
  ecdf_S_star <- ecdf(S_star)
  q_star <- quantile(ecdf_S_star, 1-alpha)
  
  ## confidence interval
  sci <- function(con, ...){
    center <- con%*%meanvec
    CI <- c(center - sqrt(q_star* t(con)%*%D%*%con/N), center + sqrt(q_star* t(con)%*%D%*%con/N))
    # test statistic for this contrast:
    Q_N_l <- N*t(center)%*%(t(con)%*%D%*%con)^(-1)%*%center
    p_val <- 1-ecdf_S_star(Q_N_l)
    
    result <- c(center, CI, p_val)
    names(result) <- c("Estimate", "Lower", "Upper", "p-value")
    return(result)
  }
  
  scis <- t(apply(contmat, 1, sci))
  if (contrast == "pairwise"){
    rownames(scis) <- rownames(M)
  }
  
  alpha <- 100*(1-object$input$alpha)
  cat(alpha, "% confidence intervals", "\n", sep = "")
  print(scis)
}
