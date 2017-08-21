#' Tests for Multivariate Data in Semi-Parametric Factorial Designs
#' 
#' The MANOVA function calculates the Wald-type statistic (WTS), the ANOVA-type 
#' statistic (ATS) and a modified ATS (MATS) as well as resampling versions of 
#' these test statistics for 
#' semi-parametric multivariate data.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side 
#'   contains the response variable and the right hand side contains the factor 
#'   variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. Data must be in long format.
#' @param subject The column name of the subjects in the data.
#' @param iter The number of iterations used for calculating the resampled 
#'   statistic. The default option is 10,000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "paramBS"
#'   (parametric bootstrap approach) and "WildBS" (wild bootstrap approach with
#'   Rademacher weights). The Wild Bootstrap is calculated for all test statistics.
#' @param CPU The number of cores used for parallel computing. If omitted, cores are
#'   detected via \code{\link[parallel]{detectCores}}.
#' @param seed A random seed for the resampling procedure. If omitted, no 
#'   reproducible seed is set.
#' @param nested.levels.unique A logical specifying whether the levels of the nested factor(s)
#'   are labeled uniquely or not. Default is FALSE, i.e., the levels of the nested 
#'   factor are the same for each level of the main factor. For an example and more explanations
#'   see the GFD package and the corresponding vignette.
#'  
#'   
#' @details The MANOVA() function provides the Wald-type statistic (WTS) as well as
#'   the ANOVA-type statistic (ATS) for multivariate designs with metric data as described in 
#'   Konietschke et al. (2015). Furthermore, it contains a modified ATS (MATS), which is invariant
#'   under scale transformations of the components and applicable to designs with singular covariance
#'   matrices, see Friedrich and Pauly (2017) for details.
#'   These tests are even applicable for non-normal error terms, 
#'   different sample sizes and/or heteroscedastic variances.
#'   They are implemented for designs with an arbitrary number of
#'   crossed factors or for nested designs. In addition to the asymptotic
#'   p-values, the function also provides p-values based on resampling approaches.
#'  
#' @section NOTE: The number of resampling iterations has been set to 100 in the examples due to run time 
#' restrictions on CRAN. Usually it is recommended to use at least 1000 iterations. 
#' For more information and detailed examples also refer to the package vignette.
#'     
#' @return A \code{MANOVA} object containing the following components: 
#'   \item{Descriptive}{Some descriptive statistics of the data for all factor 
#'   level combinations. Displayed are the number of individuals per factor 
#'   level combination and the vector of means (one column per dimension).}
#'   \item{Covariance}{The estimated covariance matrix.} 
#'   \item{WTS}{The value of the WTS along with degrees of freedom of the
#'   central chi-square distribution and p-value.} 
#'   \item{ATS}{The value of the
#'   ATS, degrees of freedom of the central F distribution and the corresponding
#'   p-value.}
#'   \item{MATS}{The value of the MATS.} 
#'   \item{resampling}{p-values for the test statistic based on the
#'   chosen resampling approach.}
#'  
#'  
#' @examples data(EEG)
#' EEG_mod <- MANOVA(resp ~ sex * diagnosis, 
#'                     data = EEG, subject = "id", resampling = "paramBS", 
#'                     alpha = 0.05, iter = 100, CPU = 1)
#' summary(EEG_mod)
#' 
#' @seealso \code{\link{RM}}
#'   
#' @references Konietschke, F., Bathke, A. C., Harrar, S. W. and Pauly, M. (2015). 
#'   Parametric and nonparametric bootstrap methods for general MANOVA. Journal 
#'   of Multivariate Analysis, 140, 291-301.
#'   
#'   Friedrich, S., Brunner, E. and Pauly, M. (2017). Permuting longitudinal data
#'   in spite of the dependencies. Journal of Multivariate Analysis, 153, 255-265.
#'   
#'    Bathke, A., Friedrich, S., Konietschke, F., Pauly, M., Staffen, W., Strobl, N. and Hoeller, Y. (2016).
#'    Using EEG, SPECT, and Multivariate Resampling Methods
#'    to Differentiate Between Alzheimer's and other Cognitive Impairments. arXiv preprint arXiv:1606.09004.
#'   
#'   Friedrich, S., Konietschke, F., Pauly, M. (2016). GFD - An 
#'   R-package for the Analysis of General Factorial Designs. Accepted for publication in 
#'   Journal of Statistical Software.
#'   
#'   Friedrich, S., and Pauly, M. (2017). MATS: Inference for potentially singular and
#'   heteroscedastic MANOVA. arXiv preprint arXiv:1704.03731.
#'    
#'   
#' @importFrom graphics axis legend par plot title abline points
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var cov rbinom quantile
#' @importFrom utils read.table
#' @importFrom methods hasArg
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster parSapply detectCores
#' @importFrom ellipse ellipse
#'   
#' @export

MANOVA <- function(formula, data, subject,
                   iter = 10000, alpha = 0.05, resampling = "paramBS", CPU,
                   seed, nested.levels.unique = FALSE){
  
  if (!(resampling %in% c("paramBS", "WildBS"))){
    stop("Resampling must be one of 'paramBS' and 'WildBS'!")
  }
  
  input_list <- list(formula = formula, data = data,
                     subject = subject, 
                     iter = iter, alpha = alpha, resampling = resampling)
  
  test1 <- hasArg(CPU)
  if(!test1){
    CPU <- parallel::detectCores()
  }
  
  test2 <- hasArg(seed)
  if(!test2){
    seed <- 0
  }
  
  dat <- model.frame(formula, data)
  if (!(subject %in% names(data))){
    stop("The subject variable is not found!")
  }
  subject <- data[, subject]
  # no. of dimensions 
  p <- length(subject)/length(unique(subject))
  dat2 <- data.frame(dat, subject = subject)
  nf <- ncol(dat) - 1
  nadat <- names(dat)
  nadat2 <- nadat[-1]
  fl <- NA
  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, (aa + 1)]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, (jj + 1)]))
  }
  lev_names <- expand.grid(levels)
  if (nf == 1) {
    # one-way layout
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    dat2 <- dat2[order(dat2[, 2], dat2$subject), ]
    response <- dat2[, 1]    
    # contrast matrix
    hypo <- (diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)) %x% diag(p)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(unique(subject)),
                     .drop = F)$Measure
    WTS_out <- matrix(NA, ncol = 3, nrow = 1)
    ATS_out <- matrix(NA, ncol = 4, nrow = 1)
    MATS_out <- NA
    WTPS_out <- rep(NA, 3)
    quantiles <- matrix(NA, 2, 1)
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    names(WTPS_out) <- fac_names
    results <- MANOVA.Stat(data = response, n = n, hypo, iter = iter, alpha, resampling, n.groups = fl, p, CPU, seed, nf)    
    WTS_out <- results$WTS
    ATS_out <- results$ATS
    MATS_out <- results$MATS
    WTPS_out <- results$WTPS
    quantiles <- results$quantiles
    names(quantiles) <- c("WTS_resampling", "MATS_resampling")
    mean_out <- matrix(results$Mean, ncol = p, byrow = TRUE)
    Var_out <- results$Cov
    #    CI <- results$CI
    #    colnames(CI) <- c("lower", "upper")
    descriptive <- cbind(lev_names, n, mean_out)
    colnames(descriptive) <- c(nadat2, "n", rep("Means", p))   
    names(WTS_out) <- cbind ("Test statistic", "df",
                             "p-value")
    names(ATS_out) <- cbind("Test statistic", "df1", "df2", "p-value")
    names(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(ATS)"), paste(resampling, "(MATS)"))
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
    #    output$CI <- CI
    output$Covariance <- Var_out
    output$WTS <- WTS_out
    output$ATS <- ATS_out
    output$MATS <- MATS_out
    output$resampling <- WTPS_out
    output$quantile <- quantiles
    output$nf <- nf
    output$factors <- fac_names
    output$H <- hypo
    output$p <- p
    output$fl <- fl
    output$Means <- mean_out
    # end one-way layout ------------------------------------------------------
  } else {
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    gr <- nadat2[1]
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)
    n <- n$Measure/p
    
    if (length(fac_names) == nf) {
      # nested
      
      # if nested factor is named uniquely
      if (nested.levels.unique){
        # delete factorcombinations which don't exist
        n <- n[n != 0]
        # create correct level combinations
        blev <- list()
        lev_names <- list()
        for (ii in 1:length(levels[[1]])) {
          blev[[ii]] <- levels(as.factor(dat[, 3][dat[, 2] == levels[[1]][ii]]))
          lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
        }
        if (nf == 2) {
          lev_names <- as.factor(unlist(lev_names))
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev)
        } else {
          lev_names <- lapply(lev_names, rep,
                              length(levels[[3]]) / length(levels[[2]]))
          lev_names <- lapply(lev_names, sort)
          lev_names <- as.factor(unlist(lev_names))
          blev <- lapply(blev, rep, length(levels[[3]]) / length(levels[[2]]))
          blev <- lapply(blev, sort)
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev, as.factor(levels[[3]]))
        }
        # correct for wrong counting of nested factors
        if (nf == 2) {
          fl[2] <- fl[2] / fl[1]
        } else if (nf == 3) {
          fl[3] <- fl[3] / fl[2]
          fl[2] <- fl[2] / fl[1]
        }
      }
      hypo_matrices <- HN_MANOVA(fl, p)
    } else {
      # crossed
      hypo_matrices <- HC_MANOVA(fl, perm_names, fac_names, p)[[1]]
      fac_names <- HC_MANOVA(fl, perm_names, fac_names, p)[[2]]
    }
    
    # ---------------------- error detection ------------------------------------
    
    # mixture of nested and crossed designs is not possible
    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is
           not implemented!")
    }
    # only 3-way nested designs are possible
    if (length(fac_names) == nf && nf >= 4) {
      stop("Four- and higher way nested designs are
           not implemented!")
    }
    # no factor combinations with less than 2 observations
    if (0 %in% n || 1 %in% n) {
      stop("There is at least one factor-level combination
           with less than 2 observations!")
    }
    
    if (length(fac_names) != length(hypo_matrices)) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }
    
    #--------------------------------------------------------------------------#
    
    
    n.groups <- prod(fl)
    WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
    ATS_out <- matrix(NA, ncol = 4, nrow = length(hypo_matrices))
    WTPS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 3)
    MATS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 1)
    quantiles <- matrix(NA, ncol = 2, nrow = length(hypo_matrices))
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    rownames(WTPS_out) <- fac_names
    rownames(MATS_out) <- fac_names
    rownames(quantiles) <- fac_names
    colnames(ATS_out) <- c("Test statistic", "df1", "df2", "p-value")
    colnames(MATS_out) <- "Test statistic"
    colnames(quantiles) <- c("WTS_resampling", "MATS_resampling")
    # calculate results
    for (i in 1:length(hypo_matrices)) {
      results <- MANOVA.Stat(data = response, n, hypo_matrices[[i]],
                             iter, alpha, resampling, n.groups, p, CPU, seed, nf)
      WTS_out[i, ] <- results$WTS
      ATS_out[i, ] <- results$ATS
      WTPS_out[i, ] <- results$WTPS
      MATS_out[i] <- results$MATS
      quantiles[i, ] <- results$quantiles
    }
    # time needed for resampling calculations
    time <- results$time
    mean_out <- matrix(results$Mean, ncol = p, byrow = TRUE)
    Var_out <- results$Cov
    descriptive <- cbind(lev_names, n, mean_out)
    colnames(descriptive) <- c(nadat2, "n", paste(rep("Mean", p), 1:p))
    
    # Output ------------------------------------------------------
    colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
    colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(ATS)"), paste(resampling, "(MATS)"))
    output <- list()
    output$time <- time
    output$input <- input_list
    output$Descriptive <- descriptive
    output$Covariance <- Var_out
    output$Means <- mean_out
    output$MATS <- MATS_out
    output$WTS <- WTS_out
    output$ATS <- ATS_out
    output$resampling <- WTPS_out
    output$quantile <- quantiles
    output$nf <- nf
    output$H <- hypo_matrices
    output$factors <- fac_names
    output$p <- p
    output$fl <- fl
  }
  class(output) <- "MANOVA"
  return(output)
}
