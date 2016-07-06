#' Tests for Multivariate Data in Semi-Parametric Factorial Designs
#' 
#' The MANOVA function calculates the Wald-type statistic (WTS), the ANOVA-type 
#' statistic (ATS) as well as resampling versions of these test statistics for 
#' semi-parametric multivariate data.
#' 
#' @param formula A model \code{\link{formula}} object. The left hand side 
#'   contains the response variable and the right hand side contains the factor 
#'   variables of interest. An interaction term must be specified.
#' @param data A data.frame, list or environment containing the variables in 
#'   \code{formula}. 
#' @param subject the column name of the subjects in the data.
#' @param iter The number of iterations used for calculating the resampled 
#'   statistic. The default option is 10000.
#' @param alpha A number specifying the significance level; the default is 0.05.
#' @param resampling The resampling method to be used, one of "paramBS"
#'   (parametric bootstrap approach) and "WildBS" (wild bootstrap approach with
#'   Rademacher weights). The Wild Bootstrap is applied to both test statistics.
#' @param CPU The number of cores used for parallel computing. If omitted, cores are
#'   detected via \code{\link[parallel]{detectCores}}.
#'  
#'   
#' @details The MANOVA() function provides the Wald-type statistic as well as
#'   the ANOVA-type statistic for multivariate designs with metric data as described in 
#'   Konietschke et al. (2015). These tests are even applicable for non-normal error terms, 
#'   different sample sizes and/or heteroscedastic variances.
#'   They are implemented for designs with an arbitrary number of
#'   crossed factors or for nested designs. In addition to the asymptotic
#'   p-values, it also provides p-values based on resampling approaches.
#'  
#' @section NOTE: The number of resampling iterations has been set to 100 in the examples due to run time 
#' restrictions on CRAN. Usually it is recommended to use at least 1000 iterations.
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
#'   Friedrich, S, Brunner, E and Pauly, M (2016). Permuting longitudinal data 
#'   despite all the dependencies. arXiv preprint arXiv:1509.05570v2
#'   
#'    Bathke, A., Friedrich, S., Konietschke, F., Pauly, M., Staffen, W., Strobl, N. and Hoeller, Y. (2016). Using EEG, SPECT, and Multivariate Resampling Methods
#' to Differentiate Between Alzheimer's and other Cognitive Impairments. arXiv preprint arXiv:1606.09004.
#'   
#'   Friedrich, S., Konietschke, F., Pauly, M.(2016). GFD - An 
#'   R-package for the Analysis of General Factorial Designs. Accepted for publication in 
#'   Journal of Statistical Software.
#'    
#'   
#' @importFrom graphics axis legend par plot title
#' @importFrom stats ecdf formula model.frame pchisq pf qt terms var cov rbinom
#' @importFrom utils read.table
#' @importFrom methods hasArg
#' @importFrom MASS mvrnorm
#' @importFrom parallel makeCluster parSapply detectCores
#'   
#' @export

MANOVA <- function(formula, data, subject,
                   iter = 10000, alpha = 0.05, resampling = "paramBS", CPU){
  
  if (!(resampling %in% c("paramBS", "WildBS"))){
    stop("Resampling must be one of 'paramBS' and 'WildBS'!")
  }
  
  input_list <- list(formula = formula, data = data,
                     subject = subject, 
                     iter = iter, alpha = alpha, resampling = resampling)
  
  test <- hasArg(CPU)
  if(!test){
    CPU <- parallel::detectCores()
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
    dat2 <- dat2[order(dat2$subject), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    # contrast matrix
    hypo <- (diag(fl) - matrix(1 / fl, ncol = fl, nrow = fl)) %x% diag(p)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(unique(subject)),
                     .drop = F)$Measure
    WTS_out <- matrix(NA, ncol = 3, nrow = 1)
    ATS_out <- matrix(NA, ncol = 4, nrow = 1)
    WTPS_out <- rep(NA, 2)
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    names(WTPS_out) <- fac_names
    results <- MANOVA.Stat(data = response, n = n, hypo, iter = iter, alpha, resampling, n.groups = fl, p, CPU)    
    WTS_out <- results$WTS
    ATS_out <- results$ATS
    WTPS_out <- results$WTPS
    mean_out <- matrix(results$Mean, ncol = p, byrow = TRUE)
    Var_out <- results$Cov
#    CI <- results$CI
#    colnames(CI) <- c("lower", "upper")
    descriptive <- cbind(lev_names, n, mean_out)
    colnames(descriptive) <- c(nadat2, "n", rep("Means", p))    
    names(WTS_out) <- cbind ("Test statistic", "df",
                                "p-value")
    names(ATS_out) <- cbind("Test statistic", "df1", "df2", "p-value")
    names(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(ATS)"))
    output <- list()
    output$input <- input_list
    output$Descriptive <- descriptive
#    output$CI <- CI
    output$Covariance <- Var_out
    output$WTS <- WTS_out
    output$ATS <- ATS_out
    output$resampling <- WTPS_out
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
      # delete factorcombinations which don't exist
      n <- n[n != 0]
    }
    
# ---------------------- error detection ------------------------------------
    
    # mixture of nested and crossed designs is not possible
    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is
           not impemented!")
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

    # correct labeling of factors in nested design
    if (length(fac_names) == nf) {
      if (nf == 2) {
        if (all(levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][1]]))
                == levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][2]])))) {
          stop("The levels of the nested factor must be
               named without repetitions!")
        }
        } else if (nf == 3) {
          if (all(levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][1]]))
                  == levels(as.factor(dat2[, 3][dat2[, 2] == levels[[1]][2]]))) ||
                all(levels(as.factor(dat2[, 4][dat2[, 3] == levels[[2]][1]]))
                    == levels(as.factor(dat2[, 4][dat2[, 3] == levels[[2]][fl[2] / fl[1] + 1]])))) {
            stop("The levels of the nested factor must be
                 named without repetitions!")
          }
      }
    }

#--------------------------------------------------------------------------#

    if (length(fac_names) == nf) {
      # nested
      hypo_matrices <- HN_MANOVA(fl, p)
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
    } else {
      # crossed
      hypo_matrices <- HC_MANOVA(fl, perm_names, fac_names, p)[[1]]
      fac_names <- HC_MANOVA(fl, perm_names, fac_names, p)[[2]]
    }
    
    if (length(fac_names) != length(hypo_matrices)) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }

    n.groups <- prod(fl)
    WTS_out <- matrix(NA, ncol = 3, nrow = length(hypo_matrices))
    ATS_out <- matrix(NA, ncol = 4, nrow = length(hypo_matrices))
    WTPS_out <- matrix(NA, nrow = length(hypo_matrices), ncol = 2)
    rownames(WTS_out) <- fac_names
    rownames(ATS_out) <- fac_names
    rownames(WTPS_out) <- fac_names
    colnames(ATS_out) <- c("Test statistic", "df1", "df2", "p-value")
    # calculate results
    for (i in 1:length(hypo_matrices)) {
        results <- MANOVA.Stat(data = response, n, hypo_matrices[[i]],
                           iter, alpha, resampling, n.groups, p, CPU)
        WTS_out[i, ] <- results$WTS
        ATS_out[i, ] <- results$ATS
        WTPS_out[i, ] <- results$WTPS
    }
    # time needed for resampling calculations
    time <- results$time
    mean_out <- matrix(results$Mean, ncol = p, byrow = TRUE)
    Var_out <- results$Cov
#    lev_names2 <- matrix(rep(as.matrix(lev_names), each = p), nrow = n.groups * p)
#    CI <- data.frame(lev_names2, rep(1:p, n.groups), results$CI)
#    colnames(CI) <- c(nadat2, "Dimension", "lower", "upper")
    descriptive <- cbind(lev_names, n, mean_out)
    colnames(descriptive) <- c(nadat2, "n", paste(rep("Mean", p), 1:p))
    
    # Output ------------------------------------------------------
    colnames(WTS_out) <- cbind ("Test statistic", "df", "p-value")
    colnames(WTPS_out) <- cbind(paste(resampling, "(WTS)"), paste(resampling, "(ATS)"))
    output <- list()
    output$time <- time
    output$input <- input_list
    output$Descriptive <- descriptive
    output$Covariance <- Var_out
#    output$CI <- CI
    output$WTS <- WTS_out
    output$ATS <- ATS_out
    output$resampling <- WTPS_out
#    output$plotting <- list(levels, fac_names, nf, mu, lower, upper, fac_names_original, dat2, fl, alpha, nadat2, lev_names)
#    names(output$plotting) <- c("levels", "fac_names", "nf", "mu", "lower", "upper", "fac_names_original", "dat2", "fl", "alpha", "nadat2", "lev_names")
  }
  class(output) <- "MANOVA"
  return(output)
}
