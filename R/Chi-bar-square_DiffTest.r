#' Chi-bar-square difference test
#'
#' This function calculates the Chi-bar-square difference test of the RI-CLPM versus the CLPM or a general variance test. The first is discussed in Hamaker, Kuiper, and Grasman (2015) and the second in Stoel et al. (2006). There is also an interactive web application on my website: Chi-bar-square difference test (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param q The total number of variances (q = n + k), that is, the sum of the number of nuisance variances (n) and the number of constrained variances (k) whose values may be on the boundary of the parameter space under the null and the alternative hypotheses, that is, the number of random intercepts. When testing RI-CLPM vs CLPM (see Hamaker et al., 2015), q = k; in case of another test of variances (see Stoel et al., 2006), q >= k.
#' Note: u = k*n + k*(k-1)/2 is the number of unconstrained variances and unconstrained covariances.
#' @param S The k times k covariance matrix of the k constrained variances.
#' @param Chi2_clpm The Chi-square value for the CLPM. In case of another test of variances, this refers to the Chi-square of the model with the k constrained variances (or, in general, the more constrained model). By default, Chi2_clpm = NULL; in that case, only the Chi-bar-square weighst will be given and not the Chi-bar-square statistic with correspoding p-value.
#' @param Chi2_riclpm The Chi-square value for the RI-CLPM. In case of another test of variances, this refers to the Chi-square of the model without the k constrained variances (or, in general, the less constrained model). By default, Chi2_riclpm = NULL; in that case, only the Chi-bar-square weighst will be given and not the Chi-bar-square statistic with correspoding p-value.
#' @param df_clpm The degrees of freedom of the Chi-square test in the CLPM (i.e., the model with the k constrained variances or the more constrained model). By default, df_clpm = NULL; in that case, only the Chi-bar-square weighst will be given and not the Chi-bar-square statistic with correspoding p-value.
#' @param df_riclpm The degrees of freedom of the Chi-square test in the RI-CLPM (i.e., the model without the k constrained variances or the less constrained model). By default, df_riclpm = NULL; in that case, only the Chi-bar-square weighst will be given and not the Chi-bar-square statistic with correspoding p-value.
#' @param alpha The alpha level in determining the significance of the Chi-bar-square difference test. By default, alpha = 0.05.
#' @param bootstrap Indicator (FALSE/TRUE) to determine the Chi-bar-square weigths based on bootstrap (instead of using the package ic.infer). By default, bootstrap = FALSE.
#' @param seed The seed number, used in case bootstrap = TRUE. By changing this number, the sensitivity of the results can be inspected. By default, seed = 123.
#' @param iter The number of iterations, used in case bootstrap = TRUE. By changing this number, the sensitivity/precision of the results can be inspected. By default, iter = 100000.
#' @param u The number of unconstrained parameters of interest. In case of a more general test than a variance test, 'u' can easily be specified and there is no 'q'. For more details, see Stoel et al. (2006). By default, a variance test is assumed and thus the argument u = NULL implying that u = k*n + k*(k-1)/2 is used in the calculation. If 'u' is specified, then 'q' is discarded and can be set to NULL.
#' @param PrintPlot Optional. Indicator whether Chi-bar-square distribution plot should be printed (TRUE; default) or not (FALSE); together with the critical value.
#' @param Min Optional. Minimum value used in the Chi-bar-square distribution plot. By default, Min = 0.
#' @param Max Optional. Maximum time used in the Chi-bar-square distribution plot. By default, Max = 20. If Max is lower than critical value c2, then it is changed to c2+1.
#' @param Step Optional. The step-size taken in the Chi-bar-square distribution plot. By default, Step = 1.
#'
#' @return The output comprises, among others, the Chi-bar-square weigths and critical value for the Chi-bar-square statistic (if asked for) the Chi-bar-square statistic and corresponding p-value.
#' @importFrom quadprog solve.QP
#' @importFrom mvtnorm rmvnorm
#' @importFrom ic.infer ic.weights
#' @importFrom nleqslv nleqslv
#' @export
#' @examples
#'
#' # Compare fit CLPM vs RI-CLPM
#'
#' # Input needed in examples below:
#' # There are 2 random intercepts in the RI-CLPM (omega and kappa): q = k = 2
#' q <- 2
#' # The full covariance matrix of the random intercepts is:
#' S <- matrix(c(1.232, -0.118, -0.118, 0.539), byrow = T, ncol = q)
#' # Chi-square values
#' Chi2_clpm <- 20.6779  #The Chi-square value of CLPM is 20.6779
#' Chi2_riclpm <- 3.2127 #The Chi-square value of RI-CLPM is 3.2127
#' #Degrees of freedom (df)
#' df_clpm <- 4   #The df in CLPM is 4
#' df_riclpm <- 1 #The df in RI-CLPM is 1
#'
#' # Run function to obtain Chi-bar-square weigths and critical value for the Chi-bar-square statistic
#' ChiBarSq.DiffTest(q, S)
#'
#' # Run function to obtain Chi-bar-square weigths based on bootstrap
#' ChiBarSq.DiffTest(q, S, bootstrap = T, seed = 123, iter = 100000)
#'
#' # Run function to do Chi-bar-square test (and also obtain Chi-bar-square weigths and critical value)
#' ChiBarSq.DiffTest(q, S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm)
#'
#' # Run function based on using 'u' as input to do Chi-bar-square test (and also obtain Chi-bar-square weigths and critical value)
#' # For simplicity, we use the same example as above and calculate u based on input above - normally you would know u (and not q).
#' q <- 2
#' k <- 2
#' n <- q-k
#' u <- k*n + k*(k-1)/2 # This expression holds in case of k constrained variances and thus also in case of k random intercepts.
#' #
#' ChiBarSq.DiffTest(q, S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm, u = u)
#' # Or
#' ChiBarSq.DiffTest(NULL, S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm, u = u)
#' # But NOT:
#' #ChiBarSq.DiffTest(S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm, u = u) # Note: This does not work (properly) and gives an error.
#' #
#' # Note: This is now again based on testing the CLPM versus the RI-CLPM, but this code is most helpful in case of a more general test than a 'k constrained variance test'. For more details, see Stoel et al. (2006).
#'

ChiBarSq.DiffTest <- function(q, S, Chi2_clpm = NULL, Chi2_riclpm = NULL, df_clpm = NULL, df_riclpm = NULL, alpha = 0.05, bootstrap = FALSE, seed = 123, iter = 100000, u = NULL, PrintPlot = TRUE, Min = 0, Max = 20, Step = 1) {

  # Checks:
  if(is.null(u)){
    if(length(q) != 1){
      print(paste("The argument q should be an integer (i.e., a scalar which is an integer and not multiple (integer) values."))
      stop()
    }
    if(q %% 1 != 0){
      print(paste("The argument q should be an integer."))
      stop()
    }
  }
  #
  # Check on S
  if(length(S) != 1){
    if(is.null(dim(S))){
      if(!is.null(length(S))){
        print(paste("The argument S is not a matrix of size k times k."))
        stop()
      }else{
        print(paste("The argument S is not found: The k times k covariance matrix of the k constrained variances is unknown, but should be part of the input."))
        stop()
      }
    }
    k <- dim(S)[1]
    #
    if(dim(S)[1] != dim(S)[2]){
      print(paste("The covariance matrix of the k constrained variances S should be a square matrix of size k times k, with k = ", k, "."))
      stop()
    }
    if(length(dim(S)) > 2){
      print(paste("The covariance matrix of the k constrained variances S should be an k times k matrix, with k = ", k, "."))
      stop()
    }
  } else{
    k <- 1
  }

  if(is.null(u)){
    n <- q-k
    if(n < 0){
      print(paste("The input for k (i.e., the number of constrained variances) cannot exceed q (i.e., the number of latent variables, that is, the total number of variances). Either q or the dimension of S is incorrect, since n = q - k < 0, with q = ", q, "and k = ", k, "."))
      stop()
    }
    u <- k*n + k*(k-1)/2
  }else{
    if(length(u) != 1){
      print(paste("The argument u should be a scalar, that is, one number, that is, a vector with one element."))
      stop()
    }
  }
  #u = unconstrained parameter of interest = number of unconstrained variances and unconstrained covariances; u <- k*n + k*(k-1)/2, n = number of nuisance parameters
  #k = number of constrained variances
  #The other parameters are called nuisance parameters and are not needed in the calculation; the number of nuisance parameters is denoted by n.
  #S = the k-times-k covariance matrix of the (k) constrained variances.

  if(!is.null(Chi2_clpm) & length(Chi2_clpm) != 1){
    print(paste("The argument Chi2_clpm should be a scalar, that is, one number, that is, a vector with one element (or NULL)."))
    stop()
  }
  if(!is.null(Chi2_riclpm) & length(Chi2_riclpm) != 1){
    print(paste("The argument Chi2_riclpm should be a scalar, that is, one number, that is, a vector with one element (or NULL)."))
    stop()
  }
  if(!is.null(df_clpm)){
    if(length(df_clpm) != 1){
      print(paste("The argument df_clpm should be an integer (i.e., a scalar which is an integer and not multiple (integer) values."))
      stop()
    }
    if(df_clpm %% 1 != 0){
      print(paste("The argument df_clpm should be an integer."))
      stop()
    }
  }
  if(!is.null(df_riclpm)){
    if(length(df_riclpm) != 1){
      print(paste("The argument df_riclpm should be an integer (i.e., a scalar which is an integer and not multiple (integer) values."))
      stop()
    }
    if(df_riclpm %% 1 != 0){
      print(paste("The argument df_riclpm should be an integer."))
      stop()
    }
  }
  #
  if(length(alpha) != 1){
    print(paste("The argument alpha should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  #
  if(!is.logical(bootstrap)){
    print(paste("The argument bootstrap should be logical, that is, have the value T(RUE) or F(ALSE)."))
    stop()
  }
  if(bootstrap == TRUE){
    if(length(seed) != 1){
      print(paste("The argument seed should be a scalar, that is, one number, that is, a vector with one element."))
      stop()
    }
    if(length(iter) != 1){
      print(paste("The argument iter should be a scalar, that is, one number, that is, a vector with one element."))
      stop()
    }
  }
  #
  if(PrintPlot != T & PrintPlot != F){
    print(paste("The argument 'PrintPlot' should be TRUE or FALSE, not ", PrintPlot, "."))
    stop()
  }
  if(length(Min) != 1){
    print(paste("The argument Min should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  if(length(Max) != 1){
    print(paste("The argument Max should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }
  if(length(Step) != 1){
    print(paste("The argument Step should be a scalar, that is, one number, that is, a vector with one element."))
    stop()
  }


  ##########################################################################

  # min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.

  #if (!require("quadprog")) install.packages("quadprog")
  #if (!require("mvtnorm")) install.packages("mvtnorm")
  #if (!require("nleqslv")) install.packages("nleqslv")
  #library(quadprog)
  #library(mvtnorm)
  #library(nleqslv)
  #
  #if (!require("ic.infer")) install.packages("ic.infer")
  #library(ic.infer)

  ##########################################################################
  ##########################################################################


  # LPs/Weights
  if(bootstrap == TRUE){
    set.seed(seed)
    Z <- rmvnorm(n=iter, mean=rep(0, k), sigma=S)
    invS = solve(S)
    nact <- apply(Z, 1, function(z){
            dvec <- 2*(z %*% invS)
            QP <- solve.QP(Dmat=2*invS, dvec, Amat=diag(k), bvec=matrix(0, ncol=k, nrow=1), meq=0)
    	if(QP$iact[1] != 0){
    	length(QP$iact)}else{
    	0}
    })
    dimsol <- k - nact
    weight <- sapply(1:(k+1), function(x) sum(x == (dimsol+1)))/iter
  } else{
    wt.bar <- ic.weights(S)
    weight <- wt.bar[length(wt.bar):1]
  }


  # We want to find c^2 such that sum_{i=1}^{k} [weight[i] pr(Chi^2(i) >= c^2)] = alpha
  FindC2_fie <- function(weight, k, u, alpha) {
    # x = c2
    function(x) {
      y <- numeric(1)
      p = 0
      for(teller in 1:(k+1)){
        p = p + weight[teller] * (1-pchisq(x, (u+teller-1))) #df_Chi2 = u until u+k
      }
      y <- p - alpha
      y
    }
  }
  FindC2 <- FindC2_fie(weight, k, u, alpha)
  xstart <- 5
  #FindC2(xstart)
  sol <- nleqslv(xstart, FindC2, control=list(btol=.001, allowSingular = TRUE))
  if(sol$fvec > 0.01){
    message = "For these data, the critical value (c^2) cannot be calculated."
    final <- list(message = message, CovMx_of_k_constrained_variances = S, ChiBar2_weights = LP)
    return(final)
  }else{
    c2 <- sol$x


    # Plot
    if(PrintPlot == T){
      if(Max < c2){Max <- c2 + 1}
      X <- seq(Min,Max,by=Step)
      ChiMix <- 0
      for(i in 1:k){
        ChiMix <- ChiMix + weight[i+1]*dchisq(X, (u+i))
      }
      plot(X, ChiMix, xlab = "x", ylab = "Chi-bar-square(x)", main = "Chi-bar-square distribution")
      lines(X, ChiMix)
      abline(v=c2, col = "red")
      #
      legend("topright",
             legend = c("critical value"),
             lty = 1,
             col = c("red")
      )
    }


    #if( (Chi2_clpm == 0) & (Chi2_riclpm == 0) ){
    if( is.null(Chi2_clpm) & is.null(Chi2_riclpm) ){
      message = paste0("The p-value of the Chi2-bar difference test will be calculated after filling in the Chi2's of the CLPM and RI-CLPM (and their degrees of freedom).")
      final <- list(message = message,
                    k=k, u=u, S = S, ChiBar2_weights = weight,
                    critical_value = c2)
      return(final)
    }else{
      # -- Determine p-value --
      names(weight) <- NULL
      c2_compare = Chi2_clpm - Chi2_riclpm
      CritValSmallerDiffChi2 <- (c2 < c2_compare)
      p = 0
      for(teller in 1:(k+1)){
        p = p + weight[teller] * (1-pchisq(c2_compare, (u+teller-1))) #df_Chi2 = u until u+k
      }
      pSmallerAlpha <- (p < alpha)
      #if(p < 0.001){
      #  p <- paste(p, " < .001")
      #}


      # Plot
      if(PrintPlot == T){
        if(Max < c2_compare){
          Max <- c2_compare + 1
          X <- seq(Min,Max,by=Step)
          ChiMix <- 0
          for(i in 1:k){
            ChiMix <- ChiMix + weight[i+1]*dchisq(X, (u+i))
          }
        }
        plot(X, ChiMix, xlab = "x", ylab = "Chi-bar-square(x)", main = "Chi-bar-square distribution")
        lines(X, ChiMix)
        abline(v=c2, col = "red")
        abline(v=c2_compare, col = "blue")
        #
        legend("topright",
               legend = c("critical value", "observed value"),
               lty = 1,
               col = c("red", "blue")
        )
      }


      diff_df = df_clpm - df_riclpm
      if(diff_df != k*(k+1)/2){
        #Could also calculate p-values based on diff_df, but that is weird, most probably df_diff consists of k and u (so, diff_df is not equal to k)
        #p_diffdf = 0
        #for(teller in 1:(k+1)){
        #  p_diffdf = p_diffdf + weight[teller] * (1-pchisq(c2_compare, (u+teller-1))) #df_Chi2 = u until u+k
        #}
        message = paste0("The difference in degrees of freedom, ", diff_df, ", does not equal k*(k+1)/2 = ", (k*(k+1)/2), ", please check your input. The p-value below is determined based on k = ", k, " and u = ", u, ".")
        final <- list(message = message,
                      k=k, u=u, S = S, ChiBar2_weights = weight,
                      critical_value = c2, DiffChi2 = c2_compare, CritValSmallerDiffChi2 = CritValSmallerDiffChi2, p_value=p, pSmallerAlpha = pSmallerAlpha )
        #k=k, u=u, CovMx_of_k_constrained_variances = S, ChiBar2_weights = weight,
        #critical_value = c2, DiffChi2 = c2_compare, CritValSmallerDiffChi2 = CritValSmallerDiffChi2, p_value=p, pSmallerAlpha = pSmallerAlpha )
        #p_value_diffdf=p_diffdf,
        return(final)
      }else if(c2_compare < 0){
        message = paste0("The difference in Chi2's (Chi2_clpm - Chi2_riclpm = ", c2_compare, ") is negative. You probably have switched them around.")
        final <- list(message = message,
        k=k, u=u, S = S, ChiBar2_weights = weight,
        critical_value = c2, DiffChi2 = c2_compare, CritValSmallerDiffChi2 = CritValSmallerDiffChi2, p_value=p, pSmallerAlpha = pSmallerAlpha )
      }else{
        #final <- list(k=k, u=u, CovMx_of_k_constrained_variances = S, ChiBar2_weights = weight,
        #              critical_value = c2, DiffChi2 = c2_compare, CritValSmallerDiffChi2 = CritValSmallerDiffChi2, p_value=p, pSmallerAlpha = pSmallerAlpha)
        final <- list(k=k, u=u, S = S, ChiBar2_weights = weight,
                      critical_value = c2, DiffChi2 = c2_compare, CritValSmallerDiffChi2 = CritValSmallerDiffChi2, p_value=p, pSmallerAlpha = pSmallerAlpha)
        return(final)
      }

    }
  }


} # end of function

