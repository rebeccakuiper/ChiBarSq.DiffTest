#' Chi-bar-square difference test
#'
#' This function calculates the Chi-bar-square difference test of the RI-CLPM versus the CLPM or a general variance test. The first is discussed in Hamaker, Kuiper, and Grasman (2015) and the second in Stoel et al. (2006). There is also an interactive web application on my website: Chi-bar-square difference test (\url{https://www.uu.nl/staff/RMKuiper/Websites\%20\%2F\%20Shiny\%20apps}).
#'
#' @param q The total number of variances (q = n + k), that is, the sum of the number of nuisance variances (n) and the number of constrained variances (k) whose values may be on the boundary of the parameter space under the null and the alternative hypotheses, that is, the number of random intercepts. When testing RI-CLPM vs CLPM (see Hamaker et al., 2015), q = k; in case of another test of variances (see Stoel et al., 2006), q >= k.
#' Note: u = k*n + k*(k-1)/2 is the number of unconstrained variances and unconstrained covariances.
#' @param S The k times k covariance matrix of the k constrained variances.
#' @param Chi2_clpm The Chi-square value for the CLPM. In case of another test of variances, this refers to the Chi-square of the model with the k constrained variances (or, in general, the more constrained model). By default, Chi2_clpm = NULL; in that case, only the Chi-bar-square weights will be given and not the Chi-bar-square statistic with corresponding p-value.
#' @param Chi2_riclpm The Chi-square value for the RI-CLPM. In case of another test of variances, this refers to the Chi-square of the model without the k constrained variances (or, in general, the less constrained model). By default, Chi2_riclpm = NULL; in that case, only the Chi-bar-square weights will be given and not the Chi-bar-square statistic with corresponding p-value.
#' @param df_clpm The degrees of freedom of the Chi-square test in the CLPM (i.e., the model with the k constrained variances or the more constrained model). By default, df_clpm = NULL; in that case, only the Chi-bar-square weights will be given and not the Chi-bar-square statistic with corresponding p-value.
#' @param df_riclpm The degrees of freedom of the Chi-square test in the RI-CLPM (i.e., the model without the k constrained variances or the less constrained model). By default, df_riclpm = NULL; in that case, only the Chi-bar-square weights will be given and not the Chi-bar-square statistic with corresponding p-value.
#' @param SB_clpm Optional. The Satorra-Bentler scaled Chi-square value for the CLPM (i.e., the model with the k constrained variances or the more constrained model). By default, SB_clpm = NULL; in that case, the difference test will be solely based on the unadjusted/regular Chi-square value.
#' @param SB_riclpm Optional. The Satorra-Bentler scaled Chi-square value for the RI-CLPM (i.e., the model without the k constrained variances or the less constrained model). By default, SB_riclpm = NULL; in that case, the difference test will be solely based on the unadjusted/regular Chi-square value.
#' @param alpha The alpha level in determining the significance of the Chi-bar-square difference test. By default, alpha = 0.05.
#' @param bootstrap Indicator (TRUE/FALSE or 1/0) whether the Chi-bar-square weights are determined using bootstrap (TRUE or 1) or using the package ic.infer (FALSE or 0; default). By default, bootstrap = FALSE (and thus ic.infer is used).
#' @param seed The seed number, used in case bootstrap = TRUE. By changing this number, the sensitivity of the results can be inspected. By default, seed = 123.
#' @param iter The number of iterations, used in case bootstrap = TRUE. By changing this number, the sensitivity/precision of the results can be inspected. By default, iter = 100000.
#' @param u The number of unconstrained parameters of interest. In case of a more general test than a variance test, 'u' can easily be specified and there is no 'q'. For more details, see Stoel et al. (2006). By default, a variance test is assumed and thus the argument u = NULL implying that u = k*n + k*(k-1)/2 is used in the calculation. If 'u' is specified, then 'q' is discarded and can be set to NULL.
#' @param PrintPlot Optional. Indicator (TRUE/FALSE or 1/0) whether Chi-bar-square distribution plot should be printed (TRUE or 1) or not (FALSE or 0); together with the critical value. By default, PrintPlot = TRUE (and thus a plot is rendered).
#' @param Min Optional. Minimum value used in the Chi-bar-square distribution plot. By default, Min = 0.
#' @param Max Optional. Maximum time used in the Chi-bar-square distribution plot. By default, Max = 20. If Max is lower than critical value c2, then it is changed to c2+1.
#' @param Step Optional. The step-size taken in the Chi-bar-square distribution plot. By default, Step = 1.
#'
#' @return The output comprises, among other things, the Chi-bar-square weights and critical value for the Chi-bar-square statistic (if asked for) the Chi-bar-square statistic and corresponding p-value.
#' @importFrom quadprog solve.QP
#' @importFrom mvtnorm rmvnorm
#' @importFrom ic.infer ic.weights
#' @importFrom nleqslv nleqslv
#' @importFrom ggplot2 ggplot
#' @export print.ChiBar2
#' @export summary.ChiBar2
#' @export
#' @examples
#'
#' # library(ChiBarSq.DiffTest)
#'
#' # Compare fit CLPM vs RI-CLPM
#' # For an elaborate example, see https://ellenhamaker.github.io/RI-CLPM/lavaan.html#(bar{chi}^{2})-test.
#'
#' # Specification of input needed in the examples below:
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
#' # Run function to obtain Chi-bar-square weights and critical value for the Chi-bar-square statistic
#' ChiBarSq.DiffTest(q, S)
#'
#' # Run function to obtain Chi-bar-square weights based on bootstrap
#' ChiBarSq.DiffTest(q, S, bootstrap = T, seed = 123, iter = 100000)
#'
#' # Run function to do Chi-bar-square test (and also obtain Chi-bar-square weights and critical value)
#' ChiBarSq.DiffTest(q, S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm)
#'
#' # Different types of output options are possible:
#' ChiBar2 <- ChiBarSq.DiffTest(q, S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm)
#' ChiBar2 # same as print(ChiBar2)
#' summary(ChiBar2)
#' print(ChiBar2, digits = 4)
#' summary(ChiBar2, digits = 4)
#' # In Rstudio, use 'ChiBar2$' to see what output there is:
#' # ChiBar2$message; ChiBar2$k; ChiBar2$u; ChiBar2$S; ChiBar2$ChiBar2_weights; ChiBar2$critical_value; ChiBar2$DiffChi2; ChiBar2$CritValSmallerDiffChi2; ChiBar2$p_value; ChiBar2$pSmallerAlpha
#' # and possibly: ChiBar2$ChiBar2_plot
#'
#'
#' # Run function based on using 'u' as input to do Chi-bar-square test (and also obtain Chi-bar-square weights and critical value)
#' # For simplicity, we use the same example as above and calculate u based on input above - normally you would know u (and not q).
#' q <- 2
#' k <- 2
#' n <- q-k
#' u <- k*n + k*(k-1)/2 # This expression holds in case of k constrained variances and thus also in case of k random intercepts.
#' #
#' ChiBarSq.DiffTest(q, S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm, u = u)
#' # Or
#' ChiBarSq.DiffTest(NULL, S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm, u = u)
#' # But do NOT leave the q-argument out like this:
#' #ChiBarSq.DiffTest(S, Chi2_clpm, Chi2_riclpm, df_clpm, df_riclpm, u = u) # Note: This does not work (properly) and gives an error.
#' # It is possible to use: ChiBarSq.DiffTest(S = S, Chi2_clpm = Chi2_clpm, Chi2_riclpm = Chi2_riclpm, df_clpm = df_clpm, df_riclpm = df_riclpm, u = u)
#' #
#' # Note: This is now again based on testing the CLPM versus the RI-CLPM, but this code is most helpful in case of a more general test than a 'k constrained variance test'. For more details, see Stoel et al. (2006).
#'

ChiBarSq.DiffTest <- function(q, S, Chi2_clpm = NULL, Chi2_riclpm = NULL, df_clpm = NULL, df_riclpm = NULL, SB_clpm = NULL, SB_riclpm = NULL, alpha = 0.05, bootstrap = FALSE, seed = 123, iter = 100000, u = NULL, PrintPlot = TRUE, Min = 0, Max = 20, Step = 1) {

  # TO DO RMK: q can be NULL, then u is needed, make sure that that works!
  # Check whether SAs already changed that!

  # Checks:
  if(length(q) != 1){
    print(paste0("The argument q should be an integer (i.e., a scalar which is an integer and not multiple (integer) values. Currently, q = ", q))
    stop()
  }
  if(q %% 1 != 0){
    print(paste0("The argument q should be an integer. Currently, q = ", q))
    stop()
  }
  #
  # Check on S
  if(length(S) != 1){
    if(is.null(dim(S))){
      if(!is.null(length(S))){
        print(paste0("The argument S is not a matrix of size k times k. Currently, it is of size ", dim(S)[1], " times ", dim(S)[2], "."))
        stop()
      }else{
        print(paste0("The argument S is not found: The k times k covariance matrix of the k constrained variances S is unknown, but should be part of the input."))
        stop()
      }
    }
    k <- dim(S)[1]
    #
    if(length(dim(S)) < 2){
      print(paste0("The covariance matrix of the k constrained variances S should be an k times k matrix."))
      stop()
    }
    if(length(dim(S)) > 2){
      print(paste0("The covariance matrix of the k constrained variances S should be an k times k matrix. Currently, it is of size ", dim(S)))
      stop()
    }
    if(dim(S)[1] != dim(S)[2]){
      print(paste0("The covariance matrix of the k constrained variances S should be a square matrix. It should be a matrix of size k times k. Currently, it is of size ", dim(S)[1], " times ", dim(S)[2], "."))
      stop()
    }
  } else{
    k <- 1
    S <- as.matrix(S)
  }

  if(is.null(u)){
    n <- q-k
    if(n < 0){
      print(paste0("The input for k (i.e., the number of constrained variances, that is, the dimension of S) cannot exceed q (i.e., the number of latent variables, that is, the total number of variances). Either q or k is incorrect, since n = q - k < 0, with q = ", q, "and k = ", k, "."))
      stop()
    }
    u <- k*n + k*(k-1)/2
  }else{
    if(length(u) != 1){
      print(paste0("The argument u should be a scalar, that is, one number, that is, a vector with one element. Currently, u = ", u))
      stop()
    }
  }
  #u = unconstrained parameter of interest = number of unconstrained variances and unconstrained covariances; u <- k*n + k*(k-1)/2, n = number of nuisance parameters
  #k = number of constrained variances
  #The other parameters are called nuisance parameters and are not needed in the calculation; the number of nuisance parameters is denoted by n.
  #S = the k-times-k covariance matrix of the (k) constrained variances.

  if(!is.null(Chi2_clpm) & length(Chi2_clpm) != 1){
    print(paste0("The argument Chi2_clpm should be a scalar, that is, one number, that is, a vector with one element (or NULL). Currently, Chi2_clpm = ", Chi2_clpm))
    stop()
  }
  if(!is.null(Chi2_riclpm) & length(Chi2_riclpm) != 1){
    print(paste0("The argument Chi2_riclpm should be a scalar, that is, one number, that is, a vector with one element (or NULL). Currently, Chi2_riclpm = ", Chi2_riclpm))
    stop()
  }
  if(!is.null(df_clpm)){
    if(length(df_clpm) != 1){
      print(paste0("The argument df_clpm should be an integer (i.e., a scalar which is an integer and not multiple (integer) values. Currently, df_clpm = ", df_clpm))
      stop()
    }
    if(df_clpm %% 1 != 0){
      print(paste0("The argument df_clpm should be an integer. Currently, df_clpm = ", df_clpm))
      stop()
    }
  }
  if(!is.null(df_riclpm)){
    if(length(df_riclpm) != 1){
      print(paste0("The argument df_riclpm should be an integer (i.e., a scalar which is an integer and not multiple (integer) values. Currently, df_riclpm = ", df_riclpm))
      stop()
    }
    if(df_riclpm %% 1 != 0){
      print(paste0("The argument df_riclpm should be an integer. Currently, df_riclpm = ", df_riclpm))
      stop()
    }
  }
  #
  if(!is.null(SB_clpm) & length(SB_clpm) != 1){
    print(paste0("The argument SB_clpm should be a scalar, that is, one number, that is, a vector with one element (or NULL). Currently, SB_clpm = ", SB_clpm))
    stop()
  }
  if(!is.null(SB_riclpm) & length(SB_riclpm) != 1){
    print(paste0("The argument SB_riclpm should be a scalar, that is, one number, that is, a vector with one element (or NULL). Currently, SB_riclpm = ", SB_riclpm))
    stop()
  }
  #
  if(length(alpha) != 1){
    print(paste0("The argument alpha should be a scalar, that is, one number, that is, a vector with one element. Currently, alpha = ", alpha))
    stop()
  }
  #
  if(!is.logical(bootstrap) & bootstrap != FALSE & bootstrap != TRUE){
    print(paste0("The argument bootstrap should be logical, that is, have the value T(RUE) or F(ALSE); or 1 or 0; not ", bootstrap))
    stop()
  }
  if(bootstrap == TRUE){
    if(length(seed) != 1){
      print(paste0("The argument seed should be a scalar, that is, one number, that is, a vector with one element. Currently, seed = ", seed))
      stop()
    }
    if(length(iter) != 1){
      print(paste0("The argument iter should be a scalar, that is, one number, that is, a vector with one element. Currently, iter = ", iter))
      stop()
    }
  }
  #
  if(!is.logical(PrintPlot) & PrintPlot != FALSE & PrintPlot != TRUE){
    print(paste0("The argument 'PrintPlot' should be T(RUE) or F(ALSE); or 1 or 0; not ", PrintPlot))
    stop()
  }
  if(length(Min) != 1){
    print(paste0("The argument Min should be a scalar, that is, one number, that is, a vector with one element. Currently, Min = ", Min))
    stop()
  }
  if(length(Max) != 1){
    print(paste0("The argument Max should be a scalar, that is, one number, that is, a vector with one element. Currently, Max = ", Max))
    stop()
  }
  if(length(Step) != 1){
    print(paste0("The argument Step should be a scalar, that is, one number, that is, a vector with one element. Currently, Step = ", Step))
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
    message = "For these data, the critical value for the Chi-bar-square difference test (c^2) cannot be calculated."
    final <- list(message = message,
                  k=k, u=u, S = S, ChiBar2_weights = weight)
  }else{
    c2 <- sol$x

    #if( (Chi2_clpm == 0) & (Chi2_riclpm == 0) ){
    if( is.null(Chi2_clpm) & is.null(Chi2_riclpm) ){

      # Plot
      if(PrintPlot == T){
        if(Max < c2){Max <- c2 + 1}
        X <- seq(Min, Max, by=Step)
        ChiMix <- 0
        for(i in 0:k){
          ChiMix <- ChiMix + weight[i+1]*dchisq(X, (u+i))
        }
        df <- data.frame(
          X = X,
          ChiMix = ChiMix
        )
        #
        Xlab <- "x"
        Ylab <- expression({bar(Chi)^2} (x)) #"Chi-bar-square(x)"
        Title <- expression({bar(Chi)^2}~('Chi'-bar-square)~distribution) #"Chi-bar-square distribution"
        Lty <- 1
        Col <- c("red")
        legendT <- c("critical value")
        Labels <- rep(as.character(legendT), length(X))
        #
        #library(ggplot2)
        ChiBar2_plot <- ggplot(df, aes(X, ChiMix)) +
          geom_line(lwd = 0.75, color = "black") +
          geom_vline(aes(xintercept = c2, color = "critical_value")) +
          scale_linetype_manual(name = " ", values = Lty, labels = legendT) +
          scale_color_manual(name = " ", values = c(critical_value = Col[1]), labels = legendT) +
          ylab(Ylab) +
          xlab(Xlab) +
          ggtitle(Title) +
          theme_classic() +
          theme(plot.title = element_text(margin = margin(t = 20))) +
          ylim(0,1) +
          theme(
            legend.key.width = unit(1, "lines"),
            legend.spacing.x = unit(1, "lines"),
            legend.text = element_text(size = 12)
          ) #; ChiBar2_plot
        #
        #plot(X, ChiMix, xlab = "x", ylab = "Chi-bar-square(x)", main = "Chi-bar-square distribution")
        #lines(X, ChiMix)
        #abline(v=c2, col = "red")
        ##
        #legend("topright",
        #       legend = c("critical value"),
        #       lty = 1,
        #       col = c("red")
        #)
      }

      message = paste0("The observed value and corresponding p-value of the Chi-bar-square difference test will be calculated after filling in the Chi2's of the CLPM and RI-CLPM (and their degrees of freedom).")
      if(PrintPlot == T){
        final <- list(message = message,
                      k=k, u=u, S = S, ChiBar2_weights = weight,
                      critical_value = c2,
                      ChiBar2_plot = ChiBar2_plot)
        print(ChiBar2_plot)
      }else{
        final <- list(message = message,
                      k=k, u=u, S = S, ChiBar2_weights = weight,
                      critical_value = c2)
      }
    }else{
      # -- Determine p-value --
      names(weight) <- NULL
      if(!is.null(SB_clpm) & !is.null(SB_clpm)){
      scf_clpm = Chi2_clpm/SB_clpm
      scf_riclpm = Chi2_riclpm/SB_riclpm
      cd = (df_clpm * scf_clpm - df_riclpm * scf_riclpm)/(df_clpm - df_riclpm)
      c2_compare = (SB_clpm * scf_clpm - SB_riclpm * scf_riclpm)/cd
      } else{
        c2_compare <- Chi2_clpm - Chi2_riclpm
      }
      CritValSmallerDiffChi2 <- (c2 < c2_compare)
      p <- 0
      for(teller in 1:(k+1)){
        p <- p + weight[teller] * (1-pchisq(c2_compare, (u+teller-1))) #df_Chi2 = u until u+k
      }
      pSmallerAlpha <- (p < alpha)
      #if(p < 0.001){
      #  p <- paste0(p, " < .001")
      #}


      # Plot
      if(PrintPlot == T){
        if(Max < c2_compare){
          Max <- c2_compare + 1
        }
        X <- seq(Min, Max, by=Step)
        ChiMix <- 0
        for(i in 0:k){
          ChiMix <- ChiMix + weight[i+1]*dchisq(X, (u+i))
        }
        df <- data.frame(
          X = X,
          ChiMix = ChiMix
        )
        #
        Xlab <- "x"
        Ylab <- expression({bar(Chi)^2} (x)) #"Chi-bar-square(x)"
        Title <- expression({bar(Chi)^2}~('Chi'-bar-square)~distribution) #"Chi-bar-square distribution"
        Lty <- 1
        Col <- c("red", "blue")
        legendT <- c("critical value", "observed value")
        Labels <- rep(as.character(legendT), length(X))
        #
        #library(ggplot2)
        ChiBar2_plot <- ggplot(df, aes(X, ChiMix)) +
          geom_line(lwd = 0.75, color = "black") +
          geom_vline(aes(xintercept = c2, color = "critical")) +
          geom_vline(aes(xintercept = c2_compare, color = "observed")) +
          scale_linetype_manual(name = " ", values = Lty, labels = legendT) +
          scale_color_manual(name = " ", values = c(critical = Col[1], observed = Col[2]), labels = legendT) +
          ylab(Ylab) +
          xlab(Xlab) +
          ggtitle(Title) +
          theme_classic() +
          theme(plot.title = element_text(margin = margin(t = 20))) +
          ylim(0,1) +
          theme(
            legend.key.width = unit(1, "lines"),
            legend.spacing.x = unit(1, "lines"),
            legend.text = element_text(size = 12)
          ) ; ChiBar2_plot
        #
        #plot(X, ChiMix, xlab = "x", ylab = "Chi-bar-square(x)", main = "Chi-bar-square distribution")
        #lines(X, ChiMix)
        #abline(v=c2, col = "red")
        #abline(v=c2_compare, col = "blue")
        ##
        #legend("topright",
        #       legend = c("critical value", "observed value"),
        #       lty = 1,
        #       col = c("red", "blue")
        #)
      }


      diff_df = df_clpm - df_riclpm
      # TO DO RMK: Why do I check this at the end? If I do it before, I can also adjust the values based on the errors below...
      if(diff_df != (u + k)){ # NOT k*(k+1)/2 (= number of elements in block matrix 2,2, so neglecting block matrix 1,2; while u + k = nr of elt's of both those block matrices )
        #Could also calculate p-values based on diff_df, but that is weird, most probably df_diff consists of k and u (so, diff_df is not equal to k - it is equal to u+k)
        #p_diffdf = 0
        #for(teller in 1:(k+1)){
        #  p_diffdf = p_diffdf + weight[teller] * (1-pchisq(c2_compare, (u+teller-1))) #df_Chi2 = u until u+k
        #}
        message = paste0("The difference in degrees of freedom, ", diff_df, ", does not equal u + k = ", (u + k), ", please check your input. The rendered p-value is determined based on k = ", k, " and u = ", u, ".")
        #TO DO
      }else if(c2_compare < 0){
        message = paste0("The difference in Chi-square's (Chi2_clpm - Chi2_riclpm = ", c2_compare, ") is negative. The rendered p-value is determined based on switching them around.")
        # TO DO RMK: I do not actually do this right... I am also not sure whether I should do that (and does it matter?)
        # TO DO RMK: Btw, I can/should also do this when diff_df < 0, this makes more sense - then start with this check.
        # BUT: I do not use diff_df (or the df) in my calculations, so not helpful. I only use it as check on the plausability of the input here.
      }else{
        message = NULL
      }

      if(PrintPlot == T){
        final <- list(message = message,
                      k=k, u=u, S = S, ChiBar2_weights = weight,
                      critical_value = c2, DiffChi2 = c2_compare, CritValSmallerDiffChi2 = CritValSmallerDiffChi2, p_value=p, pSmallerAlpha = pSmallerAlpha,
                      ChiBar2_plot = ChiBar2_plot)
        print(ChiBar2_plot)
      }else{
        final <- list(message = message,
                      k=k, u=u, S = S, ChiBar2_weights = weight,
                      critical_value = c2, DiffChi2 = c2_compare, CritValSmallerDiffChi2 = CritValSmallerDiffChi2, p_value=p, pSmallerAlpha = pSmallerAlpha)
      }

    }
  }

  class(final) <- c("ChiBar2", "list")
  final

} # end of function

