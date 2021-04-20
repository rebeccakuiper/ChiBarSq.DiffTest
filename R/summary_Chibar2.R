#' @importFrom jtools md_table
#' @S3method summary ChiBar2
#' @export summary.ChiBar2
#' @export


summary.ChiBar2 <- function(x, digits = NULL)
{
  #library(jtools)

  x <- as.list(x)

  if(is.null(digits)){
    NrDigits <- 3
    sig.digits <- TRUE
    align <- NULL
  }else{
    NrDigits <- digits
    sig.digits <- FALSE
    align <- 'c'
  }

  if(is.null(x$critical_value)) {
    cat("\n")
    cat(x$message)
    cat("\n")
  }else{
    w <- x$ChiBar2_weights
    names(w) <- 0:(length(w)-1)
    DF1 <- data.frame(ChiBar2_weights = w)
    if(is.null(x$DiffChi2)) {
      DF2 <- data.frame(critical_value = x$critical_value)
    }else{
      DF2 <- data.frame(critical_value = x$critical_value,
                        DiffChi2 = x$DiffChi2,
                        p_value = x$p_value)
    }
    #
    cat("\n")
    cat("The Chi-bar-square difference test\n")
    print(md_table(DF2, digits = NrDigits, sig.digits = sig.digits, align = align, row.names = FALSE))
    cat("\n")
    cat("The Chi-bar-square weights used in the test\n")
    print(md_table(DF1, digits = NrDigits, sig.digits = sig.digits, align = align))
    cat("\n")
    #
    if(!is.null(x$message)){
      cat(paste0("Note: ", x$message))
      cat("\n")
      cat("\n")
    }
  }
}
