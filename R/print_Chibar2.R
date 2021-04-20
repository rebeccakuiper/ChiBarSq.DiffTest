#' @S3method print ChiBar2
#' @export print.ChiBar2
#' @export

print.ChiBar2 <- function(x, digits = NULL)
{

  x <- as.list(x)

  if(is.null(digits)){
    NrDigits <- 3
  }else{
    NrDigits <- digits
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
      DF2 <- data.frame(critical_value = x$critical_value,
                        row.names = "")
    }else{
      DF2 <- data.frame(critical_value = x$critical_value,
                        DiffChi2 = x$DiffChi2,
                        p_value = x$p_value,
                        row.names = "")
    }
    #
    #
    cat("\n")
    cat("The Chi-bar-square test results\n")
    #cat("\n")
    print(DF2, digits = NrDigits, right = F)
    #cat("\n")
    cat("\n")
    cat("The Chi-bar-square weights used in the test\n")
    #cat("\n")
    print(DF1, digits = NrDigits, right = F)
    cat("\n")
    #
    #
    if(!is.null(x$message)){
      #cat("\n")
      cat(paste0("Note: ", x$message))
      cat("\n")
      cat("\n")
    }
  }

  return(invisible(x))

}

