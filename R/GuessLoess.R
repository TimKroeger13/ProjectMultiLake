#'GuessLoess
#'@description
#'Guess the best span for any function. Resolution for span = 0.01
#'@param intervall Intervall for the x-Achxes.
#'@param values values for they-achxes.
#'@param overspan Percent added to the minmal possible span.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger.
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

GuessLoess = function(intervall, values, overspan = 20){

  runningSpan = 0
  SpanNotFound = TRUE

  while (SpanNotFound) {

    runningSpan = runningSpan + 0.01

    loessMod = try(suppressWarnings(loess(values ~ intervall, span=runningSpan)),silent = TRUE)

    if(class(loessMod)!="try-error"){
      if(sum(is.na(loessMod$s)) == 0){

        SpanNotFound = FALSE

      }
    }

    if (runningSpan>1){

      SpanNotFound = FALSE
      warning("Span could not converge")

    }
  }

  return(runningSpan + (runningSpan * overspan /100))

}
