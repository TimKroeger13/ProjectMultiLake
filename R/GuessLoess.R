#'GuessLoess
#'@description
#'Guess the best span for any function. Resolution for span = 0.001
#'@param intervall Intervall for the x-Achxes.
#'@param values values for they-achxes.
#'@param overspan Percent added to the minmal possible span.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger.
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

GuessLoess = function(intervall, values, overspan = 30){

  runningSpan = 0
  SpanNotFound = TRUE

  while (SpanNotFound) {

    runningSpan = runningSpan + 0.001

    loessMod = try(suppressWarnings(loess(values ~ intervall, span=runningSpan)),silent = TRUE)

    CatchedHotPotato <- tryCatch(loess(values ~ intervall, span=runningSpan), error = function(e) e, warning = function(w) w)

    if(class(loessMod)!="try-error"){
      if(!any(class(CatchedHotPotato) == "warning")){
        if(loessMod$s < Inf){

          SpanNotFound = FALSE

        }
      }
    }

    if (runningSpan>1){

      SpanNotFound = FALSE
      warning("Span could not converge")

    }
  }

  runningSpan = runningSpan + (runningSpan * overspan /100)

  if (runningSpan>1){

    runningSpan = 1

  }

  return(runningSpan)

}
