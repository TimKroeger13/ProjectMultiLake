#'TOCAsyncTabel
#'@description
#'Calculates the TOC as a Table given by the TOC function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

TOCAsyncTabel = function(data, intervallBy = 100){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating TOC Table ...",sep="")

  TOCNames = ls(data$Carbon)

  MinEv = NULL
  MaxEv = NULL

  for (main in 1:length(TOCNames)){

    TOCData = data$Carbon[[TOCNames[main]]][["TOC"]]

    if(!is.null(TOCData)){

      MinEv = c(MinEv,min(TOCData[,1]))
      MaxEv = c(MaxEv,max(TOCData[,1]))

    }
  }

  #Create TOCAllInOneTabel

  TOCAllInOneTabel = matrix(NA, ncol = length(TOCNames), nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(TOCAllInOneTabel) = TOCNames
  rownames(TOCAllInOneTabel) = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (main in 1:length(TOCNames)){

    TOCData = data$Carbon[[TOCNames[main]]][["TOC"]]

    if(!is.null(TOCData)){

      for (k in 1:dim(TOCData)[1]){

        TOCAllInOneTabel[((TOCData[k,1] - min(MinEv)) / intervallBy)+1,main] = TOCData[k,2]

      }
    }
  }

  #Create RTOCAsyncTabel

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (i in 1:(((max(MaxEv)-min(MinEv))/intervallBy)+1)){

    if(sum(!is.na(TOCAllInOneTabel[i,]))==1){

      AsyncTabel[i,2] = as.numeric(na.omit(TOCAllInOneTabel[i,]))
      AsyncTabel[i,3] = as.numeric(na.omit(TOCAllInOneTabel[i,]))
      AsyncTabel[i,4] = as.numeric(na.omit(TOCAllInOneTabel[i,]))
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = 1

    }else{

      TOCtest = t.test(TOCAllInOneTabel[i,])

      AsyncTabel[i,2] = TOCtest$estimate
      AsyncTabel[i,3] = TOCtest$conf.int[1]
      AsyncTabel[i,4] = TOCtest$conf.int[2]
      AsyncTabel[i,5] = TOCtest$p.value
      AsyncTabel[i,6] = sum(!is.na(TOCAllInOneTabel[i,]))

    }
  }

  #AsyncTabel[which(AsyncTabel[,3]<0),3] = 0 #<---------------------------------- Not for TOC

  data[["TOCMatrix"]] = AsyncTabel

  return(data)

}
