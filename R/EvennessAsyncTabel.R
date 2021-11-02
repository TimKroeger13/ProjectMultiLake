#'EvennessAsyncTabel
#'@description
#'Calculates the evennes given by the evennes function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

EvennessAsyncTabel = function(data, intervallBy = 100){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating evenness Table ...",sep="")

  DiatomNames = ls(data$Diatom)

  MinEv = NULL
  MaxEv = NULL

  for (main in 1:length(DiatomNames)){

    evennessData = data$Diatom[[DiatomNames[main]]][["evenness"]]

    if(!is.null(evennessData)){

    MinEv = c(MinEv,min(evennessData[,1]))
    MaxEv = c(MaxEv,max(evennessData[,1]))

    }
  }

  #Create EvennessAllInOneTabel

  EvennessAllInOneTabel = matrix(NA, ncol = length(DiatomNames), nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(EvennessAllInOneTabel) = DiatomNames
  rownames(EvennessAllInOneTabel) = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (main in 1:length(DiatomNames)){

    evennessData = data$Diatom[[DiatomNames[main]]][["evenness"]]

    if(!is.null(evennessData)){

      for (k in 1:dim(evennessData)[1]){

        EvennessAllInOneTabel[((evennessData[k,1] - min(MinEv)) / intervallBy)+1,main] = evennessData[k,2]

      }
    }
  }

  #Create RateOfChangeAsyncTabel

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (i in 1:(((max(MaxEv)-min(MinEv))/intervallBy)+1)){

    if(sum(!is.na(EvennessAllInOneTabel[i,]))==1){

      AsyncTabel[i,2] = as.numeric(na.omit(EvennessAllInOneTabel[i,]))
      AsyncTabel[i,3] = as.numeric(na.omit(EvennessAllInOneTabel[i,]))
      AsyncTabel[i,4] = as.numeric(na.omit(EvennessAllInOneTabel[i,]))
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = 1

    }else{

      Evennesstest = t.test(EvennessAllInOneTabel[i,])

      AsyncTabel[i,2] = Evennesstest$estimate
      AsyncTabel[i,3] = Evennesstest$conf.int[1]
      AsyncTabel[i,4] = Evennesstest$conf.int[2]
      AsyncTabel[i,5] = Evennesstest$p.value
      AsyncTabel[i,6] = sum(!is.na(EvennessAllInOneTabel[i,]))

    }
  }

  AsyncTabel[which(AsyncTabel[,3]<0),3] = 0

  data[["EvennessMatrix"]] = AsyncTabel

  return(data)

}







