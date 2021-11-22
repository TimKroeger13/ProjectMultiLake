#'InverseSimpsionAsyncTabel
#'@description
#'Calculates the inverse simpsion given by the inverseSimpsion function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

InverseSimpsionAsyncTabel = function(data, intervallBy = 100){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating inverse Simpsion Table ...",sep="")

  DiatomNames = ls(data$Diatom)

  MinEv = NULL
  MaxEv = NULL

  for (main in 1:length(DiatomNames)){

    inverseSimpsionData = data$Diatom[[DiatomNames[main]]][["inverseSimpsion"]]

    if(!is.null(inverseSimpsionData)){

      MinEv = c(MinEv,min(inverseSimpsionData[,1]))
      MaxEv = c(MaxEv,max(inverseSimpsionData[,1]))

    }
  }

  #Create inverseSimpsionAllInOneTabel

  inverseSimpsionAllInOneTabel = matrix(NA, ncol = length(DiatomNames), nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(inverseSimpsionAllInOneTabel) = DiatomNames
  rownames(inverseSimpsionAllInOneTabel) = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (main in 1:length(DiatomNames)){

    inverseSimpsionData = data$Diatom[[DiatomNames[main]]][["inverseSimpsion"]]

    if(!is.null(inverseSimpsionData)){

      for (k in 1:dim(inverseSimpsionData)[1]){

        inverseSimpsionAllInOneTabel[((inverseSimpsionData[k,1] - min(MinEv)) / intervallBy)+1,main] = inverseSimpsionData[k,2]

      }
    }
  }

  #Create RinverseSimpsionAsyncTabel

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (i in 1:(((max(MaxEv)-min(MinEv))/intervallBy)+1)){

    if(sum(!is.na(inverseSimpsionAllInOneTabel[i,]))==1){

      AsyncTabel[i,2] = as.numeric(na.omit(inverseSimpsionAllInOneTabel[i,]))
      AsyncTabel[i,3] = as.numeric(na.omit(inverseSimpsionAllInOneTabel[i,]))
      AsyncTabel[i,4] = as.numeric(na.omit(inverseSimpsionAllInOneTabel[i,]))
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = 1

    }else{

      inverseSimpsiontest = t.test(inverseSimpsionAllInOneTabel[i,])

      AsyncTabel[i,2] = inverseSimpsiontest$estimate
      AsyncTabel[i,3] = inverseSimpsiontest$conf.int[1]
      AsyncTabel[i,4] = inverseSimpsiontest$conf.int[2]
      AsyncTabel[i,5] = inverseSimpsiontest$p.value
      AsyncTabel[i,6] = sum(!is.na(inverseSimpsionAllInOneTabel[i,]))

    }
  }

  AsyncTabel[which(AsyncTabel[,3]<0),3] = 0

  data[["InverseSimpsionMatrix"]] = AsyncTabel

  return(data)

}







