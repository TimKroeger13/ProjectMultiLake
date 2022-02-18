#'RateOfChangeAsyncTabel
#'@description
#'Calculates the Rate of change given by the rateofChange function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param Importname importname after data$Diatom$DiatomNames$
#'@param Exportname data$Diatom$DiatomNames$
#'@param CreateVector Shuld the vector plot build on this dataset.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

RateOfChangeAsyncTabel = function(data, intervallBy = 100, Importname = "", Exportname = "", CreateVector = FALSE){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating rate of change Table ...",sep="")


  DiatomNames = ls(data$Diatom)

  MinRoC = NULL
  MaxRoC = NULL

  for (main in 1:length(DiatomNames)){

    RocData = data$Diatom[[DiatomNames[main]]][[Importname]]

    if(!is.null(RocData)){

      MinRoC = c(MinRoC,min(RocData[,1]))
      MaxRoC = c(MaxRoC,max(RocData[,1]))

    }
  }

  #Create RocAllInOneTabel

  RocAllInOneTabel = matrix(NA, ncol = length(DiatomNames), nrow = ((max(MaxRoC)-min(MinRoC))/intervallBy)+1)
  colnames(RocAllInOneTabel) = DiatomNames
  rownames(RocAllInOneTabel) = seq(from = min(MinRoC), to = max(MaxRoC), by = intervallBy)

  for (main in 1:length(DiatomNames)){

    RocData = data$Diatom[[DiatomNames[main]]][[Importname]]

    if(!is.null(RocData)){

      for (k in 1:dim(RocData)[1]){

        RocAllInOneTabel[((RocData[k,1] - min(MinRoC)) / intervallBy)+1,main] = RocData[k,2]

      }
    }
  }

  #Create RateOfChangeAsyncTabel

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxRoC)-min(MinRoC))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinRoC), to = max(MaxRoC), by = intervallBy)

  for (i in 1:(((max(MaxRoC)-min(MinRoC))/intervallBy)+1)){

    if(sum(!is.na(RocAllInOneTabel[i,]))==1){

      AsyncTabel[i,2] = as.numeric(na.omit(RocAllInOneTabel[i,]))
      AsyncTabel[i,3] = as.numeric(na.omit(RocAllInOneTabel[i,]))
      AsyncTabel[i,4] = as.numeric(na.omit(RocAllInOneTabel[i,]))
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = 1

    }else{

      RoCTtest = t.test(RocAllInOneTabel[i,])

      AsyncTabel[i,2] = RoCTtest$estimate
      AsyncTabel[i,3] = RoCTtest$conf.int[1]
      AsyncTabel[i,4] = RoCTtest$conf.int[2]
      AsyncTabel[i,5] = RoCTtest$p.value
      AsyncTabel[i,6] = sum(!is.na(RocAllInOneTabel[i,]))

    }
  }

  if(CreateVector){

    #VectorMatrix

    VectorTable = DataSignalAfterTable(DataAllInOneTabel = RocAllInOneTabel,BasicAsyncTable = AsyncTabel,ValueCantBeSamlerThanZero = TRUE)

    data[["Vector_RocMatrix"]] = VectorTable

  }

  AsyncTabel[which(AsyncTabel[,3]<0),3] = 0

  data[[Exportname]] = AsyncTabel

  return(data)

}


















