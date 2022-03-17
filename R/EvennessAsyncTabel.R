#'EvennessAsyncTabel
#'@description
#'Calculates the evennes given by the evennes function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param Importname importname after data$Diatom$DiatomNames$
#'@param Exportname data$Diatom$DiatomNames$
#'@param CreateVector Shuld the vector plot build on this dataset.
#'@param TransformAllData Z Transformas all Data.
#'@param vectorName Outpurname for the vector data.
#'@param ValueCantBeSamlerThanZero If true, values smaller than 0 become 0.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

EvennessAsyncTabel = function(data, intervallBy = 100, Importname = "", Exportname = "", CreateVector = FALSE, TransformAllData = FALSE, vectorName,
                              ValueCantBeSamlerThanZero = FALSE){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating evenness Table ...",sep="")

  DiatomNames = ls(data$Diatom)

  MinEv = NULL
  MaxEv = NULL

  for (main in 1:length(DiatomNames)){

    evennessData = data$Diatom[[DiatomNames[main]]][[Importname]]

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

    evennessData = data$Diatom[[DiatomNames[main]]][[Importname]]

    if(!is.null(evennessData)){

      for (k in 1:dim(evennessData)[1]){

        EvennessAllInOneTabel[((evennessData[k,1] - min(MinEv)) / intervallBy)+1,main] = evennessData[k,2]

      }
    }
  }

  UntransformedTable = EvennessAllInOneTabel

  #TransformAllData

  if(TransformAllData){

    EvennessAllInOneTabel = scale(EvennessAllInOneTabel,center = T, scale = T)

  }

  #Create RateOfChangeAsyncTabel

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (i in 1:(((max(MaxEv)-min(MinEv))/intervallBy)+1)){

    if(sum(!is.na(EvennessAllInOneTabel[i,]))<4){

      AsyncTabel[i,2] = mean(as.numeric(na.omit(EvennessAllInOneTabel[i,])))
      AsyncTabel[i,3] = NA
      AsyncTabel[i,4] = NA
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = sum(!is.na(EvennessAllInOneTabel[i,]))

    }else{

      Evennesstest = t.test(EvennessAllInOneTabel[i,])

      AsyncTabel[i,2] = Evennesstest$estimate
      AsyncTabel[i,3] = Evennesstest$conf.int[1]
      AsyncTabel[i,4] = Evennesstest$conf.int[2]
      AsyncTabel[i,5] = Evennesstest$p.value
      AsyncTabel[i,6] = sum(!is.na(EvennessAllInOneTabel[i,]))

    }
  }

  if(CreateVector){

    #VectorMatrix

    VectorTable = DataSignalAfterTable(DataAllInOneTabel = EvennessAllInOneTabel, BasicAsyncTable = AsyncTabel, ValueCantBeSamlerThanZero = ValueCantBeSamlerThanZero,
                                       UntransformedTable = UntransformedTable)

    data[[vectorName]] = VectorTable

  }

  if (ValueCantBeSamlerThanZero){

    AsyncTabel[which(AsyncTabel[,2]<0),2] = 0
    AsyncTabel[which(AsyncTabel[,3]<0),3] = 0
    AsyncTabel[which(AsyncTabel[,4]<0),4] = 0

  }

  data[[Exportname]] = AsyncTabel

  return(data)

}







