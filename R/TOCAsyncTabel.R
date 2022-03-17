#'TOCAsyncTabel
#'@description
#'Calculates the TOC as a Table given by the TOC function.
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

TOCAsyncTabel = function(data, intervallBy = 100, Importname = "", Exportname = "", CreateVector = FALSE, TransformAllData = FALSE, vectorName,
                         ValueCantBeSamlerThanZero = FALSE){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating TOC Table ...",sep="")

  TOCNames = ls(data$Carbon)

  MinEv = NULL
  MaxEv = NULL

  for (main in 1:length(TOCNames)){

    TOCData = data$Carbon[[TOCNames[main]]][[Importname]]

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

    TOCData = data$Carbon[[TOCNames[main]]][[Importname]]

    if(!is.null(TOCData)){

      for (k in 1:dim(TOCData)[1]){

        TOCAllInOneTabel[((TOCData[k,1] - min(MinEv)) / intervallBy)+1,main] = TOCData[k,2]

      }
    }
  }

  UntransformedTable = TOCAllInOneTabel

  #TransformAllData

  if(TransformAllData){

    TOCAllInOneTabel = scale(TOCAllInOneTabel,center = T, scale = T)

  }

  #Async Table

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (i in 1:(((max(MaxEv)-min(MinEv))/intervallBy)+1)){

    if(sum(!is.na(TOCAllInOneTabel[i,]))<4){

      AsyncTabel[i,2] = mean(as.numeric(na.omit(TOCAllInOneTabel[i,])))
      AsyncTabel[i,3] = NA
      AsyncTabel[i,4] = NA
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = sum(!is.na(TOCAllInOneTabel[i,]))

    }else{

      AsyncTest = t.test(TOCAllInOneTabel[i,])

      AsyncTabel[i,2] = AsyncTest$estimate
      AsyncTabel[i,3] = AsyncTest$conf.int[1]
      AsyncTabel[i,4] = AsyncTest$conf.int[2]
      AsyncTabel[i,5] = AsyncTest$p.value
      AsyncTabel[i,6] = sum(!is.na(TOCAllInOneTabel[i,]))

    }
  }

  if(CreateVector){

    #VectorMatrix

    VectorTable = DataSignalAfterTable(DataAllInOneTabel = TOCAllInOneTabel, BasicAsyncTable = AsyncTabel, ValueCantBeSamlerThanZero = ValueCantBeSamlerThanZero,
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
