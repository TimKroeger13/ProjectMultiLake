#'TRACETabel
#'@description
#'Calculates the TRACE as a Table.
#'@param data List of data generates by the Multivar function.
#'@param TraceName Name of type of Trace data. For the search in the given list.
#'@param intervallBy Intervals by to mean to.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

TRACETabel = function(data, TraceName = "UnknownTrace", intervallBy = 100){

  TRACENameList = ls(data$TRACE)

  for (main in 1:length(TRACENameList)){

    TRACEData = data$TRACE[[TRACENameList[main]]][[paste(TraceName,"_raw",sep="")]]

    lowerBoundry = ceiling(min(TRACEData[,1])/intervallBy)*intervallBy

    upperBoundry = floor(max(TRACEData[,1])/intervallBy)*intervallBy

    MeanSequnce = seq(from = lowerBoundry+(intervallBy/2), to = upperBoundry-(intervallBy/2), by = intervallBy)
    LowerSequnce = seq(from = lowerBoundry, to = upperBoundry-intervallBy, by = intervallBy)
    UpperSequnce = seq(from = lowerBoundry+intervallBy, to = upperBoundry, by = intervallBy)

    TRACEMeanData = matrix(NA, ncol = 2, nrow = length(MeanSequnce))

    for(k in 1:length(MeanSequnce)){

      LowerSequnce[k]

      UpperSequnce[k]

      TRACEMeanData[k,2] = mean(TRACEData[which(((TRACEData[,1]>= LowerSequnce[k]) + (TRACEData[,1]< UpperSequnce[k])) == 2),2])

    }

    TRACEMeanData[,1] = MeanSequnce

    data[["TRACE"]][[TRACENameList[main]]][[paste(TraceName,"_mean",sep="")]]=TRACEMeanData

    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        main,"/",length(TRACENameList)," Calculating mean for TRACE data: ",TraceName,sep="")

  }




  #Table

  TraceName = paste(TraceName,"_mean",sep="")

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating TRACE Table ",TraceName," ...",sep="")

  MinEv = NULL
  MaxEv = NULL

  TraceRowNames = data$TRACE[[TRACENameList[1]]][[TraceName]][,1]

  for (main in 1:length(TRACENameList)){

    TRACEData = data$TRACE[[TRACENameList[main]]][[TraceName]]

    if(!is.null(TRACEData)){

      MinEv = c(MinEv,min(TRACEData[,1]))
      MaxEv = c(MaxEv,max(TRACEData[,1]))

    }
  }

  #Create TRACEAllInOneTabel

  TRACEAllInOneTabel = matrix(NA, ncol = length(TRACENameList), nrow =  dim(TRACEData)[1])
  colnames(TRACEAllInOneTabel) = TRACENameList
  rownames(TRACEAllInOneTabel) = row.names(TraceRowNames)

  for (main in 1:length(TRACENameList)){

    TRACEData = data$TRACE[[TRACENameList[main]]][[TraceName]]

    if(!is.null(TRACEData)){

      for (k in 1:dim(TRACEData)[1]){

        TRACEAllInOneTabel[k,main] = TRACEData[k,2]

      }
    }

    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        main,"/",length(TRACENameList)," Calculating TRACE Table: ",TraceName,sep="")

  }

  #Create AsyncTabel

  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating TRACE Table as Async Mean: ",TraceName," ...",sep="")

  AsyncTabel = matrix(NA, ncol = 6, nrow = dim(TRACEData)[1])
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = TraceRowNames

  for (i in 1:dim(TRACEData)[1]){

    if(sum(!is.na(TRACEAllInOneTabel[i,]))==1){

      AsyncTabel[i,2] = as.numeric(na.omit(TRACEAllInOneTabel[i,]))
      AsyncTabel[i,3] = as.numeric(na.omit(TRACEAllInOneTabel[i,]))
      AsyncTabel[i,4] = as.numeric(na.omit(TRACEAllInOneTabel[i,]))
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = 1

    }else{

      TRACEtest = t.test(TRACEAllInOneTabel[i,])

      AsyncTabel[i,2] = TRACEtest$estimate
      AsyncTabel[i,3] = TRACEtest$conf.int[1]
      AsyncTabel[i,4] = TRACEtest$conf.int[2]
      AsyncTabel[i,5] = TRACEtest$p.value
      AsyncTabel[i,6] = sum(!is.na(TRACEAllInOneTabel[i,]))

    }

    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",dim(TRACEData)[1]," Calculating TRACE Table as Async Mean: ",TraceName,sep="")

  }

  #AsyncTabel[which(AsyncTabel[,3]<0),3] = 0 #<---------------------------------- Not for TRACE

  data[["TRACEMatrix"]][[TraceName]] = AsyncTabel

  return(data)

}
