#'CutOutPionierPhase
#'@description
#'Cut out the Pionier Phase of the Lake data.
#'\cr The new stating value is the first Value that passes the T.test that for two samples with a P_value <=0.05.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervals by to interpolate to.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterInterpolating Value for the minimum rows after filtering.
#'@param method The method for the dissimilarity calculation.
#'@param AgeBuffer How much age get Cut away.
#'@param AgeBuffer2 How much age get Cut away as the second value.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

CutOutPionierPhase = function(data,intervallBy = 100,allLoessSpans = 0.3, minimumRowsAfterInterpolating = 12,method = "bray", AgeBuffer = 1000, AgeBuffer2 = 500){

  deleteDoubles = function(doublData){

    counter = 0

    for (h in 1:(dim(doublData)[1]-1)){

      counter = counter+1

      if(doublData[counter,1]==doublData[(counter+1),1]){

        doublData=doublData[-counter,]

        counter = counter-1

      }
    }

    return(doublData)

  }

  DiatomNames = ls(data$Diatom)

  for (z in 1:length(DiatomNames)){

    rateOfChangeData = data$Diatom[[DiatomNames[z]]]$SRS_data

    if(!is.null(rateOfChangeData)){

      TempDatarateOfChangeData = rateOfChangeData[,4:dim(rateOfChangeData)[2]]

      RateofChangeAge = rateOfChangeData[,1]

      DistanceMatrix = matrix(NA, ncol = 2, nrow = (dim(TempDatarateOfChangeData)[1]-1))

      for (p in 1:(dim(TempDatarateOfChangeData)[1]-1)){

        distdata = vegdist(TempDatarateOfChangeData[p:(p+1),],method=method,na.rm = T)

        DistanceMatrix[p,2] = distdata
        DistanceMatrix[p,1] = rateOfChangeData[,1][p]

      }

      PionierAge = DistanceMatrix[dim(DistanceMatrix)[1],1]-AgeBuffer

      SplitStart = DistanceMatrix[DistanceMatrix[,1] <= PionierAge,]
      SplitEnd =  DistanceMatrix[DistanceMatrix[,1] > PionierAge,]

      SplitTTest = t.test(SplitStart,SplitEnd)

      p_value = 0.05

      PionierFreeSRS = data$Diatom[[DiatomNames[z]]]$SRS_data

      if(SplitTTest$p.value > p_value){

        PionierFreeSRS = PionierFreeSRS[PionierFreeSRS[,1] <= PionierAge,]

      }else{

        PionierAge = DistanceMatrix[dim(DistanceMatrix)[1],1]-AgeBuffer2

        SplitStart = DistanceMatrix[DistanceMatrix[,1] <= PionierAge,]
        SplitEnd =  DistanceMatrix[DistanceMatrix[,1] > PionierAge,]

        SplitTTest = t.test(SplitStart,SplitEnd)

        p_value = 0.05

        if(SplitTTest$p.value > p_value){

          PionierFreeSRS = PionierFreeSRS[PionierFreeSRS[,1] <= PionierAge,]

        }
      }

      data$Diatom[[DiatomNames[z]]][["Cut_SRS_data"]] = PionierFreeSRS

    }

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," Cut SRS Data",sep="")

  }

  return(data)

}


