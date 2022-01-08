#'CutOutPionierPhase
#'@description
#'Cut out the Pionier Phase of the Lake data.
#'\cr The new stating value is the first Value that passes the T.test that for two samples with a P_value <=0.05.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervals by to interpolate to.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterInterpolating Value for the minimum rows after filtering.
#'@param method The method for the dissimilarity calculation.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

CutOutPionierPhase = function(data,intervallBy = 100,allLoessSpans = 0.3, minimumRowsAfterInterpolating = 12,method = "bray"){

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

      depthVectorOfData = rateOfChangeData[,1]

      #Delete Doubles
      rateOfChangeData = deleteDoubles(rateOfChangeData)

      lowerBoundry = ceiling(min(depthVectorOfData)/intervallBy)*intervallBy

      upperBoundry = floor(max(depthVectorOfData)/intervallBy)*intervallBy

      dissimilarityMatrixRowNames =   approx (x = rateOfChangeData[,1],
                                              y = NULL,
                                              xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                              method = "linear",
                                              n = 50)[[1]]

      dissimilarityMatrix = matrix(NA, nrow = length(dissimilarityMatrixRowNames), ncol = dim(rateOfChangeData)[2]-3)

      matrixCounter = 0

      for (i in 4:dim(rateOfChangeData)[2]){

        matrixCounter = matrixCounter+1

        dissimilarityMatrix[,matrixCounter] =  approx (x = rateOfChangeData[,1],
                                                       y = rateOfChangeData[,i],
                                                       xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                                       method = "linear",
                                                       n = 50)[[2]]

      }

      colnames(dissimilarityMatrix) = colnames(rateOfChangeData)[4:dim(rateOfChangeData)[2]]
      rownames(dissimilarityMatrix) = dissimilarityMatrixRowNames

      #check this minRows of the interpolated data

      if(dim(dissimilarityMatrix)[1]>minimumRowsAfterInterpolating){

        dissimilarityMatrixLoess = dissimilarityMatrix
        dissimilarityMatrixLoess[]=NA

        for (i in 1:dim(dissimilarityMatrix)[2]){

          dissimilarityMatrixLoess[,i] = predict(loess(dissimilarityMatrix[,i] ~ dissimilarityMatrixRowNames, span = allLoessSpans))

        }

        dissimilarityMatrixLoess[dissimilarityMatrixLoess<0]=0

        #Sqrt transform
        dissimilarityMatrixLoess = sqrt(dissimilarityMatrixLoess)

        #Distances

        DistanceMatrix = matrix(NA, ncol = 2, nrow = (dim(dissimilarityMatrixLoess)[1]-1))

        for (p in 1:(dim(dissimilarityMatrixLoess)[1]-1)){

          distdata = vegdist(dissimilarityMatrixLoess[p:(p+1),],method=method,na.rm = T)

          DistanceMatrix[p,2] = distdata
          DistanceMatrix[p,1] = as.numeric(rownames(dissimilarityMatrixLoess)[p])

        }
      }

      #Second Loess for best distribution to Cut. To make the Pionier Phase visible

      ROC_Loess_after_interplating = 0.4 #<- Fix value. Seems to be the best fit after some testing for 0-20.000 years

      #Loess
      ValuesLoess=loess(DistanceMatrix[,2] ~ DistanceMatrix[,1], span=ROC_Loess_after_interplating)
      DistanceMatrixLoess = cbind(DistanceMatrix[,1],ValuesLoess$fitted)

      #Check for Pionier Phase

      p_value = 0.05

      SpliValue = dim(DistanceMatrixLoess)[1]
      SplitStart = DistanceMatrixLoess[(SpliValue-1):SpliValue,2]
      SplitEnd = DistanceMatrixLoess[1:(SpliValue-1),2]
      SplitTTest = t.test(SplitStart,SplitEnd)

      #correction
      #Find Pionier Phase in the interpolation

      if(SplitTTest$p.value > p_value){

        SplitvalueNotFound = TRUE

        while (SplitvalueNotFound){

          SpliValue = SpliValue-1
          SplitStart = DistanceMatrixLoess[(SpliValue-1):SpliValue,2]
          SplitEnd = DistanceMatrixLoess[1:(SpliValue-1),2]

          SplitTTest = t.test(SplitStart,SplitEnd)

          if(SplitTTest$p.value <= p_value){

            SplitvalueNotFound = FALSE

          }
        }
      }

      #Find the next value of the raw data

      PionierLoess = DistanceMatrixLoess[SpliValue,1]

      PionierFreeSRS = data$Diatom[[DiatomNames[z]]]$SRS_data

      PionierFreeSRS = PionierFreeSRS[PionierFreeSRS[,1] <= PionierLoess,]

      data$Diatom[[DiatomNames[z]]][["Cut_SRS_data"]] = PionierFreeSRS

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," Cut SRS Data",sep="")
    }
  }

  return(data)

}
































