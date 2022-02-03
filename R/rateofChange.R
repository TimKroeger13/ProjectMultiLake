#'rateofChange
#'@description
#'Calculates the Rate of change for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterInterpolating alue for the minimum rows after filtering.
#'@param method Method for calculation Dissimilarity Indices for Community Ecologists.
#'@param Importname importname after data$Diatom$DiatomNames$
#'@param Exportname data$Diatom$DiatomNames$
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

rateofChange = function(data, intervallBy = 100, allLoessSpans = 0.8, minimumRowsAfterInterpolating = 0, method = "bray", Importname = "", Exportname = ""){

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

    rateOfChangeData = data$Diatom[[DiatomNames[z]]][[Importname]]

    if(!is.null(rateOfChangeData)){

      #Delete Doubles
      rateOfChangeData = deleteDoubles(rateOfChangeData)

      depthVectorOfData = rateOfChangeData[,1]


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

        #Sqrt transform
        # Got cut. Because of the interpolation values can be between 0 and 1.
        # This causes some values to grow instead of shrinking.

        #Distances

        DistanceMatrix = matrix(NA, ncol = 2, nrow = (dim(dissimilarityMatrix)[1]-1))

        for (p in 1:(dim(dissimilarityMatrix)[1]-1)){

          distdata = vegdist(dissimilarityMatrix[p:(p+1),],method=method,na.rm = T) / intervallBy

          DistanceMatrix[p,2] = distdata
          DistanceMatrix[p,1] = (as.numeric(rownames(dissimilarityMatrix)[p]) + as.numeric(rownames(dissimilarityMatrix)[p+1])) /2

        }
      }

      data[["Diatom"]][[DiatomNames[z]]][[Exportname]] = DistanceMatrix

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," calculating ",Importname," of change",sep="")

    }
  }

  return(data)

}






