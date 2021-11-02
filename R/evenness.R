#'evenness
#'@description
#'Calculates the evenness for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterInterpolating alue for the minimum rows after filtering.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

evenness = function(data, intervallBy = 100, allLoessSpans = 0.8, minimumRowsAfterInterpolating = 0){

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

    evennessData = data$Diatom[[DiatomNames[z]]]$SRS_data

    if(!is.null(evennessData)){

      depthVectorOfData = evennessData[,1]

      #Delete Doubles
      evennessData = deleteDoubles(evennessData)

      lowerBoundry = ceiling(min(depthVectorOfData)/intervallBy)*intervallBy

      upperBoundry = floor(max(depthVectorOfData)/intervallBy)*intervallBy

      InterpolationMatrixRowNames =   approx (x = evennessData[,1],
                                              y = NULL,
                                              xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                              method = "linear",
                                              n = 50)[[1]]

      InterpolationMatrix = matrix(NA, nrow = length(InterpolationMatrixRowNames), ncol = dim(evennessData)[2]-3)

      matrixCounter = 0

      for (i in 4:dim(evennessData)[2]){

        matrixCounter = matrixCounter+1

        InterpolationMatrix[,matrixCounter] =  approx (x = evennessData[,1],
                                                       y = evennessData[,i],
                                                       xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                                       method = "linear",
                                                       n = 50)[[2]]

      }

      colnames(InterpolationMatrix) = colnames(evennessData)[4:dim(evennessData)[2]]
      rownames(InterpolationMatrix) = InterpolationMatrixRowNames



    #check this minRows of the interpolated data

    if(dim(InterpolationMatrix)[1]>minimumRowsAfterInterpolating){

      InterpolationMatrixLoess = InterpolationMatrix
      InterpolationMatrixLoess[]=NA

      for (i in 1:dim(InterpolationMatrix)[2]){

        InterpolationMatrixLoess[,i] = predict(loess(InterpolationMatrix[,i] ~ InterpolationMatrixRowNames, span = allLoessSpans))

      }

      InterpolationMatrixLoess[InterpolationMatrixLoess<0]=0

      #Sqrt transform
      #InterpolationMatrixLoess = sqrt(InterpolationMatrixLoess)

      #evenness

      evennessMatrix = matrix(NA, ncol = 2, nrow = dim(InterpolationMatrixLoess)[1])

      for (p in 1:dim(InterpolationMatrixLoess)[1]){

        evennessVector = InterpolationMatrixLoess[p,]

        evennessMatrix[p,1] = as.numeric(rownames(InterpolationMatrixLoess)[p])
        evennessMatrix[p,2] = diversity(evennessVector)/log(specnumber(evennessVector))

      }
    }

      data[["Diatom"]][[DiatomNames[z]]][["evenness"]] = evennessMatrix

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," calculating evenness",sep="")


    }
  }

  return(data)

}
