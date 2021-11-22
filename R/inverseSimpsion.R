#'inverseSimpsion
#'@description
#'Calculates the inverse Simpsion for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterInterpolating Value for the minimum rows after filtering.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

inverseSimpsion = function(data, intervallBy = 100, allLoessSpans = 0.8, minimumRowsAfterInterpolating = 0){

  deleteDoubles = function(doublData){

    counter = 0

    for (h in 1:(dim(doublData)[1]-1)){
    #for (h in 1:39){

      counter = counter+1

      if(as.numeric(row.names(doublData))[counter]==as.numeric(row.names(doublData))[counter+1]){

        RNames = row.names(doublData)

        doublData=matrix(doublData[-counter,])
        row.names(doublData) = RNames[-counter]

        counter = counter-1

      }
    }

    return(doublData)

  }

  DiatomNames = ls(data$Diatom)

  for (z in 1:length(DiatomNames)){

    inverseSimpsionData = data$Diatom[[DiatomNames[z]]]$Species_richness$invsimpson

    if(!is.null(inverseSimpsionData)){

      #Delete Doubles
      inverseSimpsionData = deleteDoubles(inverseSimpsionData)

      depthVectorOfData = as.numeric(row.names(inverseSimpsionData))

      lowerBoundry = ceiling(min(depthVectorOfData)/intervallBy)*intervallBy

      upperBoundry = floor(max(depthVectorOfData)/intervallBy)*intervallBy

      InterpolationMatrixRowNames =   approx (x = inverseSimpsionData[,1],
                                              y = NULL,
                                              xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                              method = "linear",
                                              n = 50)[[1]]

      InterpolationMatrix = matrix(NA, nrow = length(InterpolationMatrixRowNames), ncol = dim(inverseSimpsionData)[2])


      InterpolationMatrix[,1] =  approx (x = depthVectorOfData,
                                                       y = inverseSimpsionData,
                                                       xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                                       method = "linear",
                                                       n = 50)[[2]]

      colnames(InterpolationMatrix) = "InverseSimpsion"
      rownames(InterpolationMatrix) = InterpolationMatrixRowNames

      #check this minRows of the interpolated data

      if(dim(InterpolationMatrix)[1]>minimumRowsAfterInterpolating){

        InterpolationMatrixLoess = InterpolationMatrix
        InterpolationMatrixLoess[]=NA

        InterpolationMatrixLoess[,1] = predict(loess(InterpolationMatrix ~ InterpolationMatrixRowNames, span = allLoessSpans))


        InterpolationMatrixLoess[InterpolationMatrixLoess<0]=0

        inverseSimpsionMatrix = matrix(NA, ncol = 2, nrow = dim(InterpolationMatrixLoess)[1])

        inverseSimpsionMatrix[,1] = InterpolationMatrixRowNames
        inverseSimpsionMatrix[,2] = InterpolationMatrixLoess

      }

      data[["Diatom"]][[DiatomNames[z]]][["inverseSimpsion"]] = inverseSimpsionMatrix

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," calculating inverseSimpsion",sep="")

    }
  }

  return(data)

}
