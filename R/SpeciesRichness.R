#'speciesRichness
#'@description
#'Calculates the inverse species richness for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterInterpolating Value for the minimum rows after filtering.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

speciesRichness = function(data, intervallBy = 100, allLoessSpans = 0.8, minimumRowsAfterInterpolating = 0){

  deleteDoubles = function(doublData){

    counter = 0

    for (h in 1:(dim(doublData)[1]-1)){

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

    speciesRichnessData = data$Diatom[[DiatomNames[z]]]$Species_richness$richness

    if(!is.null(speciesRichnessData)){

      #Delete Doubles
      speciesRichnessData = deleteDoubles(speciesRichnessData)

      depthVectorOfData = as.numeric(row.names(speciesRichnessData))

      lowerBoundry = ceiling(min(depthVectorOfData)/intervallBy)*intervallBy

      upperBoundry = floor(max(depthVectorOfData)/intervallBy)*intervallBy

      InterpolationMatrixRowNames =   approx (x = speciesRichnessData[,1],
                                              y = NULL,
                                              xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                              method = "linear",
                                              n = 50)[[1]]

      InterpolationMatrix = matrix(NA, nrow = length(InterpolationMatrixRowNames), ncol = dim(speciesRichnessData)[2])


      InterpolationMatrix[,1] =  approx (x = depthVectorOfData,
                                         y = speciesRichnessData,
                                         xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                         method = "linear",
                                         n = 50)[[2]]

      colnames(InterpolationMatrix) = "SpeciesRichness"
      rownames(InterpolationMatrix) = InterpolationMatrixRowNames

      #check this minRows of the interpolated data

      if(dim(InterpolationMatrix)[1]>minimumRowsAfterInterpolating){

        InterpolationMatrixLoess = InterpolationMatrix
        InterpolationMatrixLoess[]=NA

        InterpolationMatrixLoess[,1] = predict(loess(InterpolationMatrix ~ InterpolationMatrixRowNames, span = allLoessSpans))


        InterpolationMatrixLoess[InterpolationMatrixLoess<0]=0

        speciesRichnessMatrix = matrix(NA, ncol = 2, nrow = dim(InterpolationMatrixLoess)[1])

        speciesRichnessMatrix[,1] = InterpolationMatrixRowNames
        speciesRichnessMatrix[,2] = InterpolationMatrixLoess

      }

      data[["Diatom"]][[DiatomNames[z]]][["speciesRichness"]] = speciesRichnessMatrix

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," calculating SpeciesRichness",sep="")

    }
  }

  return(data)

}
