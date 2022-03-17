#'speciesRichness
#'@description
#'Calculates the inverse species richness for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param NonNegative Creates Positive Values after Loess calculation.
#'@param Importname1 importname 1.
#'@param Importname2 importname 2.
#'@param Exportname data$Diatom$DiatomNames$
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

speciesRichness = function(data, intervallBy = 100, NonNegative = TRUE, Importname1 = "", Importname2 = "", Exportname = ""){

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

    speciesRichnessData = data$Diatom[[DiatomNames[z]]][[Importname1]][[Importname2]]

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

      colnames(InterpolationMatrix) = Importname1
      rownames(InterpolationMatrix) = InterpolationMatrixRowNames

      #check this minRows of the interpolated data

        InterpolationMatrixLoess = InterpolationMatrix
        InterpolationMatrixLoess[]=NA

        #Guess Loess

        span = GuessLoess(intervall = InterpolationMatrixRowNames, values = InterpolationMatrix, overspan = 100)

        InterpolationMatrixLoess[,1] = predict(loess(InterpolationMatrix ~ InterpolationMatrixRowNames, span = span))


        InterpolationMatrixLoess[InterpolationMatrixLoess<0]=0

        speciesRichnessMatrix = matrix(NA, ncol = 2, nrow = dim(InterpolationMatrixLoess)[1])

        speciesRichnessMatrix[,1] = InterpolationMatrixRowNames
        speciesRichnessMatrix[,2] = InterpolationMatrixLoess

      data[["Diatom"]][[DiatomNames[z]]][[Exportname]] = speciesRichnessMatrix

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," calculating ",Importname2,sep="")

    }
  }

  return(data)

}
