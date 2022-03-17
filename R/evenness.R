#'evenness
#'@description
#'Calculates the evenness for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param NonNegative Creates Positive Values after Loess calculation.
#'@param Importname importname after data$Diatom$DiatomNames$
#'@param Exportname data$Diatom$DiatomNames$
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

evenness = function(data, intervallBy = 100, NonNegative = TRUE, Importname = "", Exportname = ""){

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

    evennessData = data$Diatom[[DiatomNames[z]]][[Importname]]

    if(!is.null(evennessData)){
      if(dim(evennessData)[1]>2){

        #Delete Doubles
        evennessData = deleteDoubles(evennessData)

        depthVectorOfData = evennessData[,1]

        #Calculation Data
        evennessData_calculation = evennessData[4:dim(evennessData)[2]]

        #Sqrt transformation

        evennessData_calculation = sqrt(evennessData_calculation)

        #interpolation

        lowerBoundry = ceiling(min(depthVectorOfData)/intervallBy)*intervallBy

        upperBoundry = floor(max(depthVectorOfData)/intervallBy)*intervallBy

        InterpolationMatrixRowNames =   approx (x = evennessData[,1],
                                                y = NULL,
                                                xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                                method = "linear",
                                                n = 50)[[1]]

        InterpolationMatrix = matrix(NA, nrow = length(InterpolationMatrixRowNames), ncol = dim(evennessData)[2]-3)

        matrixCounter = 0

        for (i in 1:dim(evennessData_calculation)[2]){

          matrixCounter = matrixCounter+1

          InterpolationMatrix[,matrixCounter] =  approx (x = evennessData[,1],
                                                         y = evennessData_calculation[,i],
                                                         xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                                         method = "linear",
                                                         n = 50)[[2]]

        }

        colnames(InterpolationMatrix) = colnames(evennessData)[4:dim(evennessData)[2]]
        rownames(InterpolationMatrix) = InterpolationMatrixRowNames

        #Calculations

        evennessMatrix = matrix(NA, ncol = 2, nrow = dim(InterpolationMatrix)[1])

        for (p in 1:dim(InterpolationMatrix)[1]){

          evennessVector = InterpolationMatrix[p,]

          evennessMatrix[p,1] = as.numeric(rownames(InterpolationMatrix)[p])
          evennessMatrix[p,2] = diversity(evennessVector)/log(specnumber(evennessVector))

        }

        #Get Correct Span

        GuessSpan = GuessLoess(intervall = evennessMatrix[,1],values = evennessMatrix[,2], overspan = 100)

        #Get Loess

        evennessLoess = predict(loess(evennessMatrix[,2] ~  evennessMatrix[,1], span = GuessSpan))

        #NonNegative

        if (NonNegative){

          evennessLoess[evennessLoess<0] = 0

        }

        #At Data

        evennessMatrixLoess = evennessMatrix
        evennessMatrixLoess[,2] = evennessLoess

        data[["Diatom"]][[DiatomNames[z]]][[Exportname]] = evennessMatrixLoess

      }
    }

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        z,"/",length(ls(data[["Diatom"]]))," calculating ",Exportname,sep="")


  }

  return(data)

}
