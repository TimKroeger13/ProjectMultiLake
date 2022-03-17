#'rateofChange
#'@description
#'Calculates the Rate of change for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@param minimumRowsAfterInterpolating alue for the minimum rows after filtering.
#'@param method Method for calculation Dissimilarity Indices for Community Ecologists.
#'@param MinAgeIntervall Minimal Intervall to interpolate to "Data Resolution".
#'@param NonNegative Creates Positive Values after Loess calculation.
#'@param Importname importname after data$Diatom$DiatomNames$
#'@param Exportname data$Diatom$DiatomNames$
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

rateofChange = function(data, intervallBy = 100, minimumRowsAfterInterpolating = 0, method = "bray", MinAgeIntervall = 1, NonNegative = TRUE, Importname = "", Exportname = ""){

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
      if(dim(rateOfChangeData)[1]>2){

        #Delete Doubles
        rateOfChangeData = deleteDoubles(rateOfChangeData)

        depthVectorOfData = rateOfChangeData[,1]

        #Calculation Data
        ROC_calculation = rateOfChangeData[4:dim(rateOfChangeData)[2]]

        #Sqrt transformation

        #ROC_calculation = sqrt(ROC_calculation)

        #MinAgeIntervall

        lowerBoundry = min(depthVectorOfData)

        upperBoundry = max(depthVectorOfData)

        SuperInterpolatedDataRowNames =  seq(from = lowerBoundry, to = upperBoundry, by = MinAgeIntervall)

        SuperInterpolatedMatrix = matrix(NA, nrow = length(SuperInterpolatedDataRowNames), ncol = dim(ROC_calculation)[2])

        #Minimal Interpolation

        for (l in 1:dim(ROC_calculation)[2]){

          CaclualtionValue = ROC_calculation[,l]

          #Interpolation

          SuperInterpolatedMatrix[,l] =  approx (x = depthVectorOfData,
                                                 y = CaclualtionValue,
                                                 xout = SuperInterpolatedDataRowNames,
                                                 method = "linear",
                                                 n = 50)[[2]]
        }

        rownames(SuperInterpolatedMatrix) = SuperInterpolatedDataRowNames

        #Calculate Dissimilarity Matrix

        DistanceMatrix = matrix(NA, ncol = 2, nrow = (dim(SuperInterpolatedMatrix)[1]-1))

        for (p in 1:(dim(SuperInterpolatedMatrix)[1]-1)){

          distdata = vegdist(SuperInterpolatedMatrix[p:(p+1),],method="euclidean",na.rm = T)

          DistanceMatrix[p,2] = distdata
          DistanceMatrix[p,1] = (as.numeric(SuperInterpolatedDataRowNames[p]) + as.numeric(SuperInterpolatedDataRowNames[p+1])) /2

        }

        #Cluster them Back to wanted Time intervals

        lowerBoundry = ceiling(min(depthVectorOfData)/intervallBy)*intervallBy

        upperBoundry = floor(max(depthVectorOfData)/intervallBy)*intervallBy

        ClusterRowNames =  seq(from = lowerBoundry, to = upperBoundry, by = intervallBy)

        ClusteryMatrix = matrix(NA, ncol = 2, nrow = length(ClusterRowNames)-1)

        for (p in 1:length(ClusterRowNames)-1){

          ClusterLocation = which((DistanceMatrix[,1] > ClusterRowNames[p]) + (DistanceMatrix[,1] < ClusterRowNames[p+1]) == 2)

          ClusterData = mean(DistanceMatrix[ClusterLocation,2])

          ClusteryMatrix[p,2] = ClusterData
          ClusteryMatrix[p,1] = (ClusterRowNames[p] + ClusterRowNames[p+1]) / 2

        }

        #Get Correct Span

        GuessSpan = GuessLoess(intervall = ClusteryMatrix[,1],values = ClusteryMatrix[,2], overspan = 100)

        #Get Loess

        RocLoess = predict(loess(ClusteryMatrix[,2] ~  ClusteryMatrix[,1], span = GuessSpan))

        #Ingnore negative values

        if (NonNegative){

          RocLoess[RocLoess<0] = 0

        }

        RocLoessMatrix = matrix(NA, ncol = 2, nrow = length(RocLoess))
        RocLoessMatrix[,1] = ClusteryMatrix[,1]
        RocLoessMatrix[,2] = RocLoess

        data[["Diatom"]][[DiatomNames[z]]][[Exportname]] = RocLoessMatrix

      }
    }

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        z,"/",length(ls(data[["Diatom"]]))," calculating ",Exportname,sep="")

  }

  return(data)

}






