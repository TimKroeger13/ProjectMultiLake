#'TOC
#'@description
#'Calculates the TOC for the multivar function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervals by to interpolate to.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterInterpolating Value for the minimum rows after filtering.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

TOC = function(data, intervallBy = 100, allLoessSpans = 0.8, minimumRowsAfterInterpolating = 0){

  deleteDoubles = function(doublData){

    counter = 0

    for (h in 1:(dim(doublData)[1]-1)){
      #for (h in 1:39){

      counter = counter+1

      if(!is.na(as.numeric(row.names(doublData))[counter])){

        if(!is.na(as.numeric(row.names(doublData))[counter+1])){

          if(as.numeric(row.names(doublData))[counter]==as.numeric(row.names(doublData))[counter+1]){

            RNames = row.names(doublData)

            doublData=matrix(doublData[-counter,])
            row.names(doublData) = RNames[-counter]

            counter = counter-1

          }
        }
      }else{

        RNames = row.names(doublData)

        doublData=matrix(doublData[-counter,])
        row.names(doublData) = RNames[-counter]

        counter = counter-1

      }
    }

    if(is.na(as.numeric(row.names(doublData))[dim(doublData)[1]])){

      RNames = row.names(doublData)

      doublData=matrix(doublData[-dim(doublData)[1],])
      row.names(doublData) = RNames[-length(RNames)]

    }

    return(doublData)

  }

  CarbonNames = ls(data$Carbon)

  for (z in 1:length(CarbonNames)){

    TOCData = matrix(data$Carbon[[CarbonNames[z]]]$rawData$TOC, ncol = 1)
    rownames(TOCData) = data$Carbon[[CarbonNames[z]]]$rawData$depth

    if(sum(is.na(TOCData))==0){
      if(!is.null(TOCData)){

        #Delete Doubles
        TOCData = deleteDoubles(na.omit(TOCData))

        depthVectorOfData = as.numeric(row.names(TOCData))

        lowerBoundry = ceiling(min(depthVectorOfData)/intervallBy)*intervallBy

        upperBoundry = floor(max(depthVectorOfData)/intervallBy)*intervallBy

        InterpolationMatrixRowNames =   approx (x = TOCData[,1],
                                                y = NULL,
                                                xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                                method = "linear",
                                                n = 50)[[1]]

        InterpolationMatrix = matrix(NA, nrow = length(InterpolationMatrixRowNames), ncol = dim(TOCData)[2])


        InterpolationMatrix[,1] =  approx (x = depthVectorOfData,
                                           y = TOCData,
                                           xout = seq(from = lowerBoundry, to = upperBoundry, by = intervallBy),
                                           method = "linear",
                                           n = 50)[[2]]

        colnames(InterpolationMatrix) = "TOC"
        rownames(InterpolationMatrix) = InterpolationMatrixRowNames

        #check this minRows of the interpolated data

        if(dim(InterpolationMatrix)[1]>minimumRowsAfterInterpolating){

          InterpolationMatrixLoess = InterpolationMatrix
          InterpolationMatrixLoess[]=NA

          InterpolationMatrixLoess[,1] = predict(loess(InterpolationMatrix ~ InterpolationMatrixRowNames, span = allLoessSpans))

          InterpolationMatrixLoess[InterpolationMatrixLoess<0]=0 #<--------------- Not for MDS

          TOCMatrix = matrix(NA, ncol = 2, nrow = dim(InterpolationMatrixLoess)[1])

          TOCMatrix[,1] = InterpolationMatrixRowNames
          TOCMatrix[,2] = InterpolationMatrixLoess

        }

        data[["Carbon"]][[CarbonNames[z]]][["TOC"]] = TOCMatrix

        #Printer
        cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
            z,"/",length(ls(data[["Carbon"]]))," calculating TOC",sep="")

      }
    }
  }

  return(data)

}
