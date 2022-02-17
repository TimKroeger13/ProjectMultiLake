#'CorrectionPlots
#'@description
#'Plots Ordinations Plots for the Pinoeer correction
#'@param data List of data generates by the Ordination function.
#'@param AllDiatomsNames All Diatom names from the Ordination function.
#'@param Allcolor The color from the Ordination function.
#'@param MaxAge MayAge from the Ordination function.
#'@param NonNegative Creates Positive Values after Loess calculation.
#'@param intervallBy Intervalls by to interpolate to.
#'@param MinAgeIntervall Minimal Intervall to interpolate to "Data Resolution".
#'@importFrom grDevices dev.off pdf
#'@importFrom graphics abline
#'@export
#'@return nothing.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

CorrectionPlots = function(data, AllDiatomsNames, Allcolor, MaxAge, NonNegative = TRUE, intervallBy = 100, MinAgeIntervall = 1){

  #Set folder

  orginalWorkingDirectoryPath=getwd()

  getwdTry <-  try(setwd(paste(orginalWorkingDirectoryPath,.Platform[2],"CorrectionPlots",sep="")),silent = TRUE)

  if(class(getwdTry) == "try-error"){

    dir.create(paste(orginalWorkingDirectoryPath,.Platform[2],"CorrectionPlots",sep=""))

    setwd(paste(orginalWorkingDirectoryPath,.Platform[2],"CorrectionPlots",sep=""))

  }

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

  GetOneAgeOverMaxAge = function(MaxAgeData, MaxAge){

    CutbyAge = length(MaxAgeData[,1][sort(MaxAgeData[,1])<=MaxAge]) + 1

    if(dim(MaxAgeData)[1] > CutbyAge){

      MaxAgeData = MaxAgeData[1:CutbyAge,]

    }

    return(MaxAgeData)

  }

  CorrectionPoints = data$CorrectionPoints
  CorrectionPoints[is.na(CorrectionPoints[,2]),2] = 0

  for (z in 1:length(AllDiatomsNames)){

    rateOfChangeData = data$Diatom[[AllDiatomsNames[z]]]$SRS_data
    DiatomsNames = AllDiatomsNames[z]
    RID = which(data$CoreList[,1]==DiatomsNames)

    if(!is.null(rateOfChangeData)){

      #Delete Doubles
      rateOfChangeData = deleteDoubles(rateOfChangeData)

      rateOfChangeData = GetOneAgeOverMaxAge(rateOfChangeData, MaxAge)

      depthVectorOfData = rateOfChangeData[,1]

      ROC_calculation = rateOfChangeData[4:dim(rateOfChangeData)[2]]

      #Sqrt transformation

      ROC_calculation = sqrt(ROC_calculation)

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

    }

    #Calculate ROC for Cut Data

    rateOfChangeData = data$Diatom[[AllDiatomsNames[z]]]$Cut_SRS_data
    DiatomsNames = AllDiatomsNames[z]
    RID = which(data$CoreList[,1]==DiatomsNames)


    if(!is.null(rateOfChangeData)){

      #Delete Doubles
      rateOfChangeData = deleteDoubles(rateOfChangeData)

      rateOfChangeData = GetOneAgeOverMaxAge(rateOfChangeData, MaxAge)

      cutdepthVectorOfData = rateOfChangeData[,1]

      ROC_calculation = rateOfChangeData[4:dim(rateOfChangeData)[2]]

      #Sqrt transformation

      ROC_calculation = sqrt(ROC_calculation)

      #MinAgeIntervall

      lowerBoundry = min(cutdepthVectorOfData)

      upperBoundry = max(cutdepthVectorOfData)

      SuperInterpolatedDataRowNames =  seq(from = lowerBoundry, to = upperBoundry, by = MinAgeIntervall)

      SuperInterpolatedMatrix = matrix(NA, nrow = length(SuperInterpolatedDataRowNames), ncol = dim(ROC_calculation)[2])

      #Minimal Interpolation

      for (l in 1:dim(ROC_calculation)[2]){

        CaclualtionValue = ROC_calculation[,l]

        #Interpolation

        SuperInterpolatedMatrix[,l] =  approx (x = cutdepthVectorOfData,
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

      lowerBoundry = ceiling(min(cutdepthVectorOfData)/intervallBy)*intervallBy

      upperBoundry = floor(max(cutdepthVectorOfData)/intervallBy)*intervallBy

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

      CutRocLoessMatrix = matrix(NA, ncol = 2, nrow = length(RocLoess))
      CutRocLoessMatrix[,1] = ClusteryMatrix[,1]
      CutRocLoessMatrix[,2] = RocLoess

      #Create Plots

      pdf(paste(RID,"_",DiatomsNames,".pdf",sep=""),width=15,height=10)

      plot(NA,
           ylim=c(min(min(RocLoessMatrix[,2]),min(CutRocLoessMatrix[,2])),
                  max(max(RocLoessMatrix[,2]),max(CutRocLoessMatrix[,2]))),
           xlim=c(max(depthVectorOfData),min(depthVectorOfData)),
           ylab="Value",
           xlab="Age",
           main=paste(RID," | ",DiatomsNames," | Rate of change control",sep="")
      )


      lines(CutRocLoessMatrix[,1],CutRocLoessMatrix[,2],col="red", lwd=3,type = "l", pch = 1)

      lines(RocLoessMatrix[,1],RocLoessMatrix[,2],col="black", lwd=2,type = "l", pch = 1, lty=2)

      abline(v = intersect(depthVectorOfData,cutdepthVectorOfData),lty = 3, col = "black")

      abline(v = setdiff(depthVectorOfData,cutdepthVectorOfData),lty = 1, col = "red")

      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          z,"/",length(ls(data[["Diatom"]]))," calculating Rate of change Control Plots",sep="")

      dev.off()

    }
  }

  setwd(orginalWorkingDirectoryPath)

}



