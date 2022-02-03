#'CorrectionPlots
#'@description
#'Plots Ordinations Plots for the Pinoeer correction
#'@param data List of data generates by the Ordination function.
#'@param AllDiatomsNames All Diatom names from the Ordination function.
#'@param Allcolor The color from the Ordination function.
#'@param MaxAge MayAge from the Ordination function.
#'@importFrom grDevices dev.off pdf
#'@importFrom graphics abline
#'@export
#'@return nothing.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

CorrectionPlots = function(data, AllDiatomsNames, Allcolor, MaxAge){

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

  DiatomNames = ls(data$Diatom)
  CorrectionPoints = data$CorrectionPoints
  CorrectionPoints[is.na(CorrectionPoints[,2]),2] = 0

  for (z in 1:length(DiatomNames)){

    rateOfChangeData = data$Diatom[[DiatomNames[z]]]$SRS_data
    DiatomsNames = AllDiatomsNames[z]
    RID = which(data$CoreList[,1]==DiatomsNames)

    #Loess values

    ROCvalues = data$Diatom[[DiatomsNames]]$RoC
    ROCvalues[,2] = ROCvalues[,2]# / 100

    if(!is.null(rateOfChangeData)){

        #Delete Doubles
        rateOfChangeData = deleteDoubles(rateOfChangeData)

        depthVectorOfData = rateOfChangeData[,1]

        dissimilarityMatrix = rateOfChangeData[4:dim(rateOfChangeData)[2]]

        rownames(dissimilarityMatrix) = depthVectorOfData

        #check this minRows of the interpolated data


        #Sqrt transform
        # Got cut. Because of the interpolation values can be between 0 and 1.
        # This causes some values to grow instead of shrinking.

        #Distances

        DistanceMatrix = matrix(NA, ncol = 2, nrow = (dim(dissimilarityMatrix)[1]-1))

        for (p in 1:(dim(dissimilarityMatrix)[1]-1)){

          TimeDifference = depthVectorOfData[p+1] - depthVectorOfData[p]

          distdata = vegdist(dissimilarityMatrix[p:(p+1),],method="bray",na.rm = T) / TimeDifference

          DistanceMatrix[p,2] = distdata
          DistanceMatrix[p,1] = (as.numeric(rownames(dissimilarityMatrix)[p]) + as.numeric(rownames(dissimilarityMatrix)[p+1])) /2

        }


        #Plot

        pdf(paste(RID,"_",DiatomsNames,".pdf",sep=""),width=15,height=10)

        plot(NA,
             ylim=c(min(min(DistanceMatrix[,2]),min(ROCvalues[,2])),max(max(DistanceMatrix[,2]),max(ROCvalues[,2]))),
             xlim=c(max(max(DistanceMatrix[,1]),max(max(ROCvalues[,1]))),min(min(DistanceMatrix[,1]),min(min(ROCvalues[,1])))),
             ylab="Value",
             xlab="Age",
             main=paste("RoC | ",RID," | ",DiatomsNames,sep="")
        )

        lines(DistanceMatrix[,1],DistanceMatrix[,2],col="black", lwd=1,type = "p", pch = 1)
        lines(DistanceMatrix[,1],DistanceMatrix[,2],col=Allcolor[RID], lwd=1,type = "l", pch = 1)

        distance = 120

        if(RID>=10){distance = 200}

        text(DistanceMatrix[dim(DistanceMatrix)[1],1]+distance,
             DistanceMatrix[dim(DistanceMatrix)[1],2],
             label=RID,
             col="white",
             cex=1.2)

        text(DistanceMatrix[dim(DistanceMatrix)[1],1]+distance,
             DistanceMatrix[dim(DistanceMatrix)[1],2],
             label=RID,
             col=Allcolor[RID])

        #Show CorrectionPoints

        PinoeerPoint = CorrectionPoints[which(CorrectionPoints[,1] == DiatomsNames),2]

        abline(v = MaxAge,lty = 3)

        if(!length(PinoeerPoint)==0){
          if(PinoeerPoint>0){

            points(DistanceMatrix[length(DistanceMatrix[,1]):(length(DistanceMatrix[,1])-(PinoeerPoint-1)),1],
                   DistanceMatrix[length(DistanceMatrix[,1]):(length(DistanceMatrix[,1])-(PinoeerPoint-1)),2],
                   col="red", lwd=1,type = "p", pch = 19)

          }
        }

        #Plot Real Rate of change

        lines(ROCvalues[,1],ROCvalues[,2],col="black", lwd=1,type = "l", pch = 1, lty=2)
        lines(ROCvalues[,1],ROCvalues[,2],col="black", lwd=1,type = "p", pch = 8)


        dev.off()


    }
  }

  setwd(orginalWorkingDirectoryPath)

}








