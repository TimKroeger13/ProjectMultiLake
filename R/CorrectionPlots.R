#'CorrectionPlots
#'@description
#'Plots Ordinations Plots for the Pinoeer correction
#'@param data List of data generates by the Ordination function.
#'@param AllDiatomsNames All Diatom names from the Ordination function.
#'@param MaxAge MayAge from the Ordination function.
#'@importFrom grDevices dev.off pdf
#'@importFrom graphics abline
#'@export
#'@return nothing.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

CorrectionPlots = function(data, AllDiatomsNames, MaxAge){

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

  GetOneAgeOverMaxAgeVector = function(MaxAgeData, MaxAge){

    CutbyAge = length(MaxAgeData[sort(MaxAgeData)<=MaxAge]) + 1

    if(dim(MaxAgeData)[1] > CutbyAge){

      MaxAgeData = MaxAgeData[1:CutbyAge]

    }

    return(MaxAgeData)

  }

  GetOneAgeByMaxAge = function(MaxAgeData, MaxAge){

    CutbyAge = length(MaxAgeData[,1][sort(MaxAgeData[,1])<=MaxAge])

    if(dim(MaxAgeData)[1] > CutbyAge){

      MaxAgeData = MaxAgeData[1:CutbyAge,]

    }

    return(MaxAgeData)

  }

  CorrectionPoints = data$CorrectionPoints
  CorrectionPoints[is.na(CorrectionPoints[,2]),2] = 0

  for (z in 1:length(AllDiatomsNames)){

    RocData = data$Diatom[[AllDiatomsNames[z]]]$RoC
    CutRocData = data$Diatom[[AllDiatomsNames[z]]]$Cut_RoC

    DiatomsNames = AllDiatomsNames[z]
    RID = which(data$CoreList[,1]==DiatomsNames)

    if(!is.null(RocData)){

      # Get Depth vectors

      depthVectorOfData = data$Diatom[[AllDiatomsNames[z]]]$SRS_data[,1]
      cutdepthVectorOfData = data$Diatom[[AllDiatomsNames[z]]]$Cut_SRS_data[,1]

      depthVectorOfData = GetOneAgeOverMaxAgeVector(depthVectorOfData, MaxAge)

      RocData = GetOneAgeByMaxAge(RocData, MaxAge)

      #Create Plots

      pdf(paste(RID,"_",DiatomsNames,".pdf",sep=""),width=15,height=10)

      plot(NA,
           ylim=c(min(min(RocData[,2]),min(CutRocData[,2])),
                  max(max(RocData[,2]),max(CutRocData[,2]))),
           xlim=c(max(depthVectorOfData),min(depthVectorOfData)),
           ylab="Value",
           xlab="Age",
           main=paste(RID," | ",DiatomsNames," | Rate of change control",sep="")
      )

      lines(CutRocData[,1],CutRocData[,2],col="red", lwd=3,type = "l", pch = 1)

      lines(RocData[,1],RocData[,2],col="black", lwd=2,type = "l", pch = 1, lty=2)

      abline(v = intersect(depthVectorOfData,cutdepthVectorOfData),lty = 3, col = "black")

      abline(v = setdiff(depthVectorOfData,cutdepthVectorOfData),lty = 1, col = "red")

      dev.off()

    }
  }

  setwd(orginalWorkingDirectoryPath)

}



