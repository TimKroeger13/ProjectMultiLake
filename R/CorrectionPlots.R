#'CorrectionPlots
#'@description
#'Plots Ordinations Plots for the Pinoeer correction
#'@param data List of data generates by the Ordination function.
#'@param minimumRowsAfterCutOutMaxAge minimum rows count from the Ordination function.
#'@param MaxAge MayAge from the Ordination function.
#'@param AllDiatomsNames All Diatom names from the Ordination function.
#'@param Allcolor The color from the Ordination function.
#'@importFrom grDevices dev.off pdf
#'@importFrom graphics abline
#'@export
#'@return nothing.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

CorrectionPlots = function(data, minimumRowsAfterCutOutMaxAge, MaxAge, AllDiatomsNames, Allcolor){

  CutOutMaxAgeForInterpolatedData = function(DataVector, Cutage, minimumRowsAfterCutOutMaxAge){

    FittingAges = as.numeric(DataVector[,1])<=Cutage

    if(sum(FittingAges)>minimumRowsAfterCutOutMaxAge){

      DataVector=cbind(DataVector[FittingAges,1],DataVector[FittingAges,2])

      return(DataVector)

    }else{

      return(NULL)

    }
  }

  #Set folder

  orginalWorkingDirectoryPath=getwd()

  getwdTry <-  try(setwd(paste(orginalWorkingDirectoryPath,.Platform[2],"CorrectionPlots",sep="")),silent = TRUE)

  if(class(getwdTry) == "try-error"){

    dir.create(paste(orginalWorkingDirectoryPath,.Platform[2],"CorrectionPlots",sep=""))

    setwd(paste(orginalWorkingDirectoryPath,.Platform[2],"CorrectionPlots",sep=""))

  }

  Xmax=0
  Ymax=0
  Xmin=Inf
  Ymin=Inf

  pdf("RoC Control Plot All.pdf",width=15,height=10)

  for (i in 1:length(data$Diatom)){

    DiatomsNames = AllDiatomsNames[i]
    Values = data$Diatom[[DiatomsNames]]$RoC
    RID = which(data$CoreList[,1]==DiatomsNames)

    if(!is.null(Values)){
      if(dim(Values)[1]>0){

        Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

        if(!is.null(Values)){

          if(max(Values[,2])>Xmax){
            Xmax=max(Values[,2])
          }
          if(min(Values[,2])<Xmin){
            Xmin=min(Values[,2])
          }
          if(max(Values[,1])>Ymax){
            Ymax=max(Values[,1])
          }
          if(min(Values[,1])<Ymin){
            Ymin=min(Values[,1])
          }
        }
      }
    }
  }

  plot(NA,
       ylim=c(Xmin,Xmax),
       xlim=c(Ymax,Ymin),
       ylab="Value",
       xlab="Age",
       main=paste("RoC Control Plot All",sep="")
  )

  for (i in 1:length(data$Diatom)){

    DiatomsNames = AllDiatomsNames[i]
    Values = data$Diatom[[DiatomsNames]]$RoC
    RID = which(data$CoreList[,1]==DiatomsNames)

    if(!is.null(Values)){
      if(dim(Values)[1]>0){

        Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

        if(!is.null(Values)){

          lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1,type = "p", pch = 1)
          lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1,type = "l", pch = 1)

        }
      }
    }
  }

  #Plot Numbers

  for (i in 1:length(data$Diatom)){

    DiatomsNames = AllDiatomsNames[i]
    Values = data$Diatom[[DiatomsNames]]$RoC
    RID = which(data$CoreList[,1]==DiatomsNames)

    if(!is.null(Values)){
      if(dim(Values)[1]>0){

        Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

        if(!is.null(Values)){

          distance = 120

          if(i>=10){distance = 200}

          text(Values[dim(Values)[1],1]+distance,
               Values[dim(Values)[1],2],
               label=RID,
               col="white",
               cex=1.2)

          text(Values[dim(Values)[1],1]+distance,
               Values[dim(Values)[1],2],
               label=RID,
               col=Allcolor[RID])

        }
      }
    }
  }

  dev.off()

  ##############################################################################
  ######################## Single plots with indicator #########################
  ##############################################################################

  CorrectionPoints = data$CorrectionPoints
  CorrectionPoints[is.na(CorrectionPoints[,2]),2] = 0

  Xmax=0
  Ymax=0
  Xmin=Inf
  Ymin=Inf


  for (i in 1:length(data$Diatom)){

    DiatomsNames = AllDiatomsNames[i]
    Values = data$Diatom[[DiatomsNames]]$RoC
    RID = which(data$CoreList[,1]==DiatomsNames)

    if(!is.null(Values)){
      if(dim(Values)[1]>0){

        Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

        if(!is.null(Values)){

          if(max(Values[,2])>Xmax){
            Xmax=max(Values[,2])
          }
          if(min(Values[,2])<Xmin){
            Xmin=min(Values[,2])
          }
          if(max(Values[,1])>Ymax){
            Ymax=max(Values[,1])
          }
          if(min(Values[,1])<Ymin){
            Ymin=min(Values[,1])
          }
        }
      }
    }
  }

  for (i in 1:length(data$Diatom)){

    DiatomsNames = AllDiatomsNames[i]
    Values = data$Diatom[[DiatomsNames]]$RoC
    RID = which(data$CoreList[,1]==DiatomsNames)

    if(!is.null(Values)){
      if(dim(Values)[1]>0){

        Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

        if(!is.null(Values)){

          pdf(paste(RID,"_",DiatomsNames,".pdf",sep=""),width=15,height=10)

          plot(NA,
               ylim=c(Xmin,Xmax),
               xlim=c(Ymax,Ymin),
               ylab="Value",
               xlab="Age",
               main=paste("RoC | ",RID," | ",DiatomsNames,sep="")
          )


          lines(Values[,1],Values[,2],col="black", lwd=1,type = "p", pch = 1)
          lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1,type = "l", pch = 1)

          distance = 120

          if(i>=10){distance = 200}

          text(Values[dim(Values)[1],1]+distance,
               Values[dim(Values)[1],2],
               label=RID,
               col="white",
               cex=1.2)

          text(Values[dim(Values)[1],1]+distance,
               Values[dim(Values)[1],2],
               label=RID,
               col=Allcolor[RID])

          #Show CorrectionPoints

          PinoeerPoint = CorrectionPoints[which(CorrectionPoints[,1] == DiatomsNames),2]

          if(!length(PinoeerPoint)==0){
            if(PinoeerPoint>0){

              abline(v = Values[(length(Values[,1])-(PinoeerPoint-1)),1],lty = 3)

              points(Values[length(Values[,1]):(length(Values[,1])-(PinoeerPoint-1)),1],
                     Values[length(Values[,1]):(length(Values[,1])-(PinoeerPoint-1)),2],
                     col="red", lwd=1,type = "p", pch = 19)

            }
          }

          dev.off()

        }
      }
    }
  }

  setwd(orginalWorkingDirectoryPath)


}








