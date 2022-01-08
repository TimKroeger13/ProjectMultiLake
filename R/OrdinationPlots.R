#'Ordination
#'@description
#'Plots Ordinations Plots.
#'@param data List of data generates by the Multivar function.
#'@param minimumRowsAfterCutOutMaxAge minimum rows count after filtering.
#'@param allspan The span for all loess functions.
#'@param MaxAge The max Age where to cut the Data.
#'@importFrom grDevices rainbow
#'@importFrom graphics lines points text legend barplot par
#'@importFrom stats na.omit
#'@importFrom grDevices dev.off pdf
#'@export
#'@return nothing.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

Ordination = function(data, minimumRowsAfterCutOutMaxAge = 12, allspan = 1, MaxAge = 20000){

  #new

  orginalWorkingDirectoryPath=getwd()

  getwdTry <-  try(setwd(paste(getwd(),.Platform[2],paste("Out_",data$Filter,sep = ""),sep="")),silent = TRUE)

  if(class(getwdTry) == "try-error"){

    dir.create(paste(getwd(),.Platform[2],paste("Out_",data$Filter,sep = ""),sep=""))

    setwd(paste(getwd(),.Platform[2],paste("Out_",data$Filter,sep = ""),sep=""))

  }

  CutOutMaxAge = function(DataVector, Cutage, minimumRowsAfterCutOutMaxAge){

    FittingAges = as.numeric(row.names(DataVector))<=Cutage

    if(sum(FittingAges)>minimumRowsAfterCutOutMaxAge){

      DataVector=cbind(as.numeric(row.names(DataVector))[FittingAges],DataVector[FittingAges])

      return(DataVector)

    }else{

      return(NULL)

    }
  }

  CutOutMaxAgeForInterpolatedData = function(DataVector, Cutage, minimumRowsAfterCutOutMaxAge){

    FittingAges = as.numeric(DataVector[,1])<=Cutage

    if(sum(FittingAges)>minimumRowsAfterCutOutMaxAge){

      DataVector=cbind(DataVector[FittingAges,1],DataVector[FittingAges,2])

      return(DataVector)

    }else{

      return(NULL)

    }
  }

  NormalizeRichnessDataTable = function(RichnessData){

    for (i in 1:dim(RichnessData)[2]){

      RichnessData[,i] = (RichnessData[,i]-min(RichnessData[,i]))/(max(RichnessData[,i])-min(RichnessData[,i]))

    }

    return(RichnessData)

  }

  AgeNormailzer = function(RichnessData){

    RichnessData = RichnessData[,1]

    return((RichnessData-min(RichnessData))/(max(RichnessData)-min(RichnessData)))

  }

  DeleteMeanNAS = function(MeanMatrix){

    MeanMatrixCounter = 0

    while (MeanMatrixCounter<dim(MeanMatrix)[1]) {

      MeanMatrixCounter = MeanMatrixCounter+1

      if(MeanMatrix[MeanMatrixCounter,2]==0){

        MeanMatrix = MeanMatrix[-MeanMatrixCounter,]

        MeanMatrixCounter = MeanMatrixCounter-1

      }
    }

    return(MeanMatrix)

  }




  ################################################################################
  ################################# Discription ##################################
  ################################################################################

  Allcolor = rainbow(dim(data$CoreList)[1])
  AllDiatomsNames = ls(data$Diatom)
  AllCarbonsNames = ls(data$Carbon)
  AllTRACENames = ls(data$TRACE)

  ################################################################################
  ##################################### MDS ######################################
  ################################################################################

  PlotsVariantsLoess = c("Normal","Loess")
  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){
      if(VariantsTransform == "Non Transformed"){
        pdf(paste("Z_old MDS_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        #Plot Limits

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(rownames(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

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
             main=paste("Z_old MDS | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

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
      }

      if(VariantsTransform == "Transformed"){
        pdf(paste("Z_old MDS_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(rownames(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

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
             main=paste("Z_old MDS | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

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
      }


    }
  }

  ################################################################################
  ################################### Richness ###################################
  ################################################################################

  nameLoess="Loess_invsimpson"
  nameNormal="invsimpson"
  minlegth = 10

  PlotsVariantsLoess = c("Normal","Loess")
  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){

      if(VariantsTransform == "Non Transformed"){

        Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])

        if(!is.null(Values)){

          pdf(paste("Z_old N2_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

          #Plot Limits

          Xmax=0
          Ymax=0
          Xmin=Inf
          Ymin=Inf

          for (i in 1:length(data$Diatom)){

            DiatomsNames = AllDiatomsNames[i]
            Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]]) #  $nMDS$Dim1
            RID = which(data$CoreList[,1]==DiatomsNames)

            if(!is.null(Values)){
              if(dim(Values)[1]>0){
                if(length(Values)>=minlegth){

                  valueNames = as.numeric(rownames(Values))
                  Values = as.numeric(Values)
                  Values = matrix(Values, ncol = 1)
                  rownames(Values) = valueNames

                  Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

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
          }

          plot(NA,
               ylim=c(Xmin,Xmax),
               xlim=c(Ymax,Ymin),
               ylab="Value",
               xlab="Age",
               main=paste("Z_old N2 | ",VariantsLoess, " | ",VariantsTransform,sep="")
          )

          for (i in 1:length(data$Diatom)){

            DiatomsNames = AllDiatomsNames[i]
            Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])
            RID = which(data$CoreList[,1]==DiatomsNames)

            if(!is.null(Values)){
              if(dim(Values)[1]>0){
                if(length(Values)>=minlegth){

                  Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

                  if(!is.null(Values)){

                    points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                    if(VariantsLoess == "Loess"){

                      #Loess
                      ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                      Values = cbind(Values[,1],ValuesLoess$fitted)

                    }

                    lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

                  }
                }
              }
            }
          }

          #Plot Numbers

          for (i in 1:length(data$Diatom)){

            DiatomsNames = AllDiatomsNames[i]
            Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])
            RID = which(data$CoreList[,1]==DiatomsNames)

            if(!is.null(Values)){
              if(dim(Values)[1]>0){
                if(length(Values)>=minlegth){

                  Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

                  if(!is.null(Values)){

                    if(VariantsLoess == "Loess"){

                      #Loess
                      ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                      Values = cbind(Values[,1],ValuesLoess$fitted)

                    }

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
          }
          dev.off()
        }
      }

      if(VariantsTransform == "Transformed"){

        Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])

        if(!is.null(Values)){

          pdf(paste("Z_old N2_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

          Xmax=0
          Ymax=0
          Xmin=Inf
          Ymin=Inf

          for (i in 1:length(data$Diatom)){

            DiatomsNames = AllDiatomsNames[i]
            Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])
            RID = which(data$CoreList[,1]==DiatomsNames)

            if(!is.null(Values)){
              if(dim(Values)[1]>0){
                if(length(Values)>=minlegth){

                  valueNames = as.numeric(rownames(Values))
                  Values = as.numeric(Values)
                  Values = matrix(Values, ncol = 1)
                  rownames(Values) = valueNames

                  Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

                  if(!is.null(Values)){

                    Values[,2] = scale(Values[,2],center = T, scale = T)

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
          }

          plot(NA,
               ylim=c(Xmin,Xmax),
               xlim=c(Ymax,Ymin),
               ylab="Value",
               xlab="Age",
               main=paste("Z_old N2 | ",VariantsLoess, " | ",VariantsTransform,sep="")
          )

          for (i in 1:length(data$Diatom)){

            DiatomsNames = AllDiatomsNames[i]
            Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])
            RID = which(data$CoreList[,1]==DiatomsNames)

            if(!is.null(Values)){
              if(dim(Values)[1]>0){
                if(length(Values)>=minlegth){

                  Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

                  if(!is.null(Values)){

                    Values[,2] = scale(Values[,2],center = T, scale = T)

                    points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                    if(VariantsLoess == "Loess"){

                      #Loess
                      ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                      Values = cbind(Values[,1],ValuesLoess$fitted)

                    }

                    lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

                  }
                }
              }
            }
          }

          #Plot Numbers

          for (i in 1:length(data$Diatom)){

            DiatomsNames = AllDiatomsNames[i]
            Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])
            RID = which(data$CoreList[,1]==DiatomsNames)

            if(!is.null(Values)){
              if(dim(Values)[1]>0){
                if(length(Values)>=minlegth){

                  Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

                  if(!is.null(Values)){

                    Values[,2] = scale(Values[,2],center = T, scale = T)

                    if(VariantsLoess == "Loess"){

                      #Loess
                      ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                      Values = cbind(Values[,1],ValuesLoess$fitted)

                    }

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
          }
          dev.off()
        }
      }
    }
  }



  ################################################################################
  ##################################### Clima ####################################
  ################################################################################

  PlotsVariantsLoess = c("Normal","Loess")
  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){

      if(VariantsTransform == "Non Transformed"){
        pdf(paste("PalioClima_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        #Plot Limits

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){ # length(data$Diatom)

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

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
             main=paste("Palio Clima | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=RID,
                     col="white",
                     cex=1.2)

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=RID,
                     col=Allcolor[RID])

              }
            }
          }
        }
        dev.off()
      }

      if(VariantsTransform == "Transformed"){
        pdf(paste("PalioClima_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))

              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

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
             main=paste("Palio Clima | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))

              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))

              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=RID,
                     col="white",
                     cex=1.2)

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=RID,
                     col=Allcolor[RID])

              }
            }
          }
        }
        dev.off()
      }
    }
  }



  ################################################################################
  #################################### Carbon ####################################
  ################################################################################

  nameLoess="Loess_invsimpson"
  nameNormal="invsimpson"

  PlotsVariantsLoess = c("Normal","Loess")
  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){

      if(VariantsTransform == "Non Transformed"){
        pdf(paste("Z_old TOC_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        #Plot Limits

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC) #  $Species_richness[[nameNormal]]
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(length(Values)>0){

              Values = data$Carbon[[DiatomsNames]]$rawData
              valueNames = Values$depth[!is.na(Values$TOC)]
              Values = Values$TOC[!is.na(Values$TOC)]

              #AgeError
              Values = Values[!is.na(valueNames)]
              valueNames = valueNames[!is.na(valueNames)]
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

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
             main=paste("Z_old TOC | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){ #length(data$Diatom)

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(length(Values)>0){

              Values = data$Carbon[[DiatomsNames]]$rawData
              valueNames = Values$depth[!is.na(Values$TOC)]
              Values = Values$TOC[!is.na(Values$TOC)]

              #AgeError
              Values = Values[!is.na(valueNames)]
              valueNames = valueNames[!is.na(valueNames)]
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(length(Values)>0){

              Values = data$Carbon[[DiatomsNames]]$rawData
              valueNames = Values$depth[!is.na(Values$TOC)]
              Values = Values$TOC[!is.na(Values$TOC)]

              #AgeError
              Values = Values[!is.na(valueNames)]
              valueNames = valueNames[!is.na(valueNames)]
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

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
      }

      if(VariantsTransform == "Transformed"){
        pdf(paste("Z_old TOC_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(length(Values)>0){

              Values = data$Carbon[[DiatomsNames]]$rawData
              valueNames = Values$depth[!is.na(Values$TOC)]
              Values = Values$TOC[!is.na(Values$TOC)]

              #AgeError
              Values = Values[!is.na(valueNames)]
              valueNames = valueNames[!is.na(valueNames)]
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] =  scale(Values[,2],center = T, scale = T)

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
             main=paste("Z_old TOC | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(length(Values)>0){

              Values = data$Carbon[[DiatomsNames]]$rawData
              valueNames = Values$depth[!is.na(Values$TOC)]
              Values = Values$TOC[!is.na(Values$TOC)]

              #AgeError
              Values = Values[!is.na(valueNames)]
              valueNames = valueNames[!is.na(valueNames)]
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] =  scale(Values[,2],center = T, scale = T)

                points(Values[,1],Values[,2],col=Allcolor[RID], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(length(Values)>0){

              Values = data$Carbon[[DiatomsNames]]$rawData
              valueNames = Values$depth[!is.na(Values$TOC)]
              Values = Values$TOC[!is.na(Values$TOC)]

              #AgeError
              Values = Values[!is.na(valueNames)]
              valueNames = valueNames[!is.na(valueNames)]
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,MaxAge,minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                Values[,2] =  scale(Values[,2],center = T, scale = T)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

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
      }


    }
  }



  #Legend

  pdf("Legend.pdf",width=15,height=10)

  # AllCoreList dim(data$CoreList)[1]

  LegendsSplit = length(Allcolor)/3

  SplitValue1 = ceiling(LegendsSplit)
  SplitValue2 = ceiling(LegendsSplit*2)
  SplitValue3 = length(Allcolor)

  AllNamesPart1 = data$CoreList[,1][1:SplitValue1]
  AllNamesPart2 = data$CoreList[,1][(SplitValue1+1):SplitValue2]
  AllNamesPart3 = data$CoreList[,1][(SplitValue2+1):SplitValue3]

  plot(NA,
       ylim=c(0,1),
       xlim=c(0,1),
       ylab="",
       xlab="",
       main="Legend",
       bty="n",
       xaxt='n',
       yaxt='n'
  )

  thickness = 1

  legend(x = "left", legend = paste(1:SplitValue1,"-",AllNamesPart1), col = Allcolor[1:SplitValue1], lty = 1, lwd = 3, text.width=0.3, cex = thickness, bty="n") # bty="n"

  legend(x = "center", legend = paste((SplitValue1+1):SplitValue2,"-",AllNamesPart2), col = Allcolor[(SplitValue1+1):SplitValue2], lty = 1, lwd = 3, text.width=0.3, cex = thickness, bty="n") # bty="n"

  legend(x = "right", legend = paste((SplitValue2+1):SplitValue3,"-",AllNamesPart3), col = Allcolor[(SplitValue2+1):SplitValue3], lty = 1, lwd = 3, text.width=0.3, cex = thickness, bty="n") # bty="n"

  dev.off()





































  ################################################################################
  ###################################### RoC #####################################
  ################################################################################

  pdf("RoC.pdf",width=15,height=10)

  RocMatrix = data$RocMatrix

  #Cut Data after x years
  CutValue = MaxAge

  RocMatrix=RocMatrix[RocMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  RocMatrix = DeleteMeanNAS(RocMatrix)

  Xmin = min(RocMatrix[,1])
  Xmax = max(RocMatrix[,1])
  Ymin = min(RocMatrix[,3])
  Ymax = max(RocMatrix[,4])

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(RocMatrix)[1])
  P_color[] = "red"
  P_color[RocMatrix[,5]<P_value] = "green"
  P_color[RocMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("RoC Mean",sep="")
  )

  #lines(RocMatrix[,1],RocMatrix[,2], col= "blue")
  #points(RocMatrix[,1],RocMatrix[,2], pch = 19, cex = 0.1, col= "black")

  for (i in 1:(dim(RocMatrix)[1]-1)){

    lines(RocMatrix[i:(i+1),1],RocMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(RocMatrix[,1],RocMatrix[,3], col= "black",lty=3,lwd=2)
  lines(RocMatrix[,1],RocMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()

  ################################################################################
  #################################     RoC     ##################################
  ############################### SpeciesRichness ################################
  ################################################################################

  pdf("Vector_RoC.pdf",width=15,height=10)

  RocMatrix = data$Vector_RocMatrix

  #Cut Data after x years
  CutValue = MaxAge

  RocMatrix=RocMatrix[RocMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  RocMatrix = DeleteMeanNAS(RocMatrix)

  Xmin = min(RocMatrix[,1])
  Xmax = max(RocMatrix[,1])
  Ymin = min(na.omit(RocMatrix[,3]))
  Ymax = max(na.omit(RocMatrix[,4]))

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(RocMatrix)[1])
  P_color[] = "red"
  P_color[RocMatrix[,5]<P_value] = "green"
  P_color[RocMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("RoC Mean",sep="")
  )

  #lines(RocMatrix[,1],RocMatrix[,2], col= "blue")
  #points(RocMatrix[,1],RocMatrix[,2], pch = 19, cex = 0.1, col= "black")

  for (i in 1:(dim(RocMatrix)[1]-1)){

    lines(RocMatrix[i:(i+1),1],RocMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(RocMatrix[,1],RocMatrix[,3], col= "black",lty=3,lwd=2)
  lines(RocMatrix[,1],RocMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()

  ################################################################################
  ################################### Evenness ###################################
  ################################################################################

  pdf("Evenness.pdf",width=15,height=10)

  EvennessMatrix = data$EvennessMatrix

  #Cut Data after x years
  CutValue = MaxAge

  EvennessMatrix=EvennessMatrix[EvennessMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  EvennessMatrix = DeleteMeanNAS(EvennessMatrix)


  Xmin = min(EvennessMatrix[,1])
  Xmax = max(EvennessMatrix[,1])
  Ymin = min(EvennessMatrix[,3])
  Ymax = max(EvennessMatrix[,4])

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(EvennessMatrix)[1])
  P_color[] = "red"
  P_color[EvennessMatrix[,5]<P_value] = "green"
  P_color[EvennessMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("Evenness Mean",CutValue,sep="")
  )

  #lines(EvennessMatrix[,1],EvennessMatrix[,2], col= "blue")
  #points(EvennessMatrix[,1],EvennessMatrix[,2], pch = 19, cex = 0.1, col= "black")

  for (i in 1:(dim(EvennessMatrix)[1]-1)){

    lines(EvennessMatrix[i:(i+1),1],EvennessMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(EvennessMatrix[,1],EvennessMatrix[,3], col= "black",lty=3,lwd=2)
  lines(EvennessMatrix[,1],EvennessMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()


  ################################################################################
  ############################### InverseSimpsion ################################
  ################################################################################

  pdf("InverseSimpsion.pdf",width=15,height=10)

  InverseSimpsionMatrix = data$InverseSimpsionMatrix

  #Cut Data after x years
  CutValue = MaxAge

  InverseSimpsionMatrix=InverseSimpsionMatrix[InverseSimpsionMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  InverseSimpsionMatrix = DeleteMeanNAS(InverseSimpsionMatrix)

  Xmin = min(InverseSimpsionMatrix[,1])
  Xmax = max(InverseSimpsionMatrix[,1])
  Ymin = min(InverseSimpsionMatrix[,3])
  Ymax = max(InverseSimpsionMatrix[,4])

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(InverseSimpsionMatrix)[1])
  P_color[] = "red"
  P_color[InverseSimpsionMatrix[,5]<P_value] = "green"
  P_color[InverseSimpsionMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("InverseSimpsion Mean",CutValue,sep="")
  )

  #lines(InverseSimpsionMatrix[,1],InverseSimpsionMatrix[,2], col= "blue")
  #points(InverseSimpsionMatrix[,1],InverseSimpsionMatrix[,2], pch = 19, cex = 0.1, col= "black")

  for (i in 1:(dim(InverseSimpsionMatrix)[1]-1)){

    lines(InverseSimpsionMatrix[i:(i+1),1],InverseSimpsionMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(InverseSimpsionMatrix[,1],InverseSimpsionMatrix[,3], col= "black",lty=3,lwd=2)
  lines(InverseSimpsionMatrix[,1],InverseSimpsionMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()


  ################################################################################
  ############################### SpeciesRichness ################################
  ################################################################################

  pdf("SpeciesRichness.pdf",width=15,height=10)

  SpeciesRichnessMatrix = data$SpeciesRichnessMatrix

  #Cut Data after x years
  CutValue = MaxAge

  SpeciesRichnessMatrix=SpeciesRichnessMatrix[SpeciesRichnessMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  SpeciesRichnessMatrix = DeleteMeanNAS(SpeciesRichnessMatrix)

  Xmin = min(SpeciesRichnessMatrix[,1])
  Xmax = max(SpeciesRichnessMatrix[,1])
  Ymin = min(SpeciesRichnessMatrix[,3])
  Ymax = max(SpeciesRichnessMatrix[,4])

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(SpeciesRichnessMatrix)[1])
  P_color[] = "red"
  P_color[SpeciesRichnessMatrix[,5]<P_value] = "green"
  P_color[SpeciesRichnessMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("SpeciesRichness Mean",CutValue,sep="")
  )

  for (i in 1:(dim(SpeciesRichnessMatrix)[1]-1)){

    lines(SpeciesRichnessMatrix[i:(i+1),1],SpeciesRichnessMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(SpeciesRichnessMatrix[,1],SpeciesRichnessMatrix[,3], col= "black",lty=3,lwd=2)
  lines(SpeciesRichnessMatrix[,1],SpeciesRichnessMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()


  ################################################################################
  ################################     Vector     ################################
  ############################### SpeciesRichness ################################
  ################################################################################

  pdf("Vector_SpeciesRichness.pdf",width=15,height=10)

  SpeciesRichnessMatrix = data$Vector_SpeciesRichnessMatrix

  #Cut Data after x years
  CutValue = MaxAge

  SpeciesRichnessMatrix=SpeciesRichnessMatrix[SpeciesRichnessMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  SpeciesRichnessMatrix = DeleteMeanNAS(SpeciesRichnessMatrix)

  Xmin = min(SpeciesRichnessMatrix[,1])
  Xmax = max(SpeciesRichnessMatrix[,1])
  Ymin = min(na.omit(SpeciesRichnessMatrix[,3]))
  Ymax = max(na.omit(SpeciesRichnessMatrix[,4]))

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(SpeciesRichnessMatrix)[1])
  P_color[] = "red"
  P_color[SpeciesRichnessMatrix[,5]<P_value] = "green"
  P_color[SpeciesRichnessMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("Vector SpeciesRichness Mean",CutValue,sep="")
  )

  for (i in 1:(dim(SpeciesRichnessMatrix)[1]-1)){

    lines(SpeciesRichnessMatrix[i:(i+1),1],SpeciesRichnessMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(SpeciesRichnessMatrix[,1],SpeciesRichnessMatrix[,3], col= "black",lty=3,lwd=2)
  lines(SpeciesRichnessMatrix[,1],SpeciesRichnessMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()


  ################################################################################
  ##################################### MDS ######################################
  ################################################################################

  pdf("MDS.pdf",width=15,height=10)

  MDSMatrix = data$MDSMatrix

  #Cut Data after x years
  CutValue = MaxAge

  MDSMatrix=MDSMatrix[MDSMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  MDSMatrix = DeleteMeanNAS(MDSMatrix)

  Xmin = min(MDSMatrix[,1])
  Xmax = max(MDSMatrix[,1])
  Ymin = min(MDSMatrix[,3])
  Ymax = max(MDSMatrix[,4])

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(MDSMatrix)[1])
  P_color[] = "red"
  P_color[MDSMatrix[,5]<P_value] = "green"
  P_color[MDSMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("MDS Mean",CutValue,sep="")
  )

  #lines(MDSMatrix[,1],MDSMatrix[,2], col= "blue")
  #points(MDSMatrix[,1],MDSMatrix[,2], pch = 19, cex = 0.1, col= "black")

  for (i in 1:(dim(MDSMatrix)[1]-1)){

    lines(MDSMatrix[i:(i+1),1],MDSMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(MDSMatrix[,1],MDSMatrix[,3], col= "black",lty=3,lwd=2)
  lines(MDSMatrix[,1],MDSMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()


  ################################################################################
  ##################################### TOC ######################################
  ################################################################################

  pdf("TOC.pdf",width=15,height=10)

  TOCMatrix = data$TOCMatrix

  #Cut Data after x years
  CutValue = MaxAge

  TOCMatrix=TOCMatrix[TOCMatrix[,1]<=CutValue,]

  #Plot Mean Rat of change

  TOCMatrix = DeleteMeanNAS(TOCMatrix)

  Xmin = min(TOCMatrix[,1])
  Xmax = max(TOCMatrix[,1])
  Ymin = min(TOCMatrix[,3])
  Ymax = max(TOCMatrix[,4])

  #Create Color after p Value from
  P_value = 0.05

  P_color <- vector( "character" , dim(TOCMatrix)[1])
  P_color[] = "red"
  P_color[TOCMatrix[,5]<P_value] = "green"
  P_color[TOCMatrix[,6]==1] = "black"

  plot(NA,
       ylim=c(Ymin,Ymax),
       xlim=c(Xmax,Xmin),
       ylab="Value",
       xlab="Age",
       main=paste("TOC Mean",CutValue,sep="")
  )

  #lines(TOCMatrix[,1],TOCMatrix[,2], col= "blue")
  #points(TOCMatrix[,1],TOCMatrix[,2], pch = 19, cex = 0.1, col= "black")

  for (i in 1:(dim(TOCMatrix)[1]-1)){

    lines(TOCMatrix[i:(i+1),1],TOCMatrix[i:(i+1),2], col= P_color[i],lwd=2)

  }

  lines(TOCMatrix[,1],TOCMatrix[,3], col= "black",lty=3,lwd=2)
  lines(TOCMatrix[,1],TOCMatrix[,4], col= "black",lty=3,lwd=2)

  dev.off()

  ################################################################################
  #################################### TRACE #####################################
  ################################################################################

  TRACEMeanPlot= function(TraceName = "UnknownTrace"){

    pdf(paste("TRACE_",TraceName,".pdf",sep=""),width=15,height=10)

    TRACEMatrix = data$TRACEMatrix[[TraceName]]

    #Cut Data after x years
    CutValue = MaxAge

    TRACEMatrix=TRACEMatrix[TRACEMatrix[,1]<=CutValue,]

    #Plot Mean Rat of change

    TRACEMatrix = DeleteMeanNAS(TRACEMatrix)

    Xmin = min(TRACEMatrix[,1])
    Xmax = max(TRACEMatrix[,1])
    Ymin = min(TRACEMatrix[,3])
    Ymax = max(TRACEMatrix[,4])

    #Create Color after p Value from
    P_value = 0.05

    P_color <- vector( "character" , dim(TRACEMatrix)[1])
    P_color[] = "red"
    P_color[TRACEMatrix[,5]<P_value] = "green"
    P_color[TRACEMatrix[,6]==1] = "black"

    plot(NA,
         ylim=c(Ymin,Ymax),
         xlim=c(Xmax,Xmin),
         ylab="Value",
         xlab="Age",
         main=paste("TRACE Mean ",TraceName," ",CutValue,sep="")
    )

    #lines(TRACEMatrix[,1],TRACEMatrix[,2], col= "blue")
    #points(TRACEMatrix[,1],TRACEMatrix[,2], pch = 19, cex = 0.1, col= "black")

    for (i in 1:(dim(TRACEMatrix)[1]-1)){

      lines(TRACEMatrix[i:(i+1),1],TRACEMatrix[i:(i+1),2], col= P_color[i],lwd=2)

    }

    '
    Xmin = 0

    lowerBoundry = ceiling(Xmin/ConvIntervall)*ConvIntervall
    upperBoundry = floor(Xmax/ConvIntervall)*ConvIntervall

    ConvNumbers = seq(from = lowerBoundry, to = upperBoundry, by = ConvIntervall)


    AgeMean = vector(mode="numeric", length=(length(ConvNumbers)-1))
    ValueMeanLower = vector(mode="numeric", length=(length(ConvNumbers)-1))
    ValueMeanUpper = vector(mode="numeric", length=(length(ConvNumbers)-1))


    for (k in 1:(length(ConvNumbers)-1)){

      AgeMean[k] = mean(TRACEMatrix[ConvNumbers[k]:ConvNumbers[k+1],1])
      ValueMeanLower[k] = mean(TRACEMatrix[ConvNumbers[k]:ConvNumbers[k+1],3])
      ValueMeanUpper[k] = mean(TRACEMatrix[ConvNumbers[k]:ConvNumbers[k+1],4])

    }

    '
    lines(TRACEMatrix[,1],TRACEMatrix[,3], col= "black",lty=3,lwd=2)
    lines(TRACEMatrix[,1],TRACEMatrix[,4], col= "black",lty=3,lwd=2)

    dev.off()

  }

  TRACEMeanPlot(TraceName = "JJA_mean")
  TRACEMeanPlot(TraceName = "DJF_mean")
  TRACEMeanPlot(TraceName = "hydrological_mean")

  ################################################################################
  ################################ Evenness Solo #################################
  ################################################################################

  PlotsVariantsLoess = c("Normal","Loess") # "Normal","Loess"
  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){

      Xmax=0
      Ymax=0
      Xmin=Inf
      Ymin=Inf

      pdf(paste("Evenness_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

      for (i in 1:length(data$Diatom)){

        DiatomsNames = AllDiatomsNames[i]
        Values = data$Diatom[[DiatomsNames]]$evenness
        RID = which(data$CoreList[,1]==DiatomsNames)

        if(!is.null(Values)){
          if(dim(Values)[1]>0){

            Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

            if(!is.null(Values)){

              if(VariantsTransform=="Transformed"){

                Values[,2] = scale(Values[,2],center = T, scale = T)

              }

              if(VariantsLoess == "Loess"){

                #Loess
                ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                Values = cbind(Values[,1],ValuesLoess$fitted)

              }

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
           main=paste("Evenness | ",VariantsLoess, " | ",VariantsTransform,sep="")
      )

      for (i in 1:length(data$Diatom)){

        DiatomsNames = AllDiatomsNames[i]
        Values = data$Diatom[[DiatomsNames]]$evenness
        RID = which(data$CoreList[,1]==DiatomsNames)

        if(!is.null(Values)){
          if(dim(Values)[1]>0){

            Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

            if(!is.null(Values)){

              if(VariantsTransform=="Transformed"){

                Values[,2] = scale(Values[,2],center = T, scale = T)

              }

              if(VariantsLoess == "Loess"){

                #Loess
                ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                Values = cbind(Values[,1],ValuesLoess$fitted)

              }

              lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

            }
          }
        }
      }

      #Plot Numbers

      for (i in 1:length(data$Diatom)){

        DiatomsNames = AllDiatomsNames[i]
        Values = data$Diatom[[DiatomsNames]]$evenness
        RID = which(data$CoreList[,1]==DiatomsNames)

        if(!is.null(Values)){
          if(dim(Values)[1]>0){

            Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

            if(!is.null(Values)){

              if(VariantsTransform=="Transformed"){

                Values[,2] = scale(Values[,2],center = T, scale = T)

              }

              if(VariantsLoess == "Loess"){

                #Loess
                ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                Values = cbind(Values[,1],ValuesLoess$fitted)

              }

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

    }
  }


  ################################################################################
  ################################### RoC Solo ###################################
  ################################################################################


  ImportVersions = c("RoC","Cut_RoC")

  for (IV in ImportVersions){

    PlotsVariantsLoess = c("Normal","Loess") # "Normal","Loess"
    PlotsVariantsTransform = c("Non Transformed","Transformed")

    for (VariantsLoess in PlotsVariantsLoess){
      for (VariantsTransform in PlotsVariantsTransform){

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        pdf(paste(IV,"_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = data$Diatom[[DiatomsNames]][[IV]]
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsTransform=="Transformed"){

                  Values[,2] = scale(Values[,2],center = T, scale = T)

                }

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

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
             main=paste(IV," | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = data$Diatom[[DiatomsNames]][[IV]]
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsTransform=="Transformed"){

                  Values[,2] = scale(Values[,2],center = T, scale = T)

                }

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = data$Diatom[[DiatomsNames]][[IV]]
          RID = which(data$CoreList[,1]==DiatomsNames)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsTransform=="Transformed"){

                  Values[,2] = scale(Values[,2],center = T, scale = T)

                }

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=allspan)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

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

      }
    }
  }

  ##############################################################################
  ############################ Inverse Simpson Solo ############################
  ##############################################################################

  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsTransform in PlotsVariantsTransform){

    Xmax=0
    Ymax=0
    Xmin=Inf
    Ymin=Inf

    pdf(paste("InverseSimpsion_",VariantsTransform,".pdf",sep=""),width=15,height=10)

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$inverseSimpsion
      RID = which(data$CoreList[,1]==DiatomsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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
         main=paste("InverseSimpsion | ",VariantsTransform,sep="")
    )

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$inverseSimpsion
      RID = which(data$CoreList[,1]==DiatomsNames)

      PointValues = data$Diatom[[DiatomsNames]]$Species_richness$invsimpson

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)
              PointValues = scale(PointValues,center = T, scale = T)

            }

            lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

            points(as.numeric(row.names(PointValues)),PointValues,col=Allcolor[RID], lwd=1, cex= 0.8)

          }
        }
      }
    }

    #Plot Numbers

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$inverseSimpsion
      RID = which(data$CoreList[,1]==DiatomsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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

  }

  ##############################################################################
  ########################### Species Richness Solo ############################
  ##############################################################################

  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsTransform in PlotsVariantsTransform){

    Xmax=0
    Ymax=0
    Xmin=Inf
    Ymin=Inf

    pdf(paste("SpeciesRichness_",VariantsTransform,".pdf",sep=""),width=15,height=10)

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$speciesRichness
      RID = which(data$CoreList[,1]==DiatomsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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
         main=paste("SpeciesRichness | ",VariantsTransform,sep="")
    )

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$speciesRichness
      RID = which(data$CoreList[,1]==DiatomsNames)

      PointValues = data$Diatom[[DiatomsNames]]$Species_richness$richness

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)
              PointValues = scale(PointValues,center = T, scale = T)

            }

            lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)
            points(as.numeric(row.names(PointValues)),PointValues,col=Allcolor[RID], lwd=1, cex= 0.8)

          }
        }
      }
    }

    #Plot Numbers

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$speciesRichness
      RID = which(data$CoreList[,1]==DiatomsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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

  }

  ##############################################################################
  ################################## MDS Solo ##################################
  ##############################################################################

  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsTransform in PlotsVariantsTransform){

    Xmax=0
    Ymax=0
    Xmin=Inf
    Ymin=Inf

    pdf(paste("MDS_",VariantsTransform,".pdf",sep=""),width=15,height=10)

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$MDS
      RID = which(data$CoreList[,1]==DiatomsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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
         main=paste("MDS | ",VariantsTransform,sep="")
    )

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$MDS
      RID = which(data$CoreList[,1]==DiatomsNames)

      PointValues = data$Diatom[[DiatomsNames]]$nMDS$Dim1

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)
              PointValues = scale(PointValues,center = T, scale = T)

            }

            lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

            points(as.numeric(row.names(PointValues)),PointValues,col=Allcolor[RID], lwd=1, cex= 0.8)

          }
        }
      }
    }

    #Plot Numbers

    for (i in 1:length(data$Diatom)){

      DiatomsNames = AllDiatomsNames[i]
      Values = data$Diatom[[DiatomsNames]]$MDS
      RID = which(data$CoreList[,1]==DiatomsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values, MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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

  }

  ##############################################################################
  ################################## TOC Solo ##################################
  ##############################################################################

  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsTransform in PlotsVariantsTransform){

    Xmax=0
    Ymax=0
    Xmin=Inf
    Ymin=Inf

    pdf(paste("TOC_",VariantsTransform,".pdf",sep=""),width=15,height=10)

    for (i in 1:length(data$Carbon)){

      CarbonsNames = AllCarbonsNames[i]
      Values = data$Carbon[[CarbonsNames]]$TOC
      RID = which(data$CoreList[,1]==CarbonsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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
         main=paste("TOC | ",VariantsTransform,sep="")
    )

    for (i in 1:length(data$Carbon)){

      CarbonsNames = AllCarbonsNames[i]
      Values = data$Carbon[[CarbonsNames]]$TOC
      RID = which(data$CoreList[,1]==CarbonsNames)

      PointValues = data$Carbon[[CarbonsNames]]$rawData$TOC

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)
              PointValues = scale(PointValues,center = T, scale = T)

            }

            lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

            points(as.numeric(data$Carbon[[CarbonsNames]]$rawData$depth),PointValues,col=Allcolor[RID], lwd=1, cex= 0.8)

          }
        }
      }
    }

    #Plot Numbers

    for (i in 1:length(data$Carbon)){

      CarbonsNames = AllCarbonsNames[i]
      Values = data$Carbon[[CarbonsNames]]$TOC
      RID = which(data$CoreList[,1]==CarbonsNames)

      if(!is.null(Values)){
        if(dim(Values)[1]>0){

          Values = CutOutMaxAgeForInterpolatedData(Values, MaxAge, minimumRowsAfterCutOutMaxAge)

          if(!is.null(Values)){

            if(VariantsTransform=="Transformed"){

              Values[,2] = scale(Values[,2],center = T, scale = T)

            }

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

  }

  ##############################################################################
  ################################# TRACE Solo #################################
  ##############################################################################

  TRACESoloPlot= function(TraceType = "UnknownTrace"){

    PlotsVariantsTransform = c("Non Transformed","Transformed")

    for (VariantsTransform in PlotsVariantsTransform){

      Xmax=0
      Ymax=0
      Xmin=Inf
      Ymin=Inf

      pdf(paste("TRACE_",TraceType," ",VariantsTransform,".pdf",sep=""),width=15,height=10)

      for (i in 1:length(data$TRACE)){

        TRACEsNames = AllTRACENames[i]
        Values = data$TRACE[[TRACEsNames]][[TraceType]]
        RID = which(data$CoreList[,1]==TRACEsNames)

        if(length(RID)>0){
          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsTransform=="Transformed"){

                  Values[,2] = scale(Values[,2],center = T, scale = T)

                }

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
      }

      plot(NA,
           ylim=c(Xmin,Xmax),
           xlim=c(Ymax,Ymin),
           ylab="Value",
           xlab="Age",
           main=paste("TRACE | ",TraceType," ",VariantsTransform,sep="")
      )

      for (i in 1:length(data$TRACE)){

        TRACEsNames = AllTRACENames[i]
        Values = data$TRACE[[TRACEsNames]][[TraceType]]
        RID = which(data$CoreList[,1]==TRACEsNames)

        #PointValues = data$TRACE[[TRACEsNames]]$rawData$TRACE

        if(length(RID)>0){
          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAgeForInterpolatedData(Values,MaxAge, minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsTransform=="Transformed"){

                  Values[,2] = scale(Values[,2],center = T, scale = T)
                  #PointValues = scale(PointValues,center = T, scale = T)

                }

                lines(Values[,1],Values[,2],col=Allcolor[RID], lwd=1)

                #points(as.numeric(data$TRACE[[TRACEsNames]]$rawData$depth),PointValues,col=Allcolor[RID], lwd=1, cex= 0.8)

              }
            }
          }
        }
      }

      #Plot Numbers

      for (i in 1:length(data$TRACE)){

        TRACEsNames = AllTRACENames[i]
        Values = data$TRACE[[TRACEsNames]][[TraceType]]
        RID = which(data$CoreList[,1]==TRACEsNames)

        if(length(RID)>0){
          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAgeForInterpolatedData(Values, MaxAge, minimumRowsAfterCutOutMaxAge)

              if(!is.null(Values)){

                if(VariantsTransform=="Transformed"){

                  Values[,2] = scale(Values[,2],center = T, scale = T)

                }

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
      }

      dev.off()

    }
  }

  TRACESoloPlot(TraceType = "JJA_mean")
  TRACESoloPlot(TraceType = "DJF_mean")
  TRACESoloPlot(TraceType = "hydrological_mean")

  ################################################################################
  ################################# StressPlot ###################################
  ################################################################################

  pdf("Stress.pdf",width=15,height=10)
  DiatomNames = ls(data$Diatom)

  LakeType = vector("character",length(DiatomNames))
  Stress = vector("numeric",length(DiatomNames))
  Names = DiatomNames

  counter = 0

  for (i in 1:length(DiatomNames)){

    counter = counter+1

    if (!sum(data$LakeData$CoreID==DiatomNames[counter])==0){

      StressValue = data$Diatom[[DiatomNames[counter]]]$nMDS$Stress

      if(!is.null(StressValue)){

      LakeType[counter] = data$LakeDat$LakeType[which(data$LakeData$CoreID==DiatomNames[counter])]
      Stress[counter] = StressValue

      }else{

        LakeType = LakeType[-counter]
        Stress = Stress[-counter]
        Names = Names[-counter]
        DiatomNames = DiatomNames[-counter]

        counter=counter-1

      }

    }else{

      LakeType = LakeType[-counter]
      Stress = Stress[-counter]
      Names = Names[-counter]
      DiatomNames = DiatomNames[-counter]

      counter=counter-1

    }
  }

  names(Stress) = Names
  names(LakeType) = Names

  Stress = Stress[!LakeType == "-"]
  LakeType = LakeType[!LakeType == "-"]

  sortStress = Stress[order(Stress,decreasing = F)]
  sortLakeType = LakeType[order(Stress,decreasing = F)]

  LakeFactor =  as.factor(LakeType)

  color = rainbow(length(levels(LakeFactor)))

  colorVector <- vector( "character" , length(LakeType)[1])

  for (c in 1:length(levels(LakeFactor))){

    colorVector[which(LakeFactor == levels(LakeFactor)[c])] = color[c]

  }

  par(mar = c(5, 9, 5, 5))
  par(mfrow = c(1,1))

  barplot(sortStress,horiz = TRUE, las = 1, col = colorVector, main = "Stress Plot",
          border = colorVector,
          space=0.3)

  legend(x = "bottom",
         legend = levels(LakeFactor),
         fill =color,
         cex = 1.3,
         bty = 'n',
         inset = c(-0.15,0))

  par(mar = c(5, 5, 5, 5))
  par(mfrow = c(1,1))

  dev.off()



  setwd(orginalWorkingDirectoryPath)

}








