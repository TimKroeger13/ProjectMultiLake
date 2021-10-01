#'Ordination
#'@description
#'Plots Ordinations Plots.
#'@param data List of data generates by the Multivar function.
#'@importFrom grDevices rainbow
#'@importFrom graphics lines points text legend
#'@importFrom stats na.omit
#'@importFrom grDevices dev.off pdf
#'@export
#'@return nothing.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

Ordination = function(data){

  orginalWorkingDirectoryPath=getwd()

  getwdTry <-  try(setwd(paste(getwd(),.Platform[2],"Output",sep="")),silent = TRUE)

  if(class(getwdTry) == "try-error"){

    dir.create(paste(getwd(),.Platform[2],"Output",sep=""))

    setwd(paste(getwd(),.Platform[2],"Output",sep=""))

  }

  CutOutMaxAge = function(DataVector, Cutage){

    FittingAges = as.numeric(row.names(DataVector))<Cutage

    if(sum(FittingAges)>0){

      DataVector=cbind(as.numeric(row.names(DataVector))[FittingAges],DataVector[FittingAges])

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




  ################################################################################
  ################################# Discription ##################################
  ################################################################################

  Allcolor = rainbow(length(data$Diatom))
  AllDiatomsNames = ls(data$Diatom)

  ################################################################################
  ##################################### MDS ######################################
  ################################################################################

  PlotsVariantsLoess = c("Normal","Loess")
  PlotsVariantsTransform = c("Non Transformed","Transformed")

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){
      if(VariantsTransform == "Non Transformed"){
        pdf(paste("MDS_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        #Plot Limits

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(rownames(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

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
             main=paste("MDS | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=0.5)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=0.5)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col="white",
                     cex=1.2)

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col=Allcolor[i])

              }
            }
          }
        }
        dev.off()
      }

      if(VariantsTransform == "Transformed"){
        pdf(paste("MDS_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(rownames(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

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
             main=paste("MDS | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=0.5)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$nMDS$Dim1)

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=0.5)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col="white",
                     cex=1.2)

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col=Allcolor[i])

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

  span=0.6

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){

      if(VariantsTransform == "Non Transformed"){
        pdf(paste("N2_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        #Plot Limits

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]]) #  $nMDS$Dim1

          if(!is.null(Values)){
            if(dim(Values)[1]>0){
              if(length(Values)>=minlegth){

                valueNames = as.numeric(rownames(Values))
                Values = as.numeric(Values)
                Values = matrix(Values, ncol = 1)
                rownames(Values) = valueNames

                Values = CutOutMaxAge(Values,20000)

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
             main=paste("N2 | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])

          if(!is.null(Values)){
            if(dim(Values)[1]>0){
              if(length(Values)>=minlegth){

                Values = CutOutMaxAge(Values,20000)

                if(!is.null(Values)){

                  points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                  if(VariantsLoess == "Loess"){

                    #Loess
                    ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                    Values = cbind(Values[,1],ValuesLoess$fitted)

                  }

                  lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

                }
              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])

          if(!is.null(Values)){
            if(dim(Values)[1]>0){
              if(length(Values)>=minlegth){

                Values = CutOutMaxAge(Values,20000)

                if(!is.null(Values)){

                  if(VariantsLoess == "Loess"){

                    #Loess
                    ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                    Values = cbind(Values[,1],ValuesLoess$fitted)

                  }

                  distance = 120

                  if(i>=10){distance = 200}

                  text(Values[dim(Values)[1],1]+distance,
                       Values[dim(Values)[1],2],
                       label=i,
                       col="white",
                       cex=1.2)

                  text(Values[dim(Values)[1],1]+distance,
                       Values[dim(Values)[1],2],
                       label=i,
                       col=Allcolor[i])

                }
              }
            }
          }
        }
        dev.off()
      }

      if(VariantsTransform == "Transformed"){
        pdf(paste("N2_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])

          if(!is.null(Values)){
            if(dim(Values)[1]>0){
              if(length(Values)>=minlegth){

                valueNames = as.numeric(rownames(Values))
                Values = as.numeric(Values)
                Values = matrix(Values, ncol = 1)
                rownames(Values) = valueNames

                Values = CutOutMaxAge(Values,20000)

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
             main=paste("N2 | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])

          if(!is.null(Values)){
            if(dim(Values)[1]>0){
              if(length(Values)>=minlegth){

                Values = CutOutMaxAge(Values,20000)

                if(!is.null(Values)){

                  Values[,2] = scale(Values[,2],center = T, scale = T)

                  points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                  if(VariantsLoess == "Loess"){

                    #Loess
                    ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                    Values = cbind(Values[,1],ValuesLoess$fitted)

                  }

                  lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

                }
              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Diatom[[DiatomsNames]]$Species_richness[[nameNormal]])

          if(!is.null(Values)){
            if(dim(Values)[1]>0){
              if(length(Values)>=minlegth){

                Values = CutOutMaxAge(Values,20000)

                if(!is.null(Values)){

                  Values[,2] = scale(Values[,2],center = T, scale = T)

                  if(VariantsLoess == "Loess"){

                    #Loess
                    ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                    Values = cbind(Values[,1],ValuesLoess$fitted)

                  }

                  distance = 120

                  if(i>=10){distance = 200}

                  text(Values[dim(Values)[1],1]+distance,
                       Values[dim(Values)[1],2],
                       label=i,
                       col="white",
                       cex=1.2)

                  text(Values[dim(Values)[1],1]+distance,
                       Values[dim(Values)[1],2],
                       label=i,
                       col=Allcolor[i])

                }
              }
            }
          }
        }
        dev.off()
      }


    }
  }



  ################################################################################
  ##################################### Clima ####################################
  ################################################################################

  PlotsVariantsLoess = c("Normal","Loess")
  PlotsVariantsTransform = c("Non Transformed","Transformed")

  span=0.1

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

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

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

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))
              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=i,
                     col="white",
                     cex=1.2)

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=i,
                     col=Allcolor[i])

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

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))

              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

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

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))

              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Clima$annual[[DiatomsNames]])

          if(!is.null(Values)){
            if(dim(Values)[1]>0){

              valueNames = as.numeric(names(Values))

              Values = as.numeric(Values)
              Values = matrix(Values, ncol = 1)
              rownames(Values) = valueNames

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                Values[,2] = scale(Values[,2],center = T, scale = T)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=i,
                     col="white",
                     cex=1.2)

                text(Values[1,1]+distance,
                     Values[1,2],
                     label=i,
                     col=Allcolor[i])

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

  span=0.3

  for (VariantsLoess in PlotsVariantsLoess){
    for (VariantsTransform in PlotsVariantsTransform){

      if(VariantsTransform == "Non Transformed"){
        pdf(paste("TOC_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        #Plot Limits

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC) #  $Species_richness[[nameNormal]]

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

              Values = CutOutMaxAge(Values,20000)

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
             main=paste("TOC | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){ #length(data$Diatom)

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)

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

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)

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

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col="white",
                     cex=1.2)

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col=Allcolor[i])

              }
            }
          }
        }
        dev.off()
      }

      if(VariantsTransform == "Transformed"){
        pdf(paste("TOC_",VariantsLoess, "_",VariantsTransform,".pdf",sep=""),width=15,height=10)

        Xmax=0
        Ymax=0
        Xmin=Inf
        Ymin=Inf

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)

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

              Values = CutOutMaxAge(Values,20000)

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
             main=paste("TOC | ",VariantsLoess, " | ",VariantsTransform,sep="")
        )

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)

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

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                Values[,2] =  scale(Values[,2],center = T, scale = T)

                points(Values[,1],Values[,2],col=Allcolor[i], lwd=1, cex= 0.8)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                lines(Values[,1],Values[,2],col=Allcolor[i], lwd=1)

              }
            }
          }
        }

        #Plot Numbers

        for (i in 1:length(data$Diatom)){

          DiatomsNames = AllDiatomsNames[i]
          Values = na.omit(data$Carbon[[DiatomsNames]]$rawData$TOC)

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

              Values = CutOutMaxAge(Values,20000)

              if(!is.null(Values)){

                Values[,2] =  scale(Values[,2],center = T, scale = T)

                if(VariantsLoess == "Loess"){

                  #Loess
                  ValuesLoess=loess(Values[,2] ~ Values[,1], span=span)
                  Values = cbind(Values[,1],ValuesLoess$fitted)

                }

                distance = 120

                if(i>=10){distance = 200}

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col="white",
                     cex=1.2)

                text(Values[dim(Values)[1],1]+distance,
                     Values[dim(Values)[1],2],
                     label=i,
                     col=Allcolor[i])

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

  LegendsSplit = length(AllDiatomsNames)/3

  SplitValue1 = ceiling(LegendsSplit)
  SplitValue2 = ceiling(LegendsSplit*2)
  SplitValue3 = length(AllDiatomsNames)

  AllDiatomsNamesPart1 = AllDiatomsNames[1:SplitValue1]
  AllDiatomsNamesPart2 = AllDiatomsNames[(SplitValue1+1):SplitValue2]
  AllDiatomsNamesPart3 = AllDiatomsNames[(SplitValue2+1):SplitValue3]

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

  thickness = 0.9

  legend(x = "left", legend = paste(1:SplitValue1,"-",AllDiatomsNamesPart1), col = Allcolor[1:SplitValue1], lty = 1, lwd = 3, text.width=0.3, cex = thickness, bty="n") # bty="n"

  legend(x = "center", legend = paste((SplitValue1+1):SplitValue2,"-",AllDiatomsNamesPart2), col = Allcolor[(SplitValue1+1):SplitValue2], lty = 1, lwd = 3, text.width=0.3, cex = thickness, bty="n") # bty="n"

  legend(x = "right", legend = paste((SplitValue2+1):SplitValue3,"-",AllDiatomsNamesPart3), col = Allcolor[(SplitValue2+1):SplitValue3], lty = 1, lwd = 3, text.width=0.3, cex = thickness, bty="n") # bty="n"

  dev.off()


  setwd(orginalWorkingDirectoryPath)

}








