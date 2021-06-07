#'Multivar
#'@description
#'Calculates several parameters from deatom data.
#'@param data List of data generates by the MultiExcelLoader function.
#'@param method Method for calculation Dissimilarity Indices for Community Ecologists.
#'@param standardize Method for data standardisation. Can be nothing "" ir sqaureroot transformation "sqrt".
#'@param percentFilterWeight Value how much percent a single species must declare at minimum from the dataset.
#'@import vegan SRS
#'@importFrom stats prcomp loess median predict qt quantile
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

Multivar = function(data,method="bray",standardize=c("","sqrt"),percentFilterWeight=0){

  DeleteNullCollums <- function(data){

    i=0

    while (!i==dim(data)[2]){

      i=i+1

      if(sum(data[,i])==0){

        data=data[,-i]

        i=i-1

      }

    }

    return(data)

  }

  dissimilarityIndex <- function(data,method){

    data=vegdist(data,method=method,na.rm = T)

    return(data)

  }

  SrsFilter = function(data){

    allCounts=NULL

    for (i in 1:length(ls(data[["Diatom"]]))){

      allCounts=c(allCounts,data[["Diatom"]][[i]][["rawData"]][,2])

    }

    SRS_Value = round(quantile(allCounts,probs = c(0.05)))

    for (i in 1:length(ls(data[["Diatom"]]))){

      Tempdata=data[["Diatom"]][[i]][["rawData"]]
      TempdataForCalculation=data[["Diatom"]][[i]][["rawData"]][,4:dim(data[["Diatom"]][[i]][["rawData"]])[2]]
      TempdataForCalculation=round(TempdataForCalculation)

      CountedValves=NULL

      for (k in 1:dim(TempdataForCalculation)[1]){

        CountedValves[k]=sum(TempdataForCalculation[k,])

      }

      Tempdata=Tempdata[CountedValves>=SRS_Value,]
      TempdataForCalculation=TempdataForCalculation[CountedValves>=SRS_Value,]

      TempTurn <- data.frame(t(TempdataForCalculation))
      TempTurn=SRS(TempTurn, Cmin=SRS_Value, set_seed = TRUE, seed = 1)
      SRSTempdataForCalculation <- data.frame(t(TempTurn))

      output=cbind(Tempdata[,1:3],SRSTempdataForCalculation)
      colnames(output) = colnames(Tempdata)
      output[,2]=SRS_Value

      data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("SRS_data")]]=output

    }

    return(data)

  }

  FilterPercentData <- function(PercentData,percentFilterWeight){

    TempDataForCalculation=PercentData[,4:dim(PercentData)[2]]

    BetterEqualMedian=TempDataForCalculation

    BetterEqualMedian[]=NA

    for(x in 1:dim(TempDataForCalculation)[1]){

      BetterEqualMedian[x,]=TempDataForCalculation[x,]>=median(TempDataForCalculation[x,][TempDataForCalculation[x,]>0])

    }

    counter=0

    while (dim(TempDataForCalculation)[2]>counter){

      counter=counter+1

      if(sum(BetterEqualMedian[,counter])<dim(BetterEqualMedian)[1]*percentFilterWeight/100){

        TempDataForCalculation=TempDataForCalculation[-counter]
        BetterEqualMedian=BetterEqualMedian[-counter]

        counter=counter-1

      }
    }

    output=PercentData[,1:3]

    output=cbind(output,TempDataForCalculation)

    return(output)

  }


  StandadizeData <- function(data, standardize, percentFilterWeight){

     for (k in 1:length(ls(data[["Diatom"]]))){

       TempData=data[["Diatom"]][[k]][["rawData"]]
       TempDataForCalculation=TempData[,4:dim(TempData)[2]]

       PercentData=TempData
       PercentData[,3:dim(PercentData)[2]]=NA

      for(i in 1:dim(TempData)[1]){

        PercentData[i,4:dim(PercentData)[2]]=TempDataForCalculation[i,]/sum(TempDataForCalculation[i,],na.rm = T)*100

      }

       PercentData = FilterPercentData(PercentData = PercentData,percentFilterWeight = percentFilterWeight) #Filter Collums

      if(standardize[1]=="sqrt"){

        PercentData[,4:dim(PercentData)[2]] =  sqrt(PercentData[,4:dim(PercentData)[2]])

      }

       data[["Diatom"]][[names(data[["Diatom"]])[[k]]]][[paste("StandadizedData")]]=PercentData

       #Printer
       cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
           k,"/",length(ls(data[["Diatom"]])),sep="")

    }

    return(data)
  }

  data=StandadizeData(data, standardize, percentFilterWeight)
  data=SrsFilter(data)

  #Uses SRS Data
  for (i in 1:length(ls(data[["Diatom"]]))){

    Tempdata=data[["Diatom"]][[i]][["SRS_data"]]
    TempdataForCalculation=data[["Diatom"]][[i]][["SRS_data"]][,4:dim(data[["Diatom"]][[i]][["SRS_data"]])[2]]
    TempdataForCalculation=round(TempdataForCalculation)
    SqrtTempdataForCalculation=sqrt(TempdataForCalculation)

    #Richness

    richness = rarefy(TempdataForCalculation, sample = sum(TempdataForCalculation[1,]), se=T)
    shannon = diversity(TempdataForCalculation, index = "shannon")
    invsimpson = diversity(TempdataForCalculation, index = "invsimpson")

    richness=cbind(richness[1,])
    colnames(richness)=c("richness")
    row.names(richness)=Tempdata[,1]
    shannon=cbind(shannon)
    row.names(shannon)=Tempdata[,1]
    invsimpson=cbind(invsimpson)
    row.names(invsimpson)=Tempdata[,1]

    data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("richness")]]=richness
    data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("shannon")]]=shannon
    data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("invsimpson")]]=invsimpson

    for (r in c("richness","shannon","invsimpson")){

      y=as.numeric(row.names(data[["Diatom"]][[i]][[r]]))
      x=data[["Diatom"]][[i]][[r]]

      loessValues=predict(loess(x ~ y, span=0.5), se=T,newdata = as.numeric(y))

      LoessMean = loessValues$fit
      LoessConfUp = loessValues$fit + qt(1-(0.05/2),loessValues$df)*loessValues$se
      LoessConfDown = loessValues$fit - qt(1-(0.05/2),loessValues$df)*loessValues$se

      LoessOut=cbind(LoessMean,LoessConfUp,LoessConfDown)
      row.names(LoessOut)=Tempdata[,1]

      data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("Loess_",r,sep = "")]]=LoessOut

    }


    #data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("SRS_data")]]=output

'
    #PCA

    FiltertSRSTempdataForCalculation=DeleteNullCollums(SRSTempdataForCalculation)

    PCA=prcomp(sqrt(FiltertSRSTempdataForCalculation), center = TRUE, scale = F)

    PCA1=scores(PCA,choices = 1,display=c("sites"))
    PCA2=scores(PCA,choices = 2,display=c("sites"))
    PCA3=scores(PCA,choices = 3,display=c("sites"))
    row.names(PCA1)=as.character(Tempdata[,1])
    row.names(PCA2)=as.character(Tempdata[,1])
    row.names(PCA3)=as.character(Tempdata[,1])

    data[[names(data)[["Diatom"]][i]]][["PCA_1_SiteScores"]]=PCA1
    data[[names(data)[["Diatom"]][i]]][["PCA_2_SiteScores"]]=PCA2
    data[[names(data)[["Diatom"]][i]]][["PCA_3_SiteScores"]]=PCA3

    PCA.var = PCA$sdev^2
    PCA.ve = PCA.var/sum(PCA.var) # Variance

    vari1=matrix(c(PCA.ve[1],rep(NA,length(Tempdata[,1])-1)),ncol = 1)
    colnames(vari1)="PCA_1_Variance"
    data[[names(data)[["Diatom"]][i]]][["PCA_1_Variance"]]=vari1

    vari2=matrix(c(PCA.ve[2],rep(NA,length(Tempdata[,1])-1)),ncol = 1)
    colnames(vari2)="PCA_2_Variance"
    data[[names(data)[["Diatom"]][i]]][["PCA_2_Variance"]]=vari2

    vari3=matrix(c(PCA.ve[3],rep(NA,length(Tempdata[,1])-1)),ncol = 1)
    colnames(vari3)="PCA_3_Variance"
    data[[names(data)[["Diatom"]][i]]][["PCA_3_Variance"]]=vari3

    #MDS

    TempTurn=SRSTempdataForCalculation
    TempTurn=DeleteNullCollums(TempTurn)

    dissimilarityData=dissimilarityIndex(TempTurn,method)

    fit <- cmdscale(dissimilarityData, eig = T, k = 3)





    MDS1=fit$points[, 1]
    MDS2=fit$points[, 2]
    MDS3=fit$points[, 3]

    MDS1=matrix(MDS1,ncol = 1)
    row.names(MDS1)=as.character(Tempdata[,1])
    colnames(MDS1)="MDS_1_SiteScores"
    data[[names(data)[["Diatom"]][i]]][["MDS_1_SiteScores"]]=MDS1

    MDS2=matrix(MDS2,ncol = 1)
    row.names(MDS2)=as.character(Tempdata[,1])
    colnames(MDS2)="MDS_2_SiteScores"
    data[[names(data)[["Diatom"]][i]]][["MDS_2_SiteScores"]]=MDS2

    MDS3=matrix(MDS3,ncol = 1)
    row.names(MDS3)=as.character(Tempdata[,1])
    colnames(MDS3)="MDS_3_SiteScores"
    data[[names(data)[["Diatom"]][i]]][["MDS_3_SiteScores"]]=MDS3


    vari1=matrix(c(round(fit$eig*100/sum(fit$eig),1)[1],rep(NA,length(Tempdata[,1])-1)),ncol = 1)
    colnames(vari1)="MDS_1_Variance"
    data[[names(data)[["Diatom"]][i]]][["MDS_1_Variance"]]=vari1

    vari2=matrix(c(round(fit$eig*100/sum(fit$eig),1)[2],rep(NA,length(Tempdata[,1])-1)),ncol = 1)
    colnames(vari2)="MDS_2_Variance"
    data[[names(data)[["Diatom"]][i]]][["MDS_2_Variance"]]=vari2

    vari3=matrix(c(round(fit$eig*100/sum(fit$eig),1)[3],rep(NA,length(Tempdata[,1])-1)),ncol = 1)
    colnames(vari3)="MDS_3_Variance"
    data[[names(data)[["Diatom"]][i]]][["MDS_3_Variance"]]=vari3
'




  }

  return(data)

}

















