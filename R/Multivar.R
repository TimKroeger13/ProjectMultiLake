#'Multivar
#'@description
#'Calculates several parameters from deatom data.
#'@param data List of data generates by the MultiExcelLoader function.
#'@param method Method for calculation Dissimilarity Indices for Community Ecologists.
#'@param standardize Method for data standardisation. Can be nothing "" ir sqaureroot transformation "sqrt".
#'@param percentFilterWeight Value how much percent a single species must relevant at minimum from the dataset.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@import vegan SRS
#'@importFrom stats prcomp loess median predict qt quantile approx
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

Multivar = function(data,method="bray",standardize=c("","sqrt"),percentFilterWeight=0,allLoessSpans=0.8){

  DeleteRowWithoutTimestamp <- function(data){

    return(data[!is.na(data[,1]),])

  }

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

  DeleteNulRows <- function(data,age){

    i=0

    while (!i==dim(data)[1]){

      i=i+1

      if(sum(data[i,])==0){

        data=data[-i,]
        age=age[-i]

        i=i-1

      }

    }

    output=list()
    output$data=data
    output$age=age

    return(output)

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

      Tempdata = Tempdata[CountedValves >= SRS_Value,]
      TempdataForCalculation = TempdataForCalculation[CountedValves >= SRS_Value,]

      TempTurn = data.frame(t(TempdataForCalculation))
      TempTurn = SRS(TempTurn, Cmin = SRS_Value, set_seed = TRUE, seed = 1)
      SRSTempdataForCalculation = data.frame(t(TempTurn))

      output = cbind(Tempdata[,1:3],SRSTempdataForCalculation)
      colnames(output) = colnames(Tempdata)
      output[,2] = SRS_Value

      data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("SRS_data")]] = output

      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          i,"/",length(ls(data[["Diatom"]]))," Rasampel Taxa Data",sep="")

    }

    return(data)

  }

  FilterPercentData <- function(PercentData,percentFilterWeight){

    TempDataForCalculation=PercentData[,4:dim(PercentData)[2]]

    BetterEqualMedian=TempDataForCalculation

    BetterEqualMedian[]=NA

    for(x in 1:dim(TempDataForCalculation)[1]){

      if(sum(TempDataForCalculation[x,]>0)==0){

        BetterEqualMedian[x,]=FALSE

      }else{

      BetterEqualMedian[x,]=TempDataForCalculation[x,]>=median(TempDataForCalculation[x,][TempDataForCalculation[x,]>0])

      }
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

       Tempdata=data[["Diatom"]][[k]][["rawData"]]
       TempdataForCalculation=Tempdata[,4:dim(Tempdata)[2]]

       PercentData=Tempdata
       PercentData[,4:dim(PercentData)[2]]=NA

      for(i in 1:dim(Tempdata)[1]){

        if(sum(TempdataForCalculation[i,],na.rm = T)==0){

          PercentData[i,4:dim(PercentData)[2]]=0

        }else{

        PercentData[i,4:dim(PercentData)[2]]=TempdataForCalculation[i,]/sum(TempdataForCalculation[i,],na.rm = T)*100

        }
      }

       PercentData = FilterPercentData(PercentData = PercentData,percentFilterWeight = percentFilterWeight) #Filter Collums

      if(standardize[1]=="sqrt"){

        PercentData[,4:dim(PercentData)[2]] =  sqrt(PercentData[,4:dim(PercentData)[2]])

      }

       data[["Diatom"]][[names(data[["Diatom"]])[[k]]]][[paste("StandadizedData")]]=PercentData

       #Printer
       cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
           k,"/",length(ls(data[["Diatom"]]))," Standadize Diatom Data",sep="")

    }

    return(data)
  }

  data=StandadizeData(data, standardize, percentFilterWeight)
  data=SrsFilter(data)





  #Main Loop
  for (i in 1:length(ls(data[["Diatom"]]))){
  #for (i in 1:10){

    #MDS

    MDSData=data[["Diatom"]][[i]][["StandadizedData"]]

    MDSData=DeleteRowWithoutTimestamp(MDSData)

    MDSAges=MDSData[,1]
    MDSData=MDSData[,4:dim(MDSData)[2]]
    MDSData=DeleteNullCollums(MDSData)
    RowDeletedData=DeleteNulRows(MDSData,MDSAges)

    MDSAges=RowDeletedData$age
    MDSData=RowDeletedData$data


    dissimilarityData=dissimilarityIndex(MDSData,method)

    monoMDSData=suppressWarnings(metaMDS(dissimilarityData, k=2, autotransform = FALSE))

    ExternalCalibrator = data[["GlobalInsolation"]]

    ExternalCalibrator = approx(y = as.numeric(ExternalCalibrator),
                                x = as.numeric(row.names(ExternalCalibrator)),
                                xout = MDSAges, method = "linear")

    MDSrotation = MDSrotate(monoMDSData, vec = ExternalCalibrator$y)


    #MDSrotation=monoMDSData #<--------------------------------------------------- Just for Testing delte laster!

    MDS1=MDSrotation$points[, 1]
    MDS2=MDSrotation$points[, 2]

    MDS1=matrix(MDS1,ncol = 1)
    row.names(MDS1)=as.character(MDSAges)
    colnames(MDS1)="nMDS_Dim1"
    data[["Diatom"]][[i]][["nMDS"]][["Dim1"]]=MDS1


    MDS2=matrix(MDS2,ncol = 1)
    row.names(MDS2)=as.character(MDSAges)
    colnames(MDS2)="nMDS_Dim2"
    data[["Diatom"]][[i]][["nMDS"]][["Dim2"]]=MDS2

    #Stress
    data[["Diatom"]][[i]][["nMDS"]][["Stress"]]=MDSrotation$stress

  }


  #Main Loop
  for (i in 1:length(ls(data[["Diatom"]]))){
  #for (i in 1:4){

    Tempdata=data[["Diatom"]][[i]][["SRS_data"]]
    TempdataForCalculation=data[["Diatom"]][[i]][["SRS_data"]][,4:dim(data[["Diatom"]][[i]][["SRS_data"]])[2]]
    TempdataForCalculation=round(TempdataForCalculation)

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

    data[["Diatom"]][[i]][["Species_richness"]][[paste("richness")]]=richness
    data[["Diatom"]][[i]][["Species_richness"]][[paste("shannon")]]=shannon
    data[["Diatom"]][[i]][["Species_richness"]][[paste("invsimpson")]]=invsimpson


    # Loess Predictor
    for (r in c("richness","shannon","invsimpson")){

      y=as.numeric(row.names(data[["Diatom"]][[i]][["Species_richness"]][[r]]))
      x=data[["Diatom"]][[i]][["Species_richness"]][[r]]

      loessValues=predict(loess(x ~ y, span=allLoessSpans), se=T,newdata = as.numeric(y)) #Critical Span value

      LoessMean = loessValues$fit
      LoessConfUp = loessValues$fit + qt(1-(0.05/2),loessValues$df)*loessValues$se.fit
      LoessConfDown = loessValues$fit - qt(1-(0.05/2),loessValues$df)*loessValues$se.fit

      LoessOut=cbind(LoessMean,LoessConfUp,LoessConfDown)
      row.names(LoessOut)=Tempdata[,1]

      data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][["Species_richness"]][[paste("Loess_",r,sep = "")]]=LoessOut

    }

    #DCA
    dca=suppressWarnings(decorana(data[["Diatom"]][[i]]$rawData[,4:dim(data[["Diatom"]][[i]]$rawData)[2]]))

    data[["Diatom"]][[i]][["DCA_Gradient_Length"]]=round(max(dca$rproj[,1:4]),digits = 3)


    #Beta diversity Turn Over

    Tempdata=data[["Diatom"]][[i]][["SRS_data"]]
    TempdataForCalculation=data[["Diatom"]][[i]][["SRS_data"]][,4:dim(data[["Diatom"]][[i]][["SRS_data"]])[2]]
    TempdataForCalculation=round(TempdataForCalculation)

    BetaDiversity=matrix(NA,ncol = 1,nrow = dim(TempdataForCalculation)[1]-1)
    rownames(BetaDiversity)=Tempdata[1:(dim(Tempdata)[1]-1),1]
    colnames(BetaDiversity)="Hierarchical beta diversity"

    for(k in (dim(TempdataForCalculation)[1]):2){

      BetaDiversity[k-1]=vegdist(TempdataForCalculation[(k-1):k,],method = "bray")

    }

    data[["Diatom"]][[i]][["BetaDiversity"]][["OffsetBy_1"]]=BetaDiversity

    for(x in 2:3){

      #Multipel Beta diversity Turn Over

      Tempdata=data[["Diatom"]][[i]][["SRS_data"]]
      TempdataForCalculation=data[["Diatom"]][[i]][["SRS_data"]][,4:dim(data[["Diatom"]][[i]][["SRS_data"]])[2]]
      TempdataForCalculation=round(TempdataForCalculation)

      SectionIntervall = x

      BetaDiversity=matrix(NA,ncol = 1,nrow = dim(TempdataForCalculation)[1]-SectionIntervall)
      rownames(BetaDiversity)=Tempdata[1:(dim(Tempdata)[1]-SectionIntervall),1]
      colnames(BetaDiversity)="Hierarchical beta diversity"

      comparisonIntervall=matrix(NA, nrow = 1, ncol = dim(TempdataForCalculation)[2])

      for(k in dim(TempdataForCalculation)[1]:SectionIntervall){

        for (z in 1:dim(TempdataForCalculation)[2]){

          comparisonIntervall[z]=mean(TempdataForCalculation[(k-SectionIntervall+1):k,z])

        }

        colnames(comparisonIntervall)=colnames(TempdataForCalculation)

        BetaDiversity[k-SectionIntervall]=vegdist(rbind(TempdataForCalculation[k-SectionIntervall,],comparisonIntervall),method = "bray")

      }

      data[["Diatom"]][[i]][["BetaDiversity"]][[paste("OffsetBy_",SectionIntervall,sep="")]]=BetaDiversity

    }

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",length(ls(data[["Diatom"]]))," Multivar",sep="")


  }

  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Done",sep="")

  return(data)

}

















