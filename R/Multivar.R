#'Multivar
#'@description
#'Calculates several parameters from deatom data.
#'@param data List of data generates by the MultiExcelLoader function.
#'@param method Method for calculation Dissimilarity Indices for Community Ecologists.
#'@param standardize Method for data standardisation. Can be nothing "" ir sqaureroot transformation "sqrt".
#'@param percentFilterWeight Value how much percent a single species must relevant at minimum from the dataset.
#'@param allLoessSpans span value for all Loess calculations made by Multivar.
#'@param minimumRowsAfterFiltering Value for the minimum rows after filtering.
#'@param InterpolationLoessSpans span value for all interpolation Loess calculations made by Multivar, where more values are given.
#'@param ResampleQuantileValue Value that indicates by which percent of quantiles all the Diatom counts should be cut.
#'@import vegan SRS
#'@importFrom stats prcomp loess median predict qt quantile approx
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger.
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

Multivar = function(data,method="bray",standardize=c("","sqrt"),percentFilterWeight=0,allLoessSpans=0.8,minimumRowsAfterFiltering = 0, InterpolationLoessSpans = 0.8, ResampleQuantileValue = 0.05){

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

  filterDataForMinimumRows = function(data,areas=c("Carbon","Element","Diatom")){

    for (i in areas){

      listNumber = 0

      for (k in 1:length(ls(data[[i]]))){

        listNumber = listNumber+1

        Tempdata = data[[i]][[listNumber]][["rawData"]]

        if(dim(Tempdata)[1]<=minimumRowsAfterFiltering){

          data[[i]]=data[[i]][-listNumber]

          listNumber=listNumber-1

        }
      }
    }

    return(data)

  }

  SrsFilter = function(data){

    allCounts=NULL

    for (i in 1:length(ls(data[["Diatom"]]))){

      allCounts=c(allCounts,data[["Diatom"]][[i]][["rawData"]][,2])

    }

    SRS_Value = round(quantile(allCounts,probs = c(ResampleQuantileValue)))

    for (i in 1:length(ls(data[["Diatom"]]))){

      Tempdata=data[["Diatom"]][[i]][["rawData"]]
      TempdataForCalculation=data[["Diatom"]][[i]][["rawData"]][,4:dim(data[["Diatom"]][[i]][["rawData"]])[2]]
      TempdataForCalculation=round(TempdataForCalculation)

      CountedValves=NULL

      for (k in 1:dim(TempdataForCalculation)[1]){

        CountedValves[k]=sum(TempdataForCalculation[k,])

      }

      if(sum(CountedValves >= SRS_Value)>0){

        Tempdata = Tempdata[CountedValves >= SRS_Value,]
        TempdataForCalculation = TempdataForCalculation[CountedValves >= SRS_Value,]

        TempTurn = data.frame(t(TempdataForCalculation))
        TempTurn = SRS(TempTurn, Cmin = SRS_Value, set_seed = TRUE, seed = 1)
        SRSTempdataForCalculation = data.frame(t(TempTurn))

        output = cbind(Tempdata[,1:3],SRSTempdataForCalculation)
        colnames(output) = colnames(Tempdata)
        output[,2] = SRS_Value

        if(dim(output)[1]>=minimumRowsAfterFiltering){

          data[["Diatom"]][[names(data[["Diatom"]])[[i]]]][[paste("SRS_data")]] = output

        }
      }

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

  data = filterDataForMinimumRows(data)

  data = SrsFilter(data)

  data = CutOutPionierPhase(data = data)


  ImportVersions = c("SRS_data","Cut_SRS_data")
  SRS_Names = c("Dim1","Dim2","Stress")
  Cut_SRS_Names = c("Cut_Dim1","Cut_Dim2","Cut_Stress")
  namesInUse = NULL

  for (IV in ImportVersions){

    if(IV == "SRS_data"){namesInUse = SRS_Names}
    if(IV == "Cut_SRS_data"){namesInUse = Cut_SRS_Names}

    #MDS Loop
    for (i in 1:length(ls(data[["Diatom"]]))){

      MDSData=data[["Diatom"]][[i]][[IV]]

      if(!is.null(MDSData)){
        if(dim(MDSData)[1]>2){

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


        #MDSrotation=monoMDSData #<--------------------------------------------------- Just for Testing delete later!

        MDS1=MDSrotation$points[, 1]
        MDS2=MDSrotation$points[, 2]

        MDS1=matrix(MDS1,ncol = 1)
        row.names(MDS1)=as.character(MDSAges)
        colnames(MDS1)="nMDS_Dim1"
        data[["Diatom"]][[i]][["nMDS"]][[namesInUse[1]]]=MDS1


        MDS2=matrix(MDS2,ncol = 1)
        row.names(MDS2)=as.character(MDSAges)
        colnames(MDS2)="nMDS_Dim2"
        data[["Diatom"]][[i]][["nMDS"]][[namesInUse[2]]]=MDS2

        #Stress
        data[["Diatom"]][[i]][["nMDS"]][[namesInUse[3]]]=MDSrotation$stress

      }
      }
    }
  }


  ImportVersions = c("SRS_data","Cut_SRS_data")
  SRS_Names = c("richness","shannon","invsimpson")
  Cut_SRS_Names = c("Cut_richness","Cut_shannon","Cut_invsimpson")
  namesInUse = NULL

  for (IV in ImportVersions){

    if(IV == "SRS_data"){namesInUse = SRS_Names}
    if(IV == "Cut_SRS_data"){namesInUse = Cut_SRS_Names}

    #Richness Loop
    for (i in 1:length(ls(data[["Diatom"]]))){

      Tempdata=data[["Diatom"]][[i]][[IV]]

      TempdataForCalculation=try(data[["Diatom"]][[i]][[IV]][,4:dim(data[["Diatom"]][[i]][[IV]])[2]],silent = T)

      if(!class(TempdataForCalculation)[1] == "try-error"){
        if(dim(TempdataForCalculation)[1]>2){

        TempdataForCalculation=data[["Diatom"]][[i]][[IV]][,4:dim(data[["Diatom"]][[i]][[IV]])[2]]
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

        data[["Diatom"]][[i]][["Species_richness"]][[namesInUse[1]]]=richness
        data[["Diatom"]][[i]][["Species_richness"]][[namesInUse[2]]]=shannon
        data[["Diatom"]][[i]][["Species_richness"]][[namesInUse[3]]]=invsimpson


        # Loess Predictor
        for (r in c(namesInUse[1],namesInUse[2],namesInUse[3])){

          y=as.numeric(row.names(data[["Diatom"]][[i]][["Species_richness"]][[r]]))
          x=data[["Diatom"]][[i]][["Species_richness"]][[r]]

          #GuessLoess

          span = GuessLoess(intervall = y,values = x,overspan = 100)

          loessValues=predict(loess(x ~ y, span=span), se=T,newdata = as.numeric(y)) #Critical Span value

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

      }
    }

      #Printer
      cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
          i,"/",length(ls(data[["Diatom"]]))," calculating species richness",sep="")

    }
  }

  #ROC

  data = rateofChange(data = data, intervallBy = 100, minimumRowsAfterInterpolating = minimumRowsAfterFiltering, method = method,
                      MinAgeIntervall = 1, NonNegative = TRUE, Importname = "SRS_data", Exportname = "RoC")
  data = rateofChange(data = data, intervallBy = 100, minimumRowsAfterInterpolating = minimumRowsAfterFiltering, method = method,
                      MinAgeIntervall = 1, NonNegative = TRUE, Importname = "Cut_SRS_data", Exportname = "Cut_RoC")

  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "RoC", Exportname = "RocMatrix", CreateVector = F, TransformAllData = FALSE, ValueCantBeSamlerThanZero = TRUE)
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_RoC", Exportname = "Cut_RocMatrix", CreateVector = T, TransformAllData = FALSE,
                    vectorName = "Vector_RocMatrix", ValueCantBeSamlerThanZero = TRUE)
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_RoC", Exportname = "Cut_RocMatrix_Transfrom", CreateVector = T, TransformAllData = TRUE,
                    vectorName = "Vector_RocMatrix_transform", ValueCantBeSamlerThanZero = FALSE)

  #Eveness

  data = evenness(data = data, intervallBy = 100,
                  NonNegative = TRUE, Importname = "SRS_data", Exportname = "evenness")
  data = evenness(data = data, intervallBy = 100,
                  NonNegative = TRUE, Importname = "Cut_SRS_data", Exportname = "Cut_evenness")

  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "evenness", Exportname = "EvennessMatrix", CreateVector = F, TransformAllData = FALSE, ValueCantBeSamlerThanZero = TRUE)
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_evenness", Exportname = "Cut_EvennessMatrix", CreateVector = T, TransformAllData = FALSE, ValueCantBeSamlerThanZero = TRUE,
                    vectorName = "Vector_Evenness")
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_evenness", Exportname = "Cut_EvennessMatrix_Transfrom", CreateVector = T, TransformAllData = TRUE, ValueCantBeSamlerThanZero = FALSE,
                    vectorName = "Vector_Evenness_transform")


  #speciesRichness

  data = CalculateBaseValues(data = data, intervallBy = 100,
                         NonNegative = TRUE, Importname1 = "Species_richness", Importname2 = "richness", Exportname = "speciesRichnessData")
  data = CalculateBaseValues(data = data, intervallBy = 100,
                         NonNegative = TRUE, Importname1 = "Species_richness", Importname2 = "Cut_richness", Exportname = "Cut_speciesRichnessData")

  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "speciesRichnessData", Exportname = "SpeciesRichnessMatrix", CreateVector = F, TransformAllData = FALSE, ValueCantBeSamlerThanZero = TRUE)
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_speciesRichnessData", Exportname = "Cut_SpeciesRichnessMatrix", CreateVector = T, TransformAllData = FALSE, ValueCantBeSamlerThanZero = TRUE,
                    vectorName = "Vector_SpeciesRichness")
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_speciesRichnessData", Exportname = "Cut_SpeciesRichnessMatrix_Transform", CreateVector = T, TransformAllData = TRUE, ValueCantBeSamlerThanZero = FALSE,
                    vectorName = "Vector_SpeciesRichness_transform")

  #inverseSimpson

  data = CalculateBaseValues(data = data, intervallBy = 100,
                             NonNegative = TRUE, Importname1 = "Species_richness", Importname2 = "invsimpson", Exportname = "inverseSimpsionData")
  data = CalculateBaseValues(data = data, intervallBy = 100,
                             NonNegative = TRUE, Importname1 = "Species_richness", Importname2 = "Cut_invsimpson", Exportname = "Cut_inverseSimpsionData")

  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "inverseSimpsionData", Exportname = "InverseSimpsionMatrix", CreateVector = F, TransformAllData = FALSE, ValueCantBeSamlerThanZero = TRUE)
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_inverseSimpsionData", Exportname = "Cut_InverseSimpsionMatrix", CreateVector = T, TransformAllData = FALSE, ValueCantBeSamlerThanZero = TRUE,
                    vectorName = "Vector_InverseSimpsion")
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_inverseSimpsionData", Exportname = "Cut_InverseSimpsionMatrix_Transform", CreateVector = T, TransformAllData = TRUE, ValueCantBeSamlerThanZero = FALSE,
                    vectorName = "Vector_InverseSimpsion_transform")

  #MDS

  data = CalculateBaseValues(data = data, intervallBy = 100,
                             NonNegative = FALSE, Importname1 = "nMDS", Importname2 = "Dim1", Exportname = "MDS")
  data = CalculateBaseValues(data = data, intervallBy = 100,
                             NonNegative = FALSE, Importname1 = "nMDS", Importname2 = "Cut_Dim1", Exportname = "Cut_MDS")

  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "MDS", Exportname = "MDSMatrix", CreateVector = F, TransformAllData = FALSE, ValueCantBeSamlerThanZero = FALSE)
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_MDS", Exportname = "Cut_MDSMatrix", CreateVector = T, TransformAllData = FALSE, ValueCantBeSamlerThanZero = FALSE,
                    vectorName = "Vector_MDS")
  data = AsyncTabel(data = data, intervallBy = 100,
                    Importname = "Cut_MDS", Exportname = "Cut_MDSMatrix_Transform", CreateVector = T, TransformAllData = TRUE, ValueCantBeSamlerThanZero = FALSE,
                    vectorName = "Vector_MDS_transform")

  #TOC

  data = TOC(data = data, intervallBy = 100, NonNegative = TRUE, Importname1 = "rawData", Importname2 = "TOC", Exportname = "TOC")

  data = TOCAsyncTabel(data = data, intervallBy = 100,
                    Importname = "TOC", Exportname = "TOCMatrix", CreateVector = T, TransformAllData = FALSE, ValueCantBeSamlerThanZero = FALSE,
                    vectorName = "Vector_TOC")
  data = TOCAsyncTabel(data = data, intervallBy = 100,
                       Importname = "TOC", Exportname = "TOCMatrix_Transform", CreateVector = T, TransformAllData = TRUE, ValueCantBeSamlerThanZero = FALSE,
                       vectorName = "Vector_TOC_transform")

  #TRACE

  data = TRACETabel(data = data, TraceName = "JJA")

  data = TRACETabel(data = data, TraceName = "DJF")

  data = TRACETabel(data = data, TraceName = "hydrological")

  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Done",sep="")


  ################################################################################
  ################################ Discrption ####################################
  ################################################################################

  #Diatom

  DiscriptionNames = ls(data$Description)

  DiatomDataNotFilterd = rep(FALSE, length(DiscriptionNames))

  for(d in 1:length(data$Description)){

    CurrentDiscriptionName = DiscriptionNames[d]

    CurrentDiscription = data$Diatom[[CurrentDiscriptionName]]$SRS_data

    if(!is.null(CurrentDiscription)){

      DiatomDataNotFilterd[d] = TRUE

    }
  }
  DiatomDataNotFilterd = DiscriptionNames[DiatomDataNotFilterd]

  #Carbon

  DiscriptionNames = ls(data$Description)

  CarbonDataNotFilterd = rep(FALSE, length(DiscriptionNames))

  for(d in 1:length(data$Description)){

    CurrentDiscriptionName = DiscriptionNames[d]

    CurrentDiscription = data$Carbon[[CurrentDiscriptionName]]$rawData$TOC

    if(sum(!is.na(CurrentDiscription))==0){

      CurrentDiscription = NULL

    }

    if(!is.null(CurrentDiscription)){

      CarbonDataNotFilterd[d] = TRUE

    }
  }
  CarbonDataNotFilterd = DiscriptionNames[CarbonDataNotFilterd]

  #Get SRS recounts

  SRSRecountsConditionData = NULL
  SRSRecountsCondition = TRUE
  SRSRecountsConditionCounter = 0

  while (SRSRecountsCondition) {

    SRSRecountsConditionCounter = SRSRecountsConditionCounter+1

    SRSRecountsConditionData = data$Diatom[[SRSRecountsConditionCounter]]$SRS_data$`Total numbers of counted diatom valves`[1]

    if(!is.null(SRSRecountsConditionData)){

      SRSRecountsCondition = FALSE

    }


    if(SRSRecountsConditionCounter == length(data$Description)){

      SRSRecountsCondition = FALSE

    }
  }
  extraDiscriptionAdder3 = matrix(c("Adjusted_Total_number_of_counted_diatom_valves",SRSRecountsConditionData),ncol = 2)

  #Add discription

  for (i in 1:length(data$Description)){

    #Diatom

    extraDiscriptionAdder1 = matrix(c("DiatomDataNotFilterd","FALSE"),ncol = 2)

    if(sum(DiatomDataNotFilterd == ls(data$Description[i]))==1){

      extraDiscriptionAdder1 = matrix(c("DiatomDataNotFilterd","TRUE"),ncol = 2)

    }

    #Carbon

    extraDiscriptionAdder2 = matrix(c("CarbonDataNotFilterd","FALSE"),ncol = 2)

    if(sum(CarbonDataNotFilterd == ls(data$Description[i]))==1){

      extraDiscriptionAdder2 = matrix(c("CarbonDataNotFilterd","TRUE"),ncol = 2)

    }

    data$Description[[ls(data$Description[i])]] = rbind(data$Description[[ls(data$Description[i])]],
                                                        extraDiscriptionAdder1,
                                                        extraDiscriptionAdder2,
                                                        extraDiscriptionAdder3)

  }

  #Get Number of not Filtered Data points
  #Diatom
  DiatomNotFilterdDataSheets = 0

  for (i in 1:length(data$Description)){

    if(data$Description[[i]][which(data$Description[[i]][,1] == "DiatomDataNotFilterd"),2]){

      DiatomNotFilterdDataSheets = DiatomNotFilterdDataSheets+1

    }
  }

  #Carbon
  CarbonNotFilterdDataSheets = 0

  for (i in 1:length(data$Description)){

    if(data$Description[[i]][which(data$Description[[i]][,1] == "CarbonDataNotFilterd"),2]){

      CarbonNotFilterdDataSheets = CarbonNotFilterdDataSheets+1

    }
  }

  cat("\n\n","Used Diatom Data (Not Filterd): ",DiatomNotFilterdDataSheets,
      "\n","Used Carbon Data (Not Filterd): ",CarbonNotFilterdDataSheets,
      "\n\n","Adjusted Total number of counted diatom valves: ",SRSRecountsConditionData,sep="")

  ##############################################################################
  #Return

  return(data)

}

















