#'SpeciesRichnessAsyncTabel
#'@description
#'Calculates the inverse species richness given by the SpeciesRichness function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

SpeciesRichnessAsyncTabel = function(data, intervallBy = 100){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating inverse Simpsion Table ...",sep="")

  DiatomNames = ls(data$Diatom)

  MinEv = NULL
  MaxEv = NULL

  for (main in 1:length(DiatomNames)){

    speciesRichnessData = data$Diatom[[DiatomNames[main]]][["speciesRichness"]]

    if(!is.null(speciesRichnessData)){

      MinEv = c(MinEv,min(speciesRichnessData[,1]))
      MaxEv = c(MaxEv,max(speciesRichnessData[,1]))

    }
  }

  #Create speciesRichnessAllInOneTabel

  speciesRichnessAllInOneTabel = matrix(NA, ncol = length(DiatomNames), nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(speciesRichnessAllInOneTabel) = DiatomNames
  rownames(speciesRichnessAllInOneTabel) = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (main in 1:length(DiatomNames)){

    speciesRichnessData = data$Diatom[[DiatomNames[main]]][["speciesRichness"]]

    if(!is.null(speciesRichnessData)){

      for (k in 1:dim(speciesRichnessData)[1]){

        speciesRichnessAllInOneTabel[((speciesRichnessData[k,1] - min(MinEv)) / intervallBy)+1,main] = speciesRichnessData[k,2]

      }
    }
  }

  #Create RspeciesRichnessAsyncTabel

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (i in 1:(((max(MaxEv)-min(MinEv))/intervallBy)+1)){

    if(sum(!is.na(speciesRichnessAllInOneTabel[i,]))==1){

      AsyncTabel[i,2] = as.numeric(na.omit(speciesRichnessAllInOneTabel[i,]))
      AsyncTabel[i,3] = as.numeric(na.omit(speciesRichnessAllInOneTabel[i,]))
      AsyncTabel[i,4] = as.numeric(na.omit(speciesRichnessAllInOneTabel[i,]))
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = 1

    }else{

      speciesRichnesstest = t.test(speciesRichnessAllInOneTabel[i,])

      AsyncTabel[i,2] = speciesRichnesstest$estimate
      AsyncTabel[i,3] = speciesRichnesstest$conf.int[1]
      AsyncTabel[i,4] = speciesRichnesstest$conf.int[2]
      AsyncTabel[i,5] = speciesRichnesstest$p.value
      AsyncTabel[i,6] = sum(!is.na(speciesRichnessAllInOneTabel[i,]))

    }
  }

  AsyncTabel[which(AsyncTabel[,3]<0),3] = 0

  data[["SpeciesRichnessMatrix"]] = AsyncTabel

  #VectorMatrix

  VectorTable = DataSignalAfterTable(DataAllInOneTabel = speciesRichnessAllInOneTabel,BasicAsyncTable = AsyncTabel,ValueCantBeSamlerThanZero = TRUE)

  data[["Vector_SpeciesRichnessMatrix"]] = VectorTable


  return(data)

}

