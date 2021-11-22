#'MDSAsyncTabel
#'@description
#'Calculates the MDS as a Table given by the MDS function.
#'@param data List of data generates by the Multivar function.
#'@param intervallBy Intervalls by to interpolate to.
#'@importFrom stats t.test
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

MDSAsyncTabel = function(data, intervallBy = 100){

  #Printer
  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Calculating MDS Table ...",sep="")

  DiatomNames = ls(data$Diatom)

  MinEv = NULL
  MaxEv = NULL

  for (main in 1:length(DiatomNames)){

    MDSData = data$Diatom[[DiatomNames[main]]][["MDS"]]

    if(!is.null(MDSData)){

      MinEv = c(MinEv,min(MDSData[,1]))
      MaxEv = c(MaxEv,max(MDSData[,1]))

    }
  }

  #Create MDSAllInOneTabel

  MDSAllInOneTabel = matrix(NA, ncol = length(DiatomNames), nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(MDSAllInOneTabel) = DiatomNames
  rownames(MDSAllInOneTabel) = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (main in 1:length(DiatomNames)){

    MDSData = data$Diatom[[DiatomNames[main]]][["MDS"]]

    if(!is.null(MDSData)){

      for (k in 1:dim(MDSData)[1]){

        MDSAllInOneTabel[((MDSData[k,1] - min(MinEv)) / intervallBy)+1,main] = MDSData[k,2]

      }
    }
  }

  #Create RMDSAsyncTabel

  AsyncTabel = matrix(NA, ncol = 6, nrow = ((max(MaxEv)-min(MinEv))/intervallBy)+1)
  colnames(AsyncTabel) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  AsyncTabel[,1] = seq(from = min(MinEv), to = max(MaxEv), by = intervallBy)

  for (i in 1:(((max(MaxEv)-min(MinEv))/intervallBy)+1)){

    if(sum(!is.na(MDSAllInOneTabel[i,]))==1){

      AsyncTabel[i,2] = as.numeric(na.omit(MDSAllInOneTabel[i,]))
      AsyncTabel[i,3] = as.numeric(na.omit(MDSAllInOneTabel[i,]))
      AsyncTabel[i,4] = as.numeric(na.omit(MDSAllInOneTabel[i,]))
      AsyncTabel[i,5] = 999
      AsyncTabel[i,6] = 1

    }else{

      MDStest = t.test(MDSAllInOneTabel[i,])

      AsyncTabel[i,2] = MDStest$estimate
      AsyncTabel[i,3] = MDStest$conf.int[1]
      AsyncTabel[i,4] = MDStest$conf.int[2]
      AsyncTabel[i,5] = MDStest$p.value
      AsyncTabel[i,6] = sum(!is.na(MDSAllInOneTabel[i,]))

    }
  }

  AsyncTabel[which(AsyncTabel[,3]<0),3] = 0

  data[["MDSMatrix"]] = AsyncTabel

  return(data)

}
