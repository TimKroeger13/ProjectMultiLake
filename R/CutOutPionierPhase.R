#'CutOutPionierPhase
#'@description
#'Cut out the Pionier Phase of the Lake data.
#'\cr The new stating value is the first Value that passes the T.test that for two samples with a P_value <=0.05.
#'@param data List of data generates by the Multivar function.
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

CutOutPionierPhase = function(data){

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

  AllDiatomsNames = ls(data$Diatom)
  CorrectionPoints = data$CorrectionPoints
  CorrectionPoints[is.na(CorrectionPoints[,2]),2] = 0

  for (z in 1:length(AllDiatomsNames)){

    PionierFreeSRS = data$Diatom[[AllDiatomsNames[z]]]$SRS_data
    rateOfChangeData = data$Diatom[[AllDiatomsNames[z]]]$SRS_data

    if(!is.null(rateOfChangeData)){

      DiatomNames = AllDiatomsNames[z]

      PinoeerPoint = CorrectionPoints[which(CorrectionPoints[,1] == DiatomNames),2]

      if(PinoeerPoint>0){

        PionierFreeSRS = PionierFreeSRS[-(dim(rateOfChangeData)[1]-PinoeerPoint+1):-(dim(rateOfChangeData)[1]),]

      }

      data$Diatom[[AllDiatomsNames[z]]][["Cut_SRS_data"]] = PionierFreeSRS

    }

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        z,"/",length(ls(data[["Diatom"]]))," Cut SRS Data",sep="")

  }

  return(data)

}

