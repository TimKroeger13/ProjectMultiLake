#'RowFilter
#'@description
#'Filters row by Total numbers of counted diatom valves.
#'@param data List of data generates by the MultiExcelLoader function.
#'@param NoCDV Custom number for numbers of counted diatom valves where the data will be filtered by.
#'@param FilterByAlphaDiversity If TRUE the data will be filtered by Alpha diversity instead of NoCDV.
#'@import vegan
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

RowFilter = function(data,NoCDV=NULL,FilterByAlphaDiversity=F){

  for (i in 1:length(ls(data))){

    Tempdata=data[[i]][["rawData"]]
    TempdataForCalculation=data[[i]][["rawData"]][,4:dim(data[[i]][["rawData"]])[2]]

    if(FilterByAlphaDiversity){

      AlphaDeversity <- diversity(TempdataForCalculation,
                                  MARGIN = 2,
                                  index = "invsimpson")

      AlphaDeversity=ifelse(AlphaDeversity==Inf,0,AlphaDeversity)

      MinValue=max(AlphaDeversity)

      cat(paste("Alpha diversity for",names(data)[i],"was",round(MinValue,digits = 2)),"\n")

    }else if(!is.null(NoCDV)){

      MinValue=NoCDV

    }else{

      stop("Please enter the Total numbers of counted diatom valves (NoCDV) or set FilterByAlphaDiversity = TRUE")
    }

    CountedValves=NULL

    for (k in 1:dim(TempdataForCalculation)[1]){

      CountedValves[k]=sum(TempdataForCalculation[k,])

    }

    Tempdata=Tempdata[CountedValves>MinValue,]

    data[[names(data)[i]]][["FilterdData"]]=Tempdata

  }

  return(data)

}
