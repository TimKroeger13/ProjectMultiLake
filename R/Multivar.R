#'Multivar
#'@description
#'Calculates several parameters from deatom data.
#'@param data List of data generates by the MultiExcelLoader function.
#'@import vegan
#'@export
#'@return Returns the same List but with new added parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

Multivar = function(data){

  for (i in 1:length(ls(data))){

    Tempdata=data[[i]][["rawData"]]
    TempdataForCalculation=data[[i]][["rawData"]][,4:dim(data[[i]][["rawData"]])[2]]
    TempdataForCalculation=round(TempdataForCalculation)

    min=Inf

    for (k in 1:dim(TempdataForCalculation)[1]){

      if(sum(TempdataForCalculation[k,])<min){

        min=sum(TempdataForCalculation[k,])

      }
    }

    richness = rarefy(TempdataForCalculation, sample = min, se=T)
    shannon = diversity(TempdataForCalculation, index = "shannon")
    invsimpson = diversity(TempdataForCalculation, index = "invsimpson")

    richness=cbind(richness[1,],richness[2,])
    colnames(richness)=c("richness","error")
    row.names(richness)=Tempdata[,1]
    shannon=cbind(shannon)
    row.names(shannon)=Tempdata[,1]
    invsimpson=cbind(invsimpson)
    row.names(invsimpson)=Tempdata[,1]

    data[[names(data)[i]]][["richness"]]=richness
    data[[names(data)[i]]][["shannon"]]=shannon
    data[[names(data)[i]]][["invsimpson"]]=invsimpson

  }

  return(data)

}

