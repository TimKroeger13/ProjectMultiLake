#'ShowMetaDataFilters
#'@description
#'Shows all the Metadata values from witch data can be filtered.
#'@param data List of data generates by the Multivar function.
#'@import varhandle
#'@export
#'@return Returns a console output with all the Metadata values from witch data can be filtered.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

ShowMetaDataFilters = function(data){

  MetadataNames = data$Description[[1]][,1]

  MetaDataTable = matrix(NA,ncol = length(MetadataNames), nrow = length(ls(data$Description)))
  colnames(MetaDataTable) = MetadataNames

  for(i in 1:length(ls(data$Description))){

    MetaDataTable[i,] = data$Description[[i]][,2]

  }

  row.names(MetaDataTable) = MetaDataTable[,which(colnames(MetaDataTable)=="Unique CoreID")]

  MetaDataList = list()

  for(i in 1:dim(MetaDataTable)[2]){

    if(sum(!check.numeric(na.omit(MetaDataTable[,i])))==0){

      SingleMetavalue = as.numeric(MetaDataTable[,i])
      names(SingleMetavalue) = MetaDataTable[,which(colnames(MetaDataTable)=="Unique CoreID")]
      MetaDataList[[MetadataNames[i]]] = SingleMetavalue

    }else if(sum(is.na(as.logical(na.omit(MetaDataTable[,i]))))==0){

      SingleMetavalue = as.logical(MetaDataTable[,i])
      names(SingleMetavalue) = MetaDataTable[,which(colnames(MetaDataTable)=="Unique CoreID")]
      MetaDataList[[MetadataNames[i]]] = SingleMetavalue

    }
  }

  for(i in 1:length(MetaDataList)){

    MetadataNamesFilterd = ls(MetaDataList)

    cat(paste('"',MetadataNamesFilterd[i],'"',"\n",sep = ""))
    print(summary(MetaDataList[[MetadataNamesFilterd[i]]]))
    cat("\n\n")

  }
}
