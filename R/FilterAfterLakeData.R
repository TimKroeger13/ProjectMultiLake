#'FilterAfterLakeData
#'@description
#'Filters data based on the Lakedata Table
#'@param data List of data generates by the Multivar function.
#'@param distanceToCoast The minimum wanted distance to coast in km.
#'@export
#'@return Returns the same List but with less parameters.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

FilterAfterLakeData = function(data, distanceToCoast = 10){

  DiatomNames = ls(data$Diatom)

  counter = 0

  Tempdata=data[["Diatom"]]

  for (i in 1:length(DiatomNames)){
    #for (i in 1:8){

    counter = counter+1

    LakeData = data$LakeData

    rowIndicator = which(LakeData$CoreID == DiatomNames[counter]) #Row

    distanceToCoastData = as.numeric(LakeData$`DistanceToCoast(km)`[rowIndicator])

    if(!is.na(distanceToCoastData)){
      if(distanceToCoastData < distanceToCoast){

        Tempdata = Tempdata[-which(DiatomNames == DiatomNames[counter])]
        DiatomNames = DiatomNames[-which(DiatomNames == DiatomNames[counter])]

        counter = counter-1

      }
    }
  }

  data$Diatom = Tempdata

  return(data)

}
