#'DataSignalAfterTable
#'@description
#'Calculates the vector orientation of a data Table.
#'@param DataAllInOneTabel Tabel of all data where Rows are the Time and columns are the variabels.
#'@param BasicAsyncTable An AsyncTable based on the mean value.
#'@param ValueCantBeSamlerThanZero If true, values smaler than 0 become 0.
#'@export
#'@return nothing.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

DataSignalAfterTable = function(DataAllInOneTabel, BasicAsyncTable, ValueCantBeSamlerThanZero = FALSE){

  VectorAllInOneTabel = DataAllInOneTabel
  VectorAllInOneTabel[] = NA

  for (k in 1:dim(DataAllInOneTabel)[2]){

    for (i in 2:dim(DataAllInOneTabel)[1]){

      StartPoint = DataAllInOneTabel[(i-1),k]
      EndPoint = DataAllInOneTabel[i,k]

      if(!is.na(StartPoint) && !is.na(EndPoint)){

        VectorAllInOneTabel[i,k] = EndPoint-StartPoint

      }
    }
  }

  #Calculate DirectionVectorTable

  DirectionVectorTable = matrix(NA,ncol = 6, nrow = dim(DataAllInOneTabel)[1])
  colnames(DirectionVectorTable) = c("Depth","MeanValue","ConfDown","ConfUp","P_value","n")
  DirectionVectorTable[,1] = as.numeric(rownames(DataAllInOneTabel))


  for (i in 1:dim(DirectionVectorTable)[1]){

    if(i==1){

      DirectionVectorTable[i,2] = 0

    }else{

      if(sum(!is.na(VectorAllInOneTabel[i,]))==1){

        DirectionVectorTable[i,2] = as.numeric(na.omit(VectorAllInOneTabel[i,]))
        DirectionVectorTable[i,3] = as.numeric(na.omit(VectorAllInOneTabel[i,]))
        DirectionVectorTable[i,4] = as.numeric(na.omit(VectorAllInOneTabel[i,]))
        DirectionVectorTable[i,5] = 999
        DirectionVectorTable[i,6] = 1

      }else{

        TableTtest = t.test(na.omit(VectorAllInOneTabel[i,]))

        DirectionVectorTable[i,2] = TableTtest$estimate + DirectionVectorTable[(i-1),2]
        DirectionVectorTable[i,3] = TableTtest$conf.int[1] + DirectionVectorTable[(i-1),2]
        DirectionVectorTable[i,4] = TableTtest$conf.int[2] + DirectionVectorTable[(i-1),2]
        DirectionVectorTable[i,5] = TableTtest$p.value
        DirectionVectorTable[i,6] = sum(!is.na(VectorAllInOneTabel[i,]))

      }

    }
  }

  # Least Square adjustment

  AdjustedDirectionVectorTable = DirectionVectorTable
  AdjustedDirectionVectorTable[] = NA

  Resudials = NULL

  for (i in 1:dim(AdjustedDirectionVectorTable)[1]){

    Resudials = c(Resudials,as.numeric(na.omit(DataAllInOneTabel[i,])) - DirectionVectorTable[i,2])

  }

  ResudialsDifference = mean(Resudials)

  AdjustedDirectionVectorTable[,1] = DirectionVectorTable[,1]
  AdjustedDirectionVectorTable[,2] = DirectionVectorTable[,2] + ResudialsDifference
  AdjustedDirectionVectorTable[,3] = DirectionVectorTable[,3] + ResudialsDifference
  AdjustedDirectionVectorTable[,4] = DirectionVectorTable[,4] + ResudialsDifference
  AdjustedDirectionVectorTable[,5] = DirectionVectorTable[,5]
  AdjustedDirectionVectorTable[,6] = DirectionVectorTable[,6]

  # Adjust conf intervals and P_values and N numbers

  for (i in 2:dim(AdjustedDirectionVectorTable)[1]){

    AdjustedDirectionVectorTable[i,3] = AdjustedDirectionVectorTable[i,2] + (BasicAsyncTable[i,3] - BasicAsyncTable[i,2])
    AdjustedDirectionVectorTable[i,4] = AdjustedDirectionVectorTable[i,2] + (BasicAsyncTable[i,4] - BasicAsyncTable[i,2])
    AdjustedDirectionVectorTable[i,5] = BasicAsyncTable[i,5]
    AdjustedDirectionVectorTable[i,6] = BasicAsyncTable[i,6]

  }

  #ValueCantBeSamlerThanZero

  if(ValueCantBeSamlerThanZero){

    Cutvalue = 20000

    FittingAges = AdjustedDirectionVectorTable[,1]<=Cutvalue
    NullAdder = min(AdjustedDirectionVectorTable[FittingAges,2])

    if(NullAdder <0){

      AdjustedDirectionVectorTable[,2] = AdjustedDirectionVectorTable[,2] + abs(NullAdder)
      AdjustedDirectionVectorTable[,3] = AdjustedDirectionVectorTable[,3] + abs(NullAdder)
      AdjustedDirectionVectorTable[,4] = AdjustedDirectionVectorTable[,4] + abs(NullAdder)

    }

    AdjustedDirectionVectorTable[AdjustedDirectionVectorTable[,3]<0,3] = 0
    AdjustedDirectionVectorTable[AdjustedDirectionVectorTable[,4]<0,4] = 0

  }

  return(AdjustedDirectionVectorTable)

}
