#'MetadataFilter
#'@description
#'Filters data based on the Lakedata Table and other Metadata. To see with data fits please use the ShowMetaDataFilters function.
#'@param data List of data generates by the Multivar function.
#'@param NAignore If TRUE, value that are NA don't get filtered.
#'@param Filter1
#'\cr\cr name = name of the metadata. Get this from the function ShowMetaDataFilters.
#'\cr\cr criteria = Can be "greater","smaller","TRUE" or"FALSE".
#'\cr\cr value = The filter value when greater or smaller.
#'@param Filter2 Values for Filter 2.
#'@param Filter3 Values for Filter 3.
#'@param Filter4 Values for Filter 4.
#'@param Filter5 Values for Filter 5.
#'@export
#'@return .Returns the same List but with filtert parameters
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.
#'\cr Comma numbers are rounded up.

MetadataFilter = function(data,
                          NAignore = TRUE,
                          Filter1=c(name = "metdataName", criteria = c("greater","smaller","TRUE","FALSE"), value = NULL),
                          Filter2=c(name = "metdataName", criteria = c("greater","smaller","TRUE","FALSE"), value = NULL),
                          Filter3=c(name = "metdataName", criteria = c("greater","smaller","TRUE","FALSE"), value = NULL),
                          Filter4=c(name = "metdataName", criteria = c("greater","smaller","TRUE","FALSE"), value = NULL),
                          Filter5=c(name = "metdataName", criteria = c("greater","smaller","TRUE","FALSE"), value = NULL)){

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

  ##############################################################################
  ##############################################################################
  ##############################################################################
  #Filter1
  ##############################################################################
  ##############################################################################
  ##############################################################################

  TempFilter = Filter1

  if(!TempFilter["name"]=="metdataName"){

    if(!is.na(TempFilter["criteria"])){

      if(TempFilter["criteria"]=="greater" || TempFilter["criteria"]=="smaller" || TempFilter["criteria"]=="TRUE" || TempFilter["criteria"]=="FALSE"){

        if(sum(ls(MetaDataList)==TempFilter["name"])>0){

          #Diatom

          TempdataDiatom=data[["Diatom"]]

          AllDiatomsNames = ls(TempdataDiatom)

          counter=0

          for (d in 1:length(TempdataDiatom)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataDiatom = TempdataDiatom[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Diatom"]] = TempdataDiatom

          #Carbon

          TempdataCarbon=data[["Carbon"]]

          AllDiatomsNames = ls(TempdataCarbon)

          counter=0

          for (d in 1:length(TempdataCarbon)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataCarbon = TempdataCarbon[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Carbon"]] = TempdataCarbon

        }else{

          stop(paste("No value with the name '",TempFilter["name"],"' found",sep = ""))
        }
      }else{

        stop(paste("No criteria with the name 'greater' or 'smaller' found",sep = ""))
      }
    }else{

      stop(paste("Please enter a criteria",sep = ""))
    }
  }

  ##############################################################################
  ##############################################################################
  ##############################################################################
  #Filter2
  ##############################################################################
  ##############################################################################
  ##############################################################################

  TempFilter = Filter2

  if(!TempFilter["name"]=="metdataName"){

    if(!is.na(TempFilter["criteria"])){

      if(TempFilter["criteria"]=="greater" || TempFilter["criteria"]=="smaller" || TempFilter["criteria"]=="TRUE" || TempFilter["criteria"]=="FALSE"){

        if(sum(ls(MetaDataList)==TempFilter["name"])>0){

          #Diatom

          TempdataDiatom=data[["Diatom"]]

          AllDiatomsNames = ls(TempdataDiatom)

          counter=0

          for (d in 1:length(TempdataDiatom)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataDiatom = TempdataDiatom[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Diatom"]] = TempdataDiatom

          #Carbon

          TempdataCarbon=data[["Carbon"]]

          AllDiatomsNames = ls(TempdataCarbon)

          counter=0

          for (d in 1:length(TempdataCarbon)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataCarbon = TempdataCarbon[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Carbon"]] = TempdataCarbon

        }else{

          stop(paste("No value with the name '",TempFilter["name"],"' found",sep = ""))
        }
      }else{

        stop(paste("No criteria with the name 'greater' or 'smaller' found",sep = ""))
      }
    }else{

      stop(paste("Please enter a criteria",sep = ""))
    }
  }

  ##############################################################################
  ##############################################################################
  ##############################################################################
  #Filter3
  ##############################################################################
  ##############################################################################
  ##############################################################################

  TempFilter = Filter3

  if(!TempFilter["name"]=="metdataName"){

    if(!is.na(TempFilter["criteria"])){

      if(TempFilter["criteria"]=="greater" || TempFilter["criteria"]=="smaller" || TempFilter["criteria"]=="TRUE" || TempFilter["criteria"]=="FALSE"){

        if(sum(ls(MetaDataList)==TempFilter["name"])>0){

          #Diatom

          TempdataDiatom=data[["Diatom"]]

          AllDiatomsNames = ls(TempdataDiatom)

          counter=0

          for (d in 1:length(TempdataDiatom)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataDiatom = TempdataDiatom[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Diatom"]] = TempdataDiatom

          #Carbon

          TempdataCarbon=data[["Carbon"]]

          AllDiatomsNames = ls(TempdataCarbon)

          counter=0

          for (d in 1:length(TempdataCarbon)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataCarbon = TempdataCarbon[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Carbon"]] = TempdataCarbon

        }else{

          stop(paste("No value with the name '",TempFilter["name"],"' found",sep = ""))
        }
      }else{

        stop(paste("No criteria with the name 'greater' or 'smaller' found",sep = ""))
      }
    }else{

      stop(paste("Please enter a criteria",sep = ""))
    }
  }

  ##############################################################################
  ##############################################################################
  ##############################################################################
  #Filter4
  ##############################################################################
  ##############################################################################
  ##############################################################################

  TempFilter = Filter4

  if(!TempFilter["name"]=="metdataName"){

    if(!is.na(TempFilter["criteria"])){

      if(TempFilter["criteria"]=="greater" || TempFilter["criteria"]=="smaller" || TempFilter["criteria"]=="TRUE" || TempFilter["criteria"]=="FALSE"){

        if(sum(ls(MetaDataList)==TempFilter["name"])>0){

          #Diatom

          TempdataDiatom=data[["Diatom"]]

          AllDiatomsNames = ls(TempdataDiatom)

          counter=0

          for (d in 1:length(TempdataDiatom)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataDiatom = TempdataDiatom[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Diatom"]] = TempdataDiatom

          #Carbon

          TempdataCarbon=data[["Carbon"]]

          AllDiatomsNames = ls(TempdataCarbon)

          counter=0

          for (d in 1:length(TempdataCarbon)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataCarbon = TempdataCarbon[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Carbon"]] = TempdataCarbon

        }else{

          stop(paste("No value with the name '",TempFilter["name"],"' found",sep = ""))
        }
      }else{

        stop(paste("No criteria with the name 'greater' or 'smaller' found",sep = ""))
      }
    }else{

      stop(paste("Please enter a criteria",sep = ""))
    }
  }

  ##############################################################################
  ##############################################################################
  ##############################################################################
  #Filter5
  ##############################################################################
  ##############################################################################
  ##############################################################################

  TempFilter = Filter5

  if(!TempFilter["name"]=="metdataName"){

    if(!is.na(TempFilter["criteria"])){

      if(TempFilter["criteria"]=="greater" || TempFilter["criteria"]=="smaller" || TempFilter["criteria"]=="TRUE" || TempFilter["criteria"]=="FALSE"){

        if(sum(ls(MetaDataList)==TempFilter["name"])>0){

          #Diatom

          TempdataDiatom=data[["Diatom"]]

          AllDiatomsNames = ls(TempdataDiatom)

          counter=0

          for (d in 1:length(TempdataDiatom)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataDiatom = TempdataDiatom[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataDiatom = TempdataDiatom[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Diatom"]] = TempdataDiatom

          #Carbon

          TempdataCarbon=data[["Carbon"]]

          AllDiatomsNames = ls(TempdataCarbon)

          counter=0

          for (d in 1:length(TempdataCarbon)){

            counter=counter+1

            DiatomsNames = AllDiatomsNames[d]

            KeyValue = MetaDataList[[TempFilter[["name"]]]][which(names(MetaDataList[[TempFilter[["name"]]]]) == DiatomsNames)]

            if(is.na(KeyValue)){

              if(NAignore==TRUE){

              }else{

                TempdataCarbon = TempdataCarbon[-counter]
                counter = counter-1

              }

            }else{

              if(TempFilter["criteria"]=="greater"){

                if(KeyValue<=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="smaller"){


                if(KeyValue>=as.numeric(TempFilter["value"])){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="TRUE"){

                if(KeyValue==FALSE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }

              if(TempFilter["criteria"]=="FALSE"){

                if(KeyValue==TRUE){

                  TempdataCarbon = TempdataCarbon[-counter]
                  counter = counter-1

                }
              }
            }
          }

          data[["Carbon"]] = TempdataCarbon

        }else{

          stop(paste("No value with the name '",TempFilter["name"],"' found",sep = ""))
        }
      }else{

        stop(paste("No criteria with the name 'greater' or 'smaller' found",sep = ""))
      }
    }else{

      stop(paste("Please enter a criteria",sep = ""))
    }
  }

  ##############################################################################
  #Filter Code name and exact name
  ##############################################################################

  DataWasFilterd = FALSE

  #Create filterDiscriptionDocument

  filterDiscriptionDocument = matrix(NA, nrow = 1, ncol = 3);
  filterDiscriptionDocument[,1:3] = c("name","criteria","value")

  if(!Filter1["name"]=="metdataName"){

    DataWasFilterd = TRUE

    if(length(Filter1)==2){Filter1 = c(Filter1,"NA")}

    filterDiscriptionDocument = rbind(filterDiscriptionDocument,Filter1)

  }
  if(!Filter2["name"]=="metdataName"){

    DataWasFilterd = TRUE

    if(length(Filter2)==2){Filter2 = c(Filter2,"NA")}

    filterDiscriptionDocument = rbind(filterDiscriptionDocument,Filter2)

  }
  if(!Filter3["name"]=="metdataName"){

    DataWasFilterd = TRUE

    if(length(Filter3)==2){Filter3 = c(Filter3,"NA")}

    filterDiscriptionDocument = rbind(filterDiscriptionDocument,Filter3)

  }
  if(!Filter4["name"]=="metdataName"){

    DataWasFilterd = TRUE

    if(length(Filter4)==2){Filter4 = c(Filter4,"NA")}

    filterDiscriptionDocument = rbind(filterDiscriptionDocument,Filter4)

  }
  if(!Filter5["name"]=="metdataName"){

    DataWasFilterd = TRUE

    if(length(Filter5)==2){Filter5 = c(Filter5,"NA")}

    filterDiscriptionDocument = rbind(filterDiscriptionDocument,Filter5)

  }

  if(DataWasFilterd){

    #GetFilterKey

    FilterKey = ""

    for (i in 2:dim(filterDiscriptionDocument)[1]){

      for (k in 1:dim(filterDiscriptionDocument)[2]){

        FilterKey = c(FilterKey,substr(filterDiscriptionDocument[i,k], 0, 2))

      }
    }

    #Delete the first underscor
    FilterKey = FilterKey[-1]

    FilterKey = paste(FilterKey, collapse = "_")

    data[["Filter"]] = FilterKey
    data[["FilterList"]] = filterDiscriptionDocument

  }


  return(data)

}
