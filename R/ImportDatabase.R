#'ImportDatabase
#'@description
#'Import data form the AWI server.
#'@details
#'AWI VPN connection is needed.
#'@param directory specifyed directory.
#'@export
#'@return Downloaded data is stored in the data Folder.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

ImportDatabase = function(directory=NULL){

  FilesUpdated = 0

  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Please Wait for statup",sep="")

  orginalWorkingDirectoryPath=getwd()

  getwdTry <-  try(setwd(paste(getwd(),.Platform[2],"data",sep="")),silent = TRUE)

  if(class(getwdTry) == "try-error"){

    dir.create(paste(getwd(),.Platform[2],"data",sep=""))

    setwd(orginalWorkingDirectoryPath)

  }else{

    setwd(orginalWorkingDirectoryPath)

  }

  if(is.null(directory)){

    for (PossibleDirectory in letters){

      getwdTry <-  try(setwd(paste(toupper(PossibleDirectory),":/ArcLakesDB/xxxHIGHLIGHT-LAKES-DATASHEETCOPYxxx/standardized_datasheets/",sep="")),
                       silent = TRUE)

      if(!class(getwdTry) == "try-error"){

        directory = paste(toupper(PossibleDirectory),":",sep = "")

      }
    }
  }

  if(is.null(directory)){

    getwdTry <-  try(setwd(paste("/Volumes/projects/p_arclakes/ArcLakesDB/xxxHIGHLIGHT-LAKES-DATASHEETCOPYxxx",sep="")),
                     silent = TRUE)

    if(!class(getwdTry) == "try-error"){

      directory = "/Volumes/projects/p_arclakes"

    }
  }

  if(is.null(directory)){

    setwd(orginalWorkingDirectoryPath)

    stop("\n\nPlease connect to the AWI VPN and the AWI Server")

  }

  setwd(paste(directory,"/ArcLakesDB/xxxHIGHLIGHT-LAKES-DATASHEETCOPYxxx/standardized_datasheets/",sep=""))

  FilesToCopy = list.files()

  for (i in 1:length(FilesToCopy)){

    CurrentFilenName = FilesToCopy[i]

    #Error handling for false imported data

    if(substr(CurrentFilenName, 0, 2)=="~$"){

      CurrentFilenName = substr(CurrentFilenName, 3, nchar(CurrentFilenName))

    }

    #Check if Files exits and if they are up to date

    overwrite = TRUE

    if(file.exists(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))){

      NewDataInfo =  file.info(paste(directory,"/ArcLakesDB/xxxHIGHLIGHT-LAKES-DATASHEETCOPYxxx/standardized_datasheets/",CurrentFilenName,sep=""))$mtime
      OldDataInfo = file.info(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))$mtime

      if(NewDataInfo == OldDataInfo){

        overwrite = FALSE

      }
    }

    if(overwrite){

      FilesUpdated=FilesUpdated+1

      file.copy(from = paste(directory,"/ArcLakesDB/xxxHIGHLIGHT-LAKES-DATASHEETCOPYxxx/standardized_datasheets/",CurrentFilenName,sep=""),
                to = paste(orginalWorkingDirectoryPath,.Platform[2],"data",sep=""),
                recursive = F,
                overwrite = TRUE,
                copy.mode = TRUE,
                copy.date = TRUE)

    }

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,
        "/",
        length(FilesToCopy)," Downloading/Updating Files from ",getwd(),
        " to ",
        paste(orginalWorkingDirectoryPath,.Platform[2],"data",sep=""),sep="")

  }

  #Extra Data

  setwd(paste(directory,"/ArcLakesDB/WORKING_GROUP_TRANSFER/Tim Kroeger/MultivarData/",sep=""))

  FilesToCopy = list.files()

  for (i in 1:length(FilesToCopy)){

    CurrentFilenName = FilesToCopy[i]

    #Error handling for false imported data

    if(substr(CurrentFilenName, 0, 2)=="~$"){

      CurrentFilenName = substr(CurrentFilenName, 3, nchar(CurrentFilenName))

    }

    #Check if Files exits and if they are up to date

    overwrite = TRUE

    if(file.exists(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))){

      NewDataInfo =  file.info(paste(directory,"/ArcLakesDB/WORKING_GROUP_TRANSFER/Tim Kroeger/MultivarData/",CurrentFilenName,sep=""))$mtime
      OldDataInfo = file.info(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))$mtime

      if(NewDataInfo == OldDataInfo){

        overwrite = FALSE

      }
    }

    if(overwrite){

      FilesUpdated=FilesUpdated+1

      file.copy(from = paste(directory,"/ArcLakesDB/WORKING_GROUP_TRANSFER/Tim Kroeger/MultivarData/",CurrentFilenName,sep=""),
                to = paste(orginalWorkingDirectoryPath,.Platform[2],"data",sep=""),
                recursive = F,
                overwrite = TRUE,
                copy.mode = TRUE,
                copy.date = TRUE)

    }

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,
        "/",
        length(FilesToCopy)," Downloading/Updating Files from ",getwd(),
        " to ",
        paste(orginalWorkingDirectoryPath,.Platform[2],"data",sep=""),sep="")

  }

  setwd(orginalWorkingDirectoryPath)

  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Done\n\n",
      "Files Updated: ", FilesUpdated,sep="")

}
