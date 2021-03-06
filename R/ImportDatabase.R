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

  MacPath = "/Volumes/projects/p_arclakes"
  ArlakesDirect = "//smb.isipd.dmawi.de/projects/p_arclakes/"
  PathDatabase = "/ArcLakesDB/xxxHIGHLIGHT-LAKES-DATASHEETCOPYxxx/standardized_datasheets/"
  PathMetadata = "/ArcLakesDB/LAKEDATA/00-METADATA/"
  PathJJA_Mean_Temp = "/TraCE/Data-Output/JJA_Mean_Temp/"
  PathDJF_Mean_Temp = "/TraCE/Data-Output/DJF_Mean_Temp/"
  PathHydrological_Year_Temp = "/TraCE/Data-Output/Hydrological_Year_Temp/"

  FilesUpdated = -1

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

      getwdTry <-  try(setwd(paste(toupper(PossibleDirectory),":",PathDatabase,sep="")),
                       silent = TRUE)

      if(!class(getwdTry) == "try-error"){

        directory = paste(toupper(PossibleDirectory),":",sep = "")

      }
    }
  }

  if(is.null(directory)){

    getwdTry <-  try(setwd(MacPath),
                     silent = TRUE)

    if(!class(getwdTry) == "try-error"){

      directory = MacPath

    }
  }


  if(is.null(directory)){

    getwdTry <-  try(setwd(ArlakesDirect),
                     silent = TRUE)

    if(!class(getwdTry) == "try-error"){

      directory = ArlakesDirect

    }
  }

  if(is.null(directory)){

    setwd(orginalWorkingDirectoryPath)

    stop("\n\nPlease connect to the AWI VPN and the AWI Server")

  }

  #Standadized data sheets

  setwd(paste(directory,PathDatabase,sep=""))

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

      NewDataInfo =  file.info(paste(directory,PathDatabase,CurrentFilenName,sep=""))$mtime
      OldDataInfo = file.info(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))$mtime

      if(NewDataInfo == OldDataInfo){

        overwrite = FALSE

      }
    }

    if(overwrite){

      FilesUpdated=FilesUpdated+1

      file.copy(from = paste(directory,PathDatabase,CurrentFilenName,sep=""),
                to = paste(orginalWorkingDirectoryPath,.Platform[2],"data",sep=""),
                recursive = F,
                overwrite = TRUE,
                copy.mode = TRUE,
                copy.date = TRUE)

      if(length(grep("Age_",CurrentFilenName))>0){

        suppressWarnings(file.remove(paste(orginalWorkingDirectoryPath,.Platform[2],"data/","FixedAge.csv",sep="")))

      }
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

  setwd(paste(directory,PathMetadata,sep=""))

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

      NewDataInfo =  file.info(paste(directory,PathMetadata,CurrentFilenName,sep=""))$mtime
      OldDataInfo = file.info(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))$mtime

      if(NewDataInfo == OldDataInfo){

        overwrite = FALSE

      }
    }

    if(overwrite){

      FilesUpdated=FilesUpdated+1

      file.copy(from = paste(directory,PathMetadata,CurrentFilenName,sep=""),
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

  #PathJJA_Mean_Temp

  setwd(paste(directory,PathJJA_Mean_Temp,sep=""))

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

      NewDataInfo =  file.info(paste(directory,PathJJA_Mean_Temp,CurrentFilenName,sep=""))$mtime
      OldDataInfo = file.info(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))$mtime

      if(NewDataInfo == OldDataInfo){

        overwrite = FALSE

      }
    }

    if(overwrite){

      FilesUpdated=FilesUpdated+1

      file.copy(from = paste(directory,PathJJA_Mean_Temp,CurrentFilenName,sep=""),
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

  #PathDJF_Mean_Temp

  setwd(paste(directory,PathDJF_Mean_Temp,sep=""))

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

      NewDataInfo =  file.info(paste(directory,PathDJF_Mean_Temp,CurrentFilenName,sep=""))$mtime
      OldDataInfo = file.info(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))$mtime

      if(NewDataInfo == OldDataInfo){

        overwrite = FALSE

      }
    }

    if(overwrite){

      FilesUpdated=FilesUpdated+1

      file.copy(from = paste(directory,PathDJF_Mean_Temp,CurrentFilenName,sep=""),
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

  #PathHydrological_Year_Temp

  setwd(paste(directory,PathHydrological_Year_Temp,sep=""))

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

      NewDataInfo =  file.info(paste(directory,PathHydrological_Year_Temp,CurrentFilenName,sep=""))$mtime
      OldDataInfo = file.info(paste(orginalWorkingDirectoryPath,.Platform[2],"data/",CurrentFilenName,sep=""))$mtime

      if(NewDataInfo == OldDataInfo){

        overwrite = FALSE

      }
    }

    if(overwrite){

      FilesUpdated=FilesUpdated+1

      file.copy(from = paste(directory,PathHydrological_Year_Temp,CurrentFilenName,sep=""),
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
