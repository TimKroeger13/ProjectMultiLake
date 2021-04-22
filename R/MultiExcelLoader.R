#'MultiExcelLoader
#'@description
#'Loads multiple Excels from a sub folder named data in your working directory.
#'@details
#'When no Sub directory is give, this function will create the folder.
#'\cr When age data is given, all depth will be given a corresponding age.
#'@import compositions vegan readxl
#'@importFrom utils read.csv read.table
#'@export
#'@return Retruns a List of all excel sheets.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

MultiExcelLoader = function(){

  FixExcelRowNames = function(NameRow){

    return(gsub("\r","",gsub("\n","",as.character(NameRow),ignore.case = T),ignore.case = T))

  }

  DisconnectNameAndDepth = function(FirstEntry){

    return(read.table(textConnection(toString(FirstEntry))))

  }

  FixExcelRowNames = function(NameRow){

    return(gsub("\r","",gsub("\n","",as.character(NameRow),ignore.case = T),ignore.case = T))

  }

  DeleteNaRows = function(DataFrame){

    i=0

    while (i<dim(DataFrame)[2]) {

      i=i+1

      if(sum(is.na(DataFrame[,i]))==dim(DataFrame)[1]){

        DataFrame=DataFrame[,-i]

        i=i-1

      }
    }

    return(DataFrame)

  }

  AddAgges = function(depth,AgeTxtName){

    Ageresult=array(data = NA, dim = length(depth))
    Age=read.table(AgeTxtName)

    if(!Age[3,1]==0.25){

      setwd(orginalWorkingDirectoryPath)

      stop("Depth intervals at which ages are calculated is not 0.25!
      Please Change is in the Bacon Model!
      Bacon(core=...,thick=...,d.by = 0.25)")
    }

    for (i in 1:length(depth)){

      Ageresult[i]=as.numeric(Age[match(depth[i],Age[,1]),5])

    }

    return(Ageresult)
  }

  SingelExcelLoader = function(Excelname,AgeTxtName){

    Diatom=suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = "Diatom"))
    TempColName=FixExcelRowNames(Diatom[5,])
    Diatom=Diatom[6:dim(Diatom)[1],]
    CoreName=toString(DisconnectNameAndDepth(Diatom[1,1])[1])
    Diatom[,1]=gsub(paste(CoreName," ",sep=""),"",as.matrix(Diatom[,1]))
    Diatom=suppressWarnings(as.data.frame(matrix(as.numeric(unlist(Diatom)),ncol = dim(Diatom)[2])))
    TempColName[1]="depth"
    colnames(Diatom)=TempColName
    Diatom=cbind(Diatom[,1:3],DeleteNaRows(Diatom[,4:dim(Diatom)[2]]))

    if(!is.null(AgeTxtName)){

      Diatom[,1]=AddAgges(Diatom[,1],AgeTxtName)

    }

    return(Diatom)

  }

  ####

  orginalWorkingDirectoryPath=getwd()

  getwdTry <-  try(setwd(paste(getwd(),.Platform[2],"data",sep="")),silent = TRUE)

  if(class(getwdTry) == "try-error"){

    dir.create(paste(getwd(),.Platform[2],"data",sep=""))

    setwd(orginalWorkingDirectoryPath)

    stop("There is no directory with the name data.
            A new directory has been created in the workspace you selected./n
            Please drag and drop the files from the AWI database into the data folder in your working directory.")

  }

  if(!length(list.files(pattern = '\\.xlsx$'))>0){

    setwd(orginalWorkingDirectoryPath)

    stop("Please drag and drop the files from the AWI database into the data folder in your working directory.")

  }

  FileNames=list.files()

  Folder=list()

  FileNamesXlsx=FileNames[grep(".xlsx",FileNames)]
  FileNamesTxt=FileNames[grep(".txt",FileNames)]

  for(i in 1:(length(list.files())/2)){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesXlsx)))[i,1]

    ChoosenFileNamesXlsx=FileNamesXlsx[i]
    ChoosenFileNamesTxt=FileNamesTxt[grep(FilenameKey,FileNamesTxt)]

    Folder[[FilenameKey]][["rawData"]]=SingelExcelLoader(ChoosenFileNamesXlsx,ChoosenFileNamesTxt)

  }

  setwd(orginalWorkingDirectoryPath)

  return(Folder)

}
