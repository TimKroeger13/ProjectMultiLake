#'GetDiscriptionCSV
#'@description
#'Get a short Discription from all Cores and save them as csv.
#'@import compositions readxl
#'@importFrom utils read.csv read.table write.table
#'@export
#'@return Save CSV in workspace.
#'@author Tim Kroeger
#'@note This function has only been developed for the Alfred Wegener Institute Helmholtz Centre for Polar and Marine Research and should therefore only be used in combination with their database.

GetDiscriptionCSV = function(){

  orginalWorkingDirectoryPath=getwd()

  getwdTry <-  try(setwd(paste(getwd(),.Platform[2],"data",sep="")),silent = TRUE)

  if(class(getwdTry) == "try-error"){

    dir.create(paste(getwd(),.Platform[2],"data",sep=""))

    setwd(orginalWorkingDirectoryPath)

    stop("There is no directory with the name data.
            A new directory has been created in the workspace you selected./n
            Please drag and drop the files from the AWI database into the data folder in your working directory.")

  }

  if(!length(list.files()[grep("data.xlsx",list.files())])>0){

    setwd(orginalWorkingDirectoryPath)

    stop("Please drag and drop the files from the AWI database into the data folder in your working directory.")

  }

  FileNames=list.files()

  FileNamesXlsx=FileNames[grep("data.xlsx",FileNames)]

  DescriptionOutput=matrix(NA,ncol = length(FileNamesXlsx),nrow = 2)
  rownames(DescriptionOutput)=c("Latitude","Longitude")
  Descriptionname = c(1:length(FileNamesXlsx))

  for(i in 1:length(list.files()[grep("data.xlsx",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesXlsx)))[i,1]

    ChoosenFileNamesXlsx=FileNamesXlsx[i]

    sheet=c("DatasetDescription")

    DatasetDescription=try(suppressMessages(read_excel(path = paste(getwd(),"/",ChoosenFileNamesXlsx,sep=""),sheet = sheet)),silent = TRUE)

    DatasetDescription=suppressWarnings(as.data.frame(matrix(as.character(unlist(DatasetDescription)),ncol = dim(DatasetDescription)[2])))

    DescriptionOutput[,i] = as.numeric(DatasetDescription[14:15,3])
    Descriptionname[i] = DatasetDescription[4,3]

  }

  setwd(orginalWorkingDirectoryPath)

  colnames(DescriptionOutput)=Descriptionname

  write.table(DescriptionOutput,file = paste(getwd(),"/","Description",".csv",sep=""),
              append = FALSE,na = "", sep = ", ", row.names = T, col.names = T)

  invisible(DescriptionOutput)

}
