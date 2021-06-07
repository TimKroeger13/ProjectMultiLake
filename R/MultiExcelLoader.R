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

  AddAgges = function(depth,age){

    if(!age$compositedepth[2]==0.25){

      setwd(orginalWorkingDirectoryPath)

      stop("Depth intervals at which ages are calculated is not 0.25!")

    }

    Ageresult=array(data = NA, dim = length(depth))

    for (i in 1:length(depth)){

      Ageresult[i]=as.numeric(age$modeloutput_mean[match(depth[i],age$compositedepth)])

    }

    return(Ageresult)
  }

  SingelExcelLoader = function(Folder,Excelname,age,FilenameKey){

    sheet=c("Diatom")

    Diatom=try(suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = sheet)),silent = TRUE)

    if(!class(Diatom)[1] == "try-error"){

      #Diatom Data
      TempColName=FixExcelRowNames(Diatom[5,])
      Diatom=Diatom[6:dim(Diatom)[1],]
      CoreName=toString(DisconnectNameAndDepth(Diatom[1,1])[1])
      Diatom[,1]=gsub(paste(CoreName," ",sep=""),"",as.matrix(Diatom[,1]))
      Diatom=suppressWarnings(as.data.frame(matrix(as.numeric(unlist(Diatom)),ncol = dim(Diatom)[2])))
      TempColName[1]="depth"
      colnames(Diatom)=enc2native(TempColName)
      Diatom=cbind(Diatom[,1:3],DeleteNaRows(Diatom[,4:dim(Diatom)[2]]))

      #Ages
      TempAge=as.data.frame(age)
      AgeFinder=read.table(textConnection(TempAge[,1]))[,1]==CoreName

      if(!sum(AgeFinder==TRUE)==0){

      TempAge=TempAge[which(AgeFinder),]
      TempAge=TempAge[order(TempAge$compositedepth,decreasing=F),]
      Diatom[,1]=AddAgges(Diatom[,1],TempAge)

      Folder[["Diatom"]][[FilenameKey]][["rawData"]]=Diatom
      Folder[["ages"]][[FilenameKey]]=TempAge

      }
    }

    return(Folder)

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

  if(!length(list.files()[grep("data.xlsx",list.files())])>0){

    setwd(orginalWorkingDirectoryPath)

    stop("Please drag and drop the files from the AWI database into the data folder in your working directory.")

  }

  FileNames=list.files()

  Folder=list()

  FileNamesXlsx=FileNames[grep("data.xlsx",FileNames)]
  FileNamesAge=FileNames[grep("ge.xlsx",FileNames)]

  if(identical(FileNamesAge, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No file named age.xlsx found")

  }

  Age=suppressMessages(read_excel(path = paste(getwd(),"/",FileNamesAge,sep=""),sheet = 1))

  if(dim(Age)[2]==1){

    setwd(orginalWorkingDirectoryPath)

    stop("Excel was read in incorrectly.\n
         Separate the columns and save the file again.")

  }

  for(i in 1:length(list.files()[grep("data.xlsx",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesXlsx)))[i,1]

    ChoosenFileNamesXlsx=FileNamesXlsx[i]

    Folder=SingelExcelLoader(Folder,ChoosenFileNamesXlsx,Age,FilenameKey)

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",length(list.files()[grep("data.xlsx",list.files())]),sep="")

  }

  setwd(orginalWorkingDirectoryPath)

  return(Folder)

}
