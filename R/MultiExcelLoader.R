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

  holoceneBorder = 11700

  FixExcelRowNames = function(NameRow){

    return(gsub("\r","",gsub("\n","",as.character(NameRow),ignore.case = T),ignore.case = T))

  }

  DisconnectNameAndDepth = function(FirstEntry){

    return(read.table(textConnection(toString(FirstEntry))))

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

    #MetaData

    sheet=c("DatasetDescription")

    DatasetDescription=try(suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = sheet)),silent = TRUE)

    if(!class(DatasetDescription)[1] == "try-error"){

      DatasetDescription=DatasetDescription[4:28,2:3]

      DatasetDescription = DatasetDescription[cbind(!is.na(DatasetDescription[,2]),!is.na(DatasetDescription[,2]))]

      DatasetDescription_matrix=matrix(as.character(unlist(DatasetDescription)),ncol =2)

      DatasetDescription_matrix[,1]=c("CoreID",
                                      "Project_expedition",
                                      "sampling_year_CE",
                                      "PI_first_name",
                                      "PI_last_name",
                                      "Email",
                                      "ORCID",
                                      "Latitude_decdeg",
                                      "Longitude_decdeg",
                                      "Water_Depth_m",
                                      "Core_Length_m",
                                      "Drilling_Device",
                                      "LakeID",
                                      "Site_Name",
                                      "Country",
                                      "Catchment_Area_km2",
                                      "Lake_Extent_km2",
                                      "Lake_Depth_m" ,
                                      "Vegetation_Zone",
                                      "Climate_Zone",
                                      "Lake_Type")

      DatasetDescriptionElemts = c("dating_points",
                                   "agemodel_available",
                                   "Diatoms_available",
                                   "Diatom_samples_before_11700",
                                   "Diatom_samples_after_11700",
                                   "TC_available",
                                   "TOC_available")

      DatasetDescriptionFile = matrix(0, nrow = dim(DatasetDescription_matrix)[1]+length(DatasetDescriptionElemts),ncol =2)

      DatasetDescriptionFile[1:dim(DatasetDescription_matrix)[1],1:dim(DatasetDescription_matrix)[2]]=DatasetDescription

      DatasetDescriptionFile[(dim(DatasetDescription_matrix)[1]+1):(dim(DatasetDescriptionFile)[1]),1]=DatasetDescriptionElemts

      Folder[["Description"]][[FilenameKey]]=DatasetDescriptionFile

    }

    #Date

    sheet=c("Age")

    AgeDate=try(suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = sheet)),silent = TRUE)

    if(!class(AgeDate)[1] == "try-error"){

      AgeDate=AgeDate[7:dim(AgeDate)[1],]
      AgeDate=suppressWarnings(as.data.frame(matrix(as.numeric(unlist(AgeDate)),ncol = dim(AgeDate)[2])))

      #DataDiscription

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="14C_available"),2]=
        dim(AgeDate)[1]

    }

    #Ages

    TempAge=as.data.frame(age)
    AgeFinder=read.table(textConnection(TempAge[,1]))[,1]==Folder[["Description"]][[FilenameKey]][1,2]

    if(!sum(AgeFinder==TRUE)==0){

      TempAge=TempAge[which(AgeFinder),]
      TempAge=TempAge[order(TempAge$compositedepth,decreasing=F),]
      Folder[["ages"]][[FilenameKey]]=TempAge

    }

    #Diatom

    sheet=c("Diatom")

    Diatom=try(suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = sheet)),silent = TRUE)

    if(!class(Diatom)[1] == "try-error"){

      AgeFinder=read.table(textConnection(TempAge[,1]))[,1]==Folder[["Description"]][[FilenameKey]][1,2]

      #Discription

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Diatoms_available"),2]=TRUE

      if(!sum(AgeFinder==TRUE)==0){

        #Diatom Data
        TempColName=FixExcelRowNames(Diatom[5,])
        Diatom=Diatom[6:dim(Diatom)[1],]
        CoreName=toString(DisconnectNameAndDepth(Diatom[1,1])[1])
        Diatom[,1]=gsub(paste(CoreName," ",sep=""),"",as.matrix(Diatom[,1]))
        Diatom=suppressWarnings(as.data.frame(matrix(as.numeric(unlist(Diatom)),ncol = dim(Diatom)[2])))
        TempColName[1]="depth"
        colnames(Diatom)=enc2native(TempColName)
        Diatom=cbind(Diatom[,1:3],DeleteNaRows(Diatom[,4:dim(Diatom)[2]]))

        Diatom[,1]=AddAgges(Diatom[,1],Folder[["ages"]][[Folder[["Description"]][[FilenameKey]][1,2]]])

        Folder[["Diatom"]][[FilenameKey]][["rawData"]]=Diatom

        #Discription

        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Diatom_samples_before_11700"),2]=
          as.character(sum(Diatom[,1]<holoceneBorder,na.rm = T))

        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Diatom_samples_after_11700"),2]=
          as.character(sum(Diatom[,1]>=holoceneBorder,na.rm = T))

        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="agemodel_available"),2]=TRUE

      }else{

        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="agemodel_available"),2]=FALSE

      }
    }else{

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="agemodel_available"),2]=FALSE
      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Diatoms_available"),2]=FALSE

    }

    #Carbon

    sheet=c("Organic")

    Carbon=try(suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = sheet)),silent = TRUE)

    if(!class(Carbon)[1] == "try-error"){

      TempColName=FixExcelRowNames(Carbon[6,])
      Carbon=Carbon[7:dim(Carbon)[1],]
      CoreName=toString(DisconnectNameAndDepth(Carbon[1,1])[1])
      Carbon[,1]=gsub(paste(CoreName," ",sep=""),"",as.matrix(Carbon[,1]))
      Carbon=suppressWarnings(as.data.frame(matrix(as.numeric(unlist(Carbon)),ncol = dim(Carbon)[2])))
      TempColName[1]="depth"
      colnames(Carbon)=enc2native(TempColName)
      #Carbon[,1]=AddAgges(Carbon[,1],Folder[["ages"]][[Folder[["Description"]][[FilenameKey]][1,2]]]) <- Use later when ages can be used

      #DataDiscription

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TC_available"),2]


      if(sum(!is.na(Carbon[,which(colnames(Carbon)=="Total Carbon (TC, %)")]))>0){
        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TC_available"),2]=TRUE
        }else{ Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TC_available"),2]=FALSE }

      if(sum(!is.na(Carbon[,which(colnames(Carbon)=="Total Organic Carbon (TOC, %)")]))>0){
        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TOC_available"),2]=TRUE
      }else{ Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TOC_available"),2]=FALSE }

    }else{

      #DataDiscription

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TC_available"),2]=FALSE
      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TOC_available"),2]=FALSE

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
  FileNamesInsolation=FileNames[grep("insolation.txt",FileNames)]

  if(identical(FileNamesAge, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No file named age.xlsx found")

  }

  age=suppressMessages(read_excel(path = paste(getwd(),"/",FileNamesAge,sep=""),sheet = 1))

  if(dim(age)[2]==1){

    colnames=unlist(strsplit(colnames(age),","))
    AmountofCollums=as.numeric(lengths(regmatches(age[1,1],gregexpr(",", age[1,1]))))+1
    age=data.frame(matrix(unlist(strsplit(unlist(age), ',')), ncol = AmountofCollums, byrow = TRUE))
    colnames(age)=colnames

  }

  #Insolation

  if(identical(FileNamesInsolation, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No file named insolation.txt found")

  }

  InsolationCurve=read.table(FileNamesInsolation,sep=";",header = T)
  InsolationCurveValue=matrix(InsolationCurve[,7],ncol = 1)
  rownames(InsolationCurveValue)=InsolationCurve[,1]*1000
  Folder[["GlobalInsolation"]]=InsolationCurveValue

  for(i in 1:length(list.files()[grep("data.xlsx",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesXlsx)))[i,1]

    Excelname=FileNamesXlsx[i]

    Folder=SingelExcelLoader(Folder,Excelname,age,FilenameKey)

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",length(list.files()[grep("data.xlsx",list.files())])," Loading Excels from ",getwd(),sep="")

  }

  setwd(orginalWorkingDirectoryPath)

  return(Folder)

}
