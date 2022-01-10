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

  DeleteNullRows = function(DataFrame){

    i=0

    while (i<dim(DataFrame)[1]) {

      i=i+1

      if(sum(DataFrame[i,3:dim(DataFrame)[2]])==0){

        DataFrame=DataFrame[-i,]

        i=i-1

      }
    }

    return(DataFrame)

  }

  DeleteUnclearDepth = function(DataFrame){

    i=0

    while (i<dim(DataFrame)[1]) {

      i=i+1

      if(is.na(DataFrame[i,1])){

        DataFrame=DataFrame[-i,]

        i=i-1

      }
    }

    return(DataFrame)

  }

  RoundDepth = function(roundData){

    for (i in 1:length(roundData)){

      if(roundData[i]!=0){

        if(roundData[i] %% 0.25 !=0){

          if(roundData[i] %% 0.25<0.125){

            roundData[i]=roundData[i] - roundData[i] %% 0.25

          }else{

            roundData[i]=roundData[i] + (0.25-roundData[i] %% 0.25)

          }
        }
      }
    }

    return(roundData)

  }


  AddAgges = function(depth,age){

    if(!(age$compositedepth[2]==0.25 || age$compositedepth[1]==0.25)){

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

      DatasetDescription = DatasetDescription[cbind(!DatasetDescription[,1]=="-------------------------------------------",!DatasetDescription[,1]=="-------------------------------------------")]

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

      DatasetDescriptionElemts = c("agemodel_available",
                                   "Diatoms_available",
                                   "Diatom_samples_=<11.7",
                                   "Diatom_samples_>11.7",
                                   "TC_available",
                                   "TOC_available",
                                   "LOI_availible",
                                   "Br_availiable")

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

    #Discription

    AgeFinder=read.table(textConnection(TempAge[,1]))[,1]==Folder[["Description"]][[FilenameKey]][1,2]

    if(!sum(AgeFinder==TRUE)==0){

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="agemodel_available"),2]=TRUE

    }else{

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="agemodel_available"),2]=FALSE

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

        Diatom=DeleteUnclearDepth(Diatom)

        Diatom[,1]=RoundDepth(Diatom[,1])

        Diatom[,1]=AddAgges(Diatom[,1],Folder[["ages"]][[Folder[["Description"]][[FilenameKey]][1,2]]])

        Diatom[is.na(Diatom)]=0 #<--------------------------------------------------------------------------------------------- Fix just for NA ins Diatom data

        Diatom=DeleteNullRows(Diatom)

        Folder[["Diatom"]][[FilenameKey]][["rawData"]]=Diatom

        #Discription

        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Diatom_samples_before_11700"),2]=
          as.character(sum(Diatom[,1]<holoceneBorder,na.rm = T))

        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Diatom_samples_after_11700"),2]=
          as.character(sum(Diatom[,1]>=holoceneBorder,na.rm = T))

      }
    }else{

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Diatoms_available"),2]=FALSE

    }

    #Carbon

    sheet=c("Organic")

    Carbon=try(suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = sheet)),silent = TRUE)

    if(!class(Carbon)[1] == "try-error"){

      AgeFinder=read.table(textConnection(TempAge[,1]))[,1]==Folder[["Description"]][[FilenameKey]][1,2]

      if(!sum(AgeFinder==TRUE)==0){

        TempColName=FixExcelRowNames(Carbon[6,])
        Carbon=Carbon[7:dim(Carbon)[1],]
        CoreName=toString(DisconnectNameAndDepth(Carbon[1,1])[1])
        Carbon[,1]=gsub(paste(CoreName," ",sep=""),"",as.matrix(Carbon[,1]))
        Carbon=suppressWarnings(as.data.frame(matrix(as.numeric(unlist(Carbon)),ncol = dim(Carbon)[2])))
        TempColName[1]="depth"
        colnames(Carbon)=enc2native(TempColName)
        Carbon[,1]=AddAgges(Carbon[,1],Folder[["ages"]][[Folder[["Description"]][[FilenameKey]][1,2]]])

        #AddBr

        Carbon=cbind(Carbon,matrix(NA,ncol = 1, nrow = dim(Carbon)[1]))
        colnames(Carbon)=c("depth","Nitrogen","TC","TOC","LOI","d13c","WaterContent","Br")

        Folder[["Carbon"]][[FilenameKey]][["rawData"]]=Carbon

      }

      if(sum(!is.na(Carbon[,which(colnames(Carbon)=="LOI")]))>0){
        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="LOI_availible"),2]=TRUE
      }else{ Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="LOI_availible"),2]=FALSE }


      if(sum(!is.na(Carbon[,which(colnames(Carbon)=="TC")]))>0){
        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TC_available"),2]=TRUE
        }else{ Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TC_available"),2]=FALSE }

      if(sum(!is.na(Carbon[,which(colnames(Carbon)=="TOC")]))>0){
        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TOC_available"),2]=TRUE
      }else{ Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TOC_available"),2]=FALSE }

    }else{

      #DataDiscription

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="LOI_availible"),2]=FALSE
      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TC_available"),2]=FALSE
      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="TOC_available"),2]=FALSE

    }

    #Element

    sheet=c("Element")

    Element=try(suppressMessages(read_excel(path = paste(getwd(),"/",Excelname,sep=""),sheet = sheet)),silent = TRUE)

    if(!class(Element)[1] == "try-error"){

      AgeFinder=read.table(textConnection(TempAge[,1]))[,1]==Folder[["Description"]][[FilenameKey]][1,2]

      if(!is.null(Folder[["Carbon"]][[FilenameKey]][["rawData"]])){

        if(!sum(AgeFinder==TRUE)==0){

          TempColName=FixExcelRowNames(Element[5,])
          Element=Element[6:dim(Element)[1],]
          CoreName=toString(DisconnectNameAndDepth(Element[1,1])[1])
          Element[,1]=gsub(paste(CoreName," ",sep=""),"",as.matrix(Element[,1]))
          Element=suppressWarnings(as.data.frame(matrix(as.numeric(unlist(Element)),ncol = dim(Element)[2])))
          TempColName[1]="depth"
          colnames(Element)=enc2native(TempColName)
          Element[,1]=AddAgges(Element[,1],Folder[["ages"]][[Folder[["Description"]][[FilenameKey]][1,2]]])

          #Add Brom

          Element$Br_Area

          Folder[["Element"]][[FilenameKey]][["rawData"]] = Element

        }
      }

      if(sum(!is.na(Element[,which(colnames(Element)=="Br_Area")]))>0){
        Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Br_availiable"),2]=TRUE
      }else{ Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Br_availiable"),2]=FALSE }

    }else{

      #DataDiscription

      Folder[["Description"]][[FilenameKey]][which(Folder[["Description"]][[FilenameKey]]=="Br_availiable"),2]=FALSE

    }

    return(Folder)

  }

  TraceExcelLoader = function(Folder,Excelname,age,FilenameKey,Outputname = "UnknownTrace"){

    TraceData = read.csv(file = paste(getwd(),"/",Excelname,sep=""))
    colnames(TraceData) = c("age","temp")
    rownames(TraceData) = TraceData[,1]

    Folder[["TRACE"]][[FilenameKey]][[paste(Outputname,"_raw",sep = "")]]=TraceData

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
            Please drag and drop the files from the AWI database into the data folder in your working directory or use the funktion ImportDatabase.")

  }

  if(!length(list.files()[grep("data.xlsx",list.files())])>0){

    setwd(orginalWorkingDirectoryPath)

    stop("Please drag and drop the files from the AWI database into the data folder in your working directory.")

  }

  FileNames=list.files()

  Folder=list()

  FileNamesXlsx=FileNames[grep("data.xlsx",FileNames)]
  FileNamesAge=FileNames[grep("Age_",FileNames)]
  FileNamesFixAge=FileNames[grep("FixedAge",FileNames)]
  FileNamesInsolation=FileNames[grep("insolation.txt",FileNames)]
  FileNamesClima=FileNames[grep("HAD3MCB-M",FileNames)]
  FileNamesMeta=FileNames[grep("LakeData_Alex",FileNames)]
  FileNamesJJA_temp=FileNames[grep("_JJA_temp_",FileNames)]
  FileNamesDJF_temp=FileNames[grep("_DJF_temp_",FileNames)]
  FileNameshydrological_year_temp=FileNames[grep("_hydrological_year_temp_",FileNames)]
  CorrectionPoints=FileNames[grep("CorrectionPoint",FileNames)]

  if(identical(FileNamesAge, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No file named age.csv found")

  }

  if(identical(FileNamesInsolation, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No file named insolation.txt found")

  }

  if(identical(FileNamesClima, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No Climate filenames found (Program was looking after '...HAD3MCB-M...)'")

  }

  if(identical(FileNamesMeta, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No file named LakeData_Alex.xlsx found")

  }

  if(identical(FileNamesJJA_temp, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No files named JJA_temp founded")

  }

  if(identical(FileNamesDJF_temp, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No files named DJF_temp founded")

  }

  if(identical(FileNameshydrological_year_temp, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No files named hydrological_year_temp founded")

  }

  if(identical(CorrectionPoints, character(0))){

    setwd(orginalWorkingDirectoryPath)

    stop("No files named CorrectionPoints.csv founded")

  }

  #Insolation

  InsolationCurve=read.table(FileNamesInsolation,sep=";",header = T)
  InsolationCurveValue=matrix(InsolationCurve[,7],ncol = 1)
  rownames(InsolationCurveValue)=InsolationCurve[,1]*1000
  Folder[["GlobalInsolation"]]=InsolationCurveValue

  #Age

  if(identical(FileNamesFixAge, character(0))){

    age=read.csv(file = paste(getwd(),"/",FileNamesAge,sep=""))

    if(dim(age)[2]==1){

      colnames=unlist(strsplit(colnames(age),","))
      AmountofCollums=as.numeric(lengths(regmatches(age[1,1],gregexpr(",", age[1,1]))))+1
      age=data.frame(matrix(unlist(strsplit(unlist(age), ',')), ncol = AmountofCollums, byrow = TRUE))
      colnames(age)=colnames

    }else{

      RawFixedAges = strsplit(age[,1]," ")
      FixedAges = NULL

      for (p in 1:dim(age)[1]){

        FixedAges = c(FixedAges,RawFixedAges[[p]][1])

        if(p %% 1000 == 0){

          cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
              p,"/",dim(age)[1]," Fixing age Names ",sep="")

        }
      }

      age[,1]=FixedAges

      write.table(age,file = paste(getwd(),"/","FixedAge",".csv",sep=""),
                  append = FALSE,na = "", sep = "; ", row.names = F, col.names = T,quote = F)

    }

  }else{ #FixAges

    age=read.csv(file = paste(getwd(),"/",FileNamesFixAge,sep=""),sep = ";")

  }

  #Climate

  for (i in 1:length(FileNamesClima)){

    ClimateName = FileNamesClima[i]

    Clima = read.table(ClimateName)

    Clima[Clima==-999] = NA

    #GetAges (Based on Heidruns Email)

    ClimaValueStart = 1950
    ClimaValueIntervalls = 200

    ClimaValueEnd = ClimaValueStart - (ClimaValueIntervalls/2)
    ClimaValueStart =  ClimaValueEnd - ClimaValueIntervalls*(dim(Clima)[2]-2)

    colnames(Clima) = c("names",(1950-seq(from = ClimaValueStart, to = ClimaValueEnd, by = ClimaValueIntervalls)))
    CilmaAges = (1950-seq(from = ClimaValueStart, to = ClimaValueEnd, by = ClimaValueIntervalls))

    #GetAnuel

    ClimaMonth = NA

#annual summer winter

    if(length(try(grep("annual",ClimateName)))>0){

      ClimaMonth="annual"

    }else if (length(try(grep("summer",ClimateName)))>0){

      ClimaMonth="summer"

    }else if (length(try(grep("winter",ClimateName)))>0){

      ClimaMonth="winter"

    }

    for (k in 1:dim(Clima)[1]){

      ClimaVektor = Clima[k,2:dim(Clima)[2]]

      if (!sum(!is.na(ClimaVektor))==0){

        Folder[["Clima"]][[ClimaMonth]][[Clima[k,1]]]=ClimaVektor

      }
    }
  }

  #Alex MetaData

  LakeData = try(suppressMessages(read_excel(path = paste(getwd(),"/",FileNamesMeta,sep=""),sheet = 1)),silent = TRUE)
  Headder = colnames(LakeData)
  LakeData=suppressWarnings(as.data.frame(matrix(as.character(unlist(LakeData)),ncol = dim(LakeData)[2])))
  colnames(LakeData) = Headder
  Folder[["LakeData"]]=LakeData

  #TraceData

  #FileNamesJJA_temp

  for(i in 1:length(list.files()[grep("_JJA_temp_",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesJJA_temp)))[i,1]

    Excelname=FileNamesJJA_temp[i]

    Folder=TraceExcelLoader(Folder,Excelname,age,FilenameKey,Outputname = "JJA")

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",length(list.files()[grep("data.xlsx",list.files())])," Loading CSV from ",getwd(),sep="")

  }

  # _DJF_temp_

  for(i in 1:length(list.files()[grep("_DJF_temp_",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesDJF_temp)))[i,1]

    Excelname=FileNamesDJF_temp[i]

    Folder=TraceExcelLoader(Folder,Excelname,age,FilenameKey,Outputname = "DJF")

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",length(list.files()[grep("data.xlsx",list.files())])," Loading CSV from ",getwd(),sep="")

  }

  # _hydrological_year_temp_

  for(i in 1:length(list.files()[grep("_hydrological_year_temp_",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNameshydrological_year_temp)))[i,1]

    Excelname=FileNameshydrological_year_temp[i]

    Folder=TraceExcelLoader(Folder,Excelname,age,FilenameKey,Outputname = "hydrological")

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",length(list.files()[grep("data.xlsx",list.files())])," Loading CSV from ",getwd(),sep="")

  }

  #CorrectionPoints

  CorrectionPoints=read.table(CorrectionPoints,sep=";",header = T)
  Folder[["CorrectionPoints"]]=CorrectionPoints

  #CoreList

  CoreList = matrix(NA,ncol = 2, nrow = length(list.files()[grep("data.xlsx",list.files())]))
  colnames(CoreList) = c("name","ID")

  for(i in 1:length(list.files()[grep("data.xlsx",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesXlsx)))[i,1]

    Excelname=FileNamesXlsx[i]

    Folder = SingelExcelLoader(Folder,Excelname,age,FilenameKey)

    CoreList[i,] = c(FilenameKey,i)

    #Printer
    cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
        i,"/",length(list.files()[grep("data.xlsx",list.files())])," Loading Excels from ",getwd(),sep="")

  }

  Folder[["CoreList"]] = CoreList

  #Add extra description

  extraDiscription1 = "RID"
  extraDiscription2 = colnames(LakeData)

  AlreadyUsedNames = Folder$Description[[1]][,1]

  extraDiscriptionList = NULL

  for(e in extraDiscription2){

    if(!(sum(e==AlreadyUsedNames)+sum(e=="CoreID")+sum(e=="TanDEM-X (Krieger et al. 2007)) tile number")+sum(e=="TanDEM-X (Krieger et al. 2007)) tile number")
         +sum(e=="BathymetricData")+sum(e=="ReferenceBathy"))>0){

      extraDiscriptionList=c(extraDiscriptionList,e)

    }
  }

  for(h in 1:length(list.files()[grep("data.xlsx",list.files())])){

    FilenameKey=read.table(textConnection(gsub("_", " ", FileNamesXlsx)))[h,1]

    extraDiscriptionMatrix = matrix(NA,ncol = 2, nrow = length(extraDiscriptionList)+1)
    extraDiscriptionMatrix[,1] = c("RID",extraDiscriptionList)
    extraDiscriptionMatrix[which(extraDiscriptionMatrix[,1] == "RID"),2] = which(CoreList[,1]==FilenameKey)

    for (b in extraDiscriptionList){

      if(sum(LakeData[,which(colnames(LakeData)=="CoreID")]==FilenameKey)>0){

      extraDiscriptionMatrix[which(extraDiscriptionMatrix[,1] == b),2] =
        LakeData[which(LakeData[,which(colnames(LakeData)=="CoreID")]==FilenameKey),which(colnames(LakeData)==b)]

      }
    }

    #MaxAge & Minage

    MaxMinAgeDiscriptionMatrix = matrix(NA,ncol = 2, nrow = 2)

    MaxMinAgeValue = Folder$Diatom[[FilenameKey]]$rawData$depth

    if(is.null(MaxMinAgeValue)){

      MaxMinAgeDiscriptionMatrix[,] = c("StartAge","EndAge",NA,NA)

    }else{

      MaxMinAgeDiscriptionMatrix[,] = c("StartAge","EndAge",max(MaxMinAgeValue),min(MaxMinAgeValue))

    }

    #Add Discription together
    Folder$Description[[FilenameKey]] = rbind(Folder$Description[[FilenameKey]],extraDiscriptionMatrix,MaxMinAgeDiscriptionMatrix)

  }

  Folder[["Filter"]] = "NotFiltert"
  Folder[["FilterList"]] = "NotFiltert"

  setwd(orginalWorkingDirectoryPath)

  cat("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n",
      "Done",sep="")

  return(Folder)

}
