





rateofChange = function(data){













}








'




#Beta diversity Turn Over

Tempdata=data[["Diatom"]][[i]][["SRS_data"]]
TempdataForCalculation=data[["Diatom"]][[i]][["SRS_data"]][,4:dim(data[["Diatom"]][[i]][["SRS_data"]])[2]]
TempdataForCalculation=round(TempdataForCalculation)

BetaDiversity=matrix(NA,ncol = 1,nrow = dim(TempdataForCalculation)[1]-1)
rownames(BetaDiversity)=Tempdata[1:(dim(Tempdata)[1]-1),1]
colnames(BetaDiversity)="Hierarchical beta diversity"

for(k in (dim(TempdataForCalculation)[1]):2){

  BetaDiversity[k-1]=vegdist(TempdataForCalculation[(k-1):k,],method = "bray")

}

data[["Diatom"]][[i]][["BetaDiversity"]][["OffsetBy_1"]]=BetaDiversity

for(x in 2:3){

  #Multipel Beta diversity Turn Over

  Tempdata=data[["Diatom"]][[i]][["SRS_data"]]
  TempdataForCalculation=data[["Diatom"]][[i]][["SRS_data"]][,4:dim(data[["Diatom"]][[i]][["SRS_data"]])[2]]
  TempdataForCalculation=round(TempdataForCalculation)

  SectionIntervall = x

  BetaDiversity=matrix(NA,ncol = 1,nrow = dim(TempdataForCalculation)[1]-SectionIntervall)
  rownames(BetaDiversity)=Tempdata[1:(dim(Tempdata)[1]-SectionIntervall),1]
  colnames(BetaDiversity)="Hierarchical beta diversity"

  comparisonIntervall=matrix(NA, nrow = 1, ncol = dim(TempdataForCalculation)[2])

  for(k in dim(TempdataForCalculation)[1]:SectionIntervall){

    for (z in 1:dim(TempdataForCalculation)[2]){

      comparisonIntervall[z]=mean(TempdataForCalculation[(k-SectionIntervall+1):k,z])

    }

    colnames(comparisonIntervall)=colnames(TempdataForCalculation)

    BetaDiversity[k-SectionIntervall]=vegdist(rbind(TempdataForCalculation[k-SectionIntervall,],comparisonIntervall),method = "bray")

  }

  data[["Diatom"]][[i]][["BetaDiversity"]][[paste("OffsetBy_",SectionIntervall,sep="")]]=BetaDiversity

}
'
