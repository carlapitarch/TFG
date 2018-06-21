#### PAQUETS
library(ggplot2)
library(reshape2)

library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

library(dplyr)

library(survival)
library(cmprsk)
library(riskRegression)
library(pec)
library(crrstep)

#### LECTURA DADES

dades<-read.csv('dades.csv')

#### BOXPLOTS MARCADORS

id<-seq(1,582,by=1)
marcadors<-data.frame(id=id,dades[,18:32])
marcadors.m<-melt(marcadors,id="id")
marcadors.m$id<-as.factor(marcadors.m$id)

# DESCRIPTIVA INICIAL

windows()
ggplot(marcadors.m,aes(x=variable,y=value)) + geom_boxplot(notch=FALSE, outlier.shape='*',
                                                           outlier.size=4, outlier.colour="black",outlier.alpha=0.5,
                                                           fill="red", alpha=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(axis.text=element_text(size=12))+ 
  theme(panel.background = element_rect(fill = 'gray93', colour = 'gray93')) +
  scale_y_continuous(breaks=seq(0, 14000, 1000)) +
  xlab ("") + ylab ("")


marcadors2<-dades[,c(2,18:32)]
marcadors.m2<-melt(marcadors2,id="group")

windows()
ggplot(marcadors.m2,aes(x=variable,y=value,fill=group)) + geom_boxplot(notch=FALSE, outlier.shape='*',
                                                           outlier.size=4, outlier.colour="black",outlier.alpha=0.5,
                                                           fill="red", alpha=0.1)+
  facet_grid(.~group)+
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(axis.text=element_text(size=12))+ 
  theme(panel.background = element_rect(fill = 'gray93', colour = 'gray93')) +
  scale_y_continuous(breaks=seq(0, 14000, 1000)) +
  xlab ("") + ylab ("")


#### TRANSFORMACIÓ DELS MARCADORS

bbddmarkers<-cbind(dades[,2],log(dades[,21:34]),dades[,35])
colnames(bbddmarkers)<-c("group","logTNFa","logIL_6","logIL_8","logMCP_1","logIP_10","logMIP_1b",
                         "logG_CSF","logGM_CSF","logIL_10","logIL_1ra","logIFNg","logEotaxin","logIL_17a",
                         "logIL_7","HNA2")

#### BOXPLOTS MARCADORS

id<-seq(1,582,by=1)
marcadors<-data.frame(id=id,bbddmarkers[,2:16])
marcadors.m<-melt(marcadors,id="id")
marcadors.m$id<-as.factor(marcadors.m$id)
# DESCRIPTIVA INICIAL


windows()
ggplot(marcadors.m,aes(x=variable,y=value)) + geom_boxplot(notch=FALSE, outlier.shape='*',
                                                           outlier.size=4, outlier.colour="black",outlier.alpha=0.5,
                                                           fill="red", alpha=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(axis.text=element_text(size=12))+ 
  theme(panel.background = element_rect(fill = 'gray93', colour = 'gray93')) +
  xlab ("") + ylab ("")


marcadors2<-bbddmarkers
marcadors.m2<-melt(marcadors2,id="group")

windows()
ggplot(marcadors.m2,aes(x=variable,y=value,fill=group)) + geom_boxplot(notch=FALSE, outlier.shape='*',
                                                                       outlier.size=4, outlier.colour="black",outlier.alpha=0.5,
                                                                       fill="red", alpha=0.1)+
  facet_grid(.~group)+
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(axis.text=element_text(size=12))+ 
  theme(panel.background = element_rect(fill = 'gray93', colour = 'gray93')) +
  xlab ("") + ylab ("")


#### BOXPLOTS MARCADORS

id<-seq(1,582,by=1)
marcadors<-data.frame(id=id,scale(bbddmarkers[,2:16]))
marcadors.m<-melt(marcadors,id="id")
marcadors.m$id<-as.factor(marcadors.m$id)
# DESCRIPTIVA INICIAL


windows()
ggplot(marcadors.m,aes(x=variable,y=value)) + geom_boxplot(notch=FALSE, outlier.shape='*',
                                                           outlier.size=4, outlier.colour="black",outlier.alpha=0.5,
                                                           fill="red", alpha=0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(axis.text=element_text(size=12))+ 
  theme(panel.background = element_rect(fill = 'gray93', colour = 'gray93')) +
  xlab ("") + ylab ("")


marcadors2<-data.frame(group,marcadors[,2:16])
marcadors.m2<-melt(marcadors2,id="group")

windows()
ggplot(marcadors.m2,aes(x=variable,y=value,fill=group)) + geom_boxplot(notch=FALSE, outlier.shape='*',
                                                                       outlier.size=4, outlier.colour="black",outlier.alpha=0.5,
                                                                       fill="red", alpha=0.1)+
  facet_grid(.~group)+
  theme(strip.text.x = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(axis.text=element_text(size=12))+ 
  theme(panel.background = element_rect(fill = 'gray93', colour = 'gray93')) +
  xlab ("") + ylab ("")


#### ANÀLISI DELS MARCADORS D'INFLAMACIÓ

# HEATMAP SEGONS EL GRUP DE PACIENT AMB VALORS INDIVIDUALS

bbddheatmap1<-bbddmarkers[order(bbddmarkers[,1]),]
m<-as.matrix(bbddheatmap1[,-1])

class_colors2<-rep(0,nrow(bbddheatmap1))

for(i in 1:nrow(bbddheatmap1)){
  if(bbddheatmap1[i,1]=='Healthy'){
    class_colors2[i]<-"steelblue1"
  }else if(bbddheatmap1[i,1]=='Compensated'){
    class_colors2[i]<-"steelblue2"
  }else if(bbddheatmap1[i,1]=='AD'){
    class_colors2[i]<-"steelblue3"
  }else if(bbddheatmap1[i,1]=='ACLF'){
    class_colors2[i]<-"steelblue4"
  }
}

colorm<-cbind(class_colors2)
colnames(colorm)<-c("")
windows()
heatmap.3(t(m), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row", Colv =FALSE,
          key=T, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=2, ColSideColors = colorm)
text(x=0.129858,y=0.8764035,     
     "H",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.1875589,y=0.8764035,     
     "CC",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.4735865,y=0.8764035,     
     "Acute Decompensation",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.8757181,y=0.8788514,     
     "ACLF",
     col = "white", 
     cex=.9,
     font=2
)

# HEATMAP SEGONS EL GRUP DE PACIENT AMB MEDIANES

bbddheatmap1_median<-bbddheatmap1 %>% group_by(group) %>% summarise_all(funs(median))
class_colors2<-rep(0,nrow(bbddheatmap1_median))

for(i in 1:nrow(bbddheatmap1_median)){
  if(bbddheatmap1_median[i,1]=='Healthy'){
    class_colors2[i]<-"steelblue1"
  }else if(bbddheatmap1_median[i,1]=='Compensated'){
    class_colors2[i]<-"steelblue2"
  }else if(bbddheatmap1_median[i,1]=='AD'){
    class_colors2[i]<-"steelblue3"
  }else if(bbddheatmap1_median[i,1]=='ACLF'){
    class_colors2[i]<-"steelblue4"
  }
}

colorm<-cbind(class_colors2)
colnames(colorm)<-c("")

mM<-as.matrix(bbddheatmap1_median[,-1])

windows()
heatmap.3(t(mM), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",Colv =FALSE,
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=3.5, ColSideColors = colorm)

text(x=0.2014702,y=0.8608905,     
     "Healthy",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.4244766,y=0.8608905,     
     "Compensated
Cirrhosis",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.651805,y=0.8608905,     
     "Acute 
Decompensation",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.8799546,y=0.8608905,     
     "ACLF",
     col = "white", 
     cex=.9,
     font=2
)


# ANÀLISI DE L'AGRUPACIÓ NATURAL DELS PACIENTS (grup de pacient)

windows()
heatmap.3(t(m), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=2, ColSideColors = colorm)


d<-dist(m)
h<-hclust(d,method="complete")

hclusters <- cutree(h, 4)
table(hclusters,bbddheatmap1$group)



# HEATMAP SEGONS FALLIDES ORGÀNIQUES EN AD AMB VALORS INDIVIDUALS

bbddheatmap2<-cbind(bbddmarkers,failure1c=dades$failure1c)
bbddheatmap2<-bbddheatmap2[c(which(bbddheatmap2$group=='AD'),which(bbddheatmap2$group=='ACLF')),]
bbddheatmap2$group<-droplevels(bbddheatmap2$group)
bbddheatmap2<-bbddheatmap2[order(bbddheatmap2$group,bbddheatmap2$failure1c),]
m<-as.matrix(bbddheatmap2[,-c(1,17)])

class_colors1<-rep(0,nrow(bbddheatmap2))
class_colors2<-rep(0,nrow(bbddheatmap2))

for(i in 1:nrow(bbddheatmap2)){
  if(bbddheatmap2[i,'group']=='AD'){
    class_colors1[i]<-"steelblue1"
  }else if(bbddheatmap2[i,'group']=='ACLF'){
    class_colors1[i]<-"steelblue3"
  }
}

for(i in 1:nrow(bbddheatmap2)){
  if(bbddheatmap2[i,'group']=='AD' && bbddheatmap2[i,'failure1c']=='No_OD_No_OF'){
    class_colors2[i]<-"red"
  }else if(bbddheatmap2[i,'group']=='AD' && bbddheatmap2[i,'failure1c']=='RD_CD'){
    class_colors2[i]<-"red3"
  }else if(bbddheatmap2[i,'group']=='AD' && bbddheatmap2[i,'failure1c']=='OF'){
    class_colors2[i]<-"red4"
  }else if(bbddheatmap2[i,'group']=='ACLF'){
    class_colors2[i]<-"white"
  }
}

colorm<-cbind(class_colors2,class_colors1)
colnames(colorm)<-c("","")

windows()
heatmap.3(t(m), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",Colv =FALSE,
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=3.5, ColSideColors = colorm)

text(x=0.392413,y=0.8892194,     
     "Acute Decompensation",
     col = "white", 
     cex=.9,
     font=2
)

text(x=0.8525732,y=0.8892194,     
     "ACLF",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.2246958,y=0.8307376,     
     "No OF-No OD",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.4668515,y=0.8307376,     
     "RD-CD",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.6476763,y=0.8307376,     
     "OF",
     col = "white", 
     cex=.9,
     font=2
)


# HEATMAP SEGONS FALLIDES EN AD AMB MEDIANES

bbddheatmap2_median<-bbddheatmap2 %>% group_by(group,failure1c) %>% summarise_all(funs(median))

class_colors1<-rep(0,nrow(bbddheatmap2_median))
class_colors2<-rep(0,nrow(bbddheatmap2_median))

for(i in 1:nrow(bbddheatmap2_median)){
  if(bbddheatmap2_median[i,'group']=='AD'){
    class_colors1[i]<-"steelblue1"
  }else if(bbddheatmap2_median[i,'group']=='ACLF'){
    class_colors1[i]<-"steelblue3"
  }
}

for(i in 1:nrow(bbddheatmap2_median)){
  if(bbddheatmap2_median[i,'group']=='AD' && bbddheatmap2_median[i,'failure1c']=='No_OD_No_OF'){
    class_colors2[i]<-"red"
  }else if(bbddheatmap2_median[i,'group']=='AD' && bbddheatmap2_median[i,'failure1c']=='RD_CD'){
    class_colors2[i]<-"red3"
  }else if(bbddheatmap2_median[i,'group']=='AD' && bbddheatmap2_median[i,'failure1c']=='OF'){
    class_colors2[i]<-"red4"
  }else if(bbddheatmap2_median[i,'group']=='ACLF'){
    class_colors2[i]<-"white"
  }
}

colorm<-cbind(class_colors2,class_colors1)
colnames(colorm)<-c("","")

mM<-as.matrix(bbddheatmap2_median[,-c(1,2)])

windows()
heatmap.3(t(mM), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",Colv =FALSE,
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=3.5, ColSideColors = colorm)

text(x=0.4125403,y=0.8910908,     
     "Acute Decompensation",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.8890654,y=0.8910908,     
     "ACLF",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.2084052,y=0.8334936,     
     "No OF-No OD",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.4237771,y=0.8334936,     
     "RD-CD",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.6633677,y=0.8334936,     
     "OF",
     col = "white", 
     cex=.9,
     font=2
)

# ANÀLISI DE L'AGRUPACIÓ NATURAL DELS PACIENTS (fallides AD)


windows()
heatmap.3(t(m), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=2, ColSideColors = colorm)

d<-dist(m)
h<-hclust(d,method="complete")

failures1<-rep(0,503)
for(i in 1:503){
  failures1[i]<-bbddheatmap2$failure1c[i]
  if(bbddheatmap2$group[i]=="ACLF"){
    failures1[i]<-4
  }
}

hclusters <- cutree(h, 4)
table(hclusters,failures1)


# HEATMAP SEGONS FALLIDES EN ACLF AMB VALORS INDIVIDUALS

bbddheatmap3<-cbind(bbddmarkers,failure5c=dades$failure5c)
bbddheatmap3<-bbddheatmap3[-which(is.na(bbddheatmap3$failure5c)),]
bbddheatmap3<-bbddheatmap3[c(which(bbddheatmap3$group=='AD'),which(bbddheatmap3$group=='ACLF')),]
bbddheatmap3<-bbddheatmap3[order(bbddheatmap3[,'group'],bbddheatmap3[,'failure5c']),]
m<-as.matrix(bbddheatmap3[,-c(1,17)])



class_colors1<-rep(0,nrow(bbddheatmap3))
class_colors2<-rep(0,nrow(bbddheatmap3))

for(i in 1:nrow(bbddheatmap3)){
  if(bbddheatmap3[i,'group']=='AD'){
    class_colors1[i]<-"steelblue1"
  }else if(bbddheatmap3[i,'group']=='ACLF'){
    class_colors1[i]<-"steelblue3"
  }
}

for(i in 1:nrow(bbddheatmap3)){
  if(bbddheatmap3[i,'group']=='AD' && bbddheatmap3[i,'failure5c']=='LF only'){
    class_colors2[i]<-"red4"
  }else if(bbddheatmap3[i,'group']=='ACLF' && bbddheatmap3[i,'failure5c']=='RF_No_CD'){
    class_colors2[i]<-"red1"
  }else if(bbddheatmap3[i,'group']=='ACLF' && bbddheatmap3[i,'failure5c']=='LF_RD_or_CD'){
    class_colors2[i]<-"red3"
  }else if(bbddheatmap3[i,'group']=='ACLF' && bbddheatmap3[i,'failure5c']=='RF_CD'){
    class_colors2[i]<-"red4"
  }
}

colorm<-cbind(class_colors2,class_colors1)
colnames(colorm)<-c("","")

windows()
heatmap.3(t(m), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",Colv =FALSE,
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=5.5, ColSideColors = colorm)

text(x=0.1981169,y=0.886444,     
     "Acute
     Decompensation",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.6391929,y=0.886444,     
     "ACLF",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.2066144,y=0.7935076,     
     "Liver
Failure",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.4933379,y=0.7935076,     
     "Renal Failure",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.7264566,y=0.7935076,     
     "Liver Failure
+RD/CD",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.9009839,y=0.7935076,     
     "Renal Failure
+CD",
     col = "white", 
     cex=.9,
     font=2
)

# HEATMAP FALLIDES EN ACLF AMB MEDIANES

bbddheatmap3_median<-bbddheatmap3 %>% group_by(group,failure5c) %>% summarise_all(funs(median))

class_colors1<-rep(0,nrow(bbddheatmap3_median))
class_colors2<-rep(0,nrow(bbddheatmap3_median))

for(i in 1:nrow(bbddheatmap3_median)){
  if(bbddheatmap3_median[i,'group']=='AD'){
    class_colors1[i]<-"steelblue1"
  }else if(bbddheatmap3_median[i,'group']=='ACLF'){
    class_colors1[i]<-"steelblue3"
  }
}

for(i in 1:nrow(bbddheatmap3_median)){
  if(bbddheatmap3_median[i,'group']=='AD' && bbddheatmap3_median[i,'failure5c']=='LF only'){
    class_colors2[i]<-"red4"
  }else if(bbddheatmap3_median[i,'group']=='ACLF' && bbddheatmap3_median[i,'failure5c']=='RF_No_CD'){
    class_colors2[i]<-"red1"
  }else if(bbddheatmap3_median[i,'group']=='ACLF' && bbddheatmap3_median[i,'failure5c']=='LF_RD_or_CD'){
    class_colors2[i]<-"red3"
  }else if(bbddheatmap3_median[i,'group']=='ACLF' && bbddheatmap3_median[i,'failure5c']=='RF_CD'){
    class_colors2[i]<-"red4"
  }
}

colorm<-cbind(class_colors2,class_colors1)
colnames(colorm)<-c("","")

mM<-as.matrix(bbddheatmap2_median[,-c(1,2)])

windows()
heatmap.3(t(mM), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",Colv =FALSE,
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=5.5, ColSideColors = colorm)

text(x=0.2073377,y=0.886444,     
     "Acute
Decompensation",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.638519,y=0.886444,     
     "ACLF",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.2015267,y=0.7954218,     
     "Liver
Failure",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.4384342,y=0.7954218,     
     "Renal Failure",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.6595375,y=0.7935076,     
     "Liver Failure
+RD/CD",
     col = "white", 
     cex=.9,
     font=2
)
text(x=0.8868865,y=0.7935076,     
     "Renal Failure
+CD",
     col = "white", 
     cex=.9,
     font=2
)


# ANÀLISI DE L'AGRUPACIÓ NATURAL DELS PACIENTS (grup de pacient)

windows()
heatmap.3(t(m), col=colorRampPalette(colors = c("darkturquoise","white","red"))(100), scale="row",
          key=TRUE, keysize=1,KeyValueName="",density.info="none", trace="none", cexRow =1.0,labCol  = F,
          ColSideColorsSize=2, ColSideColors = colorm)


d<-dist(m)
h<-hclust(d,method="complete")

hclusters <- cutree(h, 4)
table(hclusters,bbddheatmap3$failure5c)


### ANÀLISI DE SUPERVIVÈNCIA

bbdd_survival<-dades[,c(2,33,34,18:32,35:40)]
bbdd_survival<-subset(bbdd_survival,group %in% c('AD','ACLF'))
bbdd_survival$morttxcen_90dM<-factor(bbdd_survival$morttxcen_90dM,levels=levels(bbdd_survival$morttxcen_90dM)[c(2,3,1)])
bbdd_survival$group<-factor(bbdd_survival$group)
attach(bbdd_survival)

# MORTALITAT SEGONS EL GRUP DE PACIENT

table(bbdd_survival$group,bbdd_survival$morttxcen_90dM)
round(prop.table(table(bbdd_survival$group,bbdd_survival$morttxcen_90dM),1)*100,2)

# corbes d'incidència acumulada

cic_90d<-cuminc(timecomprisks_90days,morttxcen_90dM,bbdd_survival$group,cencode="Censured")
timepoints(cic_90d,times=c(30,60,90))

par(las = 1, font = 2, font.lab = 2, font.axis = 2,cex.axis=1,cex.lab=1)
plot.cuminc(list(cic_90d$`ACLF Dead`,cic_90d$`AD Dead`),lty=c(1,1),lwd=c(2,2),
            col=c("red","black"),xlab="Days",cex=1,curvlab=c("ACLF","AD"))
axis(1,at=seq(0,90,by=10))
text(x=74.41721,y=0.8997758,
     "Gray's Test p<0.001",
     col="Black",
     cex=1,
     font=2)

  #Test de Gray
  cic_90d$Tests 

  
# Model de Cox 

bbdd_survival$time_tr<-bbdd_survival$timecomprisks_90days
bbdd_survival$time_tr[bbdd_survival$morttxcen_90dM=="Transplanted"]<-90

Y<-Surv(time=bbdd_survival$time_tr,event=bbdd_survival$morttxcen_90dM=="Dead")
ph<-coxph(Y~strata(bbdd_survival$group))

coxm_d<-coxph(Y~bbdd_survival$group)

# Model de Fine i Gray
model_FG<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~group,cause="Dead",
              y=TRUE,data=bbdd_survival)


# Avaluació del supòsit de Riscos Proporcionals
  
  #plot
  windows()
  par(las = 1, font = 2, font.lab = 2, font.axis = 2,cex.axis=1,cex.lab=1)
  plot(survfit(ph),fun="cloglog", lty=c(1,1), col=c("red","black"),
       xlab="log Days",ylab="log-log S(t)", main="",las=1,cex=1,lwd=c(2,2))
  legend("topleft", c("ACLF", "AD"),lty=c(1,1), col=c("red","black"),lwd=c(2,2))
  text(x=5.249831,y=-0.8970935,
       "p=0.272",
       col="Black",
       cex=1,
       font=2)
  #test
  cox.zph(coxm_d)

# MORTALITAT SEGONS FALLIDES EN AD
  
bbdd_survival$fail1[bbdd_survival$failure1c=="No_OD_No_OF"]<-1
bbdd_survival$fail1[bbdd_survival$failure1c=="RD_CD"]<-2
bbdd_survival$fail1[bbdd_survival$failure1c=="OF"]<-3
bbdd_survival$fail1[bbdd_survival$group=="ACLF"]<-4
bbdd_survival$fail1<-factor(bbdd_survival$fail1,labels=c("No_OD_No_OF","RD_CD","OF","ACLF"))
attach(bbdd_survival)

table(bbdd_survival$fail1,bbdd_survival$morttxcen_90dM)
round(prop.table(table(bbdd_survival$fail1,bbdd_survival$morttxcen_90dM),1)*100,2)

# corbes d'incidència acumulada

cic_90d_f1<-cuminc(timecomprisks_90days,morttxcen_90dM,fail1,cencode="Censured")
timepoints(cic_90d_f1,times=c(30,60,90))


windows()
par(las = 1, font = 2, font.lab = 2, font.axis = 2,cex.axis=1,cex.lab=1)
plot.cuminc(list(cic_90d_f1$`No_OD_No_OF Dead`,cic_90d_f1$`RD_CD Dead`,cic_90d_f1$`OF Dead`,
                 cic_90d_f1$`ACLF Dead`),lty=c(1,1,1,1),lwd=c(2,2,2,2),col=c("black","green","blue","red"),
            xlab="Days",cex=1,curvlab=c("No OD-No OF","RD-CD","OF","ACLF"))
axis(1,at=seq(0,90,by=10))

text(x=65.24634,y=0.8764037,
     "Gray's Test p<0.001",
     col="Black",
     cex=1,
     font=2)

  #Test de Gray
  cic_90d_f1$Tests 
  
  cic_90d_f1.1<-cuminc(timecomprisks_90days,morttxcen_90dM,fail1 %in% c("No_OD_No_OF","RD_CD"),cencode="Censured")
  cic_90d_f1.1$Tests 
  cic_90d_f1.2<-cuminc(timecomprisks_90days,morttxcen_90dM,fail1 %in% c("No_OD_No_OF","OF"),cencode="Censured")
  cic_90d_f1.2$Tests 
  cic_90d_f1.3<-cuminc(timecomprisks_90days,morttxcen_90dM,fail1 %in% c("RD_CD","OF"),cencode="Censured")
  cic_90d_f1.3$Tests 
  cic_90d_f1.4<-cuminc(timecomprisks_90days,morttxcen_90dM,fail1 %in% c("RD_CD","ACLF"),cencode="Censured")
  cic_90d_f1.4$Tests
  cic_90d_f1.5<-cuminc(timecomprisks_90days,morttxcen_90dM,fail1 %in% c("OF","ACLF"),cencode="Censured")
  cic_90d_f1.5$Tests


# Model de Cox

ph_f1<-coxph(Y~strata(bbdd_survival$fail1))

coxm_f1_d<-coxph(Y~bbdd_survival$fail1)

# Model Fine i Gray

model_FG_f1<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~fail1,cause="Dead",
                 y=TRUE,data=bbdd_survival)

# Avaluació del supòsit de Riscos Proporcionals

  #plot
  windows()
  par(las = 1, font = 2, font.lab = 2, font.axis = 2,cex.axis=1,cex.lab=1)
  plot(survfit(ph_f1),fun="cloglog", lty=c(1,1,1,1), col=c("black","green","blue","red"),
       xlab="log Days",ylab="log-log S(t)", main="",las=1,cex=1,lwd=c(2,2,2,2))
  legend("topleft", c("No OD-No OF","RD-CD","OF","ACLF"),lty=c(1,1),lwd=c(2,2), col=c("black","green","blue","red"))
  text(x=12.84034,y=-0.9983021,
       "p=0.282",
       col="Black",
       cex=1,
       font=2)
  #test
  cox.zph(coxm_f1_d)

# MORTALITAT SEGONS FALLIDES EN ACLF
  
bbdd_survival2<-subset(bbdd_survival,bbdd_survival$failure5c %in% c('LF only','RF_No_CD','LF_RD_or_CD','RF_CD'))

table(bbdd_survival2$failure5c,bbdd_survival2$morttxcen_90dM)
round(prop.table(table(bbdd_survival2$failure5c,bbdd_survival2$morttxcen_90dM),1)*100,2)

# corbes d'incidència acumulada

cic_90d_f2<-cuminc(bbdd_survival2$timecomprisks_90days,bbdd_survival2$morttxcen_90dM,bbdd_survival2$failure5c,cencode="Censured")
timepoints(cic_90d_f2,times=c(30,60,90))

windows()
par(las = 1, font = 2, font.lab = 2, font.axis = 2,cex.axis=1,cex.lab=1)
plot.cuminc(list(cic_90d_f2$`LF only Dead`,cic_90d_f2$`RF_No_CD Dead`,cic_90d_f2$`LF_RD_or_CD Dead`,
                 cic_90d_f2$`RF_CD Dead`),lty=c(1,1,1,1),lwd=c(2,2,2,2),col=c("black","green","blue","red"),
            xlab="Days",cex=1,curvlab=c("LF (AD)","RF","LF+RD/CD","RF+CD"))
axis(1,at=seq(0,90,by=10))

text(x=62.77558,y=0.8881365,
     "Gray's Test p=0.143",
     col="Black",
     cex=1,
     font=2)

  #Test de Gray
  cic_90d_f2$Tests 
  
  cic_90d_f2.1<-cuminc(timecomprisks_90days,morttxcen_90dM,failure5c %in% c("LF only","RF_No_CD"),cencode="Censured")
  cic_90d_f2.2<-cuminc(timecomprisks_90days,morttxcen_90dM,failure5c %in% c("LF only","LF_RD_or_CD"),cencode="Censured")
  cic_90d_f2.3<-cuminc(timecomprisks_90days,morttxcen_90dM,failure5c %in% c("LF only","RF_CD"),cencode="Censured")
  cic_90d_f2.4<-cuminc(timecomprisks_90days,morttxcen_90dM,failure5c %in% c("RF_No_CD","LF_RD_or_CD"),cencode="Censured")
  cic_90d_f2.5<-cuminc(timecomprisks_90days,morttxcen_90dM,failure5c %in% c("RF_No_CD","RF_CD"),cencode="Censured")
  cic_90d_f2.6<-cuminc(timecomprisks_90days,morttxcen_90dM,failure5c %in% c("LF_RD_or_CD","RF_CD"),cencode="Censured")
  
  
  cic_90d_f2.1$Tests
  cic_90d_f2.2$Tests 
  cic_90d_f2.3$Tests 
  cic_90d_f2.4$Tests
  cic_90d_f2.5$Tests
  cic_90d_f2.6$Tests
  
# Model de Cox
  
Y<-Surv(time=bbdd_survival2$time_tr,event=bbdd_survival2$morttxcen_90dM=="Dead")

ph_f2<-coxph(Y~strata(bbdd_survival2$failure5c))

coxm_f2_d<-coxph(Y~bbdd_survival2$failure5c)

# Model de Fine i Gray

model_FG_f2<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~failure5c,cause="Dead",
                 y=TRUE,data=bbdd_survival2)

#Avaluació del supòsit de Riscos Proporcionals

  #plot
  windows()
  par(las = 1, font = 2, font.lab = 2, font.axis = 2,cex.axis=1,cex.lab=1)
  plot(survfit(ph_f2),fun="cloglog", lty=c(1,1,1,1), col=c("black","green","blue","red"),
       xlab="log Days",ylab="log-log S(t)", main="",las=1,cex=1,lwd=c(2,2,2,2))
  legend("topleft", c("LF (AD))","RF","LF+RD/CD","RF+CD"),lty=c(1,1),lwd=c(2,2), col=c("black","green","blue","red"))
  text(x=16.92501,y=-0.6253979,
       "p=0.689",
       col="Black",
       cex=1,
       font=2)
  #test
  cox.zph(coxm_f2_d)

  
# MORTALITAT EN FUNCIÓ DEL CONJUNT DE MARCADORS D'INFLAMACIÓ
  
bbddmarkers_survival<-cbind(bbddmarkers,timecomprisks_90days=dades$timecomprisks_90days,morttxcen_90dM=dades$morttxcen_90dM)
bbddmarkers_survival<-subset(bbddmarkers_survival,group %in% c('AD','ACLF'))
attach(bbddmarkers_survival)  

  # correlacions entre marcadors 
  round(cor(as.matrix(dades[18:32]),use="complete.obs"),2)
  
# MODELS UNIVARIATS DE FINE I GRAY
  
  model_FG_tnfa<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logTNFa,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_il6<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIL_6,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_il8<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIL_8,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_mcp1<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logMCP_1,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_ip10<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIP_10,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_mip1b<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logMIP_1b,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_gcsf<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logG_CSF,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_gmcsf<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logGM_CSF,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_il10<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIL_10,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_il1ra<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIL_1ra,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_ifng<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIFNg,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_eotaxin<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logEotaxin,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_il17a<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIL_17a,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_il7<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~logIL_7,cause="Dead",y=TRUE,data=bbddmarkers_survival)
  model_FG_hna2<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~HNA2,cause="Dead",y=TRUE,data=bbddmarkers_survival)

  # STEPWISE
  modelstep<-crrstep(timecomprisks_90days~logTNFa+logIL_6+logIL_8+logMCP_1+logIP_10+logMIP_1b+logG_CSF+logGM_CSF+
                       logIL_10+logIL_1ra+logIFNg+logEotaxin+logIL_17a+logIL_7+HNA2, ,morttxcen_90dM,failcode="Dead",data=bbddmarkers_survival,
                     criterion="AIC",direction="backward")
  
  modelfinal<-FGR(Hist(timecomprisks_90days,morttxcen_90dM,cens.code = "Censured")~
                    logIL_6+logIL_8+HNA2,
                  cause="Dead",y=TRUE,data=bbddmarkers_survival)
  
