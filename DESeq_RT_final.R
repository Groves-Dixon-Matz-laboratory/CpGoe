# source("http://bioconductor.org/biocLite.R") #To download DESeq package (you can comment these lines out, they only need to be run once ever)
# biocLite("DESeq")
# biocLite("DESeq2")

# source("http://bioconductor.org/biocLite.R")
# biocLite("genefilter")


library(DESeq) #To upload DESeq package - you need to do this every time you open R
library(gplots) # for venn diagram
library(ggplot2)
library(RColorBrewer)
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(plotrix)
library(reshape2)

counts=read.table("AllCountsBayMapAmilApr2014.tab",header=TRUE,row.names=1) #Reading in the table of counts per isogroup by sample - you must open R in the same folder or change the working directory

head(counts) 


length(counts[,1])  #45357 isogroups

counts=counts[,order(names(counts))]
names(counts)

readsleft=c()
for (column in names(counts)) {
	val=sum(counts[,column])
	readsleft=append(readsleft,val)}

RLtable=data.frame(cbind(names(counts),readsleft))
RLtable$readsleft=as.numeric(as.character(RLtable$readsleft))
write.csv(RLtable,"readsleft.csv",quote=F)

#######################Creating table of conditions for your experiment, origin and transplant site

origin=transplant=c(1:length(names(counts)))
origin[grep("KK",names(counts))]=transplant[grep("KK",names(counts))]="keppel"
origin[grep("OK",names(counts))]="orpheus"
transplant[grep("OK",names(counts))]="keppel"
origin[grep("KO",names(counts))]="keppel"
transplant[grep("KO",names(counts))]="orpheus"
origin[grep("OO",names(counts))]=transplant[grep("OO",names(counts))]="orpheus"

conditions=data.frame(cbind(origin,transplant))
head(conditions)

replicated=c("KK7","OO2","OO5","OO8")
counts.norep=counts
for (rc in replicated){
	reps=grep(paste(rc,"[AB]",sep=""),names(counts.norep)) 
	sr=rowSums(cbind(counts[,reps],counts[,rc]))
	counts.norep=counts.norep[,-reps]
	counts.norep[,rc]=sr
}
names(counts.norep)
names(counts)
head(counts.norep)
# removing tech reps from conditions list:
techreps=grep("[AB]",names(counts))
conditions.norep=conditions[-(techreps),]
# removing questionable samples - KK04 and OK33
counts.norep=counts.norep[,-c(10,38)]
conditions.norep=conditions.norep[-c(10,38),]


########################### esimating size factors

real=newCountDataSet(counts.norep,conditions.norep) 
real=estimateSizeFactors(real)

####all the data you ever wanted about quality control

cds=estimateDispersions(real,method="blind")
vsdBlind=varianceStabilizingTransformation(cds)

v="/Users/carlykenkel/Desktop/Collaborators/Line/AQM_v2"

arrayQualityMetrics(vsdBlind,outdir=v,intgroup=c("origin","transplant")) #check .html output file in new folder

##Host
counts.nobad=counts.norep[,-c(54)] #According to heatmap, Samples OO6 is an outlier, remove from counts and conditions
conditions.nobad=conditions.norep[-c(54),] 
conditions.nobad1=cbind(conditions.nobad,traitData[,2:8])
conditions.nobad1[,9]=as.factor(conditions.nobad1[,9])

#re-make count dataset
real=newCountDataSet(counts.nobad,conditions.nobad1) 
real=estimateSizeFactors(real)

real=estimateDispersions(real,method="pooled-CR",sharingMode="gene-est-only",modelFormula=count~origin+transplant)  #can use gene-est-only, have >7 reps per treatment
plotDispEsts(real)

vsd=getVarianceStabilizedData(real)
write.csv(vsd, file="VSD_allgenes_nobadsam_Aug2014.csv", quote=F)

######################determining quality filtering cutoff
fit0=fitNbinomGLMs(real, count ~ 1) # null model: expression does not depend on anything
fit1=fitNbinomGLMs(real, count ~ transplant)
fit1a=fitNbinomGLMs(real, count ~ origin)
fit2=fitNbinomGLMs(real, count ~ transplant+origin)
fit3=fitNbinomGLMs(real, count ~ transplant+origin+transplant:origin)

pvals.i<-nbinomGLMTest(fit3,fit2) #testing significance of interaction term
pvals.t<-nbinomGLMTest(fit1,fit0) #testing significance of transplant term
pvals.o<-nbinomGLMTest(fit1a,fit0) #testing significance of origin term
pvals.t.o<-nbinomGLMTest(fit2,fit1) #testing significance of o in addition to t
pvals.o.t<-nbinomGLMTest(fit2,fit1a) 


pvalue=pvals.t #change to each type of pval
theta=seq(from=0,to=0.8,by=0.02)

filterChoices=data.frame(`mean`=rowMeans(counts(real)),`median`=apply((counts(real)),1,median),`min`=rowMin(counts(real)),`max`=rowMax(counts(real)),`sd`=rowSds(counts(real)))
rejChoices=sapply(filterChoices,function(f) filtered_R(alpha=0.1,filter=f,test=pvalue,theta=theta,method="BH"))
library("RColorBrewer")
myColours=brewer.pal(ncol(filterChoices),"Set1")

#quartz()
matplot(theta,rejChoices,type="l",lty=1,col=myColours,lwd=2,xlab=expression(theta),ylab="number of rejections")
legend("bottomleft",legend=colnames(filterChoices),fill=myColours)

#look for peak in graph - corresponds to correct theta and best-fit line for which metric to use - pick best theta for all tests

#######################Quality Filtering Data - get rid of genes with low variance

#FOR host
rs=rowSds(counts(real)) #using standard deviation as quality filtering metric based on analyses above
theta=0.6 #peak rej now occur at theta of 0.4 using SD for all..
use=(rs>quantile(rs,probs=theta)) ###
table(use) 
# use
# FALSE  TRUE 
# 27215 18142 #have 18k genes that pass filter

realFilt=real[use,]
vsd=getVarianceStabilizedData(realFilt)
vsd2=varianceStabilizingTransformation(realFilt)

realFilt=real[use,]
vsd=getVarianceStabilizedData(realFilt)
vsd2=varianceStabilizingTransformation(realFilt)


########################Model Testing
fit0=fitNbinomGLMs(realFilt, count ~ 1) # null model: expression does not depend on anything
fit1=fitNbinomGLMs(realFilt, count ~ transplant)
fit1a=fitNbinomGLMs(realFilt, count ~ origin)
fit2=fitNbinomGLMs(realFilt, count ~ transplant+origin)
fit3=fitNbinomGLMs(realFilt, count ~ transplant+origin+transplant:origin)

#NAs in trait data, must remove those samples for individual trait models
useCarb=(!is.na(pData(realFilt)[,4]))
realCarb=realFilt[,useCarb]
fitC=fitNbinomGLMs(realCarb, count ~ CARB,)
useProt=(!is.na(pData(realFilt)[,5]))
realProt=realFilt[,useProt]
fitP=fitNbinomGLMs(realProt, count ~ PROTEIN)
useLip=(!is.na(pData(realFilt)[,6]))
realLip=realFilt[,useLip]
fitL=fitNbinomGLMs(realLip, count ~ LIPID)
useZoox=(!is.na(pData(realFilt)[,7]))
realZoox=realFilt[,useZoox]
fitZ=fitNbinomGLMs(realZoox, count ~ ZOOX)
fitA=fitNbinomGLMs(realFilt, count ~ AREA)
useGain=(!is.na(pData(realFilt)[,9]))
realGain=realFilt[,useGain]
fitG=fitNbinomGLMs(realGain, count ~ GAIN)

fitS=fitNbinomGLMs(realFilt, count ~ Survival)

#Lets see how the deviances look - the smaller the number, the more is explained by the model - can skip this if you just care abotu simple stats

#first total deviance
DevTot=data.frame(fit0$deviance-fit3$deviance)
names(DevTot)="DevTot"

#then, transplant
DevTrans=data.frame((fit0$deviance-fit1$deviance)/DevTot$DevTot)
names(DevTrans)="DevTrans"
hist(DevTrans$DevTrans,breaks=20)
mean(DevTrans$DevTrans) #.0.3801817
DevTrans$type<-'Transplant'

#then origin
DevOri=data.frame((fit0$deviance-fit1a$deviance)/DevTot$DevTot)
names(DevOri)="DevOri"
hist(DevOri$DevOri,breaks=20)
mean(DevOri$DevOri) #0.3790344
DevOri$type<-'Origin'

#plotting some pretty plots
DevBoth <- rbind(DevTrans, DevOri)
DevBothCor<-cbind(DevTrans,DevOri)
DevBoth$type=as.factor(DevBoth$type)
ggplot(DevBoth, aes(DevianceRatio, fill = type)) + geom_density(alpha = 0.2)

ggplot(DevBoth, aes(DevianceRatio, fill = type)) + geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')

plot(DevBothCor[,1]~DevBothCor[,3])

#taking only top significant genes from full interaction model and doing the same as above
d<-cbind(vsd,DevTrans,DevOri,DevTot, "pval.t" = pvals.t$pvals.t, "pval.o" = pvals.o$pvals.o,"pval.i" = pvals.i$pvals.i) 
d$X<-rownames(d)
write.csv(d, file="DevianceTableSigGenes.csv", quote=F)

ori=d[d$pval.o<0.01& !is.na(d$pval.o),]
tra=d[d$pval.t<0.01& !is.na(d$pval.t),]
int=d[d$pval.i<0.01& !is.na(d$pval.i),]
sig=data.frame(rbind(ori[!(ori$X %in% tra$X),],tra))
all=data.frame(rbind(sig[!(sig$X %in% int$X),],int))
length(sig[,1]) #2215 host; 1588 sym
length(all[,1])  #2602 host; 1899 sym
summary(all)

plot(DevOri~DevTrans,all)

nrow(all[all$DevTrans>0.1 & all$DevOri<0.1,]) #751 genes
nrow(all[all$DevOri>0.1 & all$DevTrans<0.1,]) #507 genes

p1 <- hist(all$DevTrans)                     # centered at 4
p2 <- hist(all$DevOri)                     # centered at 6
plot( p1, col=rgb(0,0,1,1/8), xlim=c(0,1),ylim=c(0,600))  # first histogram
plot( p2, col=rgb(0,1,0,1/8), xlim=c(0,1), ylim=c(0,600))  # second

df <- data.frame(all$DevTrans,all$DevOri)

ggplot(melt(df), aes(value, fill = variable)) + geom_histogram(position = "dodge")


# testing section
pvals.i<-nbinomGLMTest(fit3,fit2) #testing significance of interaction term
pvals.t<-nbinomGLMTest(fit1,fit0) #testing significance of transplant term
pvals.o<-nbinomGLMTest(fit1a,fit0) #testing significance of origin term
pvals.t.o<-nbinomGLMTest(fit2,fit1) #testing significance of o in addition to t
pvals.o.t<-nbinomGLMTest(fit2,fit1a) #testing significance of t in addition to o

pvals.C<-nbinomGLMTest(fitC,fit0)
pvals.P<-nbinomGLMTest(fitP,fit0)
pvals.L<-nbinomGLMTest(fitL,fit0)
pvals.Z<-nbinomGLMTest(fitZ,fit0)
pvals.A<-nbinomGLMTest(fitA,fit0)
pvals.G<-nbinomGLMTest(fitG,fit0)
pvals.S<-nbinomGLMTest(fitS,fit0)

#making non-convergent model p-values NA's
pvals.i=data.frame(pvals.i)
rownames(fit3)->rownames(pvals.i)
badmods=subset(fit3,(!fit3[,6]))
for ( i in rownames(badmods)){pvals.i[i,1]<-NA}

pvals.t=data.frame(pvals.t)
rownames(fit1)->rownames(pvals.t)
badmods=subset(fit1,(!fit1[,4]))
for ( i in rownames(badmods)){pvals.t[i,1]<-NA}

pvals.o=data.frame(pvals.o)
rownames(fit1a)->rownames(pvals.o)
badmods=subset(fit1a,(!fit1a[,4]))
for ( i in rownames(badmods)){pvals.o[i,1]<-NA}

pvals.t.o=data.frame(pvals.t.o)
rownames(fit2)->rownames(pvals.t.o)
badmods=subset(fit2,(!fit2[,5]))
for ( i in rownames(badmods)){pvals.t.o[i,1]<-NA}

pvals.o.t=data.frame(pvals.o.t)
rownames(fit2)->rownames(pvals.o.t)
badmods=subset(fit2,(!fit2[,5]))
for ( i in rownames(badmods)){pvals.o.t[i,1]<-NA}

pvals.C=data.frame(pvals.C)
rownames(fitC)->rownames(pvals.C)
badmods=subset(fitC,(!fitC[,4]))
for ( i in rownames(badmods)){pvals.C[i,1]<-NA}

pvals.P=data.frame(pvals.P)
rownames(fitP)->rownames(pvals.P)
badmods=subset(fitP,(!fitP[,4]))
for ( i in rownames(badmods)){pvals.P[i,1]<-NA}

pvals.L=data.frame(pvals.L)
rownames(fitL)->rownames(pvals.L)
badmods=subset(fitL,(!fitL[,4]))
for ( i in rownames(badmods)){pvals.L[i,1]<-NA}

pvals.Z=data.frame(pvals.Z)
rownames(fitZ)->rownames(pvals.Z)
badmods=subset(fitZ,(!fitZ[,4]))
for ( i in rownames(badmods)){pvals.Z[i,1]<-NA}

pvals.A=data.frame(pvals.A)
rownames(fitA)->rownames(pvals.A)
badmods=subset(fitA,(!fitA[,4]))
for ( i in rownames(badmods)){pvals.A[i,1]<-NA}

pvals.G=data.frame(pvals.G)
rownames(fitG)->rownames(pvals.G)
badmods=subset(fitG,(!fitG[,4]))
for ( i in rownames(badmods)){pvals.G[i,1]<-NA}

pvals.S=data.frame(pvals.S)
rownames(fitS)->rownames(pvals.S)
badmods=subset(fitS,(!fitS[,4]))
for ( i in rownames(badmods)){pvals.S[i,1]<-NA}

	
summary(pvals.i)
summary(pvals.t)
summary(pvals.o)
# summary(pvals.t.o)
# summary(pvals.o.t)
summary(pvals.C)
summary(pvals.P)
summary(pvals.L)
summary(pvals.Z)
summary(pvals.A)
summary(pvals.G)
summary(pvals.S)

#multiple test correction - adjust p-values 

adjp.i<-p.adjust(pvals.i$pvals.i,method="BH")
adjp.t<-p.adjust(pvals.t$pvals.t,method="BH")
adjp.o<-p.adjust(pvals.o$pvals.o,method="BH")
adjp.t.o<-p.adjust(pvals.t.o$pvals.t.o,method="BH")
adjp.o.t<-p.adjust(pvals.o.t$pvals.o.t,method="BH")

adjp.C<-p.adjust(pvals.C$pvals.C,method="BH")
adjp.P<-p.adjust(pvals.P$pvals.P,method="BH")
adjp.L<-p.adjust(pvals.L$pvals.L,method="BH")
adjp.Z<-p.adjust(pvals.Z$pvals.Z,method="BH")
adjp.A<-p.adjust(pvals.A$pvals.A,method="BH")
adjp.G<-p.adjust(pvals.G$pvals.G,method="BH")
adjp.S<-p.adjust(pvals.S$pvals.S,method="BH")


PooledPvals<-cbind(vsd, "adjp.t" = adjp.t, "adjp.o" = adjp.o, "adjp.i" = adjp.i,"adjp.C" = adjp.C,"adjp.P" = adjp.P,"adjp.L" = adjp.L,"adjp.Z" = adjp.Z,"adjp.A" = adjp.A,"adjp.G" = adjp.G,"adjp.S" = adjp.S) #creating table of all multiple test corrected p-values with variance stabilized count data

do<-cbind(PooledPvals, "pval.t" = pvals.t$pvals.t, "pval.o" = pvals.o$pvals.o,"pval.t.o" = pvals.t.o$pvals.t.o, "pval.o.t" = pvals.o.t$pvals.o.t,"pval.i" = pvals.i$pvals.i,"pval.C" = pvals.C$pvals.C,"pval.P" = pvals.P$pvals.P,"pval.L" = pvals.L$pvals.L,"pval.Z" = pvals.Z$pvals.Z,"pval.A" = pvals.A$pvals.A,"pval.G" = pvals.G$pvals.G,"pval.S" = pvals.S$pvals.S) #adding in original p-values

write.csv(do, file="VSDandPVALSnobadsam_SDtheta06_plusTraits_1sept.csv", quote=F) #writing an output file of vsd plus p-values

#Writing invdividual csv p-value and gene name csv files for GO analysis using gomwu script



#########Generating directional GO output for DESeq results

#transplant
fit1$direction=ifelse(fit1$transplantorpheus>0,1,-1)
fit1$pval<-(-log((pvals.t$pvals.t+0.0000000001),10))*(fit1$direction)

#origin
fit1a$direction=ifelse(fit1a$originorpheus>0,1,-1)
fit1a$pval<-(-log((pvals.o$pvals.o+0.0000000001),10))*(fit1a$direction)

#interaction
fit3$direction=ifelse(fit3[,4]>0,1,-1)
fit3$pval<-(-log((pvals.i$pvals.i+0.0000000001),10))*(fit3$direction)

#carb
fitC$direction=ifelse(fitC$CARB>0,1,-1)
fitC$pval<-(-log((pvals.C$pvals.C+0.0000000001),10))*(fitC$direction)
#prot
fitP$direction=ifelse(fitP$PROTEIN>0,1,-1)
fitP$pval<-(-log((pvals.P$pvals.P+0.0000000001),10))*(fitP$direction)
#lipid
fitL$direction=ifelse(fitL$LIPID>0,1,-1)
fitL$pval<-(-log((pvals.L$pvals.L+0.0000000001),10))*(fitL$direction)
#zoox
fitZ$direction=ifelse(fitZ$ZOOX>0,1,-1)
fitZ$pval<-(-log((pvals.Z$pvals.Z+0.0000000001),10))*(fitZ$direction)
#area
fitA$direction=ifelse(fitA$AREA>0,1,-1)
fitA$pval<-(-log((pvals.A$pvals.A+0.0000000001),10))*(fitA$direction)
#gain
fitG$direction=ifelse(fitG$GAIN>0,1,-1)
fitG$pval<-(-log((pvals.G$pvals.G+0.0000000001),10))*(fitG$direction)
#survival
fitS$direction=ifelse(fitS$Survival1>0,1,-1)
fitS$pval<-(-log((pvals.S$pvals.S+0.0000000001),10))*(fitS$direction)

interaction<-cbind("gene"=rownames(fit3),"pval"=fit3$pval)
transplant<-cbind("gene"=rownames(fit1),"pval"=fit1$pval) 
origin<-cbind("gene"=rownames(fit1a),"pval"=fit1a$pval) 
carb<-cbind("gene"=rownames(fitC),"pval"=fitC$pval) 
prot<-cbind("gene"=rownames(fitP),"pval"=fitP$pval) 
lipid<-cbind("gene"=rownames(fitL),"pval"=fitL$pval) 
zoox<-cbind("gene"=rownames(fitZ),"pval"=fitZ$pval) 
area<-cbind("gene"=rownames(fitA),"pval"=fitA$pval) 
growth<-cbind("gene"=rownames(fitG),"pval"=fitG$pval)
survival<-cbind("gene"=rownames(fitS),"pval"=fitS$pval) 
#transpluso<-cbind("pval"=pvals.t.o$pvals.t.o,"gene"=rownames(vsd)) 
#oriplust<-cbind("pval"=pvals.o.t$pvals.o.t,"gene"=rownames(vsd))



write.csv(interaction,file="GOinteractionAug14.csv",quote=F,row.names=F)
write.csv(transplant,file="GOtransplantAug14.csv",quote=F,row.names=F)
write.csv(origin,file="GOoriginAug14.csv",quote=F,row.names=F)

#write.csv(transpluso,file="GOtranspluso.csv",quote=F,row.names=F)
#write.c
write.csv(carb,file="trait_GO_CARBAug14.csv",quote=F,row.names=F)
write.csv(prot,file="trait_GO_PROTAug14.csv",quote=F,row.names=F)
write.csv(lipid,file="trait_GO_LIPIDAug14.csv",quote=F,row.names=F)
write.csv(zoox,file="trait_GO_ZOOXAug14.csv",quote=F,row.names=F)
write.csv(area,file="trait_GO_AREAAug14.csv",quote=F,row.names=F)
write.csv(growth,file="trait_GO_GAINAug14.csv",quote=F,row.names=F)
write.csv(survival,file="trait_GO_SURVIVALAug14.csv",quote=F,row.names=F)


check=read.csv("topGO_Eonly_Sig_HighExpOrpheus.csv")

########################## counting, venn diagram:

p=data.frame(do)
inter=row.names(p[p$adjp.i<=0.05 & !is.na(p$adjp.i),])
ori=row.names(p[p$adjp.o<=0.05 & !is.na(p$adjp.o),])
tran=row.names(p[p$adjp.t<=0.05 & !is.na(p$adjp.t),])
#o.t=row.names(p[p$adjp.o.t<=0.05 & !is.na(p$adjp.o.t),])
#t.o=row.names(p[p$adjp.t.o<=0.05 & !is.na(p$adjp.t.o),])

#ttran=union(tran,o.t)
#tori=union(ori,t.o)

p=data.frame(do)
p=read.csv("VSDandPVALSnobadsam_SDtheta06_plusTraits_1sept.csv")
#prot=row.names(p[p$adjp.P<=0.1 & !is.na(p$adjp.P),])
carb=row.names(p[p$adjp.C<=0.05 & !is.na(p$adjp.C),])
lip=row.names(p[p$adjp.L<=0.05 & !is.na(p$adjp.L),])
area=row.names(p[p$adjp.A<=0.05 & !is.na(p$adjp.A),])
zoox=row.names(p[p$adjp.Z<=0.05 & !is.na(p$adjp.Z),])
gain=row.names(p[p$adjp.G<=0.05 & !is.na(p$adjp.G),])
surv=row.names(p[p$adjp.S<=0.05 & !is.na(p$adjp.S),])

candidates=list("site of outplanting"=tran,"site of origin"=ori,"Interaction"=inter)
candidates2=list("carb"=carb,"lipid"=lip,"zoox"=zoox,"gain"=gain,"survival"=surv)
candidates3=list("origin"=ori,"transplant"=tran,"gain"=gain,"surv"=surv)
candidates4=list("carb"=carb,"zoox"=zoox,"gain"=gain,"origin"=ori,"transplant"=tran)
library(gplots)
venn(candidates2)

##########################################Pulling out top most differentially expressed genes for WGCNA

d=read.csv("VSDandPVALSnobadsam_SDtheta06.csv")
rownames(d)<-d$X
ori=d[d$pval.o<0.01& !is.na(d$pval.o),]
tra=d[d$pval.t<0.01& !is.na(d$pval.t),]
int=d[d$pval.i<0.01& !is.na(d$pval.i),]
sig=data.frame(rbind(ori[!(ori$X %in% tra$X),],tra))
all=data.frame(rbind(sig[!(sig$X %in% int$X),],int))
length(sig[,1]) #2774 genes
length(all[,1])  #2873 genes
summary(all)

rest=d[!(d$X %in% all$X),]

rest=rest[,1:2]
names(rest)=c("ProbeID","moduleColors")
write.csv(rest,"GenesNotinWGCNA.csv",quote=F,row.names=F)

wgcna=all[,1:57]
write.csv(wgcna,file="Wgcna_2873genes_SigOnlyRawPvals001.csv",quote=F,row.names=F)

##########################################################pulling out gene pattern by expression types

#website - refresher on 2-way anova plots http://www4.uwsp.edu/psych/stat/13/anova-2w.htm

############Genotype Effects Only

AP=cbind(fit1a,"adjp.o"=adjp.o,"adjp.t"=adjp.t,"adjp.i"=adjp.i) #17198 genes in comparisson 

AP=AP[complete.cases(AP),] #remove all rows where fit did not converge, 17195 isogroups left

AP$HighO<-(AP$originorpheus>0)
summary(AP) #7785

AP$HighOSig<-(AP$adjp.o<=0.1 & AP$adjp.t>0.1 & AP$adjp.i>0.1)

AP$moduleColor<-(AP$HighO & AP$HighOSig) 
summary(AP) #82 groups where origin pval <0.1

AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0

AP$gene=rownames(AP)

APCLO=AP[,10:11]

write.csv(APCLO,file="topGO_Gonly_Sig_HighExpOrpheus.csv",quote=F,row.names=F) #writing file for plotting these sig genes

AP$LowO=AP$originorpheus <0 #9410 genes

AP$LowOSig=(AP$adjp.o<=0.1 & AP$adjp.t>0.1 & AP$adjp.i>0.1)

AP$moduleColor<-(AP$LowO & AP$LowOSig) #83 genes where pval O < 0.1

AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0

APCLO=AP[,10:11]
write.csv(APCLO,file="topGO_Gonly_Sig_HighExpKeppels.csv",quote=F)

############Environment Effects Only

AP=cbind(fit1,"adjp.o"=adjp.o,"adjp.t"=adjp.t,"adjp.i"=adjp.i) #

AP=AP[complete.cases(AP),] #17195

AP$HighO=(AP$transplantorpheus>0)
summary(AP) #6978 isogroups
AP$HighOSig=(AP$adjp.t<=0.1 & AP$adjp.o>0.1 & AP$adjp.i>0.1)
AP$moduleColor<-(AP$HighO & AP$HighOSig) 
summary(AP) #506 genes where origin p-val is <0.1; 

AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0

AP$gene=rownames(AP)
APCLO=AP[,10:11]
write.csv(APCLO,file="topGO_Eonly_Sig_HighExpOrpheus.csv",quote=F,row.names=F) #writing file for plotting these sig genes

AP$LowO=AP$transplantorpheus<0
summary(AP) #10217 genes overall

AP$LowOSig=(AP$adjp.t<=0.1 & AP$adjp.o>0.1 & AP$adjp.i>0.1)
AP$moduleColor<-(AP$LowO & AP$LowOSig) #1215 genes where origin pval <0.1 *******************almost double the opposite pattern...?

AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0
#APCLO=p[(rownames(p) %in% rownames(LowOSig)),] #pulling out vsd data for these 476 significant genes for plotting

APCLO=AP[,10:11]
write.csv(APCLO,file="topGO_Eonly_Sig_HighExpKeppels.csv",quote=F,row.names=F)

############Genotype + Environment Effects
AdapPtrn=cbind(fit2,"adjp.o"=adjp.o,"adjp.o.t"=adjp.o.t,"adjp.i"=adjp.i) #

AdapPtrn=AdapPtrn[complete.cases(AdapPtrn),] #remove all rows where fit did not converge, 17185 isogroups left

AdapPtrn$testVal=(apply((AdapPtrn[,2:3]),1,prod)) #create product column to find direction of interaction

AdapPtrn$APC=(AdapPtrn$testVal>0)
summary(AdapPtrn) #9284 genes total - "correct adaptation pattern"

AdapPtrn$LowO=(AdapPtrn$transplantorpheus<0) #low exp offshore
summary(AdapPtrn) 
AdapPtrn$LowOSig=(AdapPtrn$adjp.o<=0.1 & AdapPtrn$adjp.o.t<=0.1 & AdapPtrn$adjp.i>0.1)
AdapPtrn$moduleColor<-(AdapPtrn$APC & AdapPtrn$LowO & AdapPtrn$LowOSig)#19 genes where both origin and origin + transplant pvals are <0.1 and interaction is NS
AdapPtrn$moduleColor[AdapPtrn$moduleColor]<-1
AdapPtrn$moduleColor[!AdapPtrn$moduleColor]<-0

AdapPtrn$gene=rownames(AdapPtrn)
APCLO=AdapPtrn[,13:14]
write.csv(APCLO,file="topGO_Correct_GplusE_Sig_LowExpOrpheus.csv",quote=F,row.names=F) #writing file for plotting these sig genes

AdapPtrn$LowK=(AdapPtrn$transplantorpheus>0) #low exp inshore
summary(AdapPtrn) #6982 genes total
AdapPtrn$LowKSig=(AdapPtrn$adjp.o<=0.1 & AdapPtrn$adjp.o.t<=0.1 & AdapPtrn$adjp.i>0.1)
AdapPtrn$moduleColor<-(AdapPtrn$APC & AdapPtrn$LowK & AdapPtrn$LowKSig)#10 genes where both origin and origin + transplant pvals are <0.1
AdapPtrn$moduleColor[AdapPtrn$moduleColor]<-1
AdapPtrn$moduleColor[!AdapPtrn$moduleColor]<-0
APCLO=AdapPtrn[,13:14]
write.csv(APCLO,file="topGO_Correct_GplusE_Sig_LowExpKeppels.csv",quote=F,row.names=F)

AdapPtrn$APW=(AdapPtrn$testVal<0)
summary(AdapPtrn) #7901 genes show "wrong pattern"

AdapPtrn$moduleColor<-(AdapPtrn$APW & AdapPtrn$LowO & AdapPtrn$LowOSig) #4 genes where both origin and origin + transplant pvals are <0.1
AdapPtrn$moduleColor[AdapPtrn$moduleColor]<-1
AdapPtrn$moduleColor[!AdapPtrn$moduleColor]<-0
APCLO=AdapPtrn[,13:14]

write.csv(APCLO,file="topGO_Wrong_GplusE_Sig_LowExpOrpheus.csv",quote=F,row.names=F)

AdapPtrn$moduleColor<-(AdapPtrn$APW & AdapPtrn$LowK & AdapPtrn$LowKSig)  #1 genes where both origin and origin + transplant pvals are <0.1
AdapPtrn$moduleColor[AdapPtrn$moduleColor]<-1
AdapPtrn$moduleColor[!AdapPtrn$moduleColor]<-0
APCLO=AdapPtrn[,13:14]

write.csv(APCLO,file="topGO_Wrong_GplusE_Sig_LowExpKeppels.csv",quote=F,row.names=F)

##################Interaction Term Effects
#GxE only - interaction term pvals weak; using unadjusted to explore trends.....

AP=cbind(fit3,"adjp.o"=adjp.o,"adjp.t"=adjp.t,"adjp.i"=pvals.i$pvals.i) #22626 genes in comparisson (those passing quality filter)
AP=AP[complete.cases(AP),] #17195
 
AP$GE_upreg=(AP[,4]>0)
summary(AP) #9924 genes total

AP$GE_upreg_sig=(AP$adjp.i<=0.005 & AP$adjp.o>0.1 & AP$adjp.t>0.1)
AP$moduleColor<-(AP$GE_upreg & AP$GE_upreg_sig) 
summary(AP)#233 where interaction unadj pval is <=0.01; 120 <=0.005; 49 <=0.01 

AP$gene=rownames(AP)
AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0
APCLO=AP[,12:13]

write.csv(APCLO,file="topGO_GE_Sig_UpregAtHome_pvalsi_005.csv",quote=F,row.names=F)

AP$GE_dwnreg=(AP[,4]<0)
 #7271 genes total

AP$moduleColor<-(AP$GE_dwnreg & AP$GE_upreg_sig) 
summary(AP) #159 where interaction pval is <=0.01; 82 <=0.005; 
AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0
APCLO=AP[,12:13]
write.csv(APCLO,file="topGO_GE_Sig_DownregAtHome_pvalsi_005.csv",quote=F,row.names=F)

#G+E+GxE all significant
head(AP)
AP$all_sig=(AP$adjp.i<=0.05 & AP$adjp.o<0.2 & AP$adjp.t<0.2) 
summary(AP)#might want to play with pval cutoffs here, this gives 25 candidates.  interaction pvals are actually raw pvals 
AP$moduleColor=AP$all_sig
AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0
APCLO=AP[,12:13]
write.csv(APCLO,file="topGO_G_E_GE_AllSig_NoParticularPattern.csv",quote=F,row.names=F)

#G+GxE; E not significant
head(AP)
AP$all_sig=(AP$adjp.i<=0.05 & AP$adjp.o<0.2 & AP$adjp.t>0.1)
summary(AP) #might want to play with pval cutoffs here, this gives 36 candidates.  
AP$moduleColor=AP$all_sig
AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0
APCLO=AP[,12:13] 
write.csv(APCLO,file="topGO_G_GE_AllSig_NoParticularPattern.csv",quote=F,row.names=F)

head(AP)
AP$all_sig=(AP$adjp.i<=0.05 & AP$adjp.o>0.1 & AP$adjp.t<0.2)
summary(AP) #might want to play with pval cutoffs here, this gives 392 candidates.  
AP$moduleColor=AP$all_sig
AP$moduleColor[AP$moduleColor]<-1
AP$moduleColor[!AP$moduleColor]<-0
APCLO=AP[,12:13]  
write.csv(APCLO,file="topGO_E_GE_AllSig_NoParticularPattern.csv",quote=F,row.names=F)


####################################
