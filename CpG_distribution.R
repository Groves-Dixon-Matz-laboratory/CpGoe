##--------------- DOUBLECHECKING AGAINST DIFFERENT DATASETS -----------------------------
#as a double check Misha wanted me to use this dataset in place of the expression datasets 
#used to build figure 4
setwd('/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/cpg_OE')
load("LBay_RT_Allgenes_Ori_Trans.RData")
head(allOrigin)
head(allTransplant)
write.table(allOrigin, "originOrph.txt", quote = F)
write.table(allTransplant, "transplantOrph.txt", quote = F)
###save both of these as tables so that you can change the isogroup notation to isogroup=12345 to match with this R script
##then save the revised data tables as originOrphIn.txt and transplantOrphIn.txt
originOrph = read.table("originOrphIn.txt", header = T)
transplantOrph = read.table("transplantOrphIn.txt", header = T)
originOrph$EST = rownames(originOrph) #so that it can be merged with cgm dataframe
transplantOrph$EST = rownames(transplantOrph) #so that it can be merged with cgm dataframe
head(originOrph)
head(transplantOrph)
####====================================================================================
#IMPORT DATA
# cgm = read.table('/Users/grovesdixon/Documents/lab_files/CpGoe_Project/circadian_rhythm/milleporaCpGData.txt',header=T)###new one after redoing the blasting and such DOES THE GENE BODIES
# cgm = read.table('/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/amil_CDS_CpG.txt', header = T)
cgm = read.table('/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/amil_CDS_CpG_sub1000.txt', header = T)
head(cgm)
#================================================================================
length(cgm$EST)
x = unique(cgm$EST)
length(x)
##-------------- USE ANNOTATIONS TO FILTER -------------------------
cgma = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/MergedIso2spAss5-14-14_apr2014_BioP.txt")
length(cgma[,1])
cgma = na.omit(cgma)
colnames(cgma)=c("EST","GOs")#; colnames(cgda) = c("EST")
cgm = merge(cgm, cgma, by = "EST")
#cgd = merge(cgd, cgda, by = "EST")
cgma = ''; cgda = ''
cgm = cgm[,1:length(cgm)-1] ##remove the GOs so visualization is easier
#===============================================================================
##------------- FILTER BASED ON SIZE AND CPGO/E --------------------------------
par(mar = c(5, 4, 4, 2) + 0.1)
minlength = 300
maxlength = 20000 ###original value for gene bodies was 20000
minOE = .001  ##.001
maxOE = 2
filter = function(cg, minlength, maxlength, maxOE, minOE){
  cg$cpgOE = (cg$CpG/(cg$C*cg$G))*(cg$length^2/(cg$length-1))##equation in Gavery and Roberts (oysters)
  cg$gpcOE = (cg$GpC/(cg$C*cg$G))*(cg$length^2/(cg$length-1))
  cg$tpgOE = (cg$TpG/(cg$T*cg$G))*(cg$length^2/(cg$length-1))
  too.high = length(cg[cg$cpgOE > maxOE,]$EST)
  too.low = length(cg[cg$cpgOE < minOE,]$EST)
  too.short = length(cg[cg$length < minlength,]$EST)
  print(paste("Too Short Transcript length", too.short))
  print(paste("Too high CpG gene count =",too.high))
  print(paste("Too low CpG gene count =",too.low))
  cg = cg[cg$length > minlength,]
  cg = cg[cg$length < maxlength,]
  cg = cg[cg$cpgOE<maxOE,]
  cg = cg[cg$cpgOE>minOE,] ##what are these?
  cg = na.omit(cg)
  return(cg)
}
cgm = filter(cgm, minlength, maxlength, maxOE, minOE)
#cgd = filter(cgd, minlength, maxlength, maxOE, minOE)
hist(cgm$cpgOE,breaks=30, main = 'millepora')
#hist(cgd$cpgOE,breaks=30, main = 'digitifera')
length(cgm[,1])
##==============================================================================
write.table(cgm,'/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/CGM_cpg_scores8-22-14__CDS.txt',row.names=FALSE, quote = F)
##--------------------- MIXTURE MODELING ---------------------------------------
##use this section to manipulate the model parameters and find the best ones the best parameters are hard-coded in the section above so that figures can be made with set colors.
library(mixtools)
emmix = function(cg,species,means,sigma,lambda,k){
  x = cg$cpgOE
  cpgNull <- normalmixEM(x, mu = means, sigma = sigma, lambda = lambda, k=k, arbvar=T)
  summary(cpgNull)
  plot(cpgNull, which = 2, breaks = 30, density = TRUE, cex.axis = 1.4, cex.lab = 1.5, cex.main = 1.5, main2 = species, xlab2 = "CpGo/e") 
  return(cpgNull)
}
par(mfrow=c(1,1))
##NULLs
mil2mix = emmix(cgm,"millepora",NULL,NULL,NULL,2)
#dig2mix = emmix(cgd,"digitifera",NULL,NULL,NULL,2)
model_df = function(lambda,mean,sigma){
  df = data.frame(cbind(lambda,mean,sigma))
  rownames(df) = paste("comp",rownames(df),sep="")
  return(df)
}##assigns the mixture model to a variable
mmix2 = model_df(mil2mix$lambda,mil2mix$mu,mil2mix$sigma)
##============================================================================
##------------ TESTING FOR NUMBER OF COMPONENTS --------------
par(mfrow = c(1,1))
x.norm = pnorm(cgm$cpgOE, mean = 0.704, sd = 0.2394679)
x = sample(cgm$cpgOE, 5000, replace = F)
shapiro.test(x)
#reject the null hypothesis that the distribution is normal
require(mclust)
x = cgm$cpgOE
clust = Mclust(x, G = c(1:3), modelNames = c("V"))
bic = mclustBIC(x, G = c(1:5), modelNames = c("V"))
plot(bic, G = c(1:5), modelNames = c("V"))
dev.copy2pdf(file = "/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/figure_outputs/BICfigure.pdf", width = 170/25.4, height = 170/25.4, onefile = T)
##----------- ASSIGN AND PLOT CHOSEN MIXTURE MODELS -------------------------MGP = c(2,.5,0)
MGP = c(3,1,0)
set.y = c(0,1,2,3)
plot_mix = function(dat, mixmod, species, figure, reduction.factor, At, MGP, maxY){
  x = dat$cpgOE
  N = length(dat$cpgOE)
  if (figure == FALSE){
    hist.stats = hist(x,breaks=20, main=species, xlab=expression("GpC"["O/E"]), prob = TRUE)
    make_comp = function(comp_num){
      comp = function(x){
        mixmod$lambda[comp_num]*dnorm(x, mixmod$mean[comp_num], mixmod$sigma[comp_num])
      }
      return(comp)
    }
  }
  if (figure == TRUE){
    hist.stats = hist(x, breaks=30, xlab="CpGO/E", ylab = paste(reduction.factor,"s of Genes", sep = ""), freq = F, plot = F)
    hist.stats$counts = hist.stats$counts/reduction.factor
    plot(hist.stats, axes = F, main = "", xlab=expression("CpG"["O/E"]), ylab = paste(reduction.factor,"s of Genes", sep = ""), freq = TRUE, mgp = MGP, xlim=c(0,1.5), ylim = c(0, maxY))
    axis(1, mgp = MGP)
    axis(2, at = At, labels = At, las = 1, mgp = MGP)
    NORMALIZER = max(hist.stats$density)
    print(NORMALIZER)
    BREAKS = (length(hist.stats$breaks))
    make_comp = function(comp_num){
      comp = function(x){
        mixmod$lambda[comp_num]*dnorm(x, mixmod$mean[comp_num], mixmod$sigma[comp_num])/1.32##I have no idea why this is the thing to put there but it matches the density plot
      }
      return(comp)
    }##builds a component as a function of x based on the component's parameters which are found in mixmod
  }
  colors = c("yellowgreen","red")
  bells = data.frame(x)
  for (i in 1:length(mixmod$mean)){
    comp = make_comp(i)
    curve(comp,add = TRUE,col=colors[i],lwd=2)
    bell_values = data.frame(comp(x))
    colnames(bell_values) = paste("comp",i,sep="")
    bells = cbind(bells,bell_values)
  }##function to plot the normal curves for each component
  pnormmix <- function(x,mixture){
    lambda <- mixture$lambda
    k <- length(mixture$lambda)
    pnorm.from.mix <- function(x,component){
      lambda[component]*dnorm(x,mean=mixture$mean[component],sd=mixture$sigma[component])
    }
    pnorms <- sapply(1:k,pnorm.from.mix,x=x)
    return(rowSums(pnorms))
  }
  #curve(pnormmix(x,mixmod),add=T,col="purple",lty=2,lwd=2)
  bell_values = data.frame(pnormmix(x,mixmod))
  colnames(bell_values) = "cumulative"
  bells = cbind(bells,bell_values)
  return(bells)
}##function takes input mixture model parameters, plots the components' curves and the cumulative curve, and returns a dataframe of the x and y coordiantes for each curve.
m2 = plot_mix(cgm, mmix2, "", TRUE, 1000, set.y, MGP, 1)
# z = hist(cgm$cpgOE, prob = T)
# z$counts
# t = max(z$density)
# abline(h = t)
##===========================================================================
##--------------- GET BOUNDARIES FOR COMPONENTS -----------------------------
head(m2)
m2$diffs = (m2$comp1 - m2$comp2)^2
n.r = nrow(m2)
m2.middle = m2[with(m2, order(x)),]
m2.middle = m2.middle[2000:(n.r-5000),]
separator = m2.middle[m2.middle$diffs == min(m2.middle$diffs), "x"]
left.bound = min(cgm$cpgOE)
right.bound = max(cgm$cpgOE)
mc1bounds = list(left.bound, separator)
mc2bounds = list(separator, right.bound)
mbounds = list(mc1bounds, mc2bounds)
##===========================================================================
##--------------- PLOT FIGURE 1 ----------------------------------------
par(mar = c(5, 4, 4, 2) + 0.1)
CEX = 0.9
par(mfrow=c(1,1))
par(mfrow=c(1,3))
MGP = c(3,1,0)
set.y = c(0,.5,1)
m2 = plot_mix(cgm,mmix2,"", TRUE, 1000, set.y, MGP, 1)
abline(v = separator, col = "black", lty = 2, lwd = 1.5)
head(cgm)
mtext("A", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
GpC = cgm$gpcOE[cgm$gpcOE < 1.7]
cont.hist = hist(GpC, breaks = 20, plot = F)
cont.hist$counts = cont.hist$counts/1000
MGP = c(2,.75,0)
plot(cont.hist, main = "", xlab = expression("GpC"["O/E"]) , ylab = "1000s of Genes", mgp = MGP, las = 1, axes = F, ylim = c(0,2))
axis(1, at = c(0.6, 1, 1.5), cex = CEX, mgp = MGP)
axis(2, at = c(0,1,2), las = 1, cex = CEX, mgp = MGP)
mtext("B", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
plot(tpgOE~cpgOE, data = cgm, axes = F, xlab = expression("CpG"["O/E"]), ylab = expression("TpG"["O/E"]), pch = 19, cex = 0.3, xlim = c(0,1.5), ylim = c(.7, 1.8), mgp = MGP)
lm1 = lm(tpgOE~cpgOE, data = cgm) ; summary(lm1)
abline(lm1, col = "red", lwd = 1.5)
axis(1, at = c(0, 0.75, 1.5), labels = c('0','0.75', '1.50'), mgp = MGP, cex = CEX)
axis(2, at = c(0.8,1.2, 1.6), mgp = MGP, las = 1, cex = CEX) #axis(2, at = c(1,1.4, 1.8), mgp = MGP, las = 1, cex = CEX)
mtext("C", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
dev.copy2pdf(file = "/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/figure_outputs/figure1.pdf", width = 170/25.4, height = 85/25.4, onefile = T)
#=========================================================================
#-------------------- PLOT BIO PROCESSES -----------------------------------------------
# dat = read.csv("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/mgi_GO_slim/bioPinputMergedIso2spAss5-14-14_BP_amilV2.txt", header = T)
# dat =  read.csv("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/mgi_GO_slim/bioPinputMergedIso2spAss5-14-14_amil_sep2013_iso2go.tab", header = F) #this one is my favorite
#dat = read.csv("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/bioPinputMergedIso2spAss7-27-14_BP_amilV2_NolengthFilter.txt", header = T)##THIS ONE DOES NOT FILTER ON MINIMUM LENGTH
dat = read.csv("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/bioPinputMergedIso2spAss7-27-14_BP_amilV2_1000_lengthFilter_CDSmode.txt", header = T)##THIS ONE DOES FILTER ON MINIMUM LENGTH OF 1000
colnames(dat) = c("bp","cpg")
dat$cpg = as.numeric(as.character(dat$cpg))
##GATHER THE DATA FOR PLOTTING THE ERROR BARS
ps = levels(dat$bp)
means = c()
ses = c()
n = c()
for (i in ps){
  sub = dat[dat$bp == i,]
  mn = mean(sub$cpg)
  se = 1.96 * sd(sub$cpg) / sqrt(length(sub$cpg))
  ses = append(ses, se)
  means = append(means, mn)
  n = append(n, length(sub$cpg))
}
head(dat)
##do the chi square tests
ps = levels(dat$bp)
means = c()
ses = c()
n = c()
pvalues = c()
redobs = c()
greenobs = c()
red.other = c()
green.other = c()
green = cgm[cgm$cpgOE <= separator,]
red = cgm[cgm$cpgOE > separator,]
g.tot = length(green[,1])
r.tot = length(red[,1])
for (process in ps){
  print(process)
  sub = dat[dat$bp == process,]
  mn = mean(sub$cpg)
  se = 1.96 * sd(sub$cpg) / sqrt(length(sub$cpg))
  ses = append(ses, se)
  means = append(means, mn)
  n = append(n, length(sub$cpg))
  print(length(sub[,1]))
  gsub = sub[sub$cpg <= separator,]##subet out the greens for this process
  rsub = sub[sub$cpg > separator,]##subset out the reds
  g.process.count = length(gsub[,1])
  r.process.count = length(rsub[,1])
  print(g.process.count)
  print(r.process.count)
  r.other.count = r.tot - r.process.count
  g.other.count = g.tot - g.process.count
  associated = c(g.process.count, r.process.count)
  not.associated = c(g.other.count, r.other.count)
  M = as.table(rbind(associated, not.associated))
  colnames(M) = c('green', 'red')
  rownames(M) = c(process, paste('not',process))
  #   print(M)
  p = fisher.test(M)$p.value
  pvalues = append(pvalues, p)
  greenobs = append(greenobs, g.process.count)
  redobs = append(redobs, r.process.count)
  red.other = append(red.other, r.other.count)
  green.other = append(green.other, g.other.count)
}
adj.pvals = p.adjust(pvalues, method = "BH", n = length(pvalues))
pvals = data.frame(ps, means, ses, n, greenobs, green.other, redobs, red.other, pvalues, adj.pvals)

##VERTICAL BARGRAPH
require(Hmisc)
tab <- pvals
tab = tab[tab$n>50,] #remove any bio processes with fewer than 20 genes
tab = tab[with(tab, order(means)),]
s = tab$adj.pvals
chiSqr_sigs = c()
for (x in s){
  i = as.numeric(x)
  print(i)
  if (i > 0.05){
    chiSqr_sigs = append(chiSqr_sigs, "")
    next
  }
  if (i > 0.01){
    chiSqr_sigs = append(chiSqr_sigs, '*')
    next
  }
  if (i > 0.001){
    chiSqr_sigs = append(chiSqr_sigs, "**")
    next
  }
  chiSqr_sigs = append(chiSqr_sigs, "***")
} #for loop that build the set of asterisks for chi sqaur significance
ROWNAMES = tab$ps
# bp = barplot(tab2$means,ylim=c(0.,.9),main="",names=F, plot = F)
# text(bp, par("usr")[3] - .05, srt = 30, adj = 1,labels = rownames(tab2),xpd = T)
lower = tab$means - tab$ses
upper = tab$means + tab$ses
# z = errbar(bp[,1], tab$means,upper,lower)
##========= PLOT NORMAL FIGURE ================
par(mfrow=c(1,1))
bp = barplot(tab$means,xlim=c(0.,.9),main="",names=F, horiz=T, plot = F)
library(plotrix)
par(mar=c(3,18.25,1,1))
MGP = c(1.75, .5, 0)
p = plotCI(tab$means, 1:length(tab$means), xlim = c(0.46, .8), uiw = tab$ses, err="x",yaxt="n",ylab = "",xlab = expression("CpG"["O/E"]), pch=19, cex = 0.5, lwd=2, sfrac = .0075, main="", axes = F, mgp = MGP, cex.lab = 1.1)
axis(2,at = p$y,labels = F, tick = T, mgp = MGP)
axis(1, at = c(.475, .575, .675, .775), labels = T, mgp = MGP)
y.labels = paste(ROWNAMES," (",tab$n,")",chiSqr_sigs, sep = "")
mtext(y.labels, side = 2, line = .75, at = p$y, outer=F,las=1, cex = 1.1)
dev.copy2pdf(file = "/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/figure_outputs/figure2.pdf", width = 170/25.4, height = 120/25.4, onefile = T)
#==================== PLOT STRESS FIGURE =====================================================
#=========================================================================
dat = read.csv("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/bioPinputMergedIso2spAss7-27-14_BP_amilV2_1000_lengthFilter_CDSmode_STRESS.txt", header = T)##USE THIS ONE TO BREAK DOWN STRESS
head(dat)
##set up some subsets
oxy = dat[dat$bp == "response to oxidative stress",]
wound = dat[dat$bp == "response to wounding",]
cellular = dat[dat$bp == "cellular response to stress", ]
stress = dat[dat$bp == "stress response",]
no.oxy = subset(stress, !stress$cpg %in% oxy$cpg)
##compare the means of the subsets
t.test(no.oxy$cpg, stress$cpg, alternative = 'greater')
t.test(no.oxy$cpg, oxy$cpg, alternative = 'greater') ###oxidative stress is significantly lower than stress response with oxidative stress removed.
##set up subset without oxidative stress genes or cellular response to stress genes
no.house = subset(no.oxy, !no.oxy$cpg %in% cellular$cpg)
no.house2 = subset(no.house, !no.house$cpg %in% wound$cpg)
house = rbind(oxy, cellular)
t.test(no.house$cpg, house$cpg, alternative = 'greater') #combined cellular and oxidative is WAY lower than all other stress response
no.house$bp = paste(no.house$bp, "reduced")
dat = rbind(dat, no.house) ##append this new subsetted process to the data frame with CpG values
##DO THE FISHER TESTS
ps = levels(dat$bp)
means = c()
ses = c()
n = c()
pvalues = c()
redobs = c()
greenobs = c()
red.other = c()
green.other = c()
green = cgm[cgm$cpgOE <= separator,]
red = cgm[cgm$cpgOE > separator,]
g.tot = length(green[,1])
r.tot = length(red[,1])
for (process in ps){
  print(process)
  sub = dat[dat$bp == process,]
  mn = mean(sub$cpg)
  se = 1.96 * sd(sub$cpg) / sqrt(length(sub$cpg))
  ses = append(ses, se)
  means = append(means, mn)
  n = append(n, length(sub$cpg))
  print(length(sub[,1]))
  gsub = sub[sub$cpg <= separator,]##subet out the greens for this process
  rsub = sub[sub$cpg > separator,]##subset out the reds
  g.process.count = length(gsub[,1])
  r.process.count = length(rsub[,1])
  print(g.process.count)
  print(r.process.count)
  r.other.count = r.tot - r.process.count
  g.other.count = g.tot - g.process.count
  associated = c(g.process.count, r.process.count)
  not.associated = c(g.other.count, r.other.count)
  M = as.table(rbind(associated, not.associated))
  colnames(M) = c('green', 'red')
  rownames(M) = c(process, paste('not',process))
  #   print(M)
  p = fisher.test(M)$p.value
  pvalues = append(pvalues, p)
  greenobs = append(greenobs, g.process.count)
  redobs = append(redobs, r.process.count)
  red.other = append(red.other, r.other.count)
  green.other = append(green.other, g.other.count)
}
adj.pvals = p.adjust(pvalues, method = "BH", n = length(pvalues))
pvals = data.frame(ps, means, ses, n, greenobs, green.other, redobs, red.other, pvalues, adj.pvals)

############========= SET UP THE DATA TABLE FOR PLOT ===============================
tab <- pvals
tab = tab[tab$n>=20,] #remove any bio processes with fewer than 20 genes
tab = tab[with(tab, order(means)),]
s = tab$pvalues
chiSqr_sigs = c()
for (x in s){
  i = as.numeric(x)
  print(i)
  if (i > 0.05){
    chiSqr_sigs = append(chiSqr_sigs, "")
    next
  }
  if (i > 0.01){
    chiSqr_sigs = append(chiSqr_sigs, '*')
    next
  }
  if (i > 0.001){
    chiSqr_sigs = append(chiSqr_sigs, "**")
    next
  }
  chiSqr_sigs = append(chiSqr_sigs, "***")
}
ROWNAMES = tab$ps
##----- DRAW THE PLOT ------------
par(mfrow=c(1,1))
par(mar=c(3,15,1,0))
MGP = c(1.75, .5, 0)
p = plotCI(tab$means, 1:length(tab$means), xlim = c(0.48, .80), uiw = tab$ses, err="x",yaxt="n",ylab = "",xlab = expression("CpG"["O/E"]), pch=19, cex = 0.5, lwd=2, sfrac = .0075, main="", axes = F, mgp = MGP, cex.lab = 1.1)
axis(2,at = p$y,labels = F, tick = T, mgp = MGP)
axis(1, at = c(.5, .6, .7,  .8), labels = T, mgp = MGP)
y.labels = paste(ROWNAMES, " (", tab$n, ")", chiSqr_sigs, sep = "")
mtext(y.labels, side = 2, line = .75, at = p$y, outer=F,las=1)
dev.copy2pdf(file = "/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/figure_outputs/StressBarFigure.pdf", width = 170/25.4, height = 120/25.4, onefile = T)
###----------------------- TUKEYS HSD --------------------------------
head(dat)
biop_aov <- aov(cpg ~ bp, data = dat)
anova(biop_aov)
z = TukeyHSD(biop_aov)$bp
z[z[,'p adj'] < 0.05,]
#--------------- PLOT MEAN EXPRESSION RELATIONSHIP TO CPGOE --------------
########## GRAB THE DATASETS ############
##get Carly's dataset
tdat = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/transplant/VSD_lbRT_nodup_nomiso__reduced.csv", header =T)
head(tdat)
tdat = merge(tdat, cgm, by = "EST")
length(tdat$EST)
##### SLIDING WINDOW FUNCTION #########
expression_plot = function(size, dat, cut){
  windows = seq(0,2, by = size)
  windows = quantile(dat$cpgOE, probs = seq(0, 1, by = size))
  print(windows)
  mn = c()
  x = c()
  N = c()
  sterr = c()
  cut.count = 0
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat$cpgOE>left,]
    sub = sub[sub$cpgOE<right,]
    n = length(sub[,1])
    if (n < cut){
      cut.count = cut.count + 1
      next
    }
    mn = append(mn,mean(sub$grandmeans))
    sterr = append(sterr,std.error(sub$grandmeans))
    x = append(x,windows[i])
    N = append(N, n)
  }
  print(paste(cut.count, "Windows were cut Because they had below", cut, "genes"))
  plot_dat = data.frame(mn,x,sterr, N)
  return(plot_dat)
}##function to get the window data to plot the curves size is the span of the windows, dat is the dataset, cut is the number of genes that has to be in a window for it to be plotted. It works by a 'tiling window'. The windows do not overlap. The size of the window will thus influence how many there are. Because the range of CpGoe is 0 to 2 I just use 0.02 to give 100 possible windows. Many are cut because they don't have enough genes.
anova.df = function(size, dat, cut){
	#windows = seq(0,2, by = size)
	windows = quantile(dat$cpgOE, probs = seq(0, 1, by = size))
	cut.count = 0
	expression = c()
	bins = c()
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat$cpgOE>left,]
    sub = sub[sub$cpgOE<right,]
    n = length(sub[,1])
    print(n)
    if (n < cut){
      cut.count = cut.count + 1
      next
    }
    sub$bin = i
    expression = append(expression, sub$grandmeans)
    bins = append(bins, sub$bin)
  }
  print(paste(cut.count, "Windows were cut Because they had below", cut, "genes"))
  bins = as.factor(bins)
  anova_dat = data.frame(expression, bins)
  return(anova_dat)
}
######---- PLOT FIGURE 3 tdat -------------------------  xlab=expression("CpG"["O/E"]), ylab="Gene Expression",
par(mfrow = c(1,1))
MGP = c(3, 1, 0)
#par(mfrow = c(1,1))
par(mar = c(5, 4, 4, 2) + 0.1)
MGP = c(1.45, 1, 0)
CEX = .75
Cutoff = 50
YLIM = c(4.21,4.635)
# YLIM = NULL
window_data = expression_plot(.04, tdat, Cutoff) ## data for Carly's RNA seq
window_data = na.omit(window_data)##remove the windows that had no genes in them
plot(mn~x,data=window_data, main="", pch = 1, cex = CEX, axes = F, mgp = MGP, cex.lab = CEX, xlab = NA, ylab = NA, ylim = YLIM)
plotCI(window_data$x, window_data$mn, uiw = window_data$sterr, liw = window_data$sterr, add = T, cex = CEX, scol = 'grey', lwd = 0.5)
mn = round(mean(window_data$mn), digits = 1)
axis(1, at = c(0.2, 0.6, 1.0), cex.axis = CEX, mgp = c(1, .4 , 0))
mtext(expression("CpG"["O/E"]), side=1, line=1.6, outer=F, cex= CEX, font=1, adj = .48)
AT = c(4.3, 4.4, 4.5, 4.6)
axis(2, at = AT, label = T, tick = T, las = 1, cex.axis = CEX, mgp = c(1.8, .6, 0))
#mtext("Gene Expression", side=2, line=1.7, outer=F, cex= CEX, font=1, adj = .5)
loess_fit <- loess(mn ~ x, window_data, span = .8, se = T)
lines(window_data$x, predict(loess_fit),col="red",lwd=1)
green = mmix2$mean[1]
red = mmix2$mean[2]
int = separator
xs = c(green, int, red)
Y = YLIM[1] - 0.012
ys = c(Y, Y, Y)
colors = c("yellowgreen", "black", "red")
points(xs, ys, pch = 17, col = colors, cex = 2)
#mtext("A", side=3, line=1, adj = 0, outer=F, cex=1.3, font=1)
#axis(1, at = c(0.2, 0.6, 1.0), cex.axis = CEX, mgp = c(1, .4 , 0))
#mtext(expression("CpG"["O/E"]), side=1, line=1.4, outer=F, cex= CEX, font=1, adj = .48)
dev.copy2pdf(file = "/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/figure_outputs/Figure4.pdf", width = 85/25.4, height = 95/25.4, onefile = T)
##RUN comparisons for bins
adat = anova.df(.04, tdat, Cutoff)
head(adat)
z = aov(expression~bins, data = adat)
head(adat)
anova(z)
head(tdat)
green = tdat[tdat$cpgOE < separator,]
red = tdat[tdat$cpgOE > separator,]
t.test(green$grandmeans, red$grandmeans)

#===========================================================================================
##-------------- ENRICHMENT PLASTIC GENES FROM TRANSPLANTS ---------------
tdat = read.table("/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/transplant/VSD_lbRT_nodup_nomiso__reduced.csv", header =T)
tdat = merge(tdat, cgm, by = "EST")
#function to perform fisher exact test for enrichment of differentially expressed genes in one category or another. 
#Arguments: dat = the main dataset; condition_dat = the expression dataset; column = the column name for the variable of interest in the conditional dataset; alpha = the p value cutoff; bounds1 = bounds for category1; bounds2 = bounds for 2.

###FISHER TEST FUNCTION
fisher_test = function(dat, exp.dat, p.col, cpg.col, alpha, separator){
	green = exp.dat[exp.dat[,cpg.col] < separator,]
	print("number of green genes with expression data:")
	green.count = length(green[,1])
	print(green.count)
	red = exp.dat[exp.dat[,cpg.col] > separator,]
	red.count = length(red[,1])
	print("number of red genes with expression data:")
	print(red.count)
	green.sig = green[green[,p.col] <= alpha,]
	g.sig.count = length(green.sig[,1])
	green.other = green[green[,p.col] > alpha,]
	g.other.count = length(green.other[,1])
	red.sig = red[red[,p.col] <= alpha,]
	r.sig.count = length(red.sig[,1])
	red.other = red[red[,p.col] > alpha,]
	r.other.count = length(red.other[,1])
	sig = c(g.sig.count, r.sig.count)
	other = c(g.other.count, r.other.count)
	M = as.table(rbind(sig, other))
	colnames(M) = c("green", "red")
	print(M)
	print(fisher.test(M, alternative = 'less'))
}
#TESTS USING OLD DATASET
fisher_test(cgm, tdat, 'adjp.t', 'cpgOE', 0.01, separator)
fisher_test(cgm, tdat, 'adjp.o', 'cpgOE', 0.01, separator)
######################### FISHER TEST FOR DOUBLE CHECK ########################
#origin differential expression
head(originOrph)
oot = merge(originOrph, cgm, by = "EST")
head(oot)
length(oot$EST)
fisher_test(cgm, oot, "adj.p", 'cpgOE', 0.1, separator)
fisher_test(cgm, oot, "adj.p", 'cpgOE', 0.05, separator)
fisher_test(cgm, oot, "adj.p", 'cpgOE', 0.01, separator)
fisher_test(cgm, oot, "adj.p", 'cpgOE', 0.001, separator)
#transplantation differential expression
head(transplantOrph)
ott = merge(transplantOrph, cgm, by = "EST")
head(ott)
length(ott$EST)
fisher_test(cgm, ott, "adj.p", 'cpgOE', 0.1, separator)
fisher_test(cgm, ott, "adj.p", 'cpgOE', 0.05, separator)
fisher_test(cgm, ott, "adj.p", 'cpgOE', 0.01, separator)
fisher_test(cgm, ott, "adj.p", 'cpgOE', 0.001, separator)
## so the story does not really change, but the origin genes do a little better. I also can't do the 0.05 or 0.1 genes
#====================================================================
##------------- PLOT FIGURE 4 WITH LOG TRANSFORMATION -------------------------
window_plot = function(size, dat, alpha, column, sig.cutoff, n.cutoff){
  transformed.alpha = -log(alpha, 10)
  sig.dat = dat[dat[,column] > transformed.alpha,]
  sig.tot.ratio = length(sig.dat$EST)/length(dat$EST)
  print("Ratio of significant Genes:")
  #print(sig.tot.ratio)
  print(paste("Number of Significant Genes:",length(sig.dat$EST), sep = ' '))
  windows = quantile(dat$cpgOE, probs = seq(0,1, size))
  #print(windows)
  mn = c()
  x = c()
  n = c()
  n.sig = c()
  stdv = c()
  sig.ratio = c()
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat$cpgOE>left,]
    sub = sub[sub$cpgOE<right,]##subset out the percentile
    sub.n = length(sub$cpgOE)
    sub.sig = sub[sub[,column] > transformed.alpha,]
    sub.n.sig = length(sub.sig[,1])
    n = append(n, sub.n)
    n.sig = append(n.sig, sub.n.sig)
    mn = append(mn,mean(sub[,column]))
    stdv = append(stdv, sd(sub[,column]))
    x = append(x, mean(sub$cpgOE))
    sig.ratio = append(sig.ratio, sig.tot.ratio)
  }
  plot_dat = data.frame(mn,x,stdv,n,n.sig, sig.ratio)
  plot_dat = plot_dat[plot_dat$n.sig >= sig.cutoff,] # filter quantiles with too few significant loci
  plot_dat$o.e = plot_dat$n.sig/(plot_dat$sig.ratio*plot_dat$n)
  plot_dat = plot_dat[plot_dat$n >= n.cutoff,] #filter quantiles with too few total loci
  plot_dat$log.o.e = log(plot_dat$o.e)
  return(plot_dat)
}
#============================================================
##================================================================
####getting DESeq
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq")
##----------------------------------------------------------------
load('/Users/grovesdixon/Documents/lab_files/CpGoe_Project/Data-Analysis_files/lbay_OrpheusKeppelsTransplant_filt.rdat')
write.table(fit1a, "origin_effects.tab", sep = '\t')
write.table(fit1, "transplant_effects.tab", sep = '\t')
## edit the table so that isogroup labels match how we have them here
##-------- LOOK AT TRANSPLANT EFFECT --------------------
head(cgm)
tran = read.table("transplant_effects_revised.txt")
tran$EST = rownames(tran)
colnames(tran) = c("int", "t.eff", "dev", "con", "EST")
head(tran)
t.dat = merge(cgm, tran, by = "EST")
head(t.dat)
length(t.dat$EST)
t.dat$env = abs(t.dat$t.eff)
# plot(env~cpgOE, data = t.dat, ylim = c(0,2), xlim = c(0, 1.5), pch=1,cex=.5,col="grey50", ylab = "hey")
lm1 = lm(env~cpgOE, data = t.dat); summary(lm1); abline(lm1, col = "red")
par(mar = c(5, 4, 4, 2) + 0.1)
spearmans.cor = cor(t.dat$cpgOE, t.dat$env, method = "spearman")
spear.test = cor.test(t.dat$cpgOE, t.dat$env, method = "spearman")
spear.test
#---------- LOOK AT ORIGIN EFFECT -------------------
ori = read.table("origin_effect_revised.txt", header = T)
o.dat = merge(cgm, ori, by = "EST")
colnames(o.dat)[12:15] = c("int", "o.eff", "dev", "con")
o.dat$env = abs(o.dat$o.eff)
head(o.dat)
o.dat = o.dat[o.dat$env < 20,] #remove two crazy outliers
# plot(env~cpgOE, data = o.dat, ylim = c(0,2), xlim = c(0, 1.5), pch=1,cex=.5,col="grey50", ylab = "hey")
lm1 = lm(env~cpgOE, data = o.dat); summary(lm1); abline(lm1, col = "red")
par(mar = c(5, 4, 4, 2) + 0.1)
spearmans.cor = cor(t.dat$cpgOE, o.dat$env, method = "spearman")
spear.test = cor.test(t.dat$cpgOE, o.dat$env, method = "spearman")
spear.test
#---------- SET UP PLOTTING FUNCTIONS -------------------
plot_sig_counts3 = function(dat, title, color, setY, xaxis, loess.span, MGP, YLAB, Y.AT, scale){
  plot(n.sig~x, data = dat, xlim = c(.1,1.25), ylab = YLAB, xlab = expression("CpG"["O/E"]), main = title, axes = F, pch = 1, mgp = MGP)
  axis(1, at = c(0.3, 0.6, 0.9, 1.2))
  axis(2, at = Y.AT, las = 1)
  make_comp = function(comp_num, mixmod){
    comp = function(x){
      mixmod$lambda[comp_num]*dnorm(x,mixmod$mean[comp_num],mixmod$sigma[comp_num])*(scale * max(dat$n.sig)) + (min(dat$n.sig))
    }
    return(comp)
  }
  green = make_comp(1, mmix2)
  red = make_comp(2, mmix2)
  curve(green, add = TRUE, col = "yellowgreen", lty = 2)
  curve(red, add = T, col = "red", lty = 2)
  trim.dat = dat[1: (length(dat[,1])),] ##
  loess_fit <- loess(n.sig ~ x, trim.dat, span = loess.span)
  lines(trim.dat$x, predict(loess_fit),col=color,lwd=1)
  box()
}

# ###############################plotting 4f as smaller file size
envSmall_plot = function(size, dat, cut, xaxis, yaxis, Xlim, Ylim, mnCEX, t.f, pt.size){
  P = plot(abs(dat[,yaxis])~dat[,xaxis], ylim = Ylim, xlim = Xlim, pch=t.f, cex=pt.size, xlab = XLAB, ylab = YLAB, las = 1, mgp = MGP, axes = F, col = "grey")
  windows = quantile(dat[,xaxis], probs = seq(0,1,(1/size))) #size is number of error bars you want
  mn = c()
  x = c()
  sterr = c()
  stdv = c()
  N = c()
  for (i in 1:(length(windows)-1)){
    left = windows[i]
    right = windows[i+1]
    sub = dat[dat[,xaxis]>left,]
    sub = sub[sub[,xaxis]<right,]
    n = length(sub[,1])
    if (n < cut){
      next
    }
    stdv = append(stdv, sd(sub[,yaxis]))
    mn = append(mn,mean(sub[,yaxis]))
    sterr = append(sterr, std.error(sub[,yaxis]))
    x = append(x,mean(sub[,xaxis]))
    N = append(n, length(sub[,yaxis]))
  }
  plot_dat = data.frame(mn,x,sterr,stdv, n)
  round = 2
  plot_dat2 = data.frame(x = round(plot_dat$x, round), mn = round(plot_dat$mn, round)) 
  print(head(plot_dat))
  print(head(plot_dat2))
  plot_dat <- plot_dat[!duplicated(plot_dat2),]
  CIup= plot_dat$mn + plot_dat$sterr
  CIdown = plot_dat$mn - plot_dat$sterr
  p = plotCI(plot_dat$x, plot_dat$mn, uiw = plot_dat$sterr, pch=19, cex = pt.size*1.5, lwd=1.5, sfrac = .003, main="", mgp = MGP, add = T, ylab = "ENV")
  points(plot_dat$x, plot_dat$mn, pch = 19, cex = mnCEX)
  return(plot_dat)
}###this function is exactly the same as the one above, but it removes datapoints that cover each other up so that the figure looks the same, but gives a smaller pdf file size
##------- PLOT FIGURE 4 NEW DATASETS -----------------------
ott$log2.difference = abs(ott$log2.difference)
oot$log2.difference = abs(oot$log2.difference)
par(mfrow = c(2,3))
yaxis = c(.5,1.0,1.5,2.0,2.5)
xaxis = c(.3,.6,.9,1.2)
MGP = c(2.5, .5, 0)
Size = 0.04
Span = 6
Num.Err = 12
XLAB = expression("CpG"["O/E"])
YLAB = "Environmentally Flexible Genes"
##############plot 4a
SCALE = 0.58
Y.AT = c(6, 12, 18, 24)
Y.AT = NULL
ott$log.t = -log(ott$adj.p, 10)
dat = window_plot(Size, ott, 0.01, "log.t",  0, 0)
plot_sig_counts3(dat, '', 'red', NULL, xaxis, Span, MGP, YLAB, Y.AT, SCALE)
mtext("A", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
Span = 3
mtext("A", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
##############plot 4b
YLAB = 'Environment Effect'
XLIM = c(.04, 1.5)
YLIM = ylim = c(0,2.25)
MGP = c(2.5, .5, 0)
black = envSmall_plot(Num.Err, ott, 20, "cpgOE", "log2.difference", XLIM, YLIM, 0.5, 1, .2)
axis(1, at = c(0.2, 0.7, 1.2))
axis(2, las = 1)
lm1 = lm(env~cpgOE, data = t.dat);summary(lm1); abline(lm1, col = 'red', lwd = 1)
MN = mean(t.dat$env)
mtext("B", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
##############plot 4c
Num.Err = 12
XLIM = c(0.15, 1.1)
YLIM = ylim = c(0.198, 0.36)
xs = c(XLIM[1], XLIM[1], XLIM[2], XLIM[2])
ys = c(YLIM[1], YLIM[2], YLIM[1],YLIM[2])
black = envSmall_plot(Num.Err, ott, 20, "cpgOE", "log2.difference", XLIM, YLIM, 0, 1, .2)
green = black[black[,"x"] < separator,]
#plotCI(green$x, green$mn, uiw = green$sterr, pch=".",lwd=3, col = "yellowgreen", sfrac = .004, main="", axes = F, mgp = MGP, add = T)
axis(1, at = c(.2, .6, 1))
axis(2, at = c(.2, .24, .28, .32, .36), las = 1, mgp = MGP)
green = mmix2$mean[1]
red = mmix2$mean[2]
int = separator
xs = c(green, int, red)
y = YLIM[1] - 0.0026
ys = c(y, y, y)
colors = c("yellowgreen", "black", "red")
points(xs, ys, pch = 24, bg = colors, col = 'white', cex = 2.4)
mtext("C", side=3, line=1, adj = 0, outer=F, cex= 1, font=1)
box()
#dev.copy2pdf(file = "/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/figure_outputs/Figure3.pdf", width = 170/25.4, height = 60/25.4, onefile = T)
#=================================================
XLAB = expression("CpG"["O/E"])
YLAB = 'Origin-Specific Genes'
MGP = c(2.5, .5, 0)
#############plot 4d
SCALE = 0.81
Y.AT = c(0, 2, 4, 6, 8)
oot$log.o = -log(oot$adj.p, 10)
dat = window_plot(Size, oot, 0.01, "log.o", 0, 0)
plot_sig_counts3(dat, '', 'red', NULL, xaxis, Span, MGP, YLAB, Y.AT, SCALE)
mtext("D", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
#############plot 4e
YLAB = 'Origin Effect'
XLIM = c(.04, 1.5)
YLIM = ylim = c(0,2.25)
MGP = c(2.5, .5, 0)
black = envSmall_plot(Num.Err, oot, 20, "cpgOE", "log2.difference", XLIM, YLIM, 0.5, 1, .2)
axis(1, at = c(0.2, 0.7, 1.2))
axis(2, las = 1)
lm1 = lm(env~cpgOE, data = t.dat);summary(lm1); abline(lm1, col = 'red', lwd = 1)
MN = mean(t.dat$env)
mtext("E", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
#############plot 4f
Num.Err = 12
XLIM = c(0.15, 1.1)
YLIM = ylim = c(.185, 0.32)
xs = c(XLIM[1], XLIM[1], XLIM[2], XLIM[2])
ys = c(YLIM[1], YLIM[2], YLIM[1],YLIM[2])
black = envSmall_plot(Num.Err, oot, 20, "cpgOE", "log2.difference", XLIM, YLIM, 0, 1, .2)
green = black[black[,"x"] < separator,]
#plotCI(green$x, green$mn, uiw = green$sterr, pch=".",lwd=3, col = "yellowgreen", sfrac = .004, main="", axes = F, mgp = MGP, add = T)
axis(1, at = c(.2, .4, .6, .8, 1, 1.2))
axis(2, at =c(.2, .24, .28, .32), las = 1, mgp = MGP)
green = mmix2$mean[1]
red = mmix2$mean[2]
int = separator
xs = c(green, int, red)
y = YLIM[1] - 0.0024
ys = c(y, y, y)
colors = c("yellowgreen", "black", "red")
points(xs, ys, pch = 24, bg = colors, col = 'white', cex = 2.4)
mtext("F", side=3, line=1, adj = 0, outer=F, cex=.9, font=1)
box()
dev.copy2pdf(file = "/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis/figure_outputs/Figure3.pdf", width = 170/25.4, height = 120/25.4, onefile = T)