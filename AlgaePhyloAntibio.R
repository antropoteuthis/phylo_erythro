##Microalgal Resistance Evolution##

rm(list=ls())
par(ask=F)

###Load libraries
library(reshape2)
library(dplyr)
library(tidyverse)
library(stringr)
library(vegan)
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(phylobase)
library(fields)
library(phylosignal)
library(picante)
library(geomorph)
library(FactoMineR)
library(factoextra)
library(phylopath)
library(Rphylip)
library("arbutus", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")

###Set working directory
setwd("~/Dropbox/Microalgas-antibioticos/")
getwd()

###Load raw data
read.csv("base_de_datos/erythro_data.csv") -> raw_data
read.csv("base_de_datos/controls.csv") -> controls
rownames(controls) = controls$Species
controls = controls[,-1]
names(controls) = c("cGrowth","cEQY", "cChlA", "cETR")
  #Extract the mean effect across concentrations
keyvars = c("Species","dGrowth", "dEQY", "dChlA", "dETR")
d_data = raw_data[,keyvars]
d_byspp = split(d_data, d_data$Species)
d_means = keyvars
d_SEs = keyvars
std.err <- function(a){sd(a)/sqrt(length(a))}
d_byspp$`Nannochloropsis gaditana` <- d_byspp$`Nannochloropsis gaditana`[-12,]
for(i in 1:length(d_byspp)){
  d_byspp[[i]] %>% .[,-1] %>% apply(.,2,FUN=mean) %>% c(names(d_byspp[i]),.) %>% rbind(d_means,.) -> d_means
  d_byspp[[i]] %>% .[,-1] %>% apply(.,2,FUN=std.err) %>% c(names(d_byspp[i]),.) %>% rbind(d_SEs,.) -> d_SEs
}
d_means = d_means[-1,]
d_SEs = d_SEs[-1,]
colnames(d_means)[1] = "Species"
colnames(d_SEs)[1] = "Species"
rownames(d_means) = d_means[,1]
rownames(d_SEs) = d_SEs[,1]
d_means = as.data.frame(d_means[,-1])
d_SEs = as.data.frame(d_SEs[,-1])
for(i in 1:ncol(d_means)){
  d_means[,i] <- d_means[,i] %>% as.character() %>% as.numeric()
  d_SEs[,i] <- d_SEs[,i] %>% as.character() %>% as.numeric()
}
  #Extracting the differences between minimal and maximal Erythromycin concentrations
d_change = keyvars
d_changeSEs = keyvars
for(i in 1:length(d_byspp)){
  DSP <- d_byspp[[i]]
  dSPchange = c(names(d_byspp)[i],1:4)
  for(j in 2:5){
    dSPchange[j] <- round((mean(DSP[(nrow(DSP)-3):(nrow(DSP)),j])-mean(DSP[1:3,j]))/mean(DSP[,j]),3)
  }
  d_change <- rbind(d_change, dSPchange)
}
d_change <- d_change[-1,]
colnames(d_change) = str_replace_all(keyvars, "d", "m")
d_change = d_change[,-1]
d_change = as.data.frame(apply(d_change,2, as.numeric))
rownames(d_change) = sort(unique(d_data$Species))
endosym = c(3,2,2,1,2,2,2,1)
d_all <- cbind(controls,d_means, d_change, endosym)

###Load phylogenies
Nphylo <- read.tree("Fasta/Concatenado_nuclear/extra_species_conca_nuclear")
#plot(Nphylo); nodelabels(frame = "none", font = 0.5); tiplabels(); edgelabels(frame = "none", font = 0.2, offset = 3)
Nphylo <- drop.tip(Nphylo, c(15,13,8,7,3,2,1))
Nphylo = reroot(Nphylo, 8)
Nphylo = drop.tip(Nphylo, 8)
Nphylo$tip.label <- c("Isochrysis galbana", "Dunaliella salina", "Tetraselmis suecica", "Nannochloropsis gaditana", "Cylindrotheca closterium", "Phaeodactylum tricornutum", "Chaetoceros gracilis", "Amphidinium carterae")
Nphylo$edge.length[c(5,11)] = mean(Nphylo$edge.length)/1000
Nultra = chronos(Nphylo)
d_all_N = d_all[match(Nultra$tip.label, rownames(d_all)),]
SE_N = d_SEs[match(Nultra$tip.label, rownames(d_SEs)),]

PSBphylo <- read.tree("Fasta/psbA/psbA_tree_grupo externo")
#plot(PSBphylo); nodelabels(); tiplabels()
PSBphylo = reroot(PSBphylo, 9)
PSBphylo <- drop.tip(PSBphylo, 9)
PSBphylo$tip.label = c("Isochrysis galbana", "Dunaliella salina", "Tetraselmis suecica", "Nannochloropsis gaditana", "Cylindrotheca closterium", "Phaeodactylum tricornutum", "Chaetoceros gracilis", "Amphidinium carterae")[c(7,4,6,5,3,2,1,8)]
PSBphylo$edge.length[c(3,4,6)] = mean(PSBphylo$edge.length)/1000
PSBultra = chronos(PSBphylo)
d_all_PSB = d_all[match(PSBultra$tip.label, rownames(d_all)),]
SE_PSB = d_SEs[match(PSBultra$tip.label, rownames(d_SEs)),]

Ch16Sphylo <- read.tree("Fasta/16S_cloroplastos/16S_tree_grupo_externo")
#plot(Ch16Sphylo); nodelabels(); tiplabels()
Ch16Sphylo <- reroot(Ch16Sphylo, 8)
Ch16Sphylo <- drop.tip(Ch16Sphylo, 8)
Ch16Sphylo$tip.label = c("Isochrysis galbana", "Dunaliella salina", "Tetraselmis suecica", "Nannochloropsis gaditana", "Cylindrotheca closterium", "Phaeodactylum tricornutum", "Chaetoceros gracilis", "Amphidinium carterae")[c(6,5,7,4,1,2,3,8)]
Ch16Sphylo$edge.length[c(4,9)] = mean(Ch16Sphylo$edge.length)/1000
Ch16Sultra = chronos(Ch16Sphylo)
d_all_16S = d_all[match(Ch16Sultra$tip.label, rownames(d_all)),]
SE_16S = d_SEs[match(Ch16Sultra$tip.label, rownames(d_SEs)),]

Ch23Sphylo <- read.tree("Fasta/23s/23s_tree")
#plot(Ch23Sphylo); nodelabels(); tiplabels()
Ch23Sphylo <- reroot(Ch23Sphylo, 7)
Ch23Sphylo <- drop.tip(Ch23Sphylo, 7)
Ch23Sphylo$tip.label = c("Isochrysis galbana", "Dunaliella salina", "Tetraselmis suecica", "Nannochloropsis gaditana", "Cylindrotheca closterium", "Phaeodactylum tricornutum", "Chaetoceros gracilis", "Amphidinium carterae")[c(1,4,3,6,7,5,8,2)]
Ch23Sphylo$edge.length[c(12,4,5,6)] = mean(Ch23Sphylo$edge.length)/1000
Ch23Sultra = chronos(Ch23Sphylo)
d_all_23S = d_all[match(Ch23Sultra$tip.label, rownames(d_all)),]
SE_23S = d_SEs[match(Ch23Sultra$tip.label, rownames(d_SEs)),]

### Pack-up data and trees into iterable lists
treelist <- list(Nultra, PSBultra, Ch16Sultra, Ch23Sultra)
data_list = list(d_all_N, d_all_PSB, d_all_16S, d_all_23S)
SElist = list(SE_N, SE_PSB, SE_16S, SE_23S)
names(data_list) = c("Nuclear", "PSBA", "16S", "23S")
names(SElist) = c("Nuclear", "PSBA", "16S", "23S")

### Estimate RF distance between trees
RFdiffs = matrix(nrow = 4, ncol = 4)
colnames(RFdiffs) = names(data_list)
rownames(RFdiffs) = names(data_list)
BSdiffs = RFdiffs
for(i in 1:4){
  tree_i <- treelist[[i]]
  for(j in 1:4){
    tree_j <- treelist[[j]]
    RFdiffs[i,j] = treedist(tree_i, tree_j)[1]
    BSdiffs[i,j] = treedist(tree_i, tree_j)[2]
  }
}
#write.csv(RFdiffs, "deliverables/Tree_difference/RF_tree_distances.csv")
#write.csv(BSdiffs, "deliverables/Tree_difference/BranchScore_differences.csv")

###Model testing and Phylogenetic Signals
AICc_list = list()
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  startree <- rescale(i_tree, "lambda", 0)
  class(i_tree) = "phylo"
  i_dat <- data_list[[tree]]
  i_SE <- data.frame(rep(0, nrow(i_dat)), rep(0, nrow(i_dat)),rep(0, nrow(i_dat)),rep(0, nrow(i_dat)), SElist[[tree]],SElist[[tree]], rep(0, nrow(i_dat)))
  names(i_SE)=names(i_dat)
  print(names(data_list)[tree])
  AIC_treei = matrix(nrow=ncol(i_dat), ncol = 7)
  colnames(AIC_treei) = c("Variable", "Tree", "Best_model", "2nd_Best_model", "1_2LoglikRatio_Pvalue", "K", "p_K")
  for(c in 1:ncol(i_dat)){
    C = i_dat[,c]
    names(C) = rownames(i_dat)
    Cse = i_SE[,c]
    names(Cse) = rownames(i_SE)
    PSS_c <- phylosig(i_tree, C, se=Cse, test=T)
    PSS_c <- c(PSS_c$K, PSS_c$P)
    model_matrix = matrix("NA", nrow = 10, ncol = 2)
    colnames(model_matrix) = c("aicc","lnL")
    row.names(model_matrix) = c("BM", "white", "drift", "EB", "OU", "trend", "delta", "lambda", "kappa", "starBM")
    for(j in 1:dim(model_matrix)[1]){
      if(j==nrow(model_matrix)){
        temp_model <- fitContinuous(startree, C, model="BM", SE = Cse)$opt
      }
      else{
        temp_model <- fitContinuous(i_tree, C, model=row.names(model_matrix)[j], SE = Cse)$opt
      }
      model_matrix[j, "aicc"] <- temp_model$aicc
      model_matrix[j, "lnL"] <- temp_model$lnL
    }
    model_matrix %>% as.data.frame() -> model_df
    model_df$aicc <- as.numeric(as.character(model_df$aicc))
    model_df$lnL <- as.numeric(as.character(model_df$lnL))
    model_df <- model_df[order(model_df$aicc),]
    print(model_df)
    Pchi <- (2*(model_df$lnL[1] - model_df$lnL[2])) %>% pchisq(df=1, lower.tail = F)
    print(Pchi)
    string_c <- c(names(i_dat)[c], names(data_list)[tree], rownames(model_df)[1], rownames(model_df)[2], Pchi, PSS_c[1:2])
    names(string_c) = colnames(AIC_treei)
    AIC_treei[c,] <- string_c
  }
  AICc_list[[tree]] <- as.data.frame(AIC_treei)
  print(AIC_treei)
}
fullModelTesting <- rbind(AICc_list[[1]], AICc_list[[2]], AICc_list[[3]], AICc_list[[4]])
#write.csv(fullModelTesting, "deliverables/model_support/BestModels_wSE_wLL.csv")

### Arbutus model adequacy

supportBM <- read.csv("deliverables/model_support/BestModels_wSE_wLL.csv", row.names = 1) %>% .[which(.$Best_model == "BM"),]
names(treelist) = unique(supportBM$Tree)
madlist=list()
for(tree in 1:length(treelist)){
  i_tree <- treelist[[tree]]
  class(i_tree) = "phylo"
  i_dat <- data_list[[tree]] %>% .[,which(names(.) %in% supportBM$Variable[supportBM$Tree == names(treelist)[tree]])]
  MADtable = as.data.frame(matrix(ncol=8, nrow=nrow(supportBM[which(supportBM$Tree==names(treelist)[tree]),])))
  names(MADtable) = c("Variable", "Tree", "msig", "cvar", "svar", "sasr", "shgt", "dcfd")
  for(c in 1:ncol(i_dat)){
    MADtable$Variable[c] <- names(i_dat)[c]
    MADtable$Tree[c] <- names(treelist)[tree]
    C <- i_dat[,c]
    names(C)=rownames(i_dat)
    fitC <- fitContinuous(i_tree, C, model="BM")
    UTC <- make_unit_tree(fitC)
    picstat_data <- calculate_pic_stat(UTC)
    sim <- simulate_char_unit(UTC)
    picstat_sim <- calculate_pic_stat(sim)
    compare_pic_stat(picstat_data, picstat_sim) %>% .$p.values -> MADtable[,c(3:8)]
  }
  madlist[[tree]] <- MADtable
}
madBMtable <- rbind(madlist[[1]], madlist[[2]], madlist[[3]], madlist[[4]])

  #dGrowth in 23S, EB
dGr23S = d_all_23S$dGrowth
names(dGr23S) = rownames(d_all_23S)

fC_EB_dGr23S = fitContinuous(Ch23Sultra, dGr23S, model="EB")
UTEB_dGr23S = make_unit_tree(fC_EB_dGr23S)
picstatEB_dGr23S = calculate_pic_stat(UTEB_dGr23S)
simcharEB_dGr23S = simulate_char_unit(UTEB_dGr23S)
simpicstatEB_dGr23S = calculate_pic_stat(simcharEB_dGr23S)
compare_pic_stat(picstatEB_dGr23S, simpicstatEB_dGr23S) %T>% plot() %>% .$p.values %>% as.data.frame() -> MAD_EB_dGr23S

  #mEQY in 23S, EB
mEQY23S = d_all_23S$mEQY
names(mEQY23S) = rownames(d_all_23S)

fC_EB_mEQY23S = fitContinuous(Ch23Sultra, mEQY23S, model="EB")
UTEB_mEQY23S = make_unit_tree(fC_EB_mEQY23S)
picstatEB_mEQY23S = calculate_pic_stat(UTEB_mEQY23S)
simcharEB_mEQY23S = simulate_char_unit(UTEB_mEQY23S)
simpicstatEB_mEQY23S = calculate_pic_stat(simcharEB_mEQY23S)
comparison = compare_pic_stat(picstatEB_mEQY23S, simpicstatEB_mEQY23S)
compare_pic_stat(picstatEB_mEQY23S, simpicstatEB_mEQY23S) %T>% plot() %>% .$p.values %>% as.data.frame() -> MAD_EB_mEQY23S
MADs = data.frame(MAD_EB_dGr23S, MAD_EB_mEQY23S)
names(MADs) = c("dGrowth", "mEQY")


###Contmaps
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
  evol_vars <- which(AICc_list[[tree]]$Best_model != "WN")
  i_dat = i_dat[,evol_vars]
  print(names(data_list)[tree])
for(i in 1:ncol(i_dat)){
  coloring = as.numeric(i_dat[,i])
  names(coloring) = rownames(i_dat)
  coloring = coloring[!is.na(coloring)]
  obj_i <- contMap(i_tree,coloring, plot=F)
  obj_i$cols[1:length(obj_i$cols)] = colorRampPalette(c("blue","red"), space="Lab")(length(obj_i$cols))
  pdf(paste("deliverables/contMaps/", names(data_list)[tree], colnames(i_dat)[i],".pdf",sep=""),width=6,height=10)
  plot(obj_i, fsize=0.5)
  title(names(i_dat)[i], cex.main = 0.5)
  dev.off()
}}

###Phylomorphospace scattergrams with contMaps in the diagonal
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
  evol_vars <- which(AICc_list[[tree]]$Best_model != "WN")
  i_dat = i_dat[,evol_vars]
  print(names(data_list)[tree])
    pdf(paste("deliverables/phylomorphograms/", names(data_list)[tree],".pdf",sep=""),width=24,height=16)
    fancyTree(i_tree, type="scattergram",X=as.matrix(i_dat),control=list(spin=FALSE, fsize=0.1), label = 'horizontal')
    dev.off()
}

###Phylomorphospaces in focus for 23S
vars = c("mGrowth", "mETR", "dChlA")
coloring = as.numeric(d_all_23S[,vars[3]])
names(coloring) = rownames(d_all_23S)
obj_i <- contMap(i_tree,coloring, plot=F)
obj_i$cols[1:length(obj_i$cols)] = colorRampPalette(c("blue","red"), space="Lab")(length(obj_i$cols))
pdf(paste("deliverables/chromophylomorphospaces/23S",vars[1], vars[2], vars[3], ".pdf",sep=""),width=16,height=12)
phylomorphospace(obj_i$tree, d_all_23S[,vars[1:2]], colors=obj_i$cols, label="horizontal", fsize=0.8)
title(main="X: mGrowth . Y: mETR. Color: dChla")
dev.off()


###Phylogenetic Generalized Least-Squares Regression Models
PGLSp_list = list()
logLs_list = list()
PGLStrunc_list = list()
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
  print(names(data_list)[tree])
  #evol_vars <- which(AICc_list[[tree]]$Best_model != "WN")
  #i_dat = i_dat[,evol_vars]
  PGLS_pvalues = matrix(nrow=ncol(i_dat), ncol=ncol(i_dat))
  colnames(PGLS_pvalues)=rownames(PGLS_pvalues)=colnames(i_dat)
  logLs = PGLS_pvalues
for(i in 1:ncol(i_dat)){
  CH_I=as.numeric(i_dat[,i])
  names(CH_I) = rownames(i_dat)
  CH_I = CH_I[!is.na(CH_I)]
  for(j in 1:ncol(i_dat)){
    CH_J=as.numeric(i_dat[,j])
    names(CH_J) = rownames(i_dat)
    if(!(i==j)){
      if(!identical(CH_I,CH_J)){
        print(colnames(i_dat)[c(i,j)])
        gls(CH_I ~ CH_J, correlation =  corBrownian(phy = i_tree), data = as.data.frame(cbind(CH_I,CH_J)), method = "ML") %>% summary() %>% .$tTable %>% as.data.frame() %>% .$p %>% .[2] -> PGLS_pvalues[i,j]
        pgls.Ives(i_tree, rep(CH_I,2), rep(CH_J,2)) %>% .$logL -> logLs[i,j]
      }
    }
  }}
  PGLS_pvaluesTRUNC = PGLS_pvalues
  PGLS_pvaluesTRUNC[PGLS_pvaluesTRUNC>0.05] = NA
  PGLSp_list[[tree]]<-PGLS_pvalues
  logLs_list[[tree]]<-logLs
  PGLStrunc_list[[tree]]<-PGLS_pvaluesTRUNC
}
for(tree in 1:4){
  write.csv(PGLSp_list[[tree]], paste("deliverables/PGLS/PGLSpvalues_", names(data_list)[tree], ".csv", sep=""))
  write.csv(logLs_list[[tree]], paste("deliverables/PGLS/PGLSlogLikelihoods_", names(data_list)[tree], ".csv", sep=""))
  write.csv(PGLStrunc_list[[tree]], paste("deliverables/PGLS/PGLSpvaluesTRUNC_", names(data_list)[tree], ".csv", sep=""))
}

DG23S = d_all_23S$dGrowth
names(DG23S) = rownames(d_all_23S)
MEQY23S = d_all_23S$mEQY
names(MEQY23S) = rownames(d_all_23S)
DG23S_SE = SE_23S$dGrowth
names(DG23S_SE) = rownames(d_all_23S)
MEQY23S_SE = SE_23S$dEQY
names(MEQY23S_SE) = rownames(d_all_23S)
pglsDATA = data.frame(DG23S, MEQY23S)
rownames(pglsDATA) = rownames(d_all_23S)
star23S = rescale(Ch23Sultra,"lambda",0)
Me_Dg <- pgls.SEy(model = MEQY23S~DG23S, data = pglsDATA, se=MEQY23S_SE, tree = Ch23Sultra)
Me_Dg_star <- pgls.SEy(model = MEQY23S~DG23S, data = pglsDATA, se=MEQY23S_SE, tree = star23S)
pgls23tree = Ch23Sultra
class(pgls23tree) <- "phylo"
PIC_meqy_dgrowth = Rcontrast(pgls23tree, pglsDATA, path="~/Downloads/phylip-3.695/exe")
PIC_meqy_dgrowth$Contrasts %>% cor.table()
cor.table(pglsDATA)

###Phenograms with uncertainty
Phenogram <- function(a,b,c){
  contTrait <- a
  names(contTrait) = rownames(b)
  contTrait = contTrait[!is.na(contTrait)]
  treeI = drop.tip(c,which(!(c$tip.label %in% names(contTrait))))
  A<-fastAnc(treeI,contTrait,CI=TRUE)
  paintree<-paintSubTree(treeI,node=length(c$tip)+1,"1")
  trans<-as.character(floor(0:50/2))
  trans[as.numeric(trans)<10]<- paste("0", trans[as.numeric(trans)<10],sep="")
  for(i in 0:50){
    p<-i/length(trans)
    phenogram(treeI,c(contTrait,(1-p)*A$CI95[,1]+p*A$ace), colors=setNames(paste("#0000ff",trans[i+1],sep=""),1), add=i>0, ftype="off")
    phenogram(treeI,c(contTrait,(1-p)*A$CI95[,2]+p*A$ace), colors=setNames(paste("#0000ff",trans[i+1],sep=""),1), add=TRUE, ftype="off")
  }
  phenogram(treeI,c(contTrait,A$ace),add=TRUE, colors=setNames("black",1), ftype="off")
}

for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
  evol_vars <- which(AICc_list[[tree]]$Best_model != "WN")
  i_dat = i_dat[,evol_vars]
  print(names(data_list)[tree])
  for(i in 1:ncol(i_dat)){
  pdf(paste("deliverables/phenograms/", names(data_list)[tree], names(i_dat)[i], ".pdf",sep=""),width=12,height=12)
  Phenogram(i_dat[,i], i_dat, i_tree)
  dev.off()
  }
}

### Uncertainty estimates
uncertainty = c("Variable", "Tree", "Mean_ACE", "Mean_CI95_range", "Rel_CI95_uncertainty")
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
  evol_vars <- which(AICc_list[[tree]]$Best_model != "WN")
  i_dat = i_dat[,evol_vars]
  for(c in 1:ncol(i_dat)){
    var_i <- i_dat[,c]
    names(var_i) = rownames(i_dat)
    ace_i <- fastAnc(i_tree, var_i, CI=TRUE)
    CI95r <- apply(ace_i$CI95, 1, function(x){abs(x[1]-x[2])})
    mCI95 <- mean(CI95r)
    print(mCI95)
    scaled_mCi95 <- (mCI95 - min(abs(c(var_i, CI95r))))/(max(var_i) - min(var_i))
    uncertainty <- rbind(uncertainty, c(names(i_dat)[c], names(data_list)[tree], mean(ace_i$ace), mCI95, scaled_mCi95))
  }
}
rownames(uncertainty) = 1:nrow(uncertainty)
colnames(uncertainty) = uncertainty[1,]
uncertainty = uncertainty[-1,]
uncertainty <- as.data.frame(uncertainty)
uncertainty <- sapply(uncertainty, as.character)
uncertainty[,3:5] <- sapply(uncertainty[,3:5], as.numeric)
#write.csv(uncertainty, "deliverables/model_support/uncertainty.csv")

###PCA ecotoxospace
PCPS_list = list()
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
PCA(i_dat[,-c(1:4,13)]) -> Pca
coloring = as.numeric(Pca$ind$coord[,3])
names(coloring) = rownames(Pca$ind$coord)
obj_i <- contMap(i_tree,coloring, plot=F)
obj_i$cols[1:length(obj_i$cols)] = colorRampPalette(c("yellow","purple"), space="Lab")(length(obj_i$cols))
pdf(paste("deliverables/chromophylomorphospaces/PCA_", names(data_list)[tree], ".pdf",sep=""),width=14,height=10)
phylomorphospace(obj_i$tree, Pca$ind$coord[,1:2], colors=obj_i$cols, label="horizontal", fsize=0.8, xlab = "PC1 (45%)", ylab = "PC2 (24%)")
title(paste("Color = PC3 (23%) --- Tree: ", names(data_list)[tree], sep = ""))
dev.off()
  PCdat_i <- Pca$ind$coord[,1:4]
  print(names(data_list)[tree])
  physignal(PCdat_i, i_tree) -> PCPSi 
  c(PCPSi$phy.signal, PCPSi$pvalue) -> PCPS_list[[tree]]
}
names(PCPS_list) = names(data_list)
as.data.frame(rbind(PCPS_list[[1]], PCPS_list[[2]], PCPS_list[[3]], PCPS_list[[4]])) -> PCPS_table
names(PCPS_table) = c("multivariate_K", "p_value")
rownames(PCPS_table) = names(data_list)
write.csv(PCPS_table, "deliverables/PCA/physignal.csv")

PCA(d_all[,-c(1:4,13)]) -> Pca
Pca %>% fviz_contrib(choice="var", axes=1, sort.val="desc")
Pca %>% fviz_contrib(choice="var", axes=2, sort.val="desc")
Pca %>% fviz_contrib(choice="var", axes=3, sort.val="desc")
Pca %>% fviz_pca_biplot( col.var="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
Pca %>% fviz_pca_biplot(axes=c(4,3), col.var="contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)

  #PhylPCA
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
  ppca_i <- phyl.pca(i_tree, i_dat[,-c(1:4,13)])
  ppca_i %>% biplot(main=paste(signif(summary(ppca_i)$importance[2,1]*100,3),"%"), ylab=paste(signif(summary(ppca_i)$importance[2,2]*100,3),"%"), cex = .6, expand =100)
  ppca_i$L %>% as.data.frame() %>% ggplot() + geom_point(mapping = aes(x=rownames(ppca_i$L), y=-PC1)) + theme(axis.text.x = element_text(angle = 50, hjust = 1))
  ppca_i$L %>% as.data.frame() %>% ggplot() + geom_point(mapping = aes(x=rownames(ppca_i$L), y=-PC2)) + theme(axis.text.x = element_text(angle = 50, hjust = 1))
}

###Best Model parameters
modelpars = c("Variable", "Tree", "model", "sigsq", "z0", "alpha", "slope")
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]]
  i_SE <- data.frame(rep(0, nrow(i_dat)), rep(0, nrow(i_dat)),rep(0, nrow(i_dat)),rep(0, nrow(i_dat)), SElist[[tree]],SElist[[tree]], rep(0, nrow(i_dat)))
  names(i_SE)=names(i_dat)
  evol_vars <- which(AICc_list[[tree]]$Best_model != "WN")
  i_dat = i_dat[,evol_vars]
  for(c in 1:ncol(i_dat)){
    model_c = fullModelTesting$Best_model[which(fullModelTesting$Tree == names(data_list)[tree] & fullModelTesting$Variable == names(i_dat)[c])]
    if(model_c != "WN" & names(i_dat)[c] != "endosym"){
      C = i_dat[,c]
      names(C) = rownames(i_dat)
      Cse = i_SE[,c]
      names(Cse) = rownames(i_SE)
      fit_c <- fitContinuous(i_tree, C, model=as.character(model_c), SE = Cse)
      if(as.character(model_c)=="trend"){
        SLP = fit_c$opt$slope
      }
      else{SLP = 0}
      if(as.character(model_c)=="EB"){
        alpha = fit_c$opt$a
      }
      else{alpha = 0}
      modelpars <- rbind(modelpars, c(names(i_dat)[c], names(data_list)[tree], as.character(model_c), fit_c$opt$sigsq, fit_c$opt$z0, alpha, SLP))
    }
  }
}
rownames(modelpars) = 1:nrow(modelpars)
colnames(modelpars) = modelpars[1,]
modelpars = modelpars[-1,]
modelpars = as.data.frame(modelpars)
modelpars[,1:3] <- sapply(modelpars[,1:3], as.character)
modelpars[,4:7] <- sapply(modelpars[,4:7], as.character)
modelpars[,4:7] <- sapply(modelpars[,4:7], as.numeric)
#modelpars = modelpars[,-6]     #to remove alpha if no OU
write.csv(modelpars, "deliverables/model_support/model_parameters_wSE.csv")

##Phylogenetic Path Analysis
for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]] %>% .[,c(5:8,13)]
  m <- define_model_set(
    null = c(), 
    electon_mediation = c(dGrowth~dETR, dETR~dEQY, dEQY~dChlA), 
    .common = c()
    )
  PP <- phylo_path(m, i_dat, i_tree)
  #summary(PP) %>% .$CICc %>% min() %>% print()
  print(summary(PP))
  #plot(best(PP))
  #coef_plot(best(PP), error_bar = "se", order_by = "strength", to = "a") + ggplot2::coord_flip()
  #avgPP = average(PP, avg_method = "full")
  #coef_plot(avgPP, error_bar = "se", order_by = "strength", to = "a") + ggplot2::coord_flip()
}

for(tree in 1:4){
  i_tree <- treelist[[tree]]
  i_dat <- data_list[[tree]] %>% .[,9:13]
  m <- define_model_set(
    null = c(), 
    electon_mediation = c(mGrowth~mETR, mETR~mEQY, mEQY~mChlA)) #, 
    #.common = c()
  #)
  PP <- phylo_path(m, i_dat, i_tree)
  #summary(PP) %>% .$CICc %>% min() %>% print()
  print(summary(PP))
  #plot(best(PP))
  #coef_plot(best(PP), error_bar = "se", order_by = "strength", to = "a") + ggplot2::coord_flip()
  #avgPP = average(PP, avg_method = "full")
  #coef_plot(avgPP, error_bar = "se", order_by = "strength", to = "a") + ggplot2::coord_flip()
}

