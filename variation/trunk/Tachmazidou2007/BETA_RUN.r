############################################################################

### root is the full path to the dataset, for example
root <- "/home/it103/data/"

### outRoot is the full path to the folder that contains treepic and seq2tr
### Save the perl files getUniqueM.pl and mutOrder.pl in the same folder 
### For example
outRoot <- "/home/it103/gtree90/"

### Default (FIXED=0) is to run the algorithm without fixing the number of clusters
### If you want to fix them then FIXED should be assinged to 1.
FIXED<-0

### Default (PENALIZE=0) is to run the algorithm assuming a uniform prior on the trees
### If you want to penalize large trees then PENALIZE should be assinged to 1.
PENALIZE<-0

############################################################################

### Read in the haplotypes (coded as 0,1), 
### the case/contol status (1/0 respectively) 
### and the physical positions of the SNPs

haplotypes <- read.table(paste(root,"haplotypes.txt",sep=""))
haplotypes <- as.matrix(haplotypes)
dim(haplotypes)
positions <- read.table(paste(root,"positions.txt",sep=""))
positions <- as.numeric(unlist(positions))
length(positions)
status <- read.table(paste(root,"status.txt",sep=""))
status <- as.numeric(unlist(status))
length(status)

############################################################################

### Divide the dataset into phylogenetic trees
myfindTrees <- findTrees(haplotypes)
n.genetrees <- myfindTrees$ntrees ### this is the number of trees
boundaries <- myfindTrees$boundaries

### Find the number of SNPs in each tree
treelength <- numeric(ncol(boundaries))
for(i in 1:ncol(boundaries)){
  treelength[i] <- boundaries[2,i] - boundaries[1,i] + 1
}

### Gives the first and the last SNP in each tree
Ntree <- matrix(0,ncol=2,nrow=n.genetrees)
for(i in 1:n.genetrees) Ntree[i,] <- boundaries[,i]

### Gives the haplotypes and SNP positions in each tree
subhaplo <- vector(nrow(Ntree),mode="list")
subpos <- vector(nrow(Ntree),mode="list")
for(i in 1:nrow(Ntree)){
  subhaplo[[i]] <- haplotypes[,Ntree[i,1]:Ntree[i,2]]
  subpos[[i]] <- positions[Ntree[i,1]:Ntree[i,2]]
}

### Gives the unique haplotypes in the whole dataset
mycat <- categorical(haplotypes)
mycat.haplo <- mycat$haplo ### the unique haplotypes
mycat.k <- mycat$k ### the number of unique haplotypes
mycat.x <- mycat$x ### identify
mult <- as.vector(table(mycat.x))

### Gives the unique haplotypes in each tree
mysubcat.haplo <- vector(nrow(Ntree), mode="list")### the unique haplotypes
mysubcat.k <- vector(nrow(Ntree), mode="list")### the number of unique haplotypes
mysubcat.x <- vector(nrow(Ntree), mode="list")### identify
for(i in 1:nrow(Ntree)){
  mysubcat <- categorical(subhaplo[[i]])
  mysubcat.haplo[[i]] <- mysubcat$haplo
  mysubcat.k[[i]] <- mysubcat$k
  mysubcat.x[[i]] <- mysubcat$x
}
submult <- vector(nrow(Ntree), mode="list")
for(i in 1:nrow(Ntree)) submult[[i]] <- as.vector(table(mysubcat.x[[i]]))

### Count the multiplicity of each unique haplotype, and how many of these have a case/control status 
count <- HaploCount(mycat.k, mycat.x, status)

### For each tree count the multiplicity of each unique haplotype, and how many of these have a case/control status 
subcountN <- vector(nrow(Ntree), mode="list")
subcountNDIS <- vector(nrow(Ntree), mode="list")
subcountRATIO <- vector(nrow(Ntree), mode="list")
for(i in 1:nrow(Ntree)){
  subcount <- HaploCount(mysubcat.k[[i]], mysubcat.x[[i]], status)
  subcountN[[i]] <- subcount$n.haplo
  subcountNDIS[[i]] <- subcount$n.diseased
  subcountRATIO[[i]] <- subcount$n.diseased/subcount$n.haplo
}

####################################################################################################################

fileNamesR <- NULL
for(tr in 1:nrow(Ntree)){
  fileNamesR[tr] <- paste(outRoot,sprintf("Rtree%d",tr),sep="")
}

fileNamesMULT <- NULL
for(tr in 1:nrow(Ntree)){
  fileNamesMULT[tr] <- paste(outRoot,sprintf("tree%d",tr),sep="")
}

fileNamesTR <- NULL
for(tr in 1:nrow(Ntree)){
  fileNamesTR[tr] <- paste(outRoot,sprintf("tree%d.tr",tr),sep="")
}

fileNamesPS <- NULL
for(tr in 1:nrow(Ntree)){
  fileNamesPS[tr] <- paste(outRoot,sprintf("tree%d.ps",tr),sep="")
}

if (getRPLM(outRoot, Ntree, mysubcat.haplo, submult, fileNamesR)) {
  if (getSEQ2TR(outRoot, Ntree, fileNamesMULT)) {
    if (getTREEPIC(outRoot, Ntree, fileNamesTR)) {
      myage.tree <- getMutOrder(outRoot, Ntree, fileNamesPS)
      if (myage.tree$run) {
        age.tree <- myage.tree$age.tree
      }
    }
  }
}

### Get the relative order of mutations in each tree
myAGE <- constructAGE(treelength, age.tree)
removeFiles(Ntree, fileNamesR, fileNamesMULT, fileNamesTR)

### Do not run this if you want to keep the files containing the trees (they are in the OutRoot folder)
for(tr in 1:nrow(Ntree)) unlink(fileNamesPS[tr])

#################################################################################

alpha <- 1
beta <- 1
MCMC <- 100000

if(FIXED==0){ ### NOT FIXED

indicator <- function(X){

  i <- sample(1:length(X), 1)
  pick <- X[i]

  if(pick==0){
    move <- 1
    pick<-1
  }else{
    move <- -1
    pick<-0
  }

  X[i] <- pick
  nones <- sum(X)
  if(nones>0){
    pos.ones <- which(X>0)
  }else pos.ones<- -1

  out <- list(X=X, n.centers=nones, center.pos=pos.ones, move=move)
  return(out)
  
}

select_tree<- function(Ntree){

  nt <- nrow(Ntree)
  new <- sample(1:nt,1)
  ns <- Ntree[new,2] - Ntree[new,1] + 1
  
  ind <- rbinom(ns,1,0.5)
    
  nones <- sum(ind)
  if(nones>0){
    pos.ones <- which(ind>0)
  }
  else pos.ones<- -1

  out <- list(TREE=new,X=ind, n.centers=nones, center.pos=pos.ones)
  return(out)

}

}else if(FIXED==1){ ### FIXED

indicator <- function(X){

  X <- rep(0,length(X))
  i <- sample(1:length(X), 1)
  pick <- X[i]

  if(pick==0){
    move <- 1
    pick<-1
  }else{
    move <- -1
    pick<-0
  }
  
  X[i] <- pick
  nones <- sum(X)
  if(nones>0){
    pos.ones <- which(X>0)
  }else pos.ones<- -1

  out <- list(X=X, n.centers=nones, center.pos=pos.ones, move=move)
  return(out)
  
}

select_tree<- function(Ntree){

  nt <- nrow(Ntree)
  new <- sample(1:nt,1)
  ns <- Ntree[new,2] - Ntree[new,1] + 1
  
  X <- rep(0,ns)
  s <- sample(1:ns, 1)
  X[s] <- 1
  ind <- X
    
  nones <- sum(ind)
  if(nones>0){
    pos.ones <- which(ind>0)
  }
  else pos.ones<- -1

  out <- list(TREE=new,X=ind, n.centers=nones, center.pos=pos.ones)
  return(out)

}

}else{
  cat("\nError, FIXED can be 0 or 1.\n")
}

#################################################################################

### Run the algorithm
myntrees <- multTreesPartition(alpha, beta, status, haplotypes, myAGE, Ntree, mycat.haplo, MCMC, PENALIZE)

### Posterior probability of each tree carrying the causal SNP
post.prob.tree <- numeric(nrow(Ntree))
for(i in 1:nrow(Ntree)){
  st <- myntrees$TREERES == i
  post.prob.tree[i] <- sum(st)/length(myntrees$TREERES)
}
cat(paste("\nPosterior probability of each tree carrying the causal SNP\n"))
print(post.prob.tree)
indtree <- which(post.prob.tree==max(post.prob.tree))

### Prior probability of each tree carrying the causal SNP
prior.prob.tree <- 1/n.genetrees
### Bayes factor in favour of association for each tree
TreeBayesFactor <- (post.prob.tree/(1-post.prob.tree))/(prior.prob.tree/(1-prior.prob.tree))

### Number of clusters including the null cluster 
myNCLUST <- myntrees$NCLUST
x <- cbind(as.numeric(row.names(table(myNCLUST))), as.numeric(table(myNCLUST)))
nclust.final <- x[x[,2]==max(x[,2]),1]
cat(paste("\nThe mode of the distribution of the number of clusters (including the null cluster) is",nclust.final,"\n"))

### Marginal posterior probability of each SNP being a cluster centre 
myINDICATORS <- myntrees$INDICATORS
margProbSNPbeingCLUSTCentre <- apply(myINDICATORS,2,mean)
z <- cbind(SNP=1:ncol(myINDICATORS), marg.prob=apply(myINDICATORS,2,mean))
z <- z[z[,2]>=0.005,]
if(length(z)!=2){
  w <- z[order(z[,2],z[,1],decreasing=TRUE),]
}else{
  w <- z
}
S <- order(apply(myINDICATORS,2,sum), decreasing = TRUE)
cat(paste("\nMarginal posterior probability (>=0.005) of a SNP being a cluster centre (in descending order)\n"))
print(w)
clust.centres <- order(apply(myINDICATORS,2,sum), decreasing = TRUE)[1:(nclust.final-1)]
cat(paste("\nMost likely cluster centres (marginally) in descending order excluding the null cluster\n"))
print(clust.centres)

### SNP with the highest marginal posterior probability of being cluster centre, its position and the tree it belongs to 
SNP.MSasCO <- order(apply(myINDICATORS,2,sum), decreasing = TRUE)[1]
st <- numeric(nrow(Ntree))
for(i in 1:nrow(Ntree)) st[i] <- SNP.MSasCO%in%(Ntree[i,1]:Ntree[i,2])
TR.MSasCO <- which(st==1)
MSasCO <- c(SNP=SNP.MSasCO, position=positions[SNP.MSasCO], tree=TR.MSasCO)
cat(paste("\nResults:\n"))
print(MSasCO)

if(length(indtree)>1) indtree <- TR.MSasCO
cat(paste("\nThe tree with the highest posterior probability of carrying the causal SNP is tree",indtree,"\n"))
### Construct the tree choosen in a file called "pickedtr.ps" (in the OutRoot folder)
plotTree(outRoot, indtree, mysubcat.haplo, submult)
### Corresponding mutation names in indtree
sub2all <- numeric(length(indtree))
for(i in 1:length(indtree)){
  x <- myAGE[indtree[i],]
  sub2all[i] <- list(rbind(age.tree[[indtree[i]]], x[x!=-1]))
}
cat(paste("\nCorresponding mutation names in indtree\n"))
print(sub2all)

### Bayes factor in favour of association for each SNP
my.prior <- prior(n.genetrees, treelength)
SNPBayesFactor <- as.numeric((margProbSNPbeingCLUSTCentre/(1-margProbSNPbeingCLUSTCentre))/(my.prior/(1-my.prior)))

### Optional
### Mean posterior probability of disease for each unique haplotype in the whole dataset
myPVECTOR <- myntrees$PVECTOR
meanpr <- apply(myPVECTOR, 2, mean)#for each unique haplotype in the whole dataset
stdevpr <- sqrt(apply(myPVECTOR, 2, var))
Lpr <- numeric(ncol(myPVECTOR))
Upr <- numeric(ncol(myPVECTOR))
for(i in 1:ncol(myPVECTOR)){
  Lpr[i] <- quantile(myPVECTOR[,i],prob=0.025)
  Upr[i] <- quantile(myPVECTOR[,i],prob=0.975)
}
### Mean posterior probability of disease for each unique haplotype in the identified tree
mysubresults <- subresults(indtree, meanpr, stdevpr, Lpr, Upr, mycat.haplo, Ntree)
cat(paste("\nMean posterior probability of disease in the identified tree\n"))
print(cbind(multiplicity=subcountN[[indtree]], diseased=subcountNDIS[[indtree]], mean=mysubresults$postprob, stdev=mysubresults$poststdev, unique.haplos=mysubresults$cathaplo))

#################################################################################
### Construct 95% credible interval

Pos <- numeric(nrow(myINDICATORS))
for(f in 1:nrow(myINDICATORS)){
  selSNP <- which(myINDICATORS[f,]==1)
  Pos[f] <- mean(positions[selSNP])
}
x <- is.na(Pos)
Pos[x] <- -1
l<- Pos>-1
POS <- Pos[l]

Lci <- quantile(POS,prob=0.025)
Uci <- quantile(POS,prob=0.975)
SIM.IN95CI.L <- Lci
SIM.IN95CI.U <- Uci

cat(paste("\nLci is ",Lci,"\n"))
cat(paste("\nUci is ",Uci,"\n"))
cat(paste("\nround(Uci-Lci,3) is ",round(Uci-Lci,3),"\n"))

#################################################################################
