############################################################################

### We have created two populations using the FREGENE software. 
### One is simulated with recombination hotspots, and the other with uniform recombination rate 
### (for more details on how these datasets were simulated please see our paper).
### root is the full path to the populations, for example
root <- "/home/it103/data/"

### Below we determine the simulation parameters. 
### The values of the parameters correspond to our 'default' simulation scenario (explained in the paper).
maf <- 0.01 ### minor allele frequency of the alleles to be selected
GRR <- 1.6 ### genetic relative risk
n.markers <- 1000 ### number of markers to be selected (without including the causal SNP)
N <- 1000 ### number of case/control genotypes
MAF <- 0.05 
MAF <- round(seq(MAF-0.005,MAF+0.005,length=50),4) ### this is the maf of the causal SNP to be selected

### Dummy variable that indicates which population to sample from.
### The default is to sample from the population with hotspots (RECOMB=1). Else use RECOMB <- 0.
RECOMB <- 1 ### for variable recombination rate

### Dummy variable that indicates if the genotypes are sampled according to an additive or dominant disease model.
### The default is to use an additive model (DIS=1). Else use DIS <- 0.
DIS <- 1 ### for additive disease model

############################################################################
### Extra functions needed to simulate data

select.causal.SNP <- function(H,MAF){

  pos <- -1
  pr <- as.numeric(apply(H,2,sum)/nrow(H))
  for(i in 1:length(MAF)){
    dif <- abs(pr-MAF[i])
    candidates <- which(dif==min(dif))
    pos <- c(pos,candidates)
  }
  pos <- pos[-1]
  pos <- sort(pos)
  pos <- as.numeric(row.names(table(pos)))
  if(length(pos) > 1) pos <- sample(pos,1,prob=rep(0.5,length(pos)))

  return(list(pos=pos,maf.causal=pr[pos]))
}

haplo.mut.dat <- function(data, positions, pos, maf){

  new.pos <- -1
  n <- nrow(data)
  freq <- apply(data,2,sum)/n
  index <- which(freq > maf)
  sample.dat <- data[,index]
  mut.pos.new <- positions[index]

  for(w in 1:ncol(sample.dat)){
    if(sum(sample.dat[,w]==data[,pos])==n){
      new.pos <- w
      haplo.dat <- sample.dat
    }
  }

  if(new.pos == -1){

    new.pos <- length(mut.pos.new[mut.pos.new<=positions[pos]])+1
    y <- c(mut.pos.new, positions[pos])
    mut.pos.new <- sort(y)
    temp <- matrix(0,nrow(sample.dat),ncol(sample.dat)+1)
    if(new.pos > 1){
      for(i in 1:(new.pos-1)) temp[,i] <- sample.dat[,i]
    }
    temp[,new.pos] <- data[,pos]
    if(new.pos < ncol(temp)){
      for(i in (new.pos+1):ncol(temp))temp[,i] <- sample.dat[,(i-1)]
    }
    haplo.dat <- temp
  }
  return(list(haplo.dat=haplo.dat, mut.pos.new=mut.pos.new, new.pos=new.pos))
}

diseasedGenoPPA <- function(position, N, H, maf.causal, GRR){

  r <- GRR
  K <- 0.01 ### disease prevalence 1%
  p <- maf.causal
  # For an additive disease model I calculate the penetrances
  # assuming Hardy-Weinberg equilibrium
  f0 <- K/(1 - 2*p + 2*r*p)
  f1 <- r*f0
  f2 <- 2*r*f0-f0

  ### Pick N case genotypes
  geno.cases <- matrix(0,N,ncol(H))
  geno.cases.index <- matrix(0,N,2)
  i <- 1
  while(i <= N){

    pick <- sample(1:nrow(H),2,replace=TRUE)
    geno <- H[pick[1],] + H[pick[2],]

    if(geno[position] == 0){
      geno.status <- rbinom(1,1,f0)
    }else if(geno[position] == 1){
      geno.status <- rbinom(1,1,f1)
    }else{
      geno.status <- rbinom(1,1,f2)
    }
    if(geno.status == 1){

      geno.cases[i,] <- as.numeric(geno)
      geno.cases.index[i,] <- pick
      i<-i+1

    }
  }

  ### Pick N control genotypes
  geno.controls <- matrix(0,N,ncol(H))
  geno.controls.index <- matrix(0,N,2)
  i <- 1
  while(i <= N){

    pick <- sample(1:nrow(H),2,replace=TRUE)
    geno <- H[pick[1],] + H[pick[2],]

    if(geno[position] == 0){
      geno.status <- rbinom(1,1,f0)
    }else if(geno[position] == 1){
      geno.status <- rbinom(1,1,f1)
    }else{
      geno.status <- rbinom(1,1,f2)
    }
    if(geno.status == 0){

      geno.controls[i,] <- as.numeric(geno)
      geno.controls.index[i,] <- pick
      i<-i+1

    }
  }

  return(list(f0=f0, f1=f1, f2=f2, geno.cases=geno.cases, geno.cases.index=geno.cases.index, geno.controls=geno.controls, geno.controls.index=geno.controls.index))

}

diseasedGenoPPD <- function(position, N, H, maf.causal, GRR){

  r <- GRR
  K <- 0.01 #disease prevalence 1%
  p <- maf.causal
  # For a dominant disease model I calculate the penetrances
  # assuming Hardy-Weinberg equilibrium
  f0 <- K/(1 - 2*p + 2*r*p + p^2 - r*(p^2))
  f1 <- r*f0
  f2 <- f1

  ### Pick N case genotypes
  geno.cases <- matrix(0,N,ncol(H))
  geno.cases.index <- matrix(0,N,2)
  i <- 1
  while(i <= N){

    pick <- sample(1:nrow(H),2,replace=TRUE)
    geno <- H[pick[1],] + H[pick[2],]

    if(geno[position] == 0){
      geno.status <- rbinom(1,1,f0)
    }else if(geno[position] == 1){
      geno.status <- rbinom(1,1,f1)
    }else{
      geno.status <- rbinom(1,1,f2)
    }

    if(geno.status == 1){

      geno.cases[i,] <- as.numeric(geno)
      geno.cases.index[i,] <- pick
      i<-i+1

    }
  }

  ### Pick N control genotypes
  geno.controls <- matrix(0,N,ncol(H))
  geno.controls.index <- matrix(0,N,2)
  i <- 1
  while(i <= N){

    pick <- sample(1:nrow(H),2,replace=TRUE)
    geno <- H[pick[1],] + H[pick[2],]

    if(geno[position] == 0){
      geno.status <- rbinom(1,1,f0)
    }else if(geno[position] == 1){
      geno.status <- rbinom(1,1,f1)
    }else{
      geno.status <- rbinom(1,1,f2)
    }

    if(geno.status == 0){

      geno.controls[i,] <- as.numeric(geno)
      geno.controls.index[i,] <- pick
      i<-i+1

    }
  }

  return(list(f0=f0, f1=f1, f2=f2, geno.cases=geno.cases, geno.cases.index=geno.cases.index, geno.controls=geno.controls, geno.controls.index=geno.controls.index))

}

diseasedHaplo<- function(N, index.cases, index.controls, H){

  # If a genotype is a case, then both its haplotypes are cases
  haplo.cases <- matrix(0, N*2, ncol(H))
  haplo.controls <- matrix(0, N*2, ncol(H))

  k<-0
  for(i in 1:N){
    haplo.cases[i+k,] <-  as.numeric(H[index.cases[i,][1],])
    haplo.cases[i+k+1,] <-  as.numeric(H[index.cases[i,][2],])
    haplo.controls[i+k,] <-  as.numeric(H[index.controls[i,][1],])
    haplo.controls[i+k+1,] <-  as.numeric(H[index.controls[i,][2],])
    k<-k+1
  }

  myhaplo.dat <- rbind(haplo.cases,haplo.controls)
  myhaplo.status <- matrix(c(rep(1,2*N), rep(0,2*N)),ncol=1)

  return(list(haplo.cases=haplo.cases, haplo.controls=haplo.controls, myhaplo.dat=myhaplo.dat, myhaplo.status=myhaplo.status))

}

##########################################################################################################
### FOR VARIABLE RECOMBINATION RATE

if(RECOMB==1){

fname9 <- paste(root,"haplotypes_sim3.txt", sep = "")
fname10 <- paste(root,"positions_sim3.txt", sep = "")
mut.pos <- as.numeric(as.vector(read.table(fname10, sep = "")))
haplo.seq <- as.matrix(read.table(fname9,sep=""))
colnames(haplo.seq) <- NULL

}

##########################################################################################################
### FOR UNIFORM RECOMBINATION RATE

if(RECOMB==0){

fname9 <- paste(root,"haplotypes_sim1.txt", sep = "")
fname10 <- paste(root,"positions_sim1.txt", sep = "")
mut.pos <- as.numeric(as.vector(read.table(fname10, sep = "")))
haplo.seq <- as.matrix(read.table(fname9,sep=""))
colnames(haplo.seq) <- NULL
haplo.seq <- haplo.seq[,-1]
haplo.seq <- haplo.seq[-1,]

}

############################################################################

pos.info <- select.causal.SNP(haplo.seq,MAF)
pos <- pos.info$pos
maf.causal <- pos.info$maf.causal

sample <- haplo.mut.dat(haplo.seq, mut.pos, pos, maf)
haplo.dat <- sample$haplo.dat
mut.dat <- sample$mut.pos.new
causal.pos <- sample$new.pos

if(DIS==0){
 diseased.geno <- diseasedGenoPPD(causal.pos, N, haplo.dat, maf.causal, GRR)
}
if(DIS==1){
 diseased.geno <- diseasedGenoPPA(causal.pos, N, haplo.dat, maf.causal, GRR)
}
geno.cases <- diseased.geno$geno.cases
geno.controls <- diseased.geno$geno.controls
geno.cases.index <- diseased.geno$geno.cases.index
geno.controls.index <- diseased.geno$geno.controls.index
f0 <- diseased.geno$f0
f1 <- diseased.geno$f1
f2 <- diseased.geno$f2

diseased.haplo <- diseasedHaplo(N, geno.cases.index, geno.controls.index, haplo.dat)
haplo.cases <- diseased.haplo$haplo.cases
haplo.controls <- diseased.haplo$haplo.controls
myhaplo.dat <- diseased.haplo$myhaplo.dat
myhaplo.status <- diseased.haplo$myhaplo.status

indiv <- 4*N
mn <- apply(myhaplo.dat,2,sum)/indiv
mkr <- sort(sample(1:ncol(myhaplo.dat),n.markers+1,prob=mn*(1-mn)))
if(!(causal.pos%in%mkr)){
  del <- sample(1:length(mkr),1)
  mkr <- mkr[-del]
  mkr <- sort(c(mkr,causal.pos))
}
causal.pos <- which(mkr==causal.pos)
myhaplo.dat <- myhaplo.dat[,mkr]
geno.cases <- geno.cases[,mkr]
geno.controls <- geno.controls[,mkr]
haplo.dat <-  haplo.dat[,mkr]
mut.dat <- mut.dat[mkr]
haplo.cases <- haplo.cases[,mkr]
haplo.controls <- haplo.controls[,mkr]
causal.pos.mut <- mut.dat[causal.pos]

cat(paste("\nThe causal SNP is SNP", causal.pos,"\n"))
cat(paste("\nThe physical position of the causal SNP is", causal.pos.mut,"\n"))

### Remove the causal SNP
unobs.myhaplo.dat <- myhaplo.dat[,-causal.pos]
unobs.haplo.cases <- haplo.cases[,-causal.pos]
unobs.haplo.controls <- haplo.controls[,-causal.pos]
unobs.mut.dat <- mut.dat[-causal.pos]

### Save the data in the root folder
filehaplotypes <- paste(root,"haplotypes.txt",sep = "")
write.table(unobs.myhaplo.dat,file=filehaplotypes, row.names = FALSE, col.names = FALSE)
filepositions <- paste(root,"positions.txt",sep = "")
write.table(unobs.mut.dat,file=filepositions, row.names = FALSE, col.names = FALSE)
filestatus <- paste(root,"status.txt",sep = "")
write.table(myhaplo.status,file=filestatus, row.names = FALSE, col.names = FALSE)

#################################################################################

