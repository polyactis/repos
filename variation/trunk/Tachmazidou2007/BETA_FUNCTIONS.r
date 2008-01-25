############################################################################

findTrees<-function(seq){

  a <- nrow(seq) #number of haplotype sequences
  b <- ncol(seq) #number of segregating sites
  
  ### Initialize necessary variables ###
  ns <- numeric(b)
  compat <- numeric(b)

  rooted<-1
  c<-numeric(b)

  for(j in 1:b){
    c[j]<-1
    ns[j]<-1
    compat[j]<-1
  }

  # Check for segregating columns consistent
  for(l in 1:b){
    if((ns[l]!=0)&&(c[l]!=0)){
      
      j<-numeric(2)
      m<-0
      for(i in 1:a){
        j[seq[i,l]+1] <- j[seq[i,l]+1]+1
        m<-m+1
      }

      if((j[1]==m)||(j[2]==m)){
        ns[l]<-0
        c[l]<-0
      }

      if((rooted==1)&&(j[1]==m-1)) compat[l]<-0
    }
  }

  treeInColumn <- numeric(b)
  treeInColumn[1] <- 1
  tree_count <- 1
  tree.start <- 1
  
  for (i in 2:b) {

    for (j in (i-1):tree.start) {
      if((c[i]!=0)&&(c[j]!=0)&&(compat[i]!=0)&&(compat[j]!=0)){
        u <- table(seq[,j],seq[,i])
        condition <- (u[1,2]!=0)&&(u[2,1]!=0)&&(u[2,2]!=0)
        if (condition) {
          tree.start <- i
          tree_count <- tree_count + 1
          break
        }
      }
    }
    treeInColumn[i] <- tree_count
  }

  tree_boundaries <- matrix(0,2,tree_count)
  for(i in 1:tree_count){
    tree_boundaries[,i] <- range(which(treeInColumn==i))
  }

  return(list(ntrees=tree_count, boundaries=tree_boundaries, ns=ns))

}

categorical<-function(haplo.dat){

  wind<-as.matrix(haplo.dat)

  n <- dim(wind)[1] ### number of individuals
  m <- dim(wind)[2] ### number of SNPs in the window
  x <- numeric(n)

  haplo <- t(as.matrix(wind[1,]))
  k <- 1
  x[1] <- 1

  for(i in 2:n){

    temp <- t(as.matrix(wind[i,]))
    check <- matrix(0,k,1)
    ind <- 0

    for(j in 1:k){

      tt <- temp - haplo[j,]
      if(any(tt!=0)) check[j] <- 1
      if(all(tt==0)) ind <- j
      
    }

    if(all(check!=0)){

      haplo <- rbind(haplo,temp)
      k <- k+1
      x[i] <- k

    }else x[i] <- ind

    out <- list(haplo=haplo, k=k, x=x)
  }
}

HaploCount<-function(k, unqhapind, diseased){

  n.haplo<- numeric(k)
  n.diseased<- numeric(k)

  for(i in 1:k){

    l<- unqhapind == i
    n.haplo[i] <- sum(l)
    n.diseased[i]<- sum(diseased[l])

  }

  n.nondiseased <- n.haplo - n.diseased

  out <- list(n.haplo=n.haplo,n.diseased=n.diseased,n.nondiseased=n.nondiseased)
  return(out)

}

getRPLM<-function(outRoot,N,R,M,fileNamesR){

  setwd(outRoot)
  for(tr in 1:nrow(N)){
    myfile <- fileNamesR[tr]
    w <- as.matrix(cbind(M[[tr]],R[[tr]]))
    write.table(w, file=myfile, row.names=FALSE, col.names=FALSE) 
  }

  for(tr in 1:nrow(N)){
    if (file.exists(fileNamesR[tr])) {   
      hello <- system(sprintf("perl getUniqueM.pl %sRtree%d %stree%d",outRoot,tr,outRoot,tr),intern=TRUE)
    } else {
      return(FALSE)
    }
  } 
  return(TRUE)
}

getSEQ2TR<-function(outRoot,N,fileNamesMULT){
  
  setwd(outRoot)                       
  for(tr in 1:nrow(N)){
    if (file.exists(fileNamesMULT[tr])) {  
      hello1 <- system(sprintf("./seq2tr tree%d tree%d.tr",tr,tr),intern=TRUE)
    } else {
      return(FALSE)
    }
  }
  return(TRUE)
}

getTREEPIC<-function(outRoot,N,fileNamesTR){

  setwd(outRoot)                       
  for(tr in 1:nrow(N)){
    if (file.exists(fileNamesTR[tr])) {
     hello2 <- system(sprintf("./treepic tree%d.tr tree%d.ps",tr,tr),intern=TRUE)
    } else {
      return(FALSE)
    }
  }
  return(TRUE)
}

getMutOrder<-function(outRoot,N,fileNamesPS){
  
  setwd(outRoot)                       
  age.tree <- vector(nrow(N),mode="list")
  ok<-TRUE
  for(tr in 1:nrow(N)){
    if (file.exists(fileNamesPS[tr])) {
     hello3 <- system(sprintf("perl mutOrder.pl tree%d.ps",tr),intern=TRUE)
     file.mutOrder <- paste(outRoot,"SNPsOrder.txt", sep = "")
     if (file.exists(file.mutOrder)) {
       mutOrder.dat <- try(read.table(file.mutOrder, header= FALSE, sep = ""), silent = TRUE)
       if (inherits(mutOrder.dat, "try-error")) {
         ok<-FALSE
         break
       }
     } else {
       ok<-FALSE
       break
     }
    } else {
      ok<-FALSE
      break
    }
    if (ok) {
      age.tree[[tr]] <- mutOrder.dat[order(mutOrder.dat[,2],mutOrder.dat[,1],decreasing=FALSE),][,1]
    }
  }
  return(list(run = ok, age.tree=age.tree))
}

constructAGE<-function(trl, age){
  k <- max(trl)
  AGE <- matrix(-1,nrow=length(trl),ncol=k)
  AGE[1,1:length(age[[1]])] <- age[[1]]
  for(i in 2:nrow(AGE)){
    AGE[i,1:length(age[[i]])] <- age[[i]] + max(AGE[i-1,])
  }
  return(AGE)
}

removeFiles <- function(N, fileNamesR, fileNamesMULT, fileNamesTR){
  for(tr in 1:nrow(N)){
    unlink(fileNamesR[tr])
    unlink(fileNamesMULT[tr])
    unlink(fileNamesTR[tr])
  }
}

cluster<-function(MUTIND, centers, age){

  nc <- length(centers)
  N <- nrow(MUTIND)
  ns <- ncol(MUTIND)
  CLUSTIND <- matrix(0,N,1)

  cent<-centers
  ind<- 1
  m<-length(age)

  for(j in 1:m){

    x<-age[j]

    for(i in 1:nc){
      if(centers[i]==x){
        cent[ind]<-centers[i]
        ind<-ind+1
      }
    }
  }

  for(i in 1:nc){
    c<- cent[i]
    k<- MUTIND[,c]==1
    l<- CLUSTIND==0
    a<- k & l
    CLUSTIND[a]<- i
  }

  out <- list(CLUSTIND=CLUSTIND,corresp.centers=c(0,cent))
  return(out)

}

multTreesPartition<-function(alpha, beta, DISEASED, MUTIND, AGE, Ntree, XPRED, MCMC, PENALIZE){

  skip<- 10000 ##burn-in
  every<- 10 ##store every "every" value
  iter<- (MCMC-skip)/every 

  ns <- ncol(MUTIND) ##number of SNPs
  NPEOPLE <- nrow(MUTIND) ##number of haplotypes
  ### INITIALIZE THE PARTITION: EVERYBODY is in the null cluster 

  TREE <- 1
  nint <- Ntree[TREE,2]-Ntree[TREE,1]+1
  ind <- rep(0,nint)
  NDISEASED <- sum(DISEASED)

  part1 <- lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta)
  part2 <- lgamma(alpha + NDISEASED) + lgamma(beta + NPEOPLE - NDISEASED) - lgamma(alpha + beta + NPEOPLE)
  logmargold <- part1 + part2
  logmarginitial <- logmargold

  n.clustold <- 0
  centersold <- -1
  peopleold <- NPEOPLE
  diseaseold <- NDISEASED

  NCLUST <- numeric(iter)
  INDICATORS <- matrix(0,iter,ns)

  start <- Ntree[,1]
  end <- Ntree[,2]
  NSTREE <- end - start + 1
  TREERES <- numeric(iter)##tree result

  PVECTOR <- matrix(0,iter,nrow(XPRED))
  NPRED <- nrow(XPRED)
  TREEold <- TREE

  treelength<-numeric(nrow(Ntree))
  for(m in 1:nrow(Ntree)){
    treelength[m]<-Ntree[m,2]-Ntree[m,1]+1
  }

  for(i in 1:MCMC){

######
###MOVE STEP
########

    indold <- ind
    L <- indicator(indold)
    indnew <- L$X
    n.clust <- L$n.centers
    move <- L$move
    
    if(n.clust > 0){

      start <- Ntree[TREEold,1]
      end <- Ntree[TREEold,2]
      centers <- L$center.pos + (start-1)
      m <- end - start + 1
      temp <- AGE[TREEold,1:m]
      clust.ind <- cluster(MUTIND, centers, temp)

      # Find the number of people in each cluster
      n.people <- numeric(n.clust+1) ##number of haplotypes in each cluster
      n.diseased <- numeric(n.clust+1) ##number of diseased haplotypes in each cluster

      for(j in 1:(n.clust+1)){

        l<- clust.ind$CLUSTIND == (j-1)
        n.people[j] <- sum(l)
        n.diseased[j]<- sum(DISEASED[l])

      }

      ### EVALUATE THE MARG LIKE for the clust.ind and call it logmargnew
      temp <- lgamma(alpha + n.diseased) + lgamma(beta + n.people - n.diseased) - lgamma(alpha + beta + n.people)
      logmargnew <- (n.clust+1)*part1 + sum(temp)

    }else{

      logmargnew <- part1 + part2
      n.clust <- 0
      centers <- -1
      n.people <- NPEOPLE
      n.diseased <- NDISEASED

    }

    ### EVALUATE Bayes Factor
    lBF <- logmargnew - logmargold
    
    u <- runif(1)

    if(u < exp(lBF)){

      ind <- indnew
      logmargold <- logmargnew
      n.clustold <- n.clust
      centersold <- centers
      peopleold <- n.people
      diseaseold <- n.diseased

    }

######
###Pick tree step
##################

    indold <- ind
    L <- select_tree(Ntree)
    TREE <- L$TREE
    indnew <- L$X
    n.clust <- L$n.centers
              
    if(n.clust > 0){

      start <- Ntree[TREE,1]
      end <- Ntree[TREE,2]
      centers <- L$center.pos + start - 1
      m <- end - start + 1
      temp <- AGE[TREE,1:m]
      clust.ind <- cluster(MUTIND, centers, temp)
 
      # Find the number of people in each cluster
      n.people <- numeric(n.clust+1) ##number of haplotypes in each cluster
      n.diseased <- numeric(n.clust+1) ##number of diseased haplotypes in each cluster

      for(j in 1:(n.clust+1)){

        l<- clust.ind$CLUSTIND == (j-1)
        n.people[j] <- sum(l)
        n.diseased[j]<- sum(DISEASED[l])

      }

      ### EVALUATE THE MARG LIKE for the clust.ind and call it logmargnew
      temp <- lgamma(alpha + n.diseased) + lgamma(beta + n.people - n.diseased) - lgamma(alpha + beta + n.people)
      logmargnew <- (n.clust+1)*part1 + sum(temp)

    }else{

      logmargnew <- part1 + part2
      n.clust <- 0
      centers <- -1
      n.people <- NPEOPLE
      n.diseased <- NDISEASED

    }

    ### EVALUATE Bayes Factor
    lBF <- logmargnew - logmargold
    if((treelength[TREE]-treelength[TREEold])>0){
      logalpha <- lBF+log(0.02)*abs(treelength[TREE]-treelength[TREEold])
    }
    if((treelength[TREE]-treelength[TREEold])<0){
      logalpha <- lBF-log(0.02)*abs(treelength[TREE]-treelength[TREEold])
    }else{
      logalpha <- lBF
    }
    
    u <- runif(1)

    if(PENALIZE==0){
      st <- u < exp(lBF)
    }else if(PENALIZE==1){
      st <- u < exp(logalpha)
    }

    if(st==TRUE){

      ind <- indnew
      TREEold <- TREE
      logmargold <- logmargnew
      n.clustold <- n.clust
      centersold <- centers
      peopleold <- n.people
      diseaseold <- n.diseased

    }

    if( (i > skip) && ((i-skip)%%every == 0) ){

      NCLUST[(i-skip)/every] <- n.clustold + 1
      temp <- matrix(0,1,ns)
      start <- Ntree[TREEold,1]
      end <- Ntree[TREEold,2]
      temp[start:end] <- ind
      INDICATORS[(i-skip)/every,] <- temp
      TREERES[(i-skip)/every] <- TREEold
      

      if(n.clustold > 0){

        ### PREDICTION OF DISEASE FOR EACH UNIQUE HAPLOTYPE GIVEN THE CLUSTER IT BELONGS TO
        temp <- AGE[TREEold,1:NSTREE[TREEold]]
        clustpred <- cluster(XPRED, centersold, temp)

        for(j in 1:(n.clustold+1)){

          l <- clustpred$CLUSTIND == (j-1)
          temp <- rbeta(1, alpha + diseaseold[j], beta + peopleold[j] - diseaseold[j])
          PVECTOR[(i-skip)/every,l] <- temp

        }
      }else{

        temp <- matrix(rbeta(1, alpha + NDISEASED, beta + NPEOPLE - NDISEASED), 1, NPRED)
        PVECTOR[(i-skip)/every,] <- temp

      }
      
    }##END OF IF
    
  }##END OF MCMC

  return(list(NCLUST=NCLUST, TREERES=TREERES, INDICATORS=INDICATORS, PVECTOR=PVECTOR))

}

plotTree<-function(outRoot,tr, R, M){

  setwd(outRoot)
  myfile <- paste(outRoot,"pickedtr", sep = "")                                     
  sink(myfile)
  for(i in 1:nrow(R[[tr]])){
    cat(paste(M[[tr]][i],":"))
    cat(paste(" ",R[[tr]][i,]))
    cat("\n")
  }
  sink() 
  hello1 <- system("./seq2tr pickedtr pickedtr.tr",intern=TRUE)
  hello2 <- system("./treepic pickedtr.tr pickedtr.ps",intern=TRUE)
  unlink(paste(outRoot,"pickedtr", sep = ""))   
  unlink(paste(outRoot,"pickedtr.tr", sep = ""))   
                                       
}

subresults<-function(indtree,prob,stdev,Lpr,Upr,dat,boundaries){

  subrefined <- dat[,boundaries[indtree,][1]:boundaries[indtree,][2]]
  cat <- categorical(subrefined)
  proportSubUniqueHaplos <- as.numeric(table(cat$x))/sum(as.numeric(table(cat$x))) 
  z<-cbind(cat$x,prob,stdev,Lpr,Upr)
  postprob <- numeric(cat$k)
  poststdev <- numeric(cat$k)
  postL <- numeric(cat$k)
  postU <- numeric(cat$k)
  for(i in 1:cat$k){
    y <- z[,1]==i
    w <- z[y,]
    if(sum(y) == 1){
      postprob[i] <- w[2]
      poststdev[i] <- w[3]
      postL[i] <- w[4]
      postU[i] <- w[5]
    }else{
      postprob[i] <- mean(w[,2])
      poststdev[i] <- mean(w[,3])
      postL[i] <- mean(w[,4])
      postU[i] <- mean(w[,5])     
    }
  }

  return(list(catx=cat$x, catk=cat$k, cathaplo=cat$haplo, postprob=postprob, poststdev=poststdev, postL=postL, postU=postU))
}

prior<-function(ntree,snptree){

  n<- sum(snptree)
  PRIOR<- matrix(0,1,n)
  MCMC <- 100000
  
  for(i in 1:MCMC){
    
    tree<- sample(1:ntree,1)
    indicator<- rbinom(snptree[tree], 1, .5)
    
    if(tree==1){
      start<-1
      end<- snptree[tree]
    }
    else{
      start<- sum(snptree[1:(tree-1)])+1
      end<- sum(snptree[1:(tree-1)])+snptree[tree]
    }
    PRIOR[start:end] <- PRIOR[start:end]+indicator
    
  }### END OF MCMC

  PRIOR<- PRIOR/MCMC
  return(PRIOR)

}
