### filename: utilities.R

############################################################################


rowEntropy <- function(p) rowMeans(rowSums(log2(p^p), dims=2))


###########################################################################


rowIndepChiSqTest <- function(call1,call2){
  tmp=vector("numeric",nrow(call1))
  for(i in 1:nrow(call1)){
    tmpt= cbind(table(factor(call1[i,],levels=1:3)),
      table(factor(call2[i,],levels=1:3)))
    tmpt=tmpt[!rowSums(tmpt)==0,,drop=FALSE]
    if(nrow(tmpt)>1){
      rowtot <- rowSums(tmpt)
      coltot <- colSums(tmpt)
      e=outer(rowtot,coltot)/sum(tmpt)
      stat = sum((tmpt-e)^2/e)
      tmp[i] = 1 - pchisq(stat, (nrow(tmpt) - 1) * (ncol(tmpt) - 1))
    }
    else{
      tmp[i]=0
    }
  }
  return(tmp)
}

###########################################################################
trigammaInverse <- function(x) {
#	Solve trigamma(y) = x for y
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 March 2004.

#	Non-numeric or zero length input
	if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
	if(length(x)==0) return(numeric(0))

#	Treat out-of-range values as special cases
	omit <- is.na(x)
	if(any(omit)) {
		y <- x
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 0)
	if(any(omit)) {
		y <- x
		y[omit] <- NaN
		warning("NaNs produced")
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x > 1e7)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/sqrt(x[omit])
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 1e-6)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/x[omit]
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}

#	Newton's method
#	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
#	so iteration to solve 1/x = 1/trigamma is monotonically convergent
	y <- 0.5+1/x
	iter <- 0
	repeat {
		iter <- iter+1
		tri <- trigamma(y)
		dif <- tri*(1-tri/x)/psigamma(y,deriv=2)
		y <- y+dif
		if(max(-dif/y) < 1e-8) break
		if(iter > 50) {
			warning("Iteration limit exceeded")
			break
		}
	}
	y
}











##################################################################################################object is intensity matrix with snpid
fitAffySnpMixture <- function(object, df2=5, probs=rep(1/3,3), eps=50, subSampleSize=10^4, verbose=TRUE) {

 I <- nrow(object)/4
 J <- ncol(object)
 nprobe <- nrow(object)


 set.seed(1)
 tmp <- c( (1:I),((I+1):(2*I)))
 idx <- sort(sample(tmp, subSampleSize))
 rm(tmp)


 pis <- array(0,dim=c(I,J,3,2))   
 fs <- array(0,dim=c(I,J,2))
 snr <- array(0,dim=J)


 dimnames(fs)<-list( rownames(object)[ seq(1, nrow(object), 4) ],
                    colnames(object),
                    c("antisense","sense"))
 dimnames(pis)<-list( rownames(object)[ seq(1, nrow(object), 4) ],
                      colnames(object),
                      c("AA","AB","BB"),
                      c("antisense","sense"))
 names(snr) <- colnames(object)





 if(verbose) cat("Fitting mixture model to ",J," arrays. Epsilon must reach ",eps,".\n",sep="")


 for(j in 1:J){

  Y <- c( object[,j][seq(1, nprobe, 4)]-object[,j][seq(3, nprobe, 4)], object[,j][seq(2, nprobe, 4)]- object[,j][seq(4, nprobe, 4)] )  #row1,3 antisense, 2,4 sense
  A <- c( (object[,j][seq(1, nprobe, 4)]+object[,j][seq(3, nprobe, 4)] )/2, (object[,j][seq(2, nprobe, 4)]+ object[,j][seq(4, nprobe, 4)] ) /2 )


  mus <- quantile(Y,c(1,3,5)/6); mus[2]=0  
  sigmas <- rep(mad(c(Y[Y<mus[1]]-mus[1],Y[Y>mus[3]]-mus[3])),3)  
  sigmas[2] <- sigmas[2]/2

    
  a <- A[idx]; A=A-mean(A); a=a-mean(a)
  y <- Y[idx]
    


  weights <- apply(cbind(mus, sigmas), 1, function(p) dnorm(y, p[1], p[2]))
  PreviousLogLik <- -Inf
  change <- eps+1
  itmax <- 0
  
  
  matA <- ns(a,df2)
  while (change > eps & itmax < 100){
      itmax <- itmax+1
      
      ## E
      z <- sweep(weights, 2, probs, "*")
      LogLik <- rowSums(z)
      z <- sweep(z, 1, LogLik, "/")
      LogLik <- sum(log(LogLik))
      change <- abs(LogLik-PreviousLogLik)
      
      if (verbose){
        if (itmax > 1 | j > 1) cat(del)
        message <- paste("Array ",j,": epsilon = ", signif(change,2), "  ", sep="")
        del <- paste(rep("\b", nchar(message)), collapse="")
        cat(message)
      }
      
      PreviousLogLik <- LogLik
      probs <- colMeans(z)

      ## M
      fit1 <- lm(y~matA,weights=z[,1])
      fit2 <- sum(z[,2]*y)/sum(z[,2])
      fit3 <- lm(y~matA,weights=z[,3])
      
      sigmas[1] <- sqrt(sum(z[,1]*residuals(fit1)^2)/sum(z[,1]))
      sigmas[2] <- sqrt(sum(z[,2]*(y-fit2)^2)/sum(z[,2]))
      sigmas[3] <- sqrt(sum(z[,3]*residuals(fit3)^2)/sum(z[,3]))
      
      weights[,1] <- dnorm(y,fitted(fit1),sigmas[1])
      weights[,2] <- dnorm(y,fit2,sigmas[2])
      weights[,3] <- dnorm(y,fitted(fit3),sigmas[3])
      weights[y >= 0, 1] <- 0
      weights[y <= 0, 3] <- 0
    }


    gc()
    bigX <- cbind(1,  ns(A, knots=as.numeric(attr(matA, "knots")), Boundary.knots=attr(matA, "Boundary.knots")))
    rm(matA); gc()

    pred1 <- bigX%*%coef(fit1)
    pred2 <- rep(fit2,length(Y))
    pred3 <- bigX%*%coef(fit3)
    rm(bigX); gc()

    weights <- matrix(0,length(Y),3)
    weights[,1] <- dnorm(Y,pred1,sigmas[1])
    weights[,2] <- dnorm(Y,pred2,sigmas[2])
    weights[,3] <- dnorm(Y,pred3,sigmas[3])
    weights[Y >= pred2, 1] <- 0
    weights[Y <= pred2, 3] <- 0

    z <- sweep(weights, 2, probs, "*")
    LogLik <- rowSums(z)
    z <- sweep(z, 1, LogLik, "/")
  
    fs[,j,] <- matrix((pred3-pred1)/2,ncol=2)

    for(k in 1:3){
      pis[,j,k,] <- matrix(z[,(4-k)],ncol=2) ##4-k cause 3is1,2is2 and 1is3
    }

    snr[j] <- median(fs[,j,])^2/(sigmas[1]^2+sigmas[2]^2)

  }
  
  if(verbose) cat("Done.\n")
  return(list(f0=median(fs),fs=fs, pis=pis, snr=snr))
}




###################################################################################### object is from EM
getInitialAffySnpCalls <- function(object,subset=NULL,
                                   concordanceCutoff=0.0001,
                                   cutoffs=c(0.7,0.5,0.7),
                                   returnProbs=FALSE,
                                   verbose=FALSE){
  if(is.null(subset)) subset <- 1:(dim(object$pis)[1])
  pi1 <- object$pis[subset,,,1]
  pi2 <- object$pis[subset,,,2]; rm(object); gc()


 
  if(verbose) cat("Picking good starting value: ")
  if(verbose) cat("Computing entropy, ")
  E1 <- rowEntropy(pi1)
  E2 <- rowEntropy(pi2); gc()
  
  tmpN <- dim(pi1)
  tmpcall1 <- tmpcall2 <- matrix(NA, nrow=tmpN[1], ncol=tmpN[2])
  rownames(tmpcall1) <- rownames(tmpcall2) <- dimnames(pi1)[[1]]
  colnames(tmpcall1) <- colnames(tmpcall2) <- dimnames(pi1)[[2]]
  gc()
  if(verbose) cat("calculating calls, ")
  for (i in 1:tmpN[1]){
    for (j in 1:tmpN[2]){
      tmpcall1[i, j] <- which.max(pi1[i,j,])
      tmpcall2[i, j] <- which.max(pi2[i,j,])
    }
  }
  
  if(verbose) cat("determining non-concordant calls, ")
  concordance <- rowIndepChiSqTest(tmpcall1,tmpcall2); gc()

  if(verbose) cat("deciding which strand(s) to use")
  tc1 <- tmpcall1 == 1
  tc2 <- tmpcall1 == 2
  tc3 <- tmpcall1 == 3
  noABIndex1 <- which((rowSums(tc2)<3)*(rowSums(tc3)>0)*(rowSums(tc1)>0) == 1)
  tc1 <- tmpcall2 == 1
  tc2 <- tmpcall2 == 2
  tc3 <- tmpcall2 == 3
  noABIndex2 <- which((rowSums(tc2)<3)*(rowSums(tc3)>0)*(rowSums(tc1)>0) == 1)
  rm(tc1, tc2, tc3); gc()
  
  E1[noABIndex1] <- -Inf
  E2[noABIndex2] <- -Inf

  jointprobs <- (pi1+pi2)/2

  gc()
  noInfoIndex <- intersect(noABIndex1, noABIndex2)
  
  rm(noABIndex1, noABIndex2)
  
  jointprobs[noInfoIndex,,] <- 1/3
  rm(noInfoIndex)
  
  notBoth <- concordance < concordanceCutoff
  E1bE2 <- E1 > E2
  rm(E1, E2)
  i1 <- notBoth & E1bE2
  i2 <- notBoth & !E1bE2
  rm(notBoth, E1bE2)
  jointprobs[i1,,] <- pi1[i1,,]
  jointprobs[i2,,] <- pi2[i2,,]
  rm(i1, i2, pi1, pi2); gc()

  tmpN <- dim(jointprobs)
  tmpcall <- tmpmax <- matrix(NA, nrow=tmpN[1], ncol=tmpN[2])
  rownames(tmpcall) <- rownames(tmpmax) <- dimnames(jointprobs)[[1]]
  colnames(tmpcall) <- colnames(tmpmax) <- dimnames(jointprobs)[[2]]
  gc()

  if(verbose) cat(" finalizing"); gc()

  for (i in 1:tmpN[1]){
    for (j in 1:tmpN[2]){
      tmpcall[i, j] <- which.max(jointprobs[i, j, ])
      tmpmax[i, j] <- jointprobs[i, j, tmpcall[i,j]]
    }
  }
  
  for(i in 1:3)
    tmpcall[tmpcall==i & tmpmax<cutoffs[i]] <- NA
  if(verbose) cat("\nCompleted!\n")
  if(returnProbs) return(list(calls=tmpcall,probs=jointprobs)) else return(tmpcall)
}



###############################################################################################calculate center and scale

getGenotypeRegionParams <- function(M, initialcalls, f=0, verbose=TRUE){
  if(verbose) cat("Computing centers and scales for 3 genotypes")
  tmp <- .Call("R_HuberMatrixRows2",  M+(initialcalls-2)*f, as.integer(initialcalls), 1.5)
  if(!is.matrix(M)) M <- matrix(M,ncol=2)
  centers <- scales <- N <- array(NA,dim=c(nrow(M),3))
  dimnames(centers) <- dimnames(scales) <- dimnames(N) <- list(rownames(M),c("AA","AB","BB"))
  centers[,] <- tmp[[1]]
  scales[,2] <- tmp[[2]][,2]
  N[,2] <- tmp[[3]][,2]
  N[,-2] <- rowSums(tmp[[3]][,-2], na.rm=T)
  scales[,-2] <- sqrt(rowSums(tmp[[2]][,-2]^2*tmp[[3]][,-2], na.rm=T)/N[,-2])
  if(verbose) cat(" Done\n")
  return(list(centers=centers,scales=scales,N=N))
}



getAffySnpGenotypeRegionParams<-function(object, initialcalls, f, verbose=TRUE) {#object is intensity matrix, f$fs from EM
  

  N <- scales <- centers <- array(NA, dim=c( nrow(object)/4,3,2 ) )
  featureNames <- rownames(object) [ seq(1, nrow(object), 4) ]  
  dimnames(N) <- dimnames(scales) <- dimnames(centers) <- list( featureNames, c("AA","AB","BB"), c("antisense","sense"))


  if(verbose) cat("Computing centers and scales:\n")
  for(s in 1:2){
    gc()
    if(verbose) cat(c("\tantisense","\tsense....")[s],": ",sep="")

    if(s==1) M <- object [ seq(1, nrow(object), 4), ]- object [ seq(3, nrow(object), 4), ]  #row1,3 antisense, 2,4 sense
    if(s==2) M <- object [ seq(2, nrow(object), 4), ]- object [ seq(4, nrow(object), 4), ]  #row1,3 antisense, 2,4 sense

    tmp <- getGenotypeRegionParams(M, initialcalls, f[,,s], verbose=verbose) 
                                   
    centers[,,s] <- tmp$centers
    scales[,,s] <- tmp$scales
    N[,,s] <- tmp$N
  }
  N <- apply(N, 1:2, max, na.rm=TRUE)
  return(list(centers=centers,scales=scales,N=N,f0=median(f, na.rm=TRUE) ) ) 
}






#########################################################################################object is from genotyperegionparams
getAffySnpPriors <-  function(object,minN,subset=1:(dim(object$centers)[1]),
                              verbose=TRUE){
  if(verbose) cat("Computing priors.\n")
  N <- cbind(object$N,object$N)
#  Index <- subset[ which(rowMeans(N[subset,]>= minN)==1) ]  #I changed N[subset,] >= minN
   Index <- subset[ which(rowMeans(N[subset,]) >= minN) ]  #modified by R. J. on 19 April 2008

  N <- N[Index,]
  mus <- cbind(object$centers[Index,,1],object$centers[Index,,2])
  sigmas <- cbind(object$scales[Index,,1],object$scales[Index,,2])
  maxsigmas <- c(quantile(sigmas[,c(1,3,4,6)],.99, na.rm=TRUE),
                 quantile(sigmas[,c(2,5)],.99, na.rm=TRUE))
  
  ##variance for prior for mu
  ##need to make robust
  V <- cov(mus,use="pairwise.complete")
  
  ##prior for sigma
  zgs <- log(sigmas^2); rm(sigmas); gc()
  dgs <- N-1
  egs <- zgs-digamma(dgs/2)+log(dgs/2); rm(zgs); gc()
  n <- length(Index)
  d0s <- 2*trigammaInverse( colMeans(n*(egs-colMeans(egs, na.rm=TRUE))^2/(n-1)-trigamma(dgs/2), na.rm=TRUE) )
  s20 <- exp(colMeans(egs, na.rm=TRUE)+digamma(d0s/2)-log(d0s/2))
  return(list(V=V,d0s=d0s,s20=s20,maxsigmas=maxsigmas))
}

##minN- minimum number of points in cluster required for use



######################################################################################################object is from genotyperegionparams, priors is from snpriors
updateAffySnpParams <- function(object, priors, missingStrandIndex, minN=3,
                                maxHomoSigma=priors$maxsigma[1],
                                maxHeteSigma=priors$maxsigma[2],
                                subset=1:(dim(object$centers)[1]),
                                d0s=80, verbose=FALSE){

  object$centers <- object$centers[subset,,]
  object$scales <- object$scales[subset,,]
  object$N <- object$N[subset,]
  missingStrandIndex <- missingStrandIndex[subset]
  if(verbose) cat("Updating centers and scales")

 ##variances

  for(i in 1:2){
	#modified by R. J., to consider Homozygotes ONLY
    for(j in 1){ ##1 and 3 are the same
    #for(j in 1:2){ ##1 and 3 are the same

      if(j==2) N <- object$N[,2] else N<-rowSums(object$N[,c(1,3)],na.rm=TRUE)
      s <- object$scales[,j,i]
      if (is.null(d0s))
        d0s <- priors$d0s[3*(i-1)+j]
      s20 <- priors$s20[3*(i-1)+j] ##notice the ad-hoc choice of 3
      Index <- N>minN & !is.na(s)
      N<-N[Index];s <- s[Index]
      object$scales[Index,j,i] <- sqrt (  ( (N-1)*s^2 + d0s*s20 ) / (d0s+N-1) )
      object$scales[!Index & missingStrandIndex != i,j,i] <- sqrt(s20)
    }
  }
  object$scales[,3,] <- object$scales[,1,] ##AA=BB 
  object$scales[,2,][object$scales[,2,]>maxHeteSigma] <- maxHeteSigma
  object$scales[,c(1,3),][object$scales[,c(1,3),]>maxHomoSigma] <- maxHomoSigma
  if(verbose) cat(".")

 
 ##Means

  updateMean <- function(strandIndex, strand){
    ## strandIndex is a vector where
    ## 0 - none strands are missing on the array
    ## 1 - antisense strand is missing
    ## 2 - sense strand is missing
    ## strand argument is what do you want to be updated
    
    if (strand == "both"){
      snps <- strandIndex == 0; idx <- c(1,3,4,6) #modified by R. J., i.e., no Heterozygotes
      #snps <- strandIndex == 0; idx <- 1:6
      
	#modified by R. J.,
      N <- cbind(object$N[snps, c(1,3)],object$N[snps, c(1,3)])
      mu <- cbind(object$centers[snps, c(1,3),1],object$centers[snps, c(1,3),2])
      #N <- cbind(object$N[snps,],object$N[snps,])
      #mu <- cbind(object$centers[snps,,1],object$centers[snps,,2])

    }else if(strand == "antisense"){
      snps <- strandIndex == 2; idx <- 1:3
      N <- object$N[snps,]
      mu <- object$centers[snps,,1]
    }else if(strand == "sense"){
      snps <- strandIndex == 1; idx <- 4:6
      N <- object$N[snps,]
      mu <- object$centers[snps,,2]
    }
    Vinv <- solve(priors$V[idx, idx])
    NSinv <- t(N)/priors$s20[idx]
    tmp <- t(sapply(1:nrow(mu),function(i){
      if(verbose & i%%5000==0)  cat(".")
      mus=mu[i,]; Ns=N[i,]
      mus[Ns< minN]<-0 
      return(solve(Vinv+diag(NSinv[,i]))%*%(NSinv[,i]*mus))
    }))
    return(tmp)
  }


  if (any(missingStrandIndex == 0)){
    tmp <- updateMean(missingStrandIndex, "both")
    object$centers[missingStrandIndex == 0,c(1,3),1] <- tmp[, 1:2]
    object$centers[missingStrandIndex == 0,c(1,3),2] <- tmp[, 3:4]; rm(tmp);
    #object$centers[missingStrandIndex == 0,,1] <- tmp[,1:3]
    #object$centers[missingStrandIndex == 0,,2] <- tmp[,4:6]; rm(tmp);
  }
  if (any(missingStrandIndex == 2))
    object$centers[missingStrandIndex == 2,,1] <- updateMean(missingStrandIndex, "antisense")
  if (any(missingStrandIndex == 1))
    object$centers[missingStrandIndex == 1,,2] <- updateMean(missingStrandIndex, "sense")
  if(verbose) cat("Done.\n")
  return(object)
}






########################################################################################################## object is intensity matrix
getAffySnpDistance <- function(object, params, f=0, subset=1:(nrow(object)/4),w=NULL,verbose=FALSE){


 x1 <- object[ seq(1, nrow(object), 4), ]- object[ seq(3, nrow(object), 4), ] 
 x2 <- object[ seq(2, nrow(object), 4), ]- object[ seq(4, nrow(object), 4), ] 
 x <- cbind( x1, x2)
 x <- array(x, dim=c( nrow(object)/4, ncol(object), 2) )

  rm(object)
  Dist <- array(NA,dim=c(dim(x)[1],dim(x)[2],3,2))
  if(verbose) cat("Calculating likelihood-based distances")
  for(i in 1:2){
    if(verbose) cat(".")
    for(j in 1:3){
      tmp <- x[,,i]+(j-2)*f[subset,,i]
      Dist[,,j,i] <- 2*log(params$scales[subset,j,i]) +
        ((tmp-params$centers[subset,j,i])/params$scales[subset,j,i])^2
      if(!is.null(w)) Dist[,,j,i] <-  Dist[,,j,i] - 2*log(w[subset,,j,i])
    }
  }
  dimnames(Dist) <- list(dimnames(x)[[1]],
                         dimnames(x)[[2]], 1:3,
                         dimnames(x)[[3]])
  if(verbose) cat("Done.\n")
  return(Dist)
}



################################################################################################################ dist from snpdistance
getAffySnpCalls <- function(Dist,subset=1:(dim(Dist)[1]),
                            verbose=FALSE, option = 1){
  Dist <- Dist[subset,,,, drop=FALSE]
  
  ## Here we can even put default value to 3 for the new code. /HB
  res <- array(as.integer(-1),dim=dim(Dist)[c(1,2)]); gc()
  dimnames(res) <- dimnames(Dist)[1:2]
  Dist <- rowSums(Dist, na.rm=TRUE, dims=3)
  ##  Dist[XIndex, maleIndex, 2] <- Inf
  
  ##the following is slow!
  if(verbose) cat("Making calls for ", ncol(res), " arrays");


	if (option == 1)
{
  for(j in 1:ncol(res)){
    if(verbose) cat(".");
    D1 <- Dist[,j,1];
    D2 <- Dist[,j,2];
    D3 <- Dist[,j,3];
    d12 <- (D1 < D2);
    d23 <- (D2 < D3); rm(D2);
    d13 <- (D1 < D3); rm(D3);
    d <- rep(as.integer(3), length(D1)); rm(D1)
    d[( d12 & d13)] <- as.integer(1);
    d[(!d12 & d23)] <- as.integer(2);
    res[,j] <- d;
    rm(d);
  }
}
else # call homozygotes only
{
  for(j in 1:ncol(res)){
    if(verbose) cat(".");
    D1 <- Dist[,j,1];
    D3 <- Dist[,j,3]; 
    d13 <- (D1 < D3); rm(D3);
    d <- rep(as.integer(3), length(D1)); rm(D1)
    d[d13] <- as.integer(1);
    d[(!d13)] <- as.integer(3);
    res[,j] <- d;
    rm(d);
  }
	
}

  if(verbose) cat("Done\n")
  return(res)
}




############################################################################################################################# confidence for call

getAffySnpConfidence <- function(Dist, Calls, 
                                 subset=1:nrow(Calls), verbose=TRUE){
  Dist <- Dist[subset,,,,drop=FALSE]
  Calls <- Calls[subset,,drop=FALSE]
 
  res <- array(NA,dim=dim(Dist)[c(1,2)])
  dimnames(res) <- list(dimnames(Dist)[[1]],dimnames(Dist)[[2]])
  Dist <- rowSums(Dist, dims=3, na.rm=T)
  
  cat("Computing confidence for calls on ",ncol(res)," arrays")
  ##apply is faster apply but takes too much memory
  Index <- 1:nrow(Calls)
  for(j in 1:ncol(res)){

#	commented from the source code --R.J.
#    Index2 <- Index
#    if (verbose) cat(".")
#    tmpdist <- cbind(abs(Dist[,j,1]-Dist[,j,2]),abs(Dist[,j,2]-Dist[,j,3]))
#    tmpIndex <- split(Index2, factor(Calls[Index2,j], levels=1:3), drop=FALSE)
#    if (length(tmpIndex[[1]])>0) res[tmpIndex[[1]],j] <- tmpdist[tmpIndex[[1]],1]
#    if (length(tmpIndex[[3]])>0) res[tmpIndex[[3]],j] <- tmpdist[tmpIndex[[3]],2]
#    if (length(tmpIndex[[2]])>0) res[tmpIndex[[2]],j] <- pmin(tmpdist[tmpIndex[[2]], 1],
#                                                              tmpdist[tmpIndex[[2]], 2])
#    rm(tmpIndex, tmpdist); gc()

	### added by R. J.
	Index2 <- Index
	if (verbose) cat(".")
	res[, j] <-abs(Dist[,j,1]-Dist[,j,3])

	### end

  }
  cat("Done\n")  
  return(res)
}




##################################################################################################################

replaceAffySnpParams <- function(object,value,subset){
  object$centers[subset,,] <- value$centers
  object$scales[subset,,] <- value$scales
  object$N[subset,] <- value$N
  return(object)
}


























