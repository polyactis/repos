#recalibrate the parameters for genotype clusters

library(splines)
source("utilities-jr.R");

.fp <- function(f)
{
	as.numeric(rowSums(f,  na.rm=T)) 
}

### make calls
genocall <- function(dat)
{
#dat <- gDNA.perlegen;

ff <- fitAffySnpMixture (dat, df2=5, probs=rep(1/3,3), eps=50, verbose=TRUE)

Dist <- getAffySnpDistance (dat, params, f=ff$fs, subset=1:(nrow(dat)/4),w=NULL,verbose=FALSE)
Dist[,,-2,] <- 1.5*Dist[,,-2,]
myCalls <- getAffySnpCalls (Dist,subset=1:(dim(Dist)[1]), verbose=FALSE, option=2) #call homozygotes ONLY
LLR <- getAffySnpConfidence(Dist, myCalls, subset=1:nrow(myCalls), verbose=TRUE)

return(list(myCalls = myCalls, LLR = LLR, ff = ff));
}



### logistic regression
calpacc <- function(dat, index, ff, LLR)
{
	i <- 1;
	#index <- which(true != 0);
	#correct <- as.integer(myCalls[index]) == as.integer( true[index]);
	n <- length(index);
	snr <- rep(sqrt(ff$snr), each=n);
	llr <- as.numeric(sqrt(LLR[index]));
	idx.4 <- index *4;
	idx.3 <- idx.4 - 1;
	idx.2 <- idx.3 - 1;
	idx.1 <- idx.2 - 1;
	A <- (dat[idx.1]+dat[idx.2]+dat[idx.3]+dat[idx.4])/4;
	a <- as.numeric(A);
	fp <- .fp((ff$fs)[index,1,]);

	tmp <- cbind(1, snr, ns(llr, df = 3), ns(a, df=3), ns(fp, df=3)) %*% logit.coef;
	tmp <- 1/(1+exp(-tmp));
	return(tmp);
}

calrates <- function(pacc, true, myCalls, c)
{
	index <- which(true != 0);
	n <- length(index);
	correct <- as.integer(true[index]) == as.integer(myCalls[index]);
	callrate <- length(which(pacc > c))/n;
	accuracy <- length(which(pacc > c & correct)) / length(which(pacc > c));
	cat("callrate = ", callrate, "    accuracy = ", accuracy, "\n");
}
