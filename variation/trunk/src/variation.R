#2007-02-19 snp_variation_data = read.csv('../data/justin_data.csv', header=TRUE, na.strings='0', sep='\t', row.names=1)
snp_variation_data = read.csv('../data/justin_data.csv', header=TRUE, sep='\t', row.names=1)	#'0' is NA but princomp can't handle it
sv = snp_variation_data[,seq(2,149)]
sv.category = snp_variation_data[,1]
sv.pca = princomp(sv)	#cor=T encounters error due to 0 standard deviation among some x rows(same nt across all locus)
summary(sv.pca)
plot(sv.pca)
loadings(sv.pca)
sv.pc = predict(sv.pca)
plot(sv.pc[,1:2], type="n")
text(sv.pc[,1:2], labels=as.character(sv.category), col=3+unclass(sv.category))

#2007-02-20 try MDS
snp_variation_data = read.csv('../data/justin_data.csv', header=TRUE, na.strings='0', sep='\t', row.names=1)
sv = snp_variation_data[,seq(2,149)]
sv.category = snp_variation_data[,1]
sv.dist = dist(sv, method='binary')
no_of_NAs = sum(is.na(sv.dist))
print(no_of_NAs)

sv.scal = cmdscale(sv.dist, k=2)	#due to NA in dist, this step can't go further, sammon() or isoMDS() also works
plot(sv.scal$points, type='n')
text(sv.scal$points, labels=as.character(sv.category), col=3+unclass(sv.category), cex=0.8)

#2007-02-20 stupid loop to find NA but it's too slow
no_of_entries = dim(sv)[1]
no_of_snps = dim(sv)[2]
no_of_NA_dist = 0
for (i in seq(no_of_entries))
{
	for (j in seq(i+1, no_of_entries))
	{
		no_of_NAs_in_this_pair = 0
		for (k in seq(no_of_snps))
		{
			if (is.na(sv[i,k]) && is.na(sv[j,k]))
			{
			no_of_NAs_in_this_pair = no_of_NAs_in_this_pair + 1
			}
		}
		if (no_of_NAs_in_this_pair == no_of_snps)
		{
			no_of_NA_dist = no_of_NA_dist + 1
			cat(i, '\t', j, '\n')
		}
	}
}

#2007-02-20 try daisy() instead of dist()

snp_variation_data = read.csv('../data/justin_data.csv', header=TRUE, sep='\t', row.names=1, colClasses='factor', na.strings='0')
sv = snp_variation_data[,seq(2,149)]
sv.category = snp_variation_data[,1]

library(cluster)
sv.dist = daisy(sv)

#impute the missing value with average_dist
no_of_NAs = sum(is.na(sv.dist))
average_dist = sum(sv.dist, na.rm=TRUE)/(length(sv.dist)-no_of_NAs)
for (i in seq(length(sv.dist)))
{
	if(is.na(sv.dist[i]))
	{
	no_of_NAs = no_of_NAs  +1
	sv.dist[i] = average_dist
	}
}
#MDS
sv.scal = cmdscale(sv.dist, k=2)	#NA must be filled in
plot(sv.scal, type='n')
text(sv.scal, labels=as.character(sv.category), col=3+unclass(sv.category), cex=0.8)

library(MASS)
sv.sam = sammon(sv.dist)	#can't handle zero or negative distance
plot(sv.sam, type='n')
text(sv.sam, labels=as.character(sv.category), col=3+unclass(sv.category), cex=0.8)

sv.iso = isoMDS(sv.dist)	#same error as sammon()
plot(sv.iso, type='n')
text(sv.iso, labels=as.character(sv.category), col=3+unclass(sv.category), cex=0.8)

#try hclust
h = hclust(sv.dist)
plclust(h)

#2007-02-20 try cat package
snp_variation_data = read.csv('../data/justin_data.csv', header=TRUE, na.strings='0', sep='\t', row.names=1)
sv = snp_variation_data[,seq(2,149)]
sv.category = snp_variation_data[,1]
library(cat)
s <- prelim.cat(sv)
thetahat <- em.cat(s)
theta <- da.cat(s,thetahat,1)
ximp  <- imp.cat(s,theta)

#2007-02-27
no_of_NA_per_strain = rep(0, dim(sv)[1])
for (i in seq(dim(sv)[1]))
{
	no_of_NA_per_strain[i] = sum(is.na(sv[i,]))
}
hist(no_of_NA_per_strain)
