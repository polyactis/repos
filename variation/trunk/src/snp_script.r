#-------------------------------------------------------------------------------
#for suzi's flowering times
#read in genotypes, phenotypes, loop kruskal over each snp, write p values to file
#p_ftxx.txt
#markers are in markers.txt
#genotypes in genotypes.csv
#phenotypes in phenotypes.snp
setwd('g:/')
phen=read.table('phenotypes.csv', sep=',')
names(phen)=c('acc','ft22','ft16','ft10')
snp.n=248584
acc.n=93
log.file='log.txt'
tm=Sys.time()
for (j in 1:3) {
  gene=matrix(scan('genotypes.csv', what='character', sep=','), nrow=acc.n, ncol=snp.n, byrow=TRUE)
  p=vector(mode='numeric', length=snp.n)
  write.table(paste('Starting',names(phen)[j+1]), log.file, quote=FALSE, row.names=FALSE, col.names=FALSE), append=TRUE)
  for (snp in 1:snp.n) {
    if(length(levels(factor(gene[,snp])))<2) {
      p[snp]=1
    }else {
      p[snp]=kruskal.test(x=phen[,j+1], g=factor(gene[,snp]))$p.value
    }
  if (snp%%10000==0) {write.table(
                        paste(names(phen)[j+1], snp, Sys.time()-tm), 
                        log.file, quote=FALSE, row.names=FALSE, col.names=FALSE), 
                        append=TRUE)}
  }
  write.table(p, paste(names(phen)[j+1],'.txt', sep=''), quote=FALSE, row.names=FALSE, col.names=FALSE)
}
#-------------------------------------------------------------------------------
#draw graph of bonferroni corrected values for suzi's flowering times

setwd('g:/')
p.10=read.table('ft10.txt')
p.16=read.table('ft16.txt')
p.22=read.table('ft22.txt')
p.adj10=p.adjust(p.10[,1],'bonferroni')
p.adj16=p.adjust(p.16[,1],'bonferroni')
p.adj22=p.adjust(p.22[,1],'bonferroni')
phen=read.table('phenotypes.csv', sep=',')
names(phen)=c('acc','ft22','ft16','ft10')
markers=read.table('markers.txt', sep='.')
names(markers)=c('chr','marker')

p.adj=data.frame(p10=p.adj10, p16=p.adj16, p22=p.adj22)

#results=data.frame(markers, p.adj)

end.chr=c(59634, 92014, 142728, 186008, 248584)

last.chr=119175179
png('bonp_temps.png', width=3000*17/11, height=3000)
par(mfcol = c(3, 1))                   
titles=c('10 deg C','16 deg C','22 deg C')

for (j in 1:3) {

  xlim=c(0,  last.chr)
  ylim=c(0, 3)
  m.chr=0
  plot(0,0, type='n', axes=F, xlim=xlim, ylim=ylim, xlab='', ylab='')
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  for (i in 1:5) {
    pick=markers$chr==i
    text(m.chr+5000, 3, paste('Chromosome', i),col='light grey', pos=4 )
    lines(m.chr+markers$marker[pick], -log10(p.adj[pick,j]), type='l', col='black')
    m.chr=m.chr+max(markers$marker[pick])+1000
    abline(v=m.chr, lty='dotted', col='dark grey')
  }
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  text(last.chr, -log10(0.05), '95% significant', col='light grey', adj=c(1, -0.3))
  axis(side=2, cex=0.7)
  title(main=titles[j], ylab='-log p value')

}

dev.off()
#-------------------------------------------------------------------------------
#draw graph of raw p values for suzi's flowering times
setwd('g:/')
p.10=read.table('ft10.txt')
p.16=read.table('ft16.txt')
p.22=read.table('ft22.txt')

png('rawp_temps.png', width=1000*17/11, height=1000)
par(mfcol = c(3, 1))
                   
titles=c('10 deg C','16 deg C','22 deg C')
phen=read.table('phenotypes.csv', sep=',')
names(phen)=c('acc','ft22','ft16','ft10')
markers=read.table('markers.txt', sep='.')
names(markers)=c('chr','marker')

p=data.frame(p10=p.10, p16=p.16, p22=p.22)

for (j in 1:3) {

  xlim=c(0,  last.chr)
  ylim=c(0, 10)
  m.chr=0
  plot(0,0, type='n', axes=F, xlim=xlim, ylim=ylim, xlab='', ylab='')
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  title(main=titles[j], ylab='-log raw p value')
  for (i in 1:5) {
    pick=markers$chr==i
    text(m.chr+5000, 10, paste('Chromosome', i),col='light grey', pos=4)
    lines(m.chr+markers$marker[pick], -log10(p[pick,j]), type='l', col='black')
    m.chr=m.chr+max(markers$marker[pick])+1000
    abline(v=m.chr, lty='dotted', col='dark grey')
  }
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  text(last.chr, -log10(0.05), '95% significant', col='light grey', adj=c(1, -0.3))
  axis(side=2, cex=0.7)

}
dev.off()
#-------------------------------------------------------------------------------
#write significant markers to file
temp=c(10,16,22)
for (i in 1:3) {
pick=p.adj[,i]<=0.05
write.table(markers[pick,], paste('markers',temp[i],'.txt', sep=''), quote=FALSE, row.names=FALSE, col.names=FALSE)

}

#-------------------------------------------------------------------------------
#plot ternary diagram

library(ade4)
setwd('g:/')
png('ternary.png', width=2000, height=2000)
phen=read.table('phenotypes.csv', sep=',')
names(phen)=c('acc','ft22','ft16','ft10')
triangle.plot(phen[,2:4], label = phen$acc, clabel = 0.5, 
    cpoint = 1, draw.line = TRUE, addaxes = FALSE, addmean = FALSE, 
    labeltriangle = TRUE, sub = 'Ternary plot of flowering time at three different temperatures', csub = 1, possub = "topright", 
    show.position = TRUE, scale = TRUE, min3 = NULL, max3 = NULL, 
    box = FALSE)
dev.off()


#-------------------------------------------------------------------------------
#read in hairy data
#markers are in markers.txt
#genotypes in genotypes.csv
#phenotypes in hairy.txt
setwd('g:/')
phen=read.table('hairy.txt')
names(phen)='hairy'
pick=phen$hairy=='G'
hair=vector(length=93, mode='numeric')
hair[pick]=1
hair[!pick]=0
snp.n=248584
acc.n=93
log.file='log_hairy.txt'
tm=Sys.time()
  gene=matrix(scan('genotypes.csv', what='character', sep=','), nrow=acc.n, ncol=snp.n, byrow=TRUE)
  p=vector(mode='numeric', length=snp.n)
  #write.table(paste('Starting',names(phen)[j+1]), log.file, quote=FALSE, row.names=FALSE, col.names=FALSE), append=TRUE)
  for (snp in 1:snp.n) {
    if(length(levels(factor(gene[,snp])))<2) {
      p[snp]=1
    }else {
      p[snp]=kruskal.test(x=hair, g=factor(gene[,snp]))$p.value
    }
  if (snp%%10000==0) {write.table(
                        paste(snp, Sys.time()-tm), 
                        log.file, quote=FALSE, row.names=FALSE, col.names=FALSE, 
                        append=TRUE)}
  }
  write.table(p, paste('hairyp','.txt', sep=''), quote=FALSE, row.names=FALSE, col.names=FALSE)

#-------------------------------------------------------------------------------
#draw graph of bonferroni corrected values for suzi's hairy times

setwd('g:/')
p.hairy=read.table('hairyp.txt')
p.adj=p.adjust(p.hairy,'bonferroni')
phen=read.table('hairy.txt', sep=',')
markers=read.table('markers.txt', sep='.')
names(markers)=c('chr','marker')

last.chr=119175179
#png('bonp_temps.png', width=3000*17/11, height=3000)
                   
titles='p values for hairy/non-hairy phenotypes (bonferroni)'

  xlim=c(0,  last.chr)
  ylim=c(0, 16)
  m.chr=0
  plot(0,0, type='n', axes=F, xlim=xlim, ylim=ylim, xlab='', ylab='')
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  for (i in 1:5) {
    pick=markers$chr==i
    text(m.chr+100000, 3, paste('Chromosome', i),col='light grey', pos=4 )
    lines(m.chr+markers$marker[pick], -log10(p.adj[pick,]), type='l', col='black')
    m.chr=m.chr+max(markers$marker[pick])+100000
    abline(v=m.chr, lty='dotted', col='dark grey')
  }
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  text(last.chr, -log10(0.05), '95% significant', col='light grey', adj=c(1, -0.3))
  axis(side=2, cex=0.7)
  title(main=titles, ylab='-log p value')


#dev.off()

#-------------------------------------------------------------------------------

#for someone else's flowering times
#read in genotypes, phenotypes, loop kruskal over each snp, write p values to file
#p_ftold.txt
#markers are in markers.txt
#genotypes in genotypes.csv
#phenotypes in mariadata.csv
setwd('g:/')
phen=read.table('mariadata.csv', sep=',', header=TRUE)
snp.n=248584
acc.n=93
log.file='log_ftold.txt'
tm=Sys.time()

  gene=matrix(scan('genotypes.csv', what='character', sep=','), nrow=acc.n, ncol=snp.n, byrow=TRUE)
  p=vector(mode='numeric', length=snp.n)
  for (snp in 1:snp.n) {
    if(length(levels(factor(gene[,snp])))<2) {
      p[snp]=1
    }else {
      p[snp]=kruskal.test(x=phen$FT, g=factor(gene[,snp]))$p.value
    }
  if (snp%%10000==0) {write.table(
                        paste(snp, Sys.time()-tm), 
                        log.file, quote=FALSE, row.names=FALSE, col.names=FALSE, 
                        append=TRUE)}
  }
  write.table(p, 'p_ftold.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)

#-------------------------------------------------------------------------------
#draw graph of bonferroni corrected values for whoever else's flowering times

setwd('g:/')
p.ftold=read.table('p_ftold.txt')
p.adj=p.adjust(p.ftold[,1],'bonferroni')
markers=read.table('markers.txt', sep='.')
names(markers)=c('chr','marker')

last.chr=119175179
png('bonp_temps.png', width=3000*17/11, height=3000)                 
titles=c("Adjusted p values for Maria's data")


  xlim=c(0,  last.chr)
  ylim=c(0, 3)
  m.chr=0
  plot(0,0, type='n', axes=F, xlim=xlim, ylim=ylim, xlab='', ylab='')
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  for (i in 1:5) {
    pick=markers$chr==i
    text(m.chr+100000, 3, paste('Chromosome', i),col='light grey', pos=4 )
    lines(m.chr+markers$marker[pick], -log10(p.adj[pick]), type='l', col='black')
    m.chr=m.chr+max(markers$marker[pick])+1000
    abline(v=m.chr, lty='dotted', col='dark grey')
  }
  abline(h=-log10(0.05), lty='dashed', col='light grey')
  text(last.chr, -log10(0.05), '95% significant', col='light grey', adj=c(1, -0.3))
  axis(side=2, cex=0.7)
  title(main=titles, ylab='-log p value')

dev.off()