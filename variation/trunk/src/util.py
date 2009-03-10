"""
Various useful functions.
"""
def getRanks(values,withTies=True):
	"""
	returns ranks (large values w. large ranks)
	"""
	srt_vals = zip(values[:],range(0,len(values)))
	srt_vals.sort()
	srt_ranks = []
	if withTies:	
		i = 0
		while i < len(srt_vals):
			curVal =  srt_vals[i][0]
			j = 1
			while i+j < len(srt_vals) and srt_vals[i+j][0]==curVal:
				j += 1
			rank = i+(j+1)/2.0
			max_i = i+j
			while i<max_i:
				srt_ranks.append(rank)
				i+=1
	else:
		srt_ranks = range(1,len(srt_vals))
			
	ranks = [0.0]*len(srt_vals)
	
	for i in range(0,len(srt_vals)):
		(val,val_index) = srt_vals[i]		
		ranks[val_index] = srt_ranks[i]
	return ranks
			

def kruskal_wallis(snps,phenVals,useTieCorrection=True):
	from rpy import r
	def _kw_test_(snp,ranks,tieCorrection):
		n = len(snp)
		group_vals = list(set(snp))
		g = len(group_vals)
		snp_ranks = zip(snp,ranks)
		ns = [0]*g
		rs = [0.0]*g
		for (nt,rank) in snp_ranks:
			i = 0
			while group_vals[i] != nt:
				i += 1
			ns[i]+=1
			rs[i]+=rank
		nominator = 0.0
		mean_r = (n+1)/2.0
		for i in range(0,g):
			v = (rs[i]/ns[i])-mean_r
			nominator += ns[i]*v*v

		d = ((12.0/(n*(n+1))) * nominator)*tieCorrection
		p =  r.pchisq(d,g-1,lower_tail = False)
		return (d,p)

	
	ranks = getRanks(phenVals)
	tieCorrection = 1.0
	if useTieCorrection:
		n_total = len(ranks)
		groups = list(set(ranks))
		n_groups = len(groups)
		s = 0
		for g in groups:
			t = ranks.count(g)
			s += t*(t*t-1.0)
		s = float(s)/(n_total*(n_total*n_total-1))
		tieCorrection = 1.0/(1-s)
		print "tieCorrection:", tieCorrection
		
	ds = []
	ps = []
	for snp in snps:
		(d,p) = _kw_test_(snp, ranks,tieCorrection=tieCorrection)
		ds.append(d)
		ps.append(p)
	return {"ps":ps,"ds":ds}



def plotHist(x,y):
	pass


class Queue:
	"A first-in-first-out queue."
	def __init__(self, items=None): self.start = 0; self.A = items or []
	def __len__(self):				return len(self.A) - self.start
	def append(self, item):		   self.A.append(item)
	def extend(self, items):		  self.A.extend(items)

	def pop(self):
		A = self.A
		item = A[self.start]
		self.start += 1
		if self.start > 100 and self.start > len(A)/2:
			del A[:self.start]
			self.start = 0
		return item


def valListToStrList(l):
	newl = []
	for val in l:
		newl.append(str(val))
	return newl


def calcVar(l):
	mean = sum(l)/float(len(l))
	var = 0
	for i in range(0,len(l)):
		var = var + (l[i]-mean)*(l[i]-mean)
	var = var/float(len(l)-1)
	return var

def _transposeDoubleListsOld_(l):
	newl = [[] for i in range(0,len(l[0]))]
	for i in range(0,len(l)):
		for j in range(0,len(l[0])):
			newl[j].append(l[i][j])
	return newl

def _transposeDoubleListsOld2_(l):
	import numpy
	m = numpy.matrix(l)
	m = numpy.transpose(m)
	return m.tolist()

def transposeDoubleLists(l):
	return map(list,zip(*l))
	
def calcSD(l):
	mean = sum(l)/float(len(l))
	var = 0
	for i in range(0,len(l)):
		var = var + (l[i]-mean)*(l[i]-mean)
	var = var/float(len(l)-1)
	import math
	return math.sqrt(var)

def calcQuantiles(numbers):
	"""
	Calculates the 1st quantile.
	"""
	import math
	numbers.sort()
	
	if len(numbers)%2==0:
		#Even
		ns1 = numbers[0:len(numbers)/2]
		ns2 = numbers[len(numbers)/2:len(numbers)]
		median = (numbers[len(numbers)/2-1]+numbers[len(numbers)/2])/2.0
	else:
		#Odd
		ns1 = numbers[0:(len(numbers)-1)/2]
		ns2 = numbers[(len(numbers)+1)/2:len(numbers)]		
		median = numbers[(len(numbers)-1)/2]

	if len(ns1)%2 ==0:
		#Even
		q1=(ns1[len(ns1)/2-1]+ns1[len(ns1)/2])/2.0
		q3=(ns2[len(ns2)/2-1]+ns2[len(ns2)/2])/2.0
	else:
		q1 = ns1[(len(ns1)-1)/2]
		q3 = ns2[(len(ns2)-1)/2]

	return (q1,median,q3)
	
def calcIQR(numbers):
	"""
	Calculates the inter-quantile range.
	"""
	quantiles = calcQuantiles(numbers)
	print quantiles
	return quantiles[2]-quantiles[0]
	
	
def bin_counts(values,bins):
	counts = [0]*(len(bins)-1)
	out_of_bounds_values = []
	for value in values:
		i = 0
		while i+1 < len(bins) and value > bins[i+1]:
			i += 1
		if i+1 < len(bins) and bins[i] < value <= bins[i+1]:
			counts[i] += 1
		else: 
			out_of_bounds_values.append(value)
	if len(out_of_bounds_values):
		print "These values were OOB:",out_of_bounds_values
	return counts


if __name__=='__main__':
	l = [1,2,3,3,3,2,3,4,1,0,3,3,3]
	ranks = getRanks(l)
	print ranks
	print "Done!"

	