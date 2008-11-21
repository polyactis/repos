"""
Various useful functions.
"""
class Queue:
    "A first-in-first-out queue."
    def __init__(self, items=None): self.start = 0; self.A = items or []
    def __len__(self):                return len(self.A) - self.start
    def append(self, item):           self.A.append(item)
    def extend(self, items):          self.A.extend(items)

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
    
def calcSD(self,l):
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
	
	