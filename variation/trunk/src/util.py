"""
Various useful functions.
"""

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
    return zip(*l)
    
def calcSD(self,l):
    mean = sum(l)/float(len(l))
    var = 0
    for i in range(0,len(l)):
        var = var + (l[i]-mean)*(l[i]-mean)
    var = var/float(len(l)-1)
    import math
    return math.sqrt(var)
