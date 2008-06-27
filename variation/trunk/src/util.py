"""
Various useful functions.
"""

def valListToStrList(l):
    newl = []
    for val in l:
        newl.append(str(val))
    return newl


def calcVar(self,l):
    mean = sum(l)/float(len(l))
    var = 0
    for i in range(0,len(l)):
        var = var + (l[i]-mean)*(l[i]-mean)
    var = var/float(len(l)-1)
    return var

def transposeDoubleLists(l):
    newl = [[] for i in range(0,len(l[0]))]
    for i in range(0,len(l)):
        for j in range(0,len(l[0])):
            newl[j].append(l[i][j])
    return newl
