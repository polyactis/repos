"""
Various useful functions.
"""

def valListToStrList(l):
    newl = []
    for val in l:
        newl.append(str(val))
    return newl


def calcVar(self,list):
    mean = sum(list)/float(len(list))
    var = 0
    for i in range(0,len(list)):
        var = var + (list[i]-mean)*(list[i]-mean)
    var = var/float(len(list)-1)
    return var
