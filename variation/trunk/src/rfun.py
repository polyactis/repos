"""
A file that contains various helpful functions to deal with R.

"""
import util


def plotVectorsOnSameGraph(x,vectorList, main="", xlab="", ylab="",type=None):
    pass

def plotVectors(x, vectorList, main="", xlab="", ylab="",type=None, xname="x", ynames=None):
    """
    Writes out a simple R string to plot the vectors..
    """
    if not ynames:
        ynames = ["y"]*len(vectorList)

    x = util.valListToStrList(x)
    rstr =""
    rstr = "par(mfrow=c("+str(len(vectorList))+",1));\n"
    rstr += xname+" <- c("+",".join(x)+");\n"
    for i in range(0,len(vectorList)):
        y = util.valListToStrList(vectorList[i])
        rstr += ynames[i]+" <- c("+",".join(y)+");\n"
        rstr += 'plot('+xname+','+ynames[i]+',pch=20, main="'+main+'",xlab="'+xlab+'",ylab="'+ylab+'"'
        if type:
            rstr += ', type="'+type+'"'
        rstr += ')\n'
            
    return rstr


def plotOverlayingVectors(x, vectorList, main="", xlab="", ylab="",type=None, pch='20', xname="x", ynames=None):
    """
    Writes out a simple R string to plot the vectors..
    """
    if not ynames:
        ynames = ["y"]*len(vectorList)
    maxVal=[]
    minVal=[]
    for v in vectorList:
        maxVal.append(max(v))
        minVal.append(min(v))
    maxVal = max(maxVal)
    minVal = min(minVal)

    xmax = max(x)
    xmin = min(x)

    x = util.valListToStrList(x)
    rstr =""
    #rstr = "par(mfrow=c(1,1));\n"
    rstr += xname+" <- c("+",".join(x)+");\n"
    for i in range(0,len(vectorList)):
        y = vectorList[i]
        y = util.valListToStrList(y)
        rstr += ynames[i]+" <- c("+",".join(y)+");\n"
        if i!=0:
            rstr += "par(new=T);\n"
        rstr += 'plot('+xname+','+ynames[i]+',main="'+main+'",xlab="'+xlab+'",ylab="'+ylab+'", xlim=c('+str(xmin)+','+str(xmax)+'), ylim=c('+str(minVal)+','+str(maxVal)+'), col='+str(i+2)
        if type:
            rstr +=', type="'+type+'"'
        if pch:
            rstr +=', pch='+pch
        rstr += ')\n'
            
    return rstr


import os
def runRFile(rfile, outputFile=None):
    os.system("R --vanilla --file="+rFile+" > "+outputFile)

import tempfile
tempfile.tempdir='/tmp/'
def runR(rstr, outputFile):
    i,rFile = tempfile.mkstemp()
    f = open(rFile,"w")
    f.write(rStr)
    f.close()
    os.system("R --vanilla --file="+rFile+" > "+outputFile)
