#!/usr/bin/python
"""
Examples:
	NPUTE.py -m 0 -w 10 -i sample_data.csv -o sample_data.out

Description:
	Program to impute genotypes.
	
	Input file format:
		SNP by individual. "?" is NA.

2008-04-30 yh start modifying
"""
import sys
import os
import getopt
import time
from SNPData import *
from CircularQueue import *
from numpy import *

#Option Names
MODE_TYPE = '-m'
SING_WIN = '-w'
FILE_WIN = '-W'
RANGE_WIN = '-r'
IN_FILE = '-i'
OUT_FILE = '-o'

#Mode Types
IMP = '0'
TST = '1'

'''
This is the main NPUTE class, providing a command-line interface for imputation.
'''
import sys, os
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
	
class NPUTE:
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	option_default_dict = {('mode_type', 1, ): [IMP, 'm', 1, 'specify running mode. 0=Imputation, 1=test window sizes', ],\
							('single_window_size', 0, ): ['', 'w', 1, 'specify a window size, like 10', ],\
							('window_file', 0, ): ['', 'W', 1, 'A file with each line a number for window size. To test window sizes.', ],\
							('window_size_range', 0, ): ['', 'p', 1, 'specify a window range to test, like 5:15', ],\
							('input_fname',1, ): ['', 'i', 1, 'Input file. A plain genotype matrix.'],\
							('output_fname',1, ): ['out.csv', 'o', 1, 'Output File'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-05-01
			use ProcessOptions
		"""
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def main(self):
		'''
		2008-05-01
			yh use ProcessOptions
			
		Parses arguments, loads data, calls proper functions, and outputs results.
		'''
		
		"""
		options = {MODE_TYPE: IMP, # Mode - imputation or window test
				   SING_WIN : '', # Single Window
				   RANGE_WIN : '', # Window Range - 'start:end'
				   FILE_WIN : '', # Window File
				   IN_FILE : 'in.csv', # Input File
				   OUT_FILE : 'out.csv' # Output File
				   }
		
		optlist, args = getopt.getopt(sys.argv[1:], 'm:w:W:r:i:o:')
		for opt in optlist:
			options[opt[0]] = opt[1]
	
		Get input SNPs
		inFile = options[IN_FILE]
		"""	   
		inFile = self.input_fname
		
		if not os.path.exists(inFile):
			print "Input file '%s' not found." % inFile
			sys.exit(1)
		snpData = SNPData(inFile)
		
		# Get test windows
		L = []
		if not isEmpty(self.window_size_range):
			start,stop = self.window_size_range.split(':')
			start,stop = int(start),int(stop)
			L = range(start,stop+1)
		winFile = self.window_file
		if not isEmpty(winFile):
			if not os.path.exists(winFile):
				print "Window file '%s' not found." % winFile
				sys.exit(1)
			lines = file(winFile,'r').readlines()
			L += [int(line) for line in lines]
		if not isEmpty(self.single_window_size):
			L += [int(self.single_window_size)]
		L.sort()		
	
		mode = self.mode_type
		if mode == IMP:
			if isEmpty(self.single_window_size):
				print 'Imputation window not specified.'
				sys.exit(1)
			L = int(self.single_window_size)
			imputeData(snpData, L, self.output_fname)
		elif mode == TST:
			if isEmpty(L):
				print 'Test windows not specified.'
				sys.exit(1)
			testWindows(snpData, L, self.output_fname)
		
		
def isEmpty(x):
	'''
	Helper function to test if an object is empty (has 0 length).
	'''
	return len(x) == 0

def imputeData(snpData, L, outFile):
	'''
	Main function for doing a real imputation on a SNPData object and outputting results.
	'''
	print outFile
	start = time.time()
	c = impute(snpData, L)
	t = int(time.time()-start + 0.5)
	snpData.incorporateChanges()
	print 'Imputed %d unknowns in %dm %ds.' % (c,t/60,t%60)
	snpData.outputData(outFile)

def testWindows(snpData, Ls, outFile):
	'''
	Main function for testing the imputation accuracy of multiple windows on a SNPData object
	and outputting results.
	'''
	start = time.time()
	c, corrects = testImpute(snpData, Ls)
	t = int(time.time()-start + 0.5)
	print 'Imputed %d called values over %d windows in %dm %ds.' % (c,len(Ls),t/60,t%60)
	outputWinAccs(Ls,c,corrects,outFile)
	

def impute(snpData, L):
	'''
	Function that slides the window and calls other functions to do the actual imputation.
	'''
	
	global count
	
	snpData.changes = dict()
	count = 0
	snps = snpData.snps
	vectors = snpData.vectors
	numSNPs = len(snps)

	print "Imputing with window size " +  str(L) + "...",
	
	vectorLength = len(vectors.values()[0])
	vectorQueue = CircularQueue([L], vectorLength)
	acc = zeros(vectorLength, uint16)

	# Initialize queue
	for i in xrange(L):
		snpVector = vectors[snps[i]]
		vectorQueue.queue[i] = snpVector
		add(acc,snpVector,acc)

	# Begin impute
	for i in xrange(numSNPs):
		
		if i+L < numSNPs:
			snpVector = vectors[snps[i+L]]
			vectorQueue.enqueue(snpVector)
		else:
			vectorQueue.enqueue(zeros(vectorLength,uint16))
			
		top,bottom = vectorQueue.getEnds(0)
		add(acc,top,acc)
		subtract(acc,bottom,acc)
		snp = snps[i]
		if '?' in snp:
			imputeSNP(snpData,i,acc,snp)
						 
	print "Done"
	return count

def imputeSNP(snpData,locI,mmv,snp):
	'''
	Uses the window's mismatch vector to impute each missing value in SNP.
	'''
	global count

	for samp in xrange(snpData.numSamps):   
		if snp[samp] == '?':

			score = snpData.extractRow(mmv,samp) 

			sA = argsort(score)
			impNuc = getMinImp(snp,sA[0:-1],score)

			snpData.changes[(locI,samp)] = impNuc
			count += 1


def testImpute(snpData, Ls):
	'''
	Function that slides the window(s) and calculates accuracy of imputation
	on all called values.
	'''
	
	global corrects
	global count

	L = max(Ls)
	vectors = snpData.vectors
	snps = snpData.snps
	numSNPs = len(snps)

	print 'Running imputation window test with %d window sizes...' % len(Ls),

	count = 0
	corrects = zeros(len(Ls))
	
	vectorLength = len(vectors.values()[0])
 
	vectorQueue = CircularQueue(Ls,vectorLength)
	acc = zeros((len(Ls),vectorLength),uint16)
	# Initialize queue
	for i in xrange(L):
		snpVector = vectors[snps[i]]
		vectorQueue.queue[i] = snpVector
		for j in xrange(len(Ls)):
			if i < Ls[j]:
				add(acc[j],snpVector,acc[j])

	# Begin impute
	for i in xrange(numSNPs):
		if i+L < numSNPs:
			vectorQueue.enqueue(vectors[snps[i+L]])
		else:
			vectorQueue.enqueue(zeros(vectorLength,uint16))
			
		mid = vectorQueue.getMid()
		snp = snps[i]

		for j in xrange(len(Ls)):
			top,bottom = vectorQueue.getEnds(j)
			add(acc[j],top,acc[j])
			subtract(acc[j],bottom,acc[j])
			imputeSNPT(snpData,i,acc[j]-mid,snp,j)
	
	print 'Done'

	return count, corrects	


	
def imputeSNPT(snpData,locI,mmv,snp,j):
	'''
	Uses the window's mismatch vector to test impute each known value in SNP and
	check for correctness.
	'''
	global count
	global corrects

	# This is done so that singleton values are not attempted to be imputed	
	if snp.count('1') == 1:
		checkOne = True
	else:
		checkOne = False
	for samp in xrange(snpData.numSamps):   
		if snp[samp] != '?' and not (checkOne and snp[samp] == 1):

			if j == 0:						   
				count += 1
				 
			score = snpData.extractRow(mmv,samp) 

			sA = argsort(score)
			impNuc = getMinImp(snp,sA[0:-1],score)
			
			if impNuc == snp[samp]:
				corrects[j] += 1
				
  

def getMinImp(snp, sA, score):
	'''
	Finds nearest neighbor to sample being imputed w/ a called value.  If tere is a tie,
	uses next nearest neighbor and so on.
	'''
	lastM = 0
	points = 0
	winner = '0'

	for i in sA:
		
		m = score[i]

		if snp[i] != '?':
			if m != lastM and points > 0:
					return winner
			else:
				if snp[i] == winner:
					points += 1
				else:
					if points == 0:
						winner = snp[i]
						points = 1
					else:
						points -= 1					

			lastM = m

	return winner

def outputWinAccs(Ls,count,corrects,outFile):
	'''
	Outputs the accuracy of imputation on all called values for each window size tested
	to a CSV file.
	'''
	print "Writing estimated window accuries to '%s'..." % outFile,
	accs = corrects/float(count)
	out = ''
	for i in xrange(len(Ls)):
		out += '%d,%f\n' % (Ls[i],accs[i])
	file(outFile,'w').write(out)
	print "Done"

if __name__ == "__main__":
	from pymodule import ProcessOptions
	main_class = NPUTE
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.main()
