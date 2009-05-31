import sys, os
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv
from CircularQueue import *
from sets import *
from variation.src.snpsdata import RawSnpsData
from variation.src import dataParsers
from pymodule import PassingData

#Infinite Number (for acc array)
INF = 2**16-1 # set for unsigned 16-bit integers used in acc array

#Allele Types
MAJ = '0'
MIN = '1'
I_MAJ = 'Z'
I_MIN = 'W'

class SNPData(object):
	'''
	This class provides a datastructure for SNP datasets being imputed with NPUTE. Input
	files should be either CSV or TXT files with lines as SNPs (separated by crlf) and
	columns as samples.  The alleles in the SNPs may be separated by commas or nothing.
	The SNPs must be in order by chromosome position.  The SNPs must be ternary with
	a majority allele, a minority allele (not required), and a unknowns (not required).
	Majoriy and minority alleles can be represented by any characters, but an unknown
	must be a '?' character.
	'''
	option_default_dict = {('inFile', 0, ): ['', ],\
							('snps_name_ls', 0, ): [None, ],\
							('data_matrix', 0, ): [None, ],\
							('chromosome', 0, ): [None, ],\
							('input_NA_char', 1, ): ['0', ],\
							('input_file_format', 1, int): [1, 'f', 1, 'which file format. 1= DB_250k2data.py output. \
								2=original NPUTE input. 3=Output250KSNPs.py output wtih array Id.\
								4=API for MpiQCCall.py to call'],\
							('lower_case_for_imputation', 1, int): [0, ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		'''
		05/09/08
			change interface, inFile could be a data structure
		Constructor reads in data file, generates mismatch vectors and upper
		triangular matrix extract indices, and reports stats of data.
		'''
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		start = time.time()
		#if isinstance(inFile, str) and os.path.isfile(inFile):
		if self.input_file_format == 2:
				self.readInData(self.inFile)
		elif self.input_file_format==1:
				snps, nucs, self.chosen_snps_name_ls = self.readOneChromosomeData(self.snps_name_ls, self.data_matrix, self.chromosome)
				self.snps = array(snps)
				self.sdps = Set(snps)
				self.nucs = array(nucs)
				self.numSamps = len(self.snps[0])
		#elif isinstance(inFile, RawSnpsData):	#05/07/08 it's snpsd data structure
		elif self.input_file_format==3:
			snpsd = dataParsers.parseCSVData(self.inFile, withArrayIds=True)
			self.input_NA_char = 'NA'
			self.snps, self.sdps, self.nucs, self.numSamps = self.getDataStructureFromSNPsD(snpsd[0])
		elif self.input_file_format==4:
			self.snps, self.sdps, self.nucs, self.numSamps = self.getDataStructureFromSNPsD(self.inFile)
		else:
			sys.stderr.write('unsupported inFile format: %s\n' % self.input_file_format)
			
		
		self.genMismatchVectors()
		self.genExtractIndices()
		t = int(time.time()-start+0.5)
		sys.stderr.write('Number of Samples: %d\n' % self.numSamps)
		sys.stderr.write('Number of SNPs: %d\n' % len(self.snps))
		sys.stderr.write('Number of SDPs: %d' % len(self.sdps))
		sys.stderr.write('Time to Process: %d m %d s\n' % (t/60,t%60))
	
	def getDataStructureFromSNPsD(self, snpsd):
		"""
		05/07/08
		"""
		sys.stderr.write("Reading data ...")
		no_of_rows = len(snpsd.positions)
		no_of_cols = len(snpsd.accessions)
		snps = []
		nucs = []
		for i in range(no_of_rows):
			one_snp_ls, symbol2counts = self.get_symbol2counts(snpsd.snps, fixed_index=i, no_of_rolls=no_of_cols, by_row=0)
			
			passingdata = self.get_symbol2MAJ_MIN(symbol2counts)
			if passingdata.symbol2MAJ_MIN==3:
				sys.stderr.write("Error: SNP %s (%s) has more than 2 alleles: %s.\n"%(i, snpsd.positions[i], repr(symbol2counts)))
				sys.exit(2)
			
			map_func = lambda x: passingdata.symbol2MAJ_MIN[x]
			one_snp_ls = map(map_func, one_snp_ls)
			
			snps.append(''.join(one_snp_ls))
			nucs += [(passingdata.major, passingdata.minor)]
		passingdata = PassingData()
		passingdata.snps  = array(snps)
		passingdata.sdps = Set(snps)
		passingdata.nucs = array(nucs)
		passingdata.numSamps = no_of_cols
		sys.stderr.write("Done.\n")
		return passingdata.snps, passingdata.sdps, passingdata.nucs, passingdata.numSamps
	
	def get_symbol2counts(self, data_matrix, fixed_index=0, no_of_rolls=1, by_row=1):
		#count the frequency of symbols
		one_snp_ls = []
		symbol2counts = {}
		for i in range(no_of_rolls):
			if by_row:
				symbol = data_matrix[i][fixed_index]
			else:
				symbol = data_matrix[fixed_index][i]
			one_snp_ls.append(symbol)
			if symbol != self.input_NA_char:	#don't care about NA
				if symbol not in symbol2counts:
					symbol2counts[symbol] = 0
				symbol2counts[symbol] += 1
		return one_snp_ls, symbol2counts
	
	def get_symbol2MAJ_MIN(self, symbol2counts):
		#construct a dictionary to map input symbols to MAJ, MIN or '?'
		symbol2MAJ_MIN = {self.input_NA_char:'?'}	#'NA' is always '?'
		symbols = symbol2counts.keys()
		if len(symbols) == 0:
			major = ''
			minor = ''
		elif len(symbols) == 1:
			symbol2MAJ_MIN[symbols[0]] = MAJ
			major = symbols[0]
			minor = ''
		elif len(symbols) ==2:
			major, minor = symbols
			if symbol2counts[major]<symbol2counts[minor]:
				minor, major = symbols	#reverse them
			symbol2MAJ_MIN[major] = MAJ
			symbol2MAJ_MIN[minor] = MIN
		elif len(symbols)>2:
			major, minor = None, None
			symbol2MAJ_MIN = 3
		passingdata = PassingData()
		passingdata.symbol2MAJ_MIN = symbol2MAJ_MIN
		passingdata.major = major
		passingdata.minor = minor
		return passingdata
	
	def readOneChromosomeData(self, snps_name_ls, data_matrix, chromosome):
		"""
		2008-05-19
			snps_name could be tuple or list
		05/07/08
			replace readInData()
			
			based on chromosome info extracted from snps_name_ls, only pick data from one chromosome
		"""
		sys.stderr.write("Reading chromosome %s data ..."%(chromosome))
		no_of_rows = len(data_matrix)
		no_of_cols = len(snps_name_ls)
		snps = []
		nucs = []
		chosen_snps_name_ls = []
		for i in range(no_of_cols):
			snps_name = snps_name_ls[i]
			if isinstance(snps_name, tuple) or isinstance(snps_name, list):
				tmp_ls = snps_name
			else:
				tmp_ls = snps_name.split('_')
			chr = tmp_ls[0]
			if chr!=chromosome:	#skip the SNPs from other chromosomes
				continue
			
			chosen_snps_name_ls.append(snps_name)
			
			one_snp_ls, symbol2counts = self.get_symbol2counts(data_matrix, fixed_index=i, no_of_rolls=no_of_rows, by_row=1)
			
			passingdata = self.get_symbol2MAJ_MIN(symbol2counts)
			if passingdata.symbol2MAJ_MIN==3:
				sys.stderr.write("Error: SNP %s (%s) has more than 2 alleles: %s.\n"%(i, snps_name_ls[i], repr(symbol2counts)))
				sys.exit(2)
			
			map_func = lambda x: passingdata.symbol2MAJ_MIN[x]
			one_snp_ls = map(map_func, one_snp_ls)
			
			snps.append(''.join(one_snp_ls))
			nucs += [(passingdata.major, passingdata.minor)]
		
		del data_matrix
		sys.stderr.write("Done.\n")
		return snps, nucs, chosen_snps_name_ls
	
	def readInData(self, inFile):
		'''
		Reads in data file and normalizes alleles.
		'''
		print "Reading in SNP data from '%s'..." % inFile,
		if not os.path.exists(inFile):
			print "Input file '%s' not found." % inFile
			sys.exit(1)
		snps = []
		nucs = []
		self.numSamps = -1		
		
		lines = file(inFile, 'r').readlines()
		for i in xrange(len(lines)):
			line = self.removeExtraChars(lines[i])
			if self.numSamps == -1:
				self.numSamps = len(line)
			elif self.numSamps != len(line):
				print '\nSNP %d has an inconsistent number of sample.' % i+1
				sys.exit(1)
			major, minor = self.getAlleles(line, i)
			if minor == '':
				snp = line.replace(major, MAJ)
			else:	
				snp = (line.replace(major, MAJ)).replace(minor, MIN)
			snps += [snp]
			nucs += [(major,minor)]
			
		self.snps = array(snps)
		self.sdps = Set(snps)
		self.nucs = array(nucs)
		print "Done"

	def getAlleles(self, s, i, NA_char='?'):
		'''
		05/07/08
			add NA_char
		
		Returns the majority and minority allele for the SNP.  Ties are broken by making the first
		known allele the majority.
		'''
		s = s.upper()
		q = s.count(NA_char)
		major_minor_count_threshold = (len(s)-q+1)/2
		major = ''
		minor = ''
		for c in s:
			if c == NA_char:
				pass
			elif major == c:
				pass
			elif minor == c:
				pass
			elif major == '' and s.count(c) >= major_minor_count_threshold:
				major = c
			elif minor == '' and s.count(c) <= major_minor_count_threshold:
				minor = c
			else:
				print 'SNP %d is not ternary. Extra allele: %s.'%((i+1), c)
		return major, minor
				
		
	def removeExtraChars(self,s):
		'''
		Helper function for reading in data.  Removes formatting characters.
		'''
		s = s.replace(',','')
		s = s.replace(' ','')
		s = s.replace('\n','')
		return s
		
	def genMismatchVectors(self):
		'''
		Generates a pair-wise mismatch vector for each SDP and stores it in a dictionary.
		Match = 0
		Mismatch = 2
		Unknown = 1
		'''
		sys.stderr.write("Generating pair-wise mismatch vectors...")

		n = self.numSamps
		o = [1 for i in range(n)]
		self.vectors = dict()

		for sdp in self.sdps:
			dM = []   
			p = []
			c = []
			for i in xrange(n):
				if sdp[i] == '1':
					p += [2]
					c += [0]
				elif sdp[i] == '0':
					p += [0]
					c += [2]
				else:
					c += [1]
					p += [1]

			q = 0			 
			for i in xrange(n):
				if p[i] == 0:
					dM += p[i+1:]
				elif c[i] == 0:
					dM += c[i+1:]
				else:
					dM += o[i+1:]
			
			self.vectors[sdp] = array(dM,uint16)
		sys.stderr.write("Done.\n")
	
	def incorporateChanges(self):
		'''
		Uses the changes dictionary to replace unknowns with imputed values.
		'''
		sys.stderr.write("Incorporating imputed values into the SNP data...")
		snps = self.snps
		for loc,val in self.changes.iteritems():
			locI, samp = loc
			snp = snps[locI]
			
			if val == MIN:
				iVal = I_MIN
			elif val == MAJ:
				iVal = I_MAJ
			else:
				print loc
			
			snps[locI] = snp[0:samp] + iVal + snp[samp+1:]
		sys.stderr.write("Done.\n")

	'''
	def outputData(self, outFile):
	#Outputs the SNP data to the specified file in csv format.  Lower-case values are imputed.
		print "Writing imputed data to '%s'..." % outFile,
		out = ''
		
		for i in xrange(len(self.snps)):
			snp = self.snps[i]
			for allele in snp:
				if allele == MAJ:
					nuc = self.nucs[i,0]
				elif allele == MIN:
					nuc = self.nucs[i,1]
				elif allele == I_MAJ:
					nuc = self.nucs[i,0].lower()
				elif allele == I_MIN:
					nuc = self.nucs[i,1].lower()
				out += nuc + ','
			out = out[:-1] + '\n'
		file(outFile,'w').write(out)
		print "Done"
	'''
		
	def _outputData(self, outFile): # revised to remove bottleneck: Serge Batalov
		'''
		Outputs the SNP data to the specified file in csv format.  Lower-case values are imputed.
		'''
		print "Writing imputed data to '%s'..." % outFile,
			   
		fi = file(outFile,'w')
		for i in xrange(len(self.snps)):
			snp = self.snps[i]
			out = ''
			for allele in snp:
				out += ',' + allele
			fi.write(out[1:] + '\n')
		print "Done"
	
	def allele2output_symbol(self, allele, nucs, i, NA_char='?'):
		if allele == MAJ:
			nuc = nucs[i,0]
		elif allele == MIN:
			nuc = nucs[i,1]
		elif allele == I_MAJ:
			nuc = nucs[i,0]
			if self.lower_case_for_imputation:
				nuc = nuc.lower()
		elif allele == I_MIN:
			nuc = nucs[i,1]
			if self.lower_case_for_imputation:
				nuc = nuc.lower()
		else:
			nuc = NA_char
		return nuc
	
	def outputOneChromosomeData(self, chosen_snps_name_ls, chromosome, snps, nucs, numSamps, output_fname):
		"""
		05/07/08
			change it to SNP by Strain
		05/07/08
			a matrix with only snp name as header
		"""
		sys.stderr.write("Outputting imputed chromosome %s data to %s ..."%(chromosome, output_fname))
		writer = csv.writer(open(output_fname, 'a'), delimiter='\t')	#watch, it's 'append' mode
		no_of_snps = len(chosen_snps_name_ls)
		for i in range(no_of_snps):
			data_row = [chosen_snps_name_ls[i], chosen_snps_name_ls[i]]	#2 times
			for j in range(numSamps):
				allele = snps[i][j]
				data_row.append(allele)
			writer.writerow(data_row)
		del writer
		sys.stderr.write("Done.\n")
	
	def translate2InputSymbols(self):
		sys.stderr.write("Translating data back to input symbols ...")
		no_of_snps = len(self.snps)
		data_matrix = []
		for i in range(no_of_snps):
			data_row = []
			for j in range(self.numSamps):
				allele = self.snps[i][j]
				inputsymbol = self.allele2output_symbol(allele, self.nucs, i, NA_char=self.input_NA_char)
				data_row.append(inputsymbol)
			data_matrix.append(data_row)
		self.snps = data_matrix
		sys.stderr.write("Done.\n")
	
	def outputData(self, outFile):
		"""
		05/07/08
		"""
		if self.chromosome==None:
			self._outputData(outFile)
		else:
			self.outputOneChromosomeData(self.chosen_snps_name_ls, self.chromosome, self.snps, self.nucs, self.numSamps,outFile)
	
		
	def genExtractIndices(self):
		'''
		Generate lookup table of indices for the rows of a matrix stored in upper-triangular form.
		'''
		sys.stderr.write('Generating indices for triangular matrix row extraction...')
		numSamps = self.numSamps
		
		extractIndices = zeros((numSamps,numSamps),int)
		for samp in xrange(numSamps):
			indices = []  

			j = samp-1
			for i in xrange(samp):
				indices += [j]
				j += numSamps -(i+2)

			indices += range(j,j+numSamps-samp)
			extractIndices[samp] = indices
		
		self.extractIndices = extractIndices
		sys.stderr.write("Done.\n")
	
	def extractRow(self,tM,rowNum):
		'''
		Extracts a row from an upper-triangular mismatch matrix for the dataset.  Sets value
		for imputed sample to inf so that it will not be chosen to impute itself during
		window tests.
		'''
		row = take(tM,self.extractIndices[rowNum])
		row[rowNum] = INF
		return row		
