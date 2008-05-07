import sys, os
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv
from numpy import *
from sets import *


#Infinite Number (for acc array)
INF = 2**16-1 # set for unsigned 16-bit integers used in acc array

#Allele Types
MAJ = '0'
MIN = '1'
I_MAJ = 'Z'
I_MIN = 'W'

class SNPData:
	'''
	This class provides a datastructure for SNP datasets being imputed with NPUTE. Input
	files should be either CSV or TXT files with lines as SNPs (separated by crlf) and
	columns as samples.  The alleles in the SNPs may be separated by commas or nothing.
	The SNPs must be in order by chromosome position.  The SNPs must be ternary with
	a majority allele, a minority allele (not required), and a unknowns (not required).
	Majoriy and minority alleles can be represented by any characters, but an unknown
	must be a '?' character.
	'''

	def __init__(self, inFile, chromosome=None):
		'''
		Constructor reads in data file, generates mismatch vectors and upper
		triangular matrix extract indices, and reports stats of data.
		'''
		start = time.time()
		if chromosome==None:
			self.readInData(inFile)
		else:
			snps, nucs, self.chosen_snps_name_ls = self.readOneChromosomeData(inFile, chromosome)
			self.snps = array(snps)
			self.sdps = Set(snps)
			self.nucs = array(nucs)
			self.numSamps = len(self.snps[0])
		
		self.chromosome = chromosome
		
		self.genMismatchVectors()
		self.genExtractIndices()
		t = int(time.time()-start+0.5)
		print 'Number of Samples: %d' % self.numSamps
		print 'Number of SNPs: %d' % len(self.snps)
		print 'Number of SDPs: %d' % len(self.sdps)
		print 'Time to Process: %d m %d s\n' % (t/60,t%60)
	
	def readOneChromosomeData(self, input_fname, chromosome):
		"""
		05/07/08
			replace readInData()
		"""
		sys.stderr.write("Reading chromosome %s data from %s ..."%(chromosome, input_fname))
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix.read_data(input_fname, turn_into_integer=0)
		snps_name_ls = header[2:]
		
		no_of_rows = len(strain_acc_list)
		no_of_cols = len(snps_name_ls)
		snps = []
		nucs = []
		chosen_snps_name_ls = []
		for i in range(no_of_cols):
			snps_name = snps_name_ls[i]
			tmp_ls = snps_name.split('_')
			chr = tmp_ls[0]
			if chr!=chromosome:	#skip the SNPs from other chromosomes
				continue
			
			chosen_snps_name_ls.append(snps_name)
			
			#count the frequency of symbols
			one_snp_ls = []
			symbol2counts = {}
			for j in range(no_of_rows):
				symbol = data_matrix[j][i]
				one_snp_ls.append(symbol)
				if symbol != '0':	#don't care about NA
					if symbol not in symbol2counts:
						symbol2counts[symbol] = 0
					symbol2counts[symbol] += 1
			
			#construct a dictionary to map input symbols to MAJ, MIN or '?'
			symbol2MAJ_MIN = {'0':'?'}	#'NA' is always '?'
			symbols = symbol2counts.keys()
			if len(symbols) == 0:
				maj = ''
				min = ''
			elif len(symbols) == 1:
				symbol2MAJ_MIN[symbols[0]] = MAJ
				maj = symbols[0]
				min = ''
			elif len(symbols) ==2:
				maj, min = symbols
				if symbol2counts[maj]<symbol2counts[min]:
					min, maj = symbols	#reverse them
				symbol2MAJ_MIN[maj] = MAJ
				symbol2MAJ_MIN[min] = MIN
			elif len(symbols)>2:
				sys.stderr.write("Error: SNP %s (%s) has more than 2 alleles: %s.\n"%(i, snps_name_ls[i], repr(symbols)))
				sys.exit(2)
			
			map_func = lambda x: symbol2MAJ_MIN[x]
			one_snp_ls = map(map_func, one_snp_ls)
			
			snps.append(''.join(one_snp_ls))
			nucs += [(maj,min)]
		
		del data_matrix
		sys.stderr.write("Done.\n")
		return snps, nucs, chosen_snps_name_ls
	
	def readInData(self, inFile):
		'''
		Reads in data file and normalizes alleles.
		'''
		print "Reading in SNP data from '%s'..." % inFile,
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
		print "Generating pair-wise mismatch vectors...",

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

			
		print "Done"
		
	def incorporateChanges(self):
		'''
		Uses the changes dictionary to replace unknowns with imputed values.
		'''
		print "Incorporating imputed values into the SNP data...",
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
		print "Done"			

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
				nuc = self.allele2output_symbol(allele, self.nucs, i)
				out += ',' + nuc
			fi.write(out[1:] + '\n')
		print "Done"
	
	def allele2output_symbol(self, allele, nucs, i, NA_char='?'):
		if allele == MAJ:
			nuc = nucs[i,0]
		elif allele == MIN:
			nuc = nucs[i,1]
		elif allele == I_MAJ:
			nuc = nucs[i,0].lower()
		elif allele == I_MIN:
			nuc = nucs[i,1].lower()
		else:
			nuc = NA_char
		return nuc
	
	def outputOneChromosomeData(self, chosen_snps_name_ls, chromosome, snps, nucs, numSamps, output_fname):
		"""
		05/07/08
			a matrix with only snp name as header
		"""
		sys.stderr.write("Outputting imputed chromosome %s data to %s ..."%(chromosome, output_fname))
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(chosen_snps_name_ls)
		no_of_snps = len(chosen_snps_name_ls)
		for i in range(numSamps):
			data_row = []
			for j in range(no_of_snps):
				allele = snps[j][i]
				nuc = nuc = self.allele2output_symbol(allele, nucs, j, NA_char='0')
				data_row.append(nuc)
			writer.writerow(data_row)
		del writer
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
		print 'Generating indices for triangular matrix row extraction...',
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
		print 'Done'
	
	def extractRow(self,tM,rowNum):
		'''
		Extracts a row from an upper-triangular mismatch matrix for the dataset.  Sets value
		for imputed sample to inf so that it will not be chosen to impute itself during
		window tests.
		'''
		row = take(tM,self.extractIndices[rowNum])
		row[rowNum] = INF
		return row		
