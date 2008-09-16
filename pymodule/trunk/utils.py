import os, sys, csv

def dict_map(dict, ls, type=1):
	"""
	2008-04-03 copied from annot.codense.common
	10-13-05
		add type 2 to return item itself if mapping is not available
	2006-12-21
		add type 3 to extract a smaller map
	2007-05-14
		bug, "if value" could miss 0
	"""
	if type==3:
		new_list = {}	#it's a dictionary
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list[item] = value
	else:
		new_list = []
		for item in ls:
			value = dict.get(item)
			if value is not None:
				new_list.append(value)
			elif type==2:
				new_list.append(item)
	
	return new_list

class PassingData(object):
	"""
	05/09/08
		a class to hold any data structure
	"""
	def __init__(self, **keywords):
		"""
		2008-5-12
			add keyword handling
		"""
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)

def importNumericArray():
	"""
	2008-07-09
		numarray doesn't have int128
	2008-05-18
		give same numpy types (int, int8 ...) to other numeric modules
	2008-05-18
		add "import array as num"
		should put ImportError in except. but whatever
	2008-05-11
		import whatever available array module
	"""
	try:
		import numpy as num
	except:
		numpy_type2other_ls = ['int', 'int8', 'int16', 'int32', 'int64']
		try:
			import numarray as num
		except:
			import Numeric as num
		for numpy_type in numpy_type2other_ls:	#make sure it has same type names
			numpy_type_in_other = numpy_type[0].upper() + numpy_type[1:]
			setattr(num, numpy_type, getattr(num, numpy_type_in_other))
	return num

def figureOutDelimiter(input_fname, report=0, delimiter_choice_ls = ['\t', ',']):
	"""
	2008-08-28
		nothing weird on hpc-cmb. it's a bug in other code.
		back to 'return None' if input_fname escapes all condition checking.
	2008-08-28
		try 'open(input_fname)' anyway if input_fname escapes all condition checking.
		something weird happened during a mpi job on hpc-cmb. the file is there. but escape the first condition.
	2008-05-25
		now 3 possible types of input_fname
		1. a file name (path)
		2. input_fname is a file object
		3. input_fname is input data, string
		
		for a file object or input file name:
		it could be binary file which doesn't have readline(). have to use this dumb approach due to '\n' might mess up sniff()
	2008-05-21
		csv.Sniffer is handy, use it figure out csv.Sniffer instead.
	2008-05-12
		try tab first
	"""
	if report:
		import sys
		sys.stderr.write("Figuring out delimiter for %s ..."%input_fname)
	cs = csv.Sniffer()
	if isinstance(input_fname, str) and os.path.isfile(input_fname):
		inf = open(input_fname)
	elif isinstance(input_fname, file):	#could be a file object
		inf = input_fname
	elif isinstance(input_fname, str) and not os.path.isfile(input_fname):	#it's the input
		import StringIO
		inf = StringIO.StringIO(input_fname)
	else:
		import sys
		sys.stderr.write("Error: %s is neither a file name nor a file object. But try 'open' anyway.\n"%input_fname)
		return None
	if getattr(inf, 'readline', None) is not None:	
		line = inf.readline()
		delimiter_chosen = cs.sniff(line).delimiter
	else:
		line = inf.read(200)	##binary file doesn't have readline(). have to use this dumb approach due to '\n' might mess up sniff()
		delimiter_chosen = None
		for delimiter in delimiter_choice_ls:
			delimiter_count = line.count(delimiter)
			if delimiter_count>0:
				delimiter_chosen = delimiter
				break
	del inf
	if report:
		sys.stderr.write("Done.\n")
	return delimiter_chosen

def get_gene_symbol2gene_id_set(curs, tax_id, table='genome.gene_symbol2id', upper_case_gene_symbol=0):
	"""
	2008-07-10 derived from annot.bin.codense.common.get_gene_symbol2gene_id()
	"""
	sys.stderr.write("Getting gene_symbol2gene_id_set...")
	gene_symbol2gene_id_set = {}
	from sets import Set
	curs.execute("select gene_id, gene_symbol from %s where tax_id=%s"%(table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		gene_id, gene_symbol = row
		if upper_case_gene_symbol:
			gene_symbol = gene_symbol.upper()
		if gene_symbol not in gene_symbol2gene_id_set:
			gene_symbol2gene_id_set[gene_symbol] = Set()
		gene_symbol2gene_id_set[gene_symbol].add(gene_id)
	sys.stderr.write(" %s entries. Done.\n"%len(gene_symbol2gene_id_set))
	return gene_symbol2gene_id_set

class FigureOutTaxID(object):
	__doc__ = "2008-07-29 class to figure out tax_id using postgres database taxonomy schema"
	option_default_dict = {('hostname', 1, ): ['localhost', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['graphdb', 'd', 1, 'database name', ],\
							('schema', 1, ): ['taxonomy', 'k', 1, 'database schema name', ],\
							('db_user', 0, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 0, ): [None, 'p', 1, 'database password', ],\
							}
	def __init__(self,  **keywords):
		"""
		2008-07-29
		"""
		from ProcessOptions import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)		
	
	def curs(self):
		from db import db_connect
		conn, curs =  db_connect(self.hostname, self.dbname, self.schema, user=self.db_user, password=self.db_passwd)
		return curs
	
	curs = property(curs)
	
	def scientific_name2tax_id(self):
		scientific_name2tax_id = {}
		curs = self.curs
		curs.execute("SELECT n.name_txt, n.tax_id FROM taxonomy.names n, taxonomy.nodes o where n.name_class='scientific name' and n.tax_id=o.tax_id and o.rank='species'")
		rows = curs.fetchall()
		for row in rows:
			scientific_name, tax_id = row
			scientific_name2tax_id[scientific_name] = tax_id
		return scientific_name2tax_id
	
	scientific_name2tax_id = property(scientific_name2tax_id)
	
	def returnTaXIDGivenScientificName(self, scientific_name):
		return self.scientific_name2tax_id.get(scientific_name)
	
	def returnTaxIDGivenSentence(self, sentence):
		"""
		2008-07-29
		"""
		tax_id_to_return = None
		for scientific_name, tax_id in self.scientific_name2tax_id.iteritems():
			if sentence.find(scientific_name)>=0:
				tax_id_to_return = tax_id
				break
		return tax_id_to_return

def getColName2IndexFromHeader(header):
	"""
	2008-09-16
		convenient function to read input files with flexible column order.
		One variable doesn't have to be in the same column in different files, as far as the name is same.
	"""
	col_name2index = {}
	for i in range(len(header)):
		column_name = header[i]
		col_name2index[column_name] = i
	return col_name2index
	
if __name__ == '__main__':
	FigureOutTaxID_ins = FigureOutTaxID()
	print FigureOutTaxID_ins.returnTaxIDGivenSentence('>gi|172045488|ref|NW_001867254.1| Physcomitrella patens subsp. patens PHYPAscaffold_10696, whole genome shotgun sequence')
