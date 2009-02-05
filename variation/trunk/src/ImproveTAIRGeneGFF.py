#!/usr/bin/env python
"""
Examples:
	ImproveTAIRGeneGFF.py -i /Users/dazhe/TAIR8/TAIR8_GFF3_genes.gff -o /tmp/TAIR8_GFF3_genes_more.gff -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf
		
Description:
	2009-2-4
		program to insert more info into the gene GFF file from TAIR:
			if the last column of a non-chromosome line has 'ID' in it, find corresponding gene id using its value
			add gene_symbol as 'Alias' and description as 'description'.
			
			Put value of gene_symbol under 'Alias', rather than 'gene_symbol' because the former is built into the gbrowse search engine.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/')))
import csv, string
from pymodule.utils import getGeneIDSetGivenAccVer, figureOutDelimiter
from DrawSNPRegion import DrawSNPRegion

class ImproveTAIRGeneGFF(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
						('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
						('dbname', 1, ): ['stock_250k', 'd', 1, '', ],\
						('db_user', 1, ): [None, 'u', 1, 'database username', ],\
						('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
						('input_fname', 1, ): [None, 'i', 1, 'input file'],\
						('output_fname', 1, ): [None, 'o', 1, 'Output Filename'],\
						("gene_annotation_picklef", 0, ): [None, 'j', 1, 'given the option, If the file does not exist yet, store a pickled gene_annotation into it. If the file exists, load gene_annotation out of it.'],\
						('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
						('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2009-2-4
		"""
		from pymodule import ProcessOptions
		self.ad=ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	gene_desc_names = ['gene_symbol', 'type_of_gene', 'protein_label', 'protein_comment', 'protein_text']
	
	def improveTAIRGeneGFF(self, input_fname, gene_symbol2gene_id_set, gene_annotation, output_fname):
		"""
		2009-2-5
			apply the improvement to any non-chromosome lines with 'ID' entry
			escape ';' by '%3B', which is regarded as a separator for every "name=value"
			escape ',' by '%2C', which is regarded as a separator for every "value"
			esacpe gene_symbol/Alias whose value matches individual chromosome (like gene 'CHR5' = 'Gene CHR5')
		2009-2-4
			if the last column has 'ID' in it, find corresponding gene id using its value and add gene_symbol and description
		"""
		sys.stderr.write("Improving TAIR Gene GFF with symbols and descriptions ...\n")
		import re
		p_ID_acc_ver = re.compile(r'ID=(\w+)\.(\d+);')
		p_ID_acc = re.compile(r'ID=(\w+);')
		p_ID_protein_acc = re.compile(r'ID=(\w+)\.(\d+)-Protein;')
		p_chr_name = re.compile(r'CHR\d+$')	#to esacpe gene_symbol/Alias whose value matches individual chromosome
		delimiter = figureOutDelimiter(input_fname)
		reader = csv.reader(open(input_fname), delimiter=delimiter)
		writer = csv.writer(open(output_fname,'w'), delimiter=delimiter, lineterminator='\n')	#lineterminator is important. GFF3Loader would break down if it's dos terminator('\r\n').
		counter = 0
		success_counter = 0
		for row in reader:
			last_col = row[-1]
			tair_id = None
			if p_ID_acc_ver.search(last_col):
				tair_id, version = p_ID_acc_ver.search(last_col).groups()
			if p_ID_acc.search(last_col):
				tair_id, = p_ID_acc.search(last_col).groups()
			if p_ID_protein_acc.search(last_col):
				tair_id, version = p_ID_protein_acc.search(last_col).groups()
			counter += 1
			if tair_id is not None and row[2]!='chromosome':
				gene_id_set = getGeneIDSetGivenAccVer(tair_id, gene_symbol2gene_id_set)
				gene_id = None
				
				if gene_id_set==None:
					sys.stderr.write("Linking to gene id failed for %s. No such gene_symbol, %s, in gene_symbol2gene_id_set.\n"%(last_col, tair_id))
				elif len(gene_id_set)==1:
					gene_id = list(gene_id_set)[0]
					success_counter += 1
				elif len(gene_id_set)>1:
					sys.stderr.write("Too many gene_ids: %s, %s.\n"%(tair_id, gene_id_set))
				elif len(gene_id_set)==0:
					sys.stderr.write("Linking to gene id failed for %s. There is gene_symbol, %s, in gene_symbol2gene_id_set but it's empty.\n"%(last_col, tair_id))
				else:
					sys.stderr.write("not supposed to happen: original_name=%s, gene_symbol=%s, gene_id_set=%s\n."%(last_col, tair_id, gene_id_set))
				if gene_id is not None:
					gene_model =  gene_annotation.gene_id2model.get(gene_id)
					if gene_model is not None:
						gene_commentary = gene_model.gene_commentaries[0]
						gene_desc_ls = DrawSNPRegion.returnGeneDescLs(self.gene_desc_names, gene_model, gene_commentary, cutoff_length=600,\
																	replaceNoneElemWithEmptyStr=1)
						local_gene_desc_names = map(string.upper, self.gene_desc_names)
						description = ',  '.join([': '.join(entry) for entry in zip(local_gene_desc_names, gene_desc_ls)])
						description = description.replace(';', '%3B')	#escape ';', which is regarded as a separator for every "name=value"
						description = description.replace(',', '%2C')	#escape ',', which is regarded as a separator for every "value"
						
						if last_col[-1]!=';':	#no ; delimiter at the end, append one
							last_col += ';'
						gene_symbol = gene_model.gene_symbol
						gene_symbol = gene_symbol.replace(';', '%3B')
						gene_symbol = gene_symbol.replace(',', '%2C')
						if p_chr_name.match(gene_symbol):	#match the chromosome name, change
							gene_symbol = 'Gene %s'%gene_symbol
						last_col += 'Alias=%s;'%gene_symbol
						last_col += 'description=%s'%description
						row[-1] = last_col
			if last_col[-1]==';':
				last_col = last_col[:-1]
				row[-1] = last_col
			if counter%5000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, success_counter, counter))
			writer.writerow(row)
		
		sys.stderr.write("%s/%s Done.\n"%(success_counter, counter))
		
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		import MySQLdb
		mysql_conn = MySQLdb.connect(db=self.dbname, host='banyan.usc.edu', user = self.db_user, passwd = self.db_passwd)
		mysql_curs = mysql_conn.cursor()
		from pymodule import get_gene_symbol2gene_id_set
		gene_symbol2gene_id_set = get_gene_symbol2gene_id_set(mysql_curs, 3702, table='genome.gene_symbol2id', upper_case_gene_symbol=1)	#3702 is At's tax id
		
		from variation.src.DrawSNPRegion import DrawSNPRegion
		DrawSNPRegion_ins = DrawSNPRegion(db_user=self.db_user, db_passwd=self.db_passwd, hostname=self.hostname, database=self.dbname,\
									input_fname='/tmp/dumb', output_dir='/tmp', debug=0)	#input_fname and output_dir are just random stuff
		gene_annotation = DrawSNPRegion_ins.dealWithGeneAnnotation(self.gene_annotation_picklef)
		self.improveTAIRGeneGFF(self.input_fname, gene_symbol2gene_id_set, gene_annotation, self.output_fname)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ImproveTAIRGeneGFF
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()