#!/usr/bin/env python
"""
Usage: Output2010InCertainSNPs.py [OPTIONS] -p XXX -o OUTPUT_FNAME

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-e ...,	ecotype_table, 'stock20071008.ecotype'(default)
	-p ...,	ecotype 2 accession mapping table 'at.accession2ecotype_complete'/'at.ecotype2accession'
	-s ...,	sequence table, 'at.sequence'(default)
	-a ..., alignment table, 'at.alignment'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-o ...,	output_fname
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	#to output 2010 data in 250k SNPs
	Output2010InCertainSNPs.py -n snps_250k -o /tmp/data_2010_x_250k.tsv
	
	#to output 2010 data in 149SNP
	Output2010InCertainSNPs.py -o /tmp/data_2010_x_149SNP.tsv

Description:
	program to output 2010 data from db at with columns/SNPs from either 250k or 149SNP
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import getopt, numpy
from variation.src.common import nt2number, number2nt

class Output2010InCertainSNPs:
	"""
	2007-12-18
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='dbsnp', ecotype_table='ecotype',\
		accession2ecotype_table=None, alignment_table='at.alignment', sequence_table='at.sequence',\
		snp_locus_table='stock20071008.snps', output_fname=None,  debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		
		self.ecotype_table = ecotype_table
		self.accession2ecotype_table = accession2ecotype_table
		self.alignment_table = alignment_table
		self.sequence_table = sequence_table
		self.snp_locus_table = snp_locus_table
		self.output_fname = output_fname
		
		self.debug = int(debug)
		self.report = int(report)
	
	def get_alignment_id2positions_to_be_checked_ls(self, curs, alignment_table):
		"""
		2007-11-05
		2007-12-18
			copied from CmpAccession2Ecotype.py
			change alignment_id2start to alignment_id2chr_start_end
		"""
		sys.stderr.write("Getting alignment_id2positions_to_be_checked_ls ...")
		alignment_id2positions_to_be_checked_ls = {}
		alignment_id2chr_start_end = {}
		curs.execute("select id, chromosome, start, end, target from %s"%alignment_table)
		rows = curs.fetchall()
		for row in rows:
			alignment_id, chromosome, start, end, target_seq = row
			positions_to_be_checked_ls = []
			for i in range(len(target_seq)):
				if target_seq[i]!='-':
					positions_to_be_checked_ls.append(i)
			alignment_id2positions_to_be_checked_ls[alignment_id] = positions_to_be_checked_ls
			alignment_id2chr_start_end[alignment_id] = [chromosome, start, end]
		sys.stderr.write("Done.\n")
		return alignment_id2positions_to_be_checked_ls, alignment_id2chr_start_end
	
	def get_SNPpos_snpacc_ls(self, curs, snp_locus_table):
		"""
		2007-12-18
		"""
		sys.stderr.write("Getting SNPpos_snpacc_ls ...")
		SNPpos_snpacc_ls = []
		curs.execute("select id, snpid, chromosome, position from %s order by chromosome, position"%(snp_locus_table))
		rows = curs.fetchall()
		for row in rows:
			snpid, snp_acc, chromosome, position = row
			SNPpos_snpacc_ls.append([chromosome, position, snp_acc])
		sys.stderr.write("Done.\n")
		return SNPpos_snpacc_ls
	
	def setup_SNP_dstruc(self, SNPpos_snpacc_ls, alignment_id2chr_start_end):
		"""
		2007-12-18
		"""
		sys.stderr.write("Setting up SNP data structure ...")
		selected_SNPpos_snpacc_ls = []
		for chromosome, position, snp_acc in SNPpos_snpacc_ls:
			for alignment_id, chr_start_end in alignment_id2chr_start_end.iteritems():
				alignment_chromosome, start, end = chr_start_end
				if chromosome==alignment_chromosome and position>=start and position<=end:
					selected_SNPpos_snpacc_ls.append((chromosome, position, snp_acc))
					break	#there're multiple alignments for one SNPpos_snpacc, break out to avoid repetition
		
		selected_SNPpos_snpacc_ls.sort()	#in chromosomal order
		SNPpos2col_index = {}
		snp_acc_ls = []
		for i in range(len(selected_SNPpos_snpacc_ls)):
			chromosome, position, snp_acc = selected_SNPpos_snpacc_ls[i]
			SNPpos2col_index[(chromosome, position)] = i
			snp_acc_ls.append(snp_acc)
		sys.stderr.write("Done.\n")
		return SNPpos2col_index, snp_acc_ls 
	
	def setup_accession_ecotype_dstruc(self, curs, accession2ecotype_table, ecotype_table):
		"""
		2007-10-31
		2007-11-06
			ecotype_id = (ecotype_id, duplicate)
			add accession_id2ecotype_id_ls
		2007-12-18
			no use for calls_table, removed
			accession2ecotype_table is 'ecotype2accession', the one-to-one mapping table
		"""
		sys.stderr.write("Setting up accession or ecotype data structure ...")
		ecotype_id2accession_id = {}	#not accession2ecotype because multiple ecotype ids correspond to same accession_id
		ecotype_id2row_index = {}
		ecotype_id_ls = []
		ecotype_id2info_ls = {}
		accession_id2row_index = {}
		accession_id_ls = []
		nativename_ls = []
		curs.execute("select distinct ea.ecotype_id, e.nativename, e.stockparent, ea.accession_id from %s ea, %s e where e.id=ea.ecotype_id order by accession_id"%(accession2ecotype_table, ecotype_table))
		rows = curs.fetchall()
		for row in rows:
			ecotype_id, nativename, stockparent, accession_id = row
			accession_id2row_index[accession_id] = len(accession_id2row_index)
			
			ecotype_id2accession_id[ecotype_id] = accession_id
			ecotype_id2row_index[ecotype_id] = len(ecotype_id2row_index)
			ecotype_id2info_ls[ecotype_id] = [nativename, stockparent]
			ecotype_id_ls.append(ecotype_id)
			accession_id_ls.append(accession_id)
			nativename_ls.append(nativename)
		sys.stderr.write("Done.\n")
		return ecotype_id2accession_id, ecotype_id2row_index, ecotype_id2info_ls, ecotype_id_ls, accession_id2row_index, accession_id_ls, nativename_ls
	
	def get_accession_X_snp_matrix(self, curs, accession_id2row_index, SNPpos2col_index, sequence_table, alignment_table, alignment_id2positions_to_be_checked_ls):
		"""
		2007-11-05
			add alignment_id2positions_to_be_checked_ls to skip '-' columns in reference genome to match the real genome position
		2007-11-06
			try to differentiate 'missing'(which is 'N') and 'not sequenced' (which is not in the db) using the data_matrix_touched
			0 means not sequenced, 1 means sequenced
		2007-11-07
			add snp_index2alignment_id
		2007-12-18
			copied from CmpAccession2Ecotype.py
		"""
		sys.stderr.write("Getting accession_X_snp_matrix ...\n")
		data_matrix = numpy.zeros([len(accession_id2row_index), len(SNPpos2col_index)], numpy.integer)
		data_matrix_touched = numpy.zeros([len(accession_id2row_index), len(SNPpos2col_index)], numpy.integer)
		snp_index2alignment_id = {}
		curs.execute("select s.accession, s.alignment, s.bases, a.chromosome, a.start from %s s, %s a where s.alignment=a.id"%(sequence_table, alignment_table))
		rows = curs.fetchmany(5000)
		counter = 0
		while rows:
			for row in rows:
				accession_id, alignment_id, bases, chromosome, start_pos = row
				if accession_id in accession_id2row_index:
					positions_to_be_checked_ls = alignment_id2positions_to_be_checked_ls[alignment_id]
					for i in range(len(positions_to_be_checked_ls)):	#i is the real index after reference '-' is removed.
						inflated_index = positions_to_be_checked_ls[i]
						SNP_pos_key = (chromosome, start_pos+i)
						if SNP_pos_key in SNPpos2col_index:
							if self.debug:
								import pdb
								pdb.set_trace()
							SNP_index = SNPpos2col_index[SNP_pos_key]
							data_matrix[accession_id2row_index[accession_id], SNP_index] = nt2number[bases[inflated_index]]
							data_matrix_touched[accession_id2row_index[accession_id], SNP_index] = 1
							if SNP_index not in snp_index2alignment_id:
								snp_index2alignment_id[SNP_index] = alignment_id
				counter += 1
			rows = curs.fetchmany(5000)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return data_matrix, data_matrix_touched, snp_index2alignment_id
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=dbname,host=hostname)
		curs = conn.cursor()
		if self.debug:
			import pdb
			pdb.set_trace()

		alignment_id2positions_to_be_checked_ls, alignment_id2chr_start_end = self.get_alignment_id2positions_to_be_checked_ls(curs, self.alignment_table)
		SNPpos_snpacc_ls = self.get_SNPpos_snpacc_ls(curs, self.snp_locus_table)
		SNPpos2col_index, snp_acc_ls = self.setup_SNP_dstruc(SNPpos_snpacc_ls, alignment_id2chr_start_end)

		ecotype_id2accession_id, ecotype_id2row_index, ecotype_id2info_ls, ecotype_id_ls, accession_id2row_index, accession_id_ls, nativename_ls = self.setup_accession_ecotype_dstruc(curs, self.accession2ecotype_table, self.ecotype_table)
		accession_X_snp_matrix, accession_X_snp_matrix_touched, snp_index2alignment_id = self.get_accession_X_snp_matrix(curs, accession_id2row_index, SNPpos2col_index, self.sequence_table, self.alignment_table, alignment_id2positions_to_be_checked_ls)

		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header = ['ecotype_id', 'nativename'] + snp_acc_ls
		FilterStrainSNPMatrix_instance.write_data_matrix(accession_X_snp_matrix, self.output_fname, header, ecotype_id_ls, nativename_ls)



if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:e:p:s:a:n:o:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	ecotype_table = 'stock20071008.ecotype'
	accession2ecotype_table = 'at.ecotype2accession'
	sequence_table = 'at.sequence'
	alignment_table = 'at.alignment'
	snp_locus_table = 'stock20071008.snps'
	output_fname = None
	debug = 0
	report = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-e",):
			ecotype_table = arg
		elif opt in ("-p",):
			accession2ecotype_table = arg
		elif opt in ("-s",):
			sequence_table = arg
		elif opt in ("-a",):
			alignment_table = arg
		elif opt in ("-n",):
			snp_locus_table = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if ecotype_table and accession2ecotype_table and sequence_table and alignment_table and snp_locus_table and hostname and dbname and schema:
		instance = Output2010InCertainSNPs(hostname, dbname, schema, ecotype_table,\
			accession2ecotype_table, alignment_table, sequence_table,\
			snp_locus_table, output_fname, debug, report)
		
		instance.run()
	else:
		print __doc__
		sys.exit(2)
