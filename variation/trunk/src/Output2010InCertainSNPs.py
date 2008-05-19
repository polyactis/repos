#!/usr/bin/env python
"""
Usage: Output2010InCertainSNPs.py [OPTIONS] -p XXX -o OUTPUT_FNAME

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-e ...,	ecotype_table, 'stock20071008.ecotype'(default)
	-p ...,	ecotype 2 accession mapping table 'at.accession2ecotype_complete' or 'at.ecotype2accession' or 'at.accession2tg_ecotypeid'(default)
	-s ...,	sequence table, 'at.sequence'(default)
	-a ..., alignment table, 'at.alignment'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-o ...,	output_fname
	-y ...,	processing bits to control which processing step should be turned on.
		default is 0000. for what each bit stands, see Description.
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	#to output 2010 data in 250k SNPs with ecotype id
	Output2010InCertainSNPs.py -n snps_250k -o /tmp/data_2010_ecotype_id_x_250k_y00.tsv -r -y00
	
	#ditto except with accession id
	Output2010InCertainSNPs.py -n snps_250k -o /tmp/data_2010_accession_id_x_250k_y10.tsv -r -y10
	
	#to output 2010 data in 149SNP with ecotype id
	Output2010InCertainSNPs.py -o /tmp/data_2010_ecotype_id_x_149SNP_y00.tsv -r -y00
	
	#ditto except with accession id
	Output2010InCertainSNPs.py -o /tmp/data_2010_accession_id_x_149SNP_y10.tsv -r -y10
	
	#to output perlegen data in 250k SNPs with accession id
	Output2010InCertainSNPs.py -n snps_250k -o /tmp/data_perlegen_accession_id_x_250k_y11.tsv -y11
	
	#to output perlegen data in 250k SNPs with ecotype id
	Output2010InCertainSNPs.py -n snps_250k -o /tmp/data_perlegen_ecotype_id_x_250k_y01.tsv -y01
	
	#to output perlegen data in 149SNP with accession id
	Output2010InCertainSNPs.py -n snps -o /tmp/data_perlegen_accession_id_x_149SNP_y11.tsv -y11
	
	#to output perlegen data in its own SNP set with ecotype id
	Output2010InCertainSNPs.py -o ../data/perlegen/data_perlegen_ecotype_id_y0102.tsv -y0102
	
	#to output 2010 data in its own SNP set with ecotype id
	Output2010InCertainSNPs.py -o ../data/2010/data_2010_ecotype_id_y0002.tsv -r -y0002
	
Description:
	program to output 2010/Perlegen data from db at with columns/SNPs from either 250k or 149SNP
	
	for 2010 data, alignment.version=3, locus.offset=0
		select g.accession, l.chromosome, l.position, al.base from %s g, %s al, %s l, %s an\
				where g.allele=al.id and l.id=al.locus and l.alignment=an.id and l.offset=0 and an.version=3"%\
				(genotype_table, allele_table, locus_table, alignment_table)
	
	for perlegen data, select ecotype, chromosome, position, mbml98 from chip.snp_combined_may_9_06_no_van

	definition of each bit in processing_bits (0=off, 1=on), default is 0.
	1st bit. 0: row id of the output is ecotype id(default); 1: accession id instead.
	2nd bit. 0: 2010 data(default); 1: perlegen data
	3rd bit controls the 2nd column. duplicate id or strain name. 0: duplicate(default). 1: name
	4th bit controls what SNP set gets outputted.
		0: all SNPs from snp_locus_table(default); 1: intersection between snp_locus_table and the SNPs where there's data
		2: all SNPs from the datat type itself (according to 2nd bit).
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
from sets import Set

class row_dstruc:
	def __init__(self, row_index=None, name=None, info=None, snp_ls=None, snp_touched_ls=None):
		self.row_index = row_index
		self.name = name
		self.info = info
		self.snp_ls = snp_ls
		self.snp_touched_ls = snp_touched_ls

class Output2010InCertainSNPs(object):
	"""
	2007-12-18
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema='dbsnp', ecotype_table='ecotype',\
		accession2ecotype_table=None, alignment_table='at.alignment', sequence_table='at.sequence',\
		snp_locus_table='stock20071008.snps', output_fname=None, processing_bits='000', debug=0, report=0):
		"""
		2007-12-29
			add processing_bits
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		
		self.ecotype_table = ecotype_table
		self.accession2ecotype_table = accession2ecotype_table
		self.alignment_table = alignment_table
		self.sequence_table = sequence_table
		self.snp_locus_table = snp_locus_table
		self.output_fname = output_fname

		self.processing_bits = [0]*len(processing_bits)
		for i in range(len(processing_bits)):
			self.processing_bits[i] = int(processing_bits[i])
		self.debug = int(debug)
		self.report = int(report)
		
		self.data_type2data_table = {0:'at.locus',
									1:'chip.snp_combined_may_9_06_no_van'}
	
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
		2008-05-12
			not used. deprecated
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
		2008-05-12
			not used. deprecated
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
		2008-01-07
			default, data_matrix[:] = -2
		"""
		sys.stderr.write("Getting accession_X_snp_matrix ...\n")
		data_matrix = numpy.zeros([len(accession_id2row_index), len(SNPpos2col_index)], numpy.integer)
		data_matrix[:] = -2
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
	
	def setup_SNP_dstruc2(self, curs, snp_locus_table, cross_linking_table=None):
		"""
		2008-05-12
			snp_acc = '%s_%s'%(chromosome, position)	#to standardize between different SNP dataset
			correct a bug when SNPpos_set is empty
			cross_linking_table is optional.
		2008-02-07
			2nd version
		"""
		sys.stderr.write("Setting up SNP data structure ...")
		SNPpos_set = Set()
		if cross_linking_table:
			curs.execute("select distinct chromosome, position from %s"%cross_linking_table)
			rows = curs.fetchall()
			for row in rows:
				chromosome, position = row
				SNPpos_set.add((chromosome, position))
		snp_acc_ls = []
		SNPpos2col_index = {}
		curs.execute("select distinct chromosome, position from %s order by chromosome, position"%(snp_locus_table))
		rows = curs.fetchall()
		for row in rows:
			chromosome, position = row
			key_tuple = (chromosome, position)
			snp_acc = '%s_%s'%(chromosome, position)	#to standardize between different SNP dataset
			if not SNPpos_set or key_tuple in SNPpos_set:
				snp_acc_ls.append(snp_acc)
				SNPpos2col_index[key_tuple] = len(SNPpos2col_index)
		sys.stderr.write(" %s SNPs. Done.\n"%len(snp_acc_ls))
		return SNPpos2col_index, snp_acc_ls
	
	def setup_row_dstruc(self, curs, SNPpos2col_index, accession_id2name, genotype_table='at.genotype', \
						allele_table='at.allele', locus_table='at.locus', alignment_table='at.alignment', \
						perlegen_table='chip.snp_combined_may_9_06_no_van', data_type=0, ecotype_name2accession_id=None):
		"""
		2008-02-07
			get row_dstruc
			SNPpos2col_index makes sure that SNPs are from snp_locus_table
			data_type=0  -> 2010
			data_type=1  -> perlegen
		"""
		sys.stderr.write("Getting row data structure and data ...\n")
		row_id2dstruc = {}
		if data_type==0:
			#offset=0 and version=3
			sql_string_func = lambda key_tuple: "select g.accession, al.base from %s g, %s al, %s l, %s an where g.allele=al.id \
				and l.id=al.locus and l.alignment=an.id and l.offset=0 and an.version=3 and l.chromosome=%s and l.position=%s"%\
				(genotype_table, allele_table, locus_table, alignment_table, key_tuple[0], key_tuple[1])
			"""
			#offset=0
			sql_string_func = lambda key_tuple: "select g.accession, al.base from %s g, %s al, %s l, %s an where g.allele=al.id \
				and l.id=al.locus and l.alignment=an.id and l.offset=0 and l.chromosome=%s and l.position=%s"%\
				(genotype_table, allele_table, locus_table, alignment_table, key_tuple[0], key_tuple[1])
			"""
		elif data_type==1 and ecotype_name2accession_id:
			sql_string_func = lambda key_tuple: "select ecotype, mbml98 from %s where chromosome=%s and position=%s"%\
				(perlegen_table, key_tuple[0], key_tuple[1])
		else:
			sys.stderr("Unsupported data type: %s or no ecotype_name2accession_id specified.\n"%data_type)
			sys.exit(2)
		counter = 0
		snp_counter = 0
		for SNPpos, col_index in SNPpos2col_index.iteritems():
			sql_string = sql_string_func(SNPpos)
			curs.execute(sql_string)
			rows = curs.fetchall()
			for row in rows:
				row_id, base = row
				if data_type==1:	#in perlegen table, ecotype is a name. all but one('sha') matches accession.name.
					row_id = ecotype_name2accession_id[row_id]
				if row_id not in row_id2dstruc:
					row_id2dstruc[row_id] = row_dstruc()
					row_id2dstruc[row_id].name = accession_id2name.get(row_id)
					row_id2dstruc[row_id].snp_ls = [-2]*len(SNPpos2col_index)	#-2 means untouched.
					row_id2dstruc[row_id].snp_touched_ls = [0]*len(SNPpos2col_index)
				row_id2dstruc[row_id].snp_ls[col_index] = nt2number[base]
				row_id2dstruc[row_id].snp_touched_ls[col_index] = 1
				counter += 1
			snp_counter += 1
			if self.report and snp_counter%1000==0:
				sys.stderr.write("\t%s%s\t%s"%('\x08'*80, snp_counter, counter))
		if self.report:
			sys.stderr.write("\t%s%s\t%s"%('\x08'*80, snp_counter, counter))
		sys.stderr.write("Done.\n")
		return row_id2dstruc
	
	def transform_row_id2dstruc_2_matrix(self, row_id2dstruc):
		"""
		2008-02-07
		"""
		sys.stderr.write("Getting Row data structure and data ...")
		row_id_ls = row_id2dstruc.keys()
		row_id_ls.sort()
		row_name_ls = []
		no_of_rows = len(row_id_ls)
		data_matrix = numpy.zeros([no_of_rows, len(row_id2dstruc[row_id_ls[0]].snp_ls)], numpy.integer)
		for i in range(no_of_rows):
			row_id = row_id_ls[i]
			data_matrix[i] = row_id2dstruc[row_id].snp_ls
			row_name_ls.append(row_id2dstruc[row_id].name)
		sys.stderr.write("Done.\n")
		return row_id_ls, row_name_ls, data_matrix
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=dbname,host=hostname)
		curs = conn.cursor()
		if self.debug:
			import pdb
			pdb.set_trace()
		"""
		#2008-02-08 old way to get 2010 data is from raw alignments. didn't realize all SNPs are put into db.
		alignment_id2positions_to_be_checked_ls, alignment_id2chr_start_end = self.get_alignment_id2positions_to_be_checked_ls(curs, self.alignment_table)
		SNPpos_snpacc_ls = self.get_SNPpos_snpacc_ls(curs, self.snp_locus_table)
		SNPpos2col_index, snp_acc_ls = self.setup_SNP_dstruc(SNPpos_snpacc_ls, alignment_id2chr_start_end)

		ecotype_id2accession_id, ecotype_id2row_index, ecotype_id2info_ls, ecotype_id_ls, accession_id2row_index, accession_id_ls, nativename_ls = self.setup_accession_ecotype_dstruc(curs, self.accession2ecotype_table, self.ecotype_table)
		accession_X_snp_matrix, accession_X_snp_matrix_touched, snp_index2alignment_id = self.get_accession_X_snp_matrix(curs, accession_id2row_index, SNPpos2col_index, self.sequence_table, self.alignment_table, alignment_id2positions_to_be_checked_ls)
		"""
		if self.processing_bits[3]==0:
			SNPpos2col_index, snp_acc_ls = self.setup_SNP_dstruc2(curs, self.snp_locus_table)
		elif self.processing_bits[3]==1:
			SNPpos2col_index, snp_acc_ls = self.setup_SNP_dstruc2(curs, self.snp_locus_table, cross_linking_table=self.data_type2data_table[self.processing_bits[1]])
		elif self.processing_bits[3]==2:
			SNPpos2col_index, snp_acc_ls = self.setup_SNP_dstruc2(curs, self.data_type2data_table[self.processing_bits[1]])
		else:
			sys.stderr.write("Error: unsupported 3rd bit in processing_bits %s.\n"%self.processing_bits[3])
			sys.exit(3)
		from variation.src.common import get_accession_id2name
		accession_id2name = get_accession_id2name(curs)
		if self.processing_bits[1]==0:
			row_id2dstruc = self.setup_row_dstruc(curs, SNPpos2col_index, accession_id2name)
		elif self.processing_bits[1]==1:
			from variation.src.common import map_perlegen_ecotype_name2accession_id
			ecotype_name2accession_id = map_perlegen_ecotype_name2accession_id(curs)
			row_id2dstruc = self.setup_row_dstruc(curs, SNPpos2col_index, accession_id2name, data_type=self.processing_bits[1], ecotype_name2accession_id=ecotype_name2accession_id)
		else:
			sys.stderr("Unsupported data type: %s or no ecotype_name2accession_id specified.\n"%self.processing_bits[1])
			sys.exit(2)
		accession_id_ls, accession_name_ls, data_matrix = self.transform_row_id2dstruc_2_matrix(row_id2dstruc)
		
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		
		#2008-02-08 which type of row id/1st column
		if self.processing_bits[0]==0:
			from variation.src.common import map_accession_id2ecotype_id
			accession_id2ecotype_id = map_accession_id2ecotype_id(curs, accession2ecotype_table=self.accession2ecotype_table)
			ecotype_id_ls = []
			rows_to_be_tossed_out=Set()
			for i in range(len(accession_id_ls)):
				ecotype_id = accession_id2ecotype_id.get(accession_id_ls[i])
				if not ecotype_id:	#mapping failed
					rows_to_be_tossed_out.add(i)
				ecotype_id_ls.append(ecotype_id)
			strain_acc_list = ecotype_id_ls
			header = ['ecotype_id']	#1st column in the header
		else:
			rows_to_be_tossed_out=Set()
			strain_acc_list = accession_id_ls
			header = ['accession_id']
		#2008-02-08 which type of 2nd column
		if self.processing_bits[2]==0:
			category_list = [1]*len(accession_name_ls)
			header.append('duplicate')	#2nd column in the header
		elif self.processing_bits[2]==1:
			category_list = accession_name_ls
			header.append('accession_name')
		else:
			category_list = accession_name_ls
			header.append('accession_name')
		
		header += snp_acc_ls
		FilterStrainSNPMatrix_instance.write_data_matrix(data_matrix, self.output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=rows_to_be_tossed_out)



if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:e:p:s:a:n:o:y:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	ecotype_table = 'stock20071008.ecotype'
	accession2ecotype_table = 'at.accession2tg_ecotypeid'
	sequence_table = 'at.sequence'
	alignment_table = 'at.alignment'
	snp_locus_table = 'stock20071008.snps'
	output_fname = None
	processing_bits = '0000'
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
		elif opt in ("-y",):
			processing_bits = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if ecotype_table and accession2ecotype_table and sequence_table and alignment_table and hostname and dbname and schema:
		instance = Output2010InCertainSNPs(hostname, dbname, schema, ecotype_table,\
			accession2ecotype_table, alignment_table, sequence_table,\
			snp_locus_table, output_fname, processing_bits, debug, report)
		
		instance.run()
	else:
		print __doc__
		sys.exit(2)
