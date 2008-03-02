#!/usr/bin/env python
"""
Usage: CmpTina2010VsMy2010In250kSNPs.py [OPTIONS] -i XXX -j XXX

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)
	-e ...,	ecotype_table, 'stock20071008.ecotype'(default)
	-p ...,	ecotype 2 accession mapping table 'at.accession2ecotype_complete'/'at.ecotype2accession'
	-s ...,	sequence table, 'at.sequence'(default)
	-a ..., alignment table, 'at.alignment'(default)
	-n ...,	snp_locus_table, 'snp_locus'(default)
	-l ...,	calls table, 'stock20071008.calls'(default)
	-i ...,	input_fname1, 149SNP data filename
	-j ...,	input_fname2, 2010 data filename
	-o ...,	latex_output_fname which stores both the summary and detailed difference tables
		specify this latex_output_fname if you want output
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	CmpTina2010VsMy2010In250kSNPs.py -i /Network/Data/250k/calls/Tina_120607/250K_PERL_2010.csv -j ~/script/variation/data/2010/data_2010_x_250k.tsv -o ~/script/variation/genotyping/149snp/analysis/250K_PERL_2010Vs2010_x_250k_diff_matrix.tex
	
Description:
	program to compare the 149snp calls based on the common strains inn 2010 pcr data and Justin's sequenom data.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from variation.src.QualityControl import QualityControl
from variation.src.common import number2nt, nt2number
import getopt

class CmpTina2010VsMy2010In250kSNPs(QualityControl):
	"""
	2007-12-20
		compare Tina's 2010 data with mine
	"""
	def __init__(self, curs='', input_fname1='', input_fname2='', latex_output_fname='', debug=0, report=0):
		QualityControl.__init__(self)
		self.curs = curs
		self.input_fname1 = input_fname1
		self.input_fname2 = input_fname2
		self.latex_output_fname = latex_output_fname
		self.debug = int(debug)
		self.report = int(report)
	
	def get_col_matching_dstruc(self, header1, header2):
		"""
		2007-12-20
			entries in header1 look like '1_5481983'
			entries in header2 look like '1_5481983_T_G'
		"""
		sys.stderr.write("Getting col matching dstruc ...\n")
		snpid2col_index1 = {}
		for i in range(2, len(header1)):
			snpid = header1[i]
			snpid2col_index1[snpid] = i-2
		snpid2col_index2 = {}
		for i in range(2, len(header2)):
			snpid = header2[i][:-4]	#remove _T_G
			snpid2col_index2[snpid] = i-2
		
		col_id12col_id2 = {}
		for snpid in snpid2col_index1:
			if snpid in snpid2col_index2:
				col_id12col_id2[snpid] = snpid
		sys.stderr.write("Done.\n")
		return snpid2col_index1, snpid2col_index2, col_id12col_id2
	
	def get_row_matching_dstruc(self, curs, strain_acc_list1, strain_acc_list2, ecotype_table='ecotype', batch_ecotype_table='batch_ecotype'):
		"""
		2008-02-07
			add curs as argument
		2007-12-20
		2008-01-18
			Tina changed the naming format in her new output: /Network/Data/250k/calls/Tina_011708/250K_PERL_2010.csv
		"""
		sys.stderr.write("Getting row matching dstruc ...\n")
		strain_acc2row_index1 = {}
		for i in range(len(strain_acc_list1)):
			strain_acc = strain_acc_list1[i]
			strain_acc2row_index1[strain_acc] = i
		
		strain_acc2row_index2 = {}
		for i in range(len(strain_acc_list2)):
			strain_acc = strain_acc_list2[i]
			ecotypeid = int(strain_acc)
			strain_acc2row_index2[ecotypeid] = i
		
		row_id12row_id2 = {}
		predefined_strain_name2ecotypeid = {
			'Bay0B':8260,
			'BAY0A':8260,
			'Bor4B':8268,
			'BOR4A':8268,
			'BroA':8269,
			'BR0A':8269,
			'BUR0A':8272,
			'BurOB':8272,
			'C24A':8273,
			'C24B':8273,
			'Cibc17':8276,
			'ColOA':8279,
			'COL0A':8279,
			'Col0B':8279,
			'CS2249':8429,	#'CS2249' should be 'CS22491' and in ecotype, its nativename is 'N13'
			'CviOA':8281,
			'CVI0A':8281,
			'CviOB':8281,
			'CVI-0A':8281,
			'El-2':8289,
			'El2':8289,
			'EST1A':8291,
			'Est1B':8291,
			'Fab-2':8292,
			'Fab-4':8293,
			'FeiOB':8294,
			'FEI0A':8294,
			'Got22':8298,
			'Got7A':8299,
			'GOT7B':8299,
			'HR10':8308,
			'Kas1':8315,
			'Kas-1':8315,	#'Kas-1' mapped to 'Kas-2' in 2010
			'Kno10':8317,
			'Kno18':8318,
			'Kz_9':8322,
			'Ler1A':8324,
			'LER1A':8324,
			'Ler1B':8324,
			'Lov_1':6043,
			'LOV5A':6046,
			'Lov5B':6046,
			'LY1':8279,	#Col-0
			'LY2':8279,	#Col-0
			'LY3':8324,	#Ler-1
			'LY4':8324,	#Ler-1
			'LY5':8400,	#Van-0
			'LY6':8400,	#Van-0
			'Ms_0':8340,
			'Nafa8B':8346,
			'Nfa10':8345,
			'NFA8A':8346,
			'Nok_3':8347,
			'Omo21':8349,
			'Omo2-1':8349,
			'Omo2_3':8350,
			'RmxA02':8370,
			'RmxA180':8371,
			'RmxA180_2':8371,
			'Rrs7B':8373,
			'Rrs10B':8372,
			'Pna_17':8359,
			'Pna_10':8358,
			'Pu2_7':8362,
			'Pu223':8361,
			'SahB':8248,
			'Se-0_2':8379,
			'Se0_2':8379,
			'ShaA':8248,
			'SHADARAA':8248,
			'Spr1_6':8383,
			'Sq_1':8384,
			'Uod_1':8398,
			'TAMM2A':8390,
			'Tamm2B':8390,
			'Tamm_27':8391,
			'Ts1B':8392,
			'TS1A':8392,
			'Tsu1B':8394,
			'VanOA':8400,
			'VAN0A':8400,
			'VanOB':8400,
			'Var2-1':8401,
			'Var2_6':8402,
			'Wa_1':8403
			}
		#2007-12-20 predefined_name2ecotypeid is for some names with swedish letters in. impossible to match the nativename in the db. made by eyeballing.
		#import re
		#p_strain_name = re.compile('(.*)_[P2][E0][R1][L0]')
		for strain_acc in strain_acc2row_index1:
			#p_strain_name_search = p_strain_name.search(strain_acc)
			#if p_strain_name_search:
			if strain_acc[-5]=='_':	#2008-01-17, Tina's old naming is like 'Omo2-1_2010'. New naming is like 'Omo212010'.
				strain_name = strain_acc[:-5]
			else:
				strain_name = strain_acc[:-4]
			curs.execute("select b.ecotypeid from %s b, %s e where e.id=b.ecotypeid and b.batchid=2 and e.nativename='%s'"%(batch_ecotype_table, ecotype_table, strain_name))
			rows = curs.fetchall()
			if not rows:
				strain_name_all_dash = strain_name.replace('_','-')	#backup strain_name
				curs.execute("select b.ecotypeid from %s b, %s e where e.id=b.ecotypeid and b.batchid=2 and e.nativename='%s'"%(batch_ecotype_table, ecotype_table, strain_name_all_dash))
				rows = curs.fetchall()
			if not rows and strain_name.find('-')==-1 and strain_name.find('_')==-1:	#2008-01-17 try one more type of strain_name if the strain_name doesn't have '-' and '_'
				strain_name_tmp = strain_name[:-1] + '-' + strain_name[-1]
				curs.execute("select b.ecotypeid from %s b, %s e where e.id=b.ecotypeid and b.batchid=2 and e.nativename='%s'"%(batch_ecotype_table, ecotype_table, strain_name_tmp))
				rows = curs.fetchall()
			if rows:
				ecotypeid = rows[0][0]
				if ecotypeid in strain_acc2row_index2:
					row_id12row_id2[strain_acc] = ecotypeid
				else:
					sys.stderr.write("\tMatching Failure: %s.\n"%(strain_acc))
			else:
				if strain_name in predefined_strain_name2ecotypeid:
					ecotypeid = predefined_strain_name2ecotypeid[strain_name]
					row_id12row_id2[strain_acc] = ecotypeid
				else:
					sys.stderr.write("\tMatching Failure: %s.\n"%(strain_acc))
		#2008-01-17 temporary, swap the mapping
		#row_id12row_id2['Omo2-1_2010'] = 8350
		#row_id12row_id2['Omo2_3_2010'] = 8349
		sys.stderr.write("Done.\n")
		return strain_acc2row_index1, strain_acc2row_index2, row_id12row_id2
	
	def readTina2010In250kSNPs(self, input_fname):
		"""
		2007-12-20
			/Network/Data/250k/calls/Tina_120607/250K_PERL_2010.csv
			the format is SNP x Strain (need transposation)
		"""
		sys.stderr.write("Reading data from %s ..."%os.path.basename(input_fname))
		import csv
		reader = csv.reader(open(input_fname), delimiter=',')
		strain_acc_list = reader.next()
		strain_acc_list = strain_acc_list[5:]
		category_list = strain_acc_list
		header = ['strain', 'strain']
		data_matrix = []
		for row in reader:
			chromosome = int(row[3])
			position = int(row[4])
			snpid = '%s_%s'%(chromosome, position)
			header.append(snpid)
			data_row = []
			for nt in row[5:]:
				if nt=='?':
					data_row.append(0)
				else:
					data_row.append(nt2number[nt])
			data_matrix.append(data_row)
		del reader
		sys.stderr.write("Done.\n")
		sys.stderr.write("Transposing the matrix ...")
		#transpose the matrix
		import numpy
		data_matrix = numpy.array(data_matrix)
		data_matrix = numpy.transpose(data_matrix)
		sys.stderr.write("Done.\n")
		return header, strain_acc_list, category_list, data_matrix
	
	def load_dstruc(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		QualityControl.load_dstruc(self)
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		self.header1, self.strain_acc_list1, self.category_list1, self.data_matrix1 = self.readTina2010In250kSNPs(self.input_fname1)
		self.header2, self.strain_acc_list2, self.category_list2, self.data_matrix2 = FilterStrainSNPMatrix_instance.read_data(self.input_fname2)
	 	
		self.col_id2col_index1, self.col_id2col_index2, self.col_id12col_id2 = self.get_col_matching_dstruc(self.header1, self.header2)
		self.row_id2row_index1, self.row_id2row_index2, self.row_id12row_id2 = self.get_row_matching_dstruc(self.curs, self.strain_acc_list1, self.strain_acc_list2)
	
	def get_row_id2info(self, row_id_ls, curs, calls_250k_duplicate_comment_table='calls_250k_duplicate_comment', ecotype_table='ecotype'):
		"""
		#2007-12-16
		"""
		row_id2info = {}
		for row_id in row_id_ls:
			row_id2info[row_id] = row_id
		return row_id2info


"""
Omo2\_3\_2010&1\_11160050&A&8350&1\_11160050&G\\

select id, chromosome, start, end, target from alignment where id=770;
| id  | chromosome | start    | end      | target                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
+-----+------------+----------+----------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 770 |          1 | 11159587 | 11160124 | CCATTGTTGTGTTATCTGAAAAAATTACTGTGGTACCAAGACGAAGAGGAGATA----CCTGAAAATTCCCTCAATTCGCAGTATTTTATCCCCTGTTTCGAGGAGACCACGTAAGGTGCTAATAATTGGAAGACCAGCTCCAACGGTGGCTTCATAAAAGTAATGTGTGTACGATTTTCGTTGAAGATCTCTGATCTTTAGATACTGCATTTACAGGAGTGTTACCAATTAGTCTCTCTGTTCATAAAATCACAGAATCTTAATGGATATATA-TCAAGATGAACATATATGTAGCTCCCGATCCAACAGTCATTTAAACAAATATCAAACACTTTACCTGATCAAGTGGTCCAGAGTTAGCCTTTTTGTTCGGAGTGACCACATGAATTCCTCGTAGCAACCAGTCGTAGTAACAGCTAGCGATGTCAGCATCGGCTGTACAATCAACCATAACAGAGTTTGGGATGAAATGATTTCCCTTCACATATTGGGTGAACTTCTCCATGTCAGCTTTTTCTCCCTCTTCTTTCATAAGCTCTCT | 

5 insertions

11160050-11159587 = 463

real index for 11160050 is 463+5

#check Omo2-3 (accession=22)

select accession, alignment, bases from sequence where accession=22 and alignment=770;

|        22 |       770 | NNNNNNNNNNNNNNNNNNNNNNNNTTACTGTGGTACCAAGACGAAGAGGAGATA----CCTGAAAATTCCCTCAATTCGCAGTATTTTATCCCCTGTTTCGAGGAGACCACGTAAGGTGCTAATAATTGGAAGACCAGCTCCAACGGTGGCTTCATAAAAGTAATGTGTGTACGATTTTCGTTGAAGATCTCTGATCTTTAGATACTGCATTTACAGGAGTGTTACCAATTAGTCTCTCTGTTCATAAAATCACAGAATCTTAATGGATATATACTCAAGATGAACATATATGTAGCTCCCGATCCAACAGTCATTTAAACAAATATCAAACACTTTACCTGATCAAGTGGTCCAGAGTTAGCCTTTTTGTTCGGAGTGACCACATGAATTCCTCGTAGCAACCAGTCGTAGTAACAGCTAGCGATGTCAGCATCGGCTGTACAATCAACCATAACAGAGTTTGGGATGAAATGATTTCCCTTCACATATTGGGTGAACTTCTCCATGTCAGCTTTTTCTCCCTCTTCTTTCATAAGCTCTCT | 

s = 'NNNNNNNNNNNNNNNNNNNNNNNNTTACTGTGGTACCAAGACGAAGAGGAGATA----CCTGAAAATTCCCTCAATTCGCAGTATTTTATCCCCTGTTTCGAGGAGACCACGTAAGGTGCTAATAATTGGAAGACCAGCTCCAACGGTGGCTTCATAAAAGTAATGTGTGTACGATTTTCGTTGAAGATCTCTGATCTTTAGATACTGCATTTACAGGAGTGTTACCAATTAGTCTCTCTGTTCATAAAATCACAGAATCTTAATGGATATATACTCAAGATGAACATATATGTAGCTCCCGATCCAACAGTCATTTAAACAAATATCAAACACTTTACCTGATCAAGTGGTCCAGAGTTAGCCTTTTTGTTCGGAGTGACCACATGAATTCCTCGTAGCAACCAGTCGTAGTAACAGCTAGCGATGTCAGCATCGGCTGTACAATCAACCATAACAGAGTTTGGGATGAAATGATTTCCCTTCACATATTGGGTGAACTTCTCCATGTCAGCTTTTTCTCCCTCTTCTTTCATAAGCTCTCT'

s[463+5]
'G'

check Omo2-1 (accession=21)

select accession, alignment, bases from sequence where accession=21 and alignment=770;

| accession | alignment | bases                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
+-----------+-----------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|        21 |       770 | NNNNNNNNNNNNNNNNNNNNNNNNTTACTGTGGCACCAAGACGAAGAGGAGATAGATACCTGAAAATTCCCTCAATTCGCAGTATTTTATCCCCTGTTTCGAGGAGACCACGTAAGGTGCTAATAATTGGAAGACCAGCTCCAACGGTGGCTTCATAAAAGTAATGTGTGTACGATTTTCTTTGAAGATCTCTGATCTTTAGATACTGCATTTACAGGAGTGTTACCAATTAGTCTCTCTGTTCATAAAACCACAGAATCTTAATGGATATATACTCAAGATGAACATAAAAGTAGCTTCCGATCCAACAGTCATTTAAGCAAATATCAAACACTTTACCTGATCAAGTGGTCCAGAGTTAGCCTTTTTGTTCGGAGTGACCACATGAATTCCTCGTAGCAACCAGTCGTAGTAACAGCTAGCGATGTCAGCATCGGCTGTACAATCAACCATAACAGAGTTTGGGATAAAATGATTTCCCTTCACATATTGGGTGAACTTCTCCATGTCAGCTTTTTCTCCTTCTTCTTTCATAAGTTCTCT | 


s2 = 'NNNNNNNNNNNNNNNNNNNNNNNNTTACTGTGGCACCAAGACGAAGAGGAGATAGATACCTGAAAATTCCCTCAATTCGCAGTATTTTATCCCCTGTTTCGAGGAGACCACGTAAGGTGCTAATAATTGGAAGACCAGCTCCAACGGTGGCTTCATAAAAGTAATGTGTGTACGATTTTCTTTGAAGATCTCTGATCTTTAGATACTGCATTTACAGGAGTGTTACCAATTAGTCTCTCTGTTCATAAAACCACAGAATCTTAATGGATATATACTCAAGATGAACATAAAAGTAGCTTCCGATCCAACAGTCATTTAAGCAAATATCAAACACTTTACCTGATCAAGTGGTCCAGAGTTAGCCTTTTTGTTCGGAGTGACCACATGAATTCCTCGTAGCAACCAGTCGTAGTAACAGCTAGCGATGTCAGCATCGGCTGTACAATCAACCATAACAGAGTTTGGGATAAAATGATTTCCCTTCACATATTGGGTGAACTTCTCCATGTCAGCTTTTTCTCCTTCTTCTTTCATAAGTTCTCT'

s2[463+5]
'A'

"""
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "commit", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:e:p:s:a:n:l:i:j:o:brh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	ecotype_table = 'stock20071008.ecotype'
	accession2ecotype_table = None
	sequence_table = 'at.sequence'
	alignment_table = 'at.alignment'
	snp_locus_table = 'stock20071008.snps'
	calls_table = 'stock20071008.calls'
	input_fname1 = None
	input_fname2 = None
	diff_type_to_be_outputted = 5
	latex_output_fname = None
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
		elif opt in ("-l",):
			calls_table = arg
		elif opt in ("-i",):
			input_fname1 = arg
		elif opt in ("-j",):
			input_fname2 = arg
		elif opt in ("-o",):
			latex_output_fname = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if ecotype_table and sequence_table and alignment_table and snp_locus_table and calls_table and input_fname1 and input_fname2 and hostname and dbname and schema:
		import MySQLdb
		conn = MySQLdb.connect(db=dbname,host=hostname)
		curs = conn.cursor()
		ins= CmpTina2010VsMy2010In250kSNPs(curs, input_fname1, input_fname2, latex_output_fname, debug, report)
		ins.load_dstruc()
		ins.plot_row_NA_mismatch_rate('Tina 2010 vs 2010 db strain wise')
		ins.plot_col_NA_mismatch_rate('Tina 2010 vs 2010 db snp wise')
		ins.diff_details_table = ''
		ins.output_diff_matrix()
		ins.qc_cross_match_table = ''
		ins.cal_row_id2pairwise_dist()
		of = open(os.path.expanduser('~/script/variation/genotyping/149snp/analysis/qc_cross_Tina2010'), 'w')
		import csv
		writer = csv.writer(of, delimiter='\t')
		for row_id, pairwise_dist_ls in ins.row_id2pairwise_dist.iteritems():
			for pairwise_dist in pairwise_dist_ls:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = pairwise_dist
				writer.writerow([row_id, row_id2, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs])
		del writer, of
		#new_row_id2pairwise_dist = ins.trim_row_id2pairwise_dist(ins.row_id2pairwise_dist, 10)
	else:
		print __doc__
		sys.exit(2)
