#!/usr/bin/env python
"""
Usage: Classify2010AccessionsIntoAlignmentGroup.py [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, at(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-t ...,	alignment table, 'alignment' (default)
	-s ...,	sequence table, 'sequence' (default)
	-p ...,	alignment_type2alignment table, 'alignment_type2alignment'(default)
	-q ...,	accession2alignment_type table, 'accession2alignment_type' (default)
	-c	commit the database submission
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	Classify2010AccessionsIntoAlignmentGroup.py  -h
	
	Classify2010AccessionsIntoAlignmentGroup.py  -c
	
Description:
	For all 2010 accessions, investigate which accession has what kind of alignments.
	And group them into distinct types.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script/pymodule')))
import psycopg, sys, getopt
from codense.common import db_connect, dict_map
from sets import Set
import networkx as nx

"""
2007-10-12
	following functions till the class are used to output 2010 alignments and 149 snps in lindna (EMBOSS software to draw linear maps of dna constructs) format.
"""

def get_alignment_id2pos(curs, alignment_table='at.alignment'):
	sys.stderr.write("Getting alignment_id2pos ...")
	alignment_id2pos = {}
	curs.execute("select id, chromosome, start, end from %s"%(alignment_table))
	rows = curs.fetchall()
	for row in rows:
		alignment_id, chromosome, start, end = row
		alignment_id2pos[alignment_id] = (chromosome, start, end)
	sys.stderr.write("Done.\n")
	return alignment_id2pos

def get_snp_id2pos(curs, snp_table='stock.snps'):
	sys.stderr.write("Getting snp_id2pos ...")
	snp_id2pos = {}
	curs.execute("select id, chromosome, position from %s"%(snp_table))
	rows = curs.fetchall()
	for row in rows:
		snp_id, chromosome, position = row
		snp_id2pos[snp_id] = (chromosome, position)
	sys.stderr.write("Done.\n")
	return snp_id2pos

def get_chr_id2pos_ls(curs, alignment_type, alignment_id2pos, snp_id2pos, alignment_type2alignment_table='at.alignment_type2alignment'):
	sys.stderr.write("Getting chr_id2pos_ls ...")
	chr_id2pos_ls = {}
	curs.execute("select alignment_type, alignment_id from %s where alignment_type=%s"%(alignment_type2alignment_table, alignment_type))
	rows = curs.fetchall()
	for row in rows:
		alignment_type, alignment_id = row
		chromosome, start, end = alignment_id2pos[alignment_id]
		if chromosome not in chr_id2pos_ls:
			chr_id2pos_ls[chromosome] = []
		chr_id2pos_ls[chromosome].append((start, end, alignment_id))
	
	for snp_id, pos in snp_id2pos.iteritems():
		chromosome, position = pos
		if chromosome not in chr_id2pos_ls:
			chr_id2pos_ls[chromosome] = []
		chr_id2pos_ls[chromosome].append((position, -1, snp_id))	#-1 is an indicator of snp versus alignment
	sys.stderr.write("Done.\n")
	return chr_id2pos_ls

def get_chr_id2size(curs, chromosome_table='at.chromosome'):
	sys.stderr.write("Getting chr_id2size ...")
	curs.execute("select id, size from %s"%chromosome_table)
	chr_id2size = {}
	rows = curs.fetchall()
	for row in rows:
		chr_id, size = row
		chr_id2size[chr_id] = size
	sys.stderr.write("Done.\n")
	return chr_id2size

def output_chr_id2pos_ls(chr_id2pos_ls, chr_id2size, output_fname, label_alignment=0):
	f = open(output_fname, 'w')
	max_length = max(chr_id2size.values())
	f.write('Start\t1\n')
	f.write('End\t'+repr(max_length)+'\n')
	chr_id_ls = chr_id2pos_ls.keys()
	chr_id_ls.sort()	#starting from 1 to 2,3,4...	
	for chr_id in chr_id_ls:
		f.write('group\n')
		f.write('chr %s\n'%chr_id)
		f.write('label\n')
		f.write('Block\t1\t2\t15\n')	#the initial white block in order for the EMBOSS lindna to draw the  line before the 1st block, 15 means white (background color).
		f.write('endlabel\n')
		pos_ls = chr_id2pos_ls[chr_id]
		pos_ls.sort()	#ascending order to avoid the default connecting line ran through blocks
		for pos in pos_ls:
			f.write('label\n')
			if pos[1] == -1:	#it's a tick
				f.write('Tick\t%s\t8\tH\n'%pos[0])	#the white tick in order for the EMBOSS lindna to draw the  line before the 1st block. 8 denotes brown color. H means label in vertical
			else:
				f.write('Block\t%s\t%s\t0\tH\n'%(pos[0], pos[1]))	#0 denotes black color
				if label_alignment:
					f.write('%s\n'%pos[2])
			f.write('endlabel\n')
		#the last white tick in order for the EMBOSS lindna to draw the line after the last block till chromosome end
		f.write('label\n')
		f.write('Block\t%s\t%s\t15\n'%(chr_id2size[chr_id]-1,chr_id2size[chr_id]))
		f.write('endlabel\n')
		f.write('endgroup\n')

"""
#2007-10-12
hostname='localhost'
dbname='stock20071008'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()

alignment_id2pos = get_alignment_id2pos(curs, alignment_table='at.alignment')
snp_id2pos = get_snp_id2pos(curs, snp_table='%s.snps'%dbname)
chr_id2size = get_chr_id2size(curs, chromosome_table='at.chromosome')

for alignment_type in range(1,106):
	chr_id2pos_ls = get_chr_id2pos_ls(curs, alignment_type, alignment_id2pos, snp_id2pos, alignment_type2alignment_table='at.alignment_type2alignment')
	output_fname = '/tmp/chr_2010_alignment_type%s_149snps.txt'%alignment_type
	output_chr_id2pos_ls(chr_id2pos_ls, chr_id2size, output_fname)
"""

class Classify2010AccessionsIntoAlignmentGroup:
	"""
	2007-10-12
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', \
		alignment_table='alignment', sequence_table='sequence', alignment_type2alignment_table='alignment_type2alignment', accession2alignment_type_table='accession2alignment_type',\
		commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.alignment_table = alignment_table
		self.sequence_table = sequence_table
		self.alignment_type2alignment_table = alignment_type2alignment_table
		self.accession2alignment_type_table = accession2alignment_type_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_accession_id2alignment_id_ls(self, curs, sequence_table):
		sys.stderr.write("Getting accession_id2alignment_id_ls ...\n")
		accession_id2alignment_id_ls = {}
		curs.execute("select distinct accession, alignment from %s"%sequence_table)
		rows = curs.fetchmany(5000)
		counter = 0
		while rows:
			for row in rows:
				accession_id, alignment_id = row
				if accession_id not in accession_id2alignment_id_ls:
					accession_id2alignment_id_ls[accession_id] = []
				accession_id2alignment_id_ls[accession_id].append(alignment_id)
				counter += 1
			rows = curs.fetchmany(5000)
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		return accession_id2alignment_id_ls
	
	def classify_alignment_group(self, accession_id2alignment_id_ls):
		sys.stderr.write("Classifying alignment into groups ...")
		alignment_id_tuple2alignment_type = {}
		alignment_type2alignment_id_ls = {}
		accession_id2alignment_type = {}
		for accession_id, alignment_id_ls in accession_id2alignment_id_ls.iteritems():
			alignment_id_ls.sort()
			alignment_id_tuple = tuple(alignment_id_ls)
			if  alignment_id_tuple not in alignment_id_tuple2alignment_type:
				alignment_id_tuple2alignment_type[alignment_id_tuple] = len(alignment_id_tuple2alignment_type)+1
			alignment_type = alignment_id_tuple2alignment_type[alignment_id_tuple]
			if alignment_type not in alignment_type2alignment_id_ls:
				alignment_type2alignment_id_ls[alignment_type] = []
				alignment_type2alignment_id_ls[alignment_type] = alignment_id_ls
			
			accession_id2alignment_type[accession_id] = alignment_type
		sys.stderr.write("Done.\n")
		return accession_id2alignment_type, alignment_type2alignment_id_ls
	
	def create_accession2alignment_type_table(self, curs, accession2alignment_type_table):
		sys.stderr.write("Creating %s ..."%accession2alignment_type_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			accession_id	integer	not null,\
			alignment_type	integer not null)"%accession2alignment_type_table)
		sys.stderr.write("Done.\n")
	
	def create_alignment_type2alignment_table(self, curs, alignment_type2alignment_table):
		sys.stderr.write("Creating %s ..."%alignment_type2alignment_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			alignment_type	integer	not null,\
			alignment_id	integer not null)"%alignment_type2alignment_table)
		sys.stderr.write("Done.\n")
	
	def submit_accession_id2alignment_type(self, curs, accession_id2alignment_type, accession2alignment_type_table):
		sys.stderr.write("Submitting  accession_id2alignment_type...")
		for accession_id, alignment_type in accession_id2alignment_type.iteritems():
			curs.execute("insert into %s(accession_id, alignment_type) values(%s, %s)"%\
			(accession2alignment_type_table, accession_id, alignment_type))
		sys.stderr.write("Done.\n")
	
	def submit_alignment_type2alignment_id_ls(self, curs, alignment_type2alignment_id_ls, alignment_type2alignment_table):
		sys.stderr.write("Submitting  alignment_type2alignment_id_ls...")
		for alignment_type, alignment_id_ls in alignment_type2alignment_id_ls.iteritems():
			for alignment_id in alignment_id_ls:
				curs.execute("insert into %s(alignment_type, alignment_id) values(%s, %s)"%\
				(alignment_type2alignment_table, alignment_type, alignment_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		accession_id2alignment_id_ls = self.get_accession_id2alignment_id_ls(curs, self.sequence_table)
		accession_id2alignment_type, alignment_type2alignment_id_ls = self.classify_alignment_group(accession_id2alignment_id_ls)
		if self.commit:
			self.create_accession2alignment_type_table(curs, self.accession2alignment_type_table)
			self.create_alignment_type2alignment_table(curs, self.alignment_type2alignment_table)
			self.submit_accession_id2alignment_type(curs, accession_id2alignment_type, self.accession2alignment_type_table)
			self.submit_alignment_type2alignment_id_ls(curs, alignment_type2alignment_id_ls, self.alignment_type2alignment_table)
	
if __name__ == '__main__':
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:t:s:p:q:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'at'
	schema = 'dbsnp'
	alignment_table = 'alignment'
	sequence_table = 'sequence'
	alignment_type2alignment_table = 'alignment_type2alignment'
	accession2alignment_type_table = 'accession2alignment_type'
	commit = 0
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
		elif opt in ("-s",):
			sequence_table = arg
		elif opt in ("-t",):
			alignment_table = arg
		elif opt in ("-p",):
			alignment_type2alignment_table = arg
		elif opt in ("-q",):
			accession2alignment_type_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if hostname and dbname and schema:
		instance = Classify2010AccessionsIntoAlignmentGroup(hostname, dbname, schema, alignment_table, sequence_table,\
			alignment_type2alignment_table, accession2alignment_type_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)