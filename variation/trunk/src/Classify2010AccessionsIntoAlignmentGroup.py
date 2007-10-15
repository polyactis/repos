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
	"""
	2007-10-12
	"""
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
	"""
	2007-10-12
	"""
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
	"""
	2007-10-12
	"""
	sys.stderr.write("Getting chr_id2pos_ls for alignment_type %s..."%alignment_type)
	chr_id2pos_ls = {}
	curs.execute("select alignment_type, alignment_id from %s where alignment_type=%s"%(alignment_type2alignment_table, alignment_type))
	rows = curs.fetchall()
	chr_id2alignment_pos_ls = {}
	for row in rows:
		alignment_type, alignment_id = row
		chromosome, start, end = alignment_id2pos[alignment_id]
		if chromosome not in chr_id2pos_ls:
			chr_id2pos_ls[chromosome] = []
			chr_id2alignment_pos_ls[chromosome] = []
		chr_id2pos_ls[chromosome].append((start, end, alignment_id))
		chr_id2alignment_pos_ls[chromosome].append((start, end, alignment_id))
	
	chr_id2snp_pos_ls = {}
	for snp_id, pos in snp_id2pos.iteritems():
		chromosome, position = pos
		if chromosome not in chr_id2pos_ls:
			chr_id2pos_ls[chromosome] = []
		if chromosome not in chr_id2snp_pos_ls:
			chr_id2snp_pos_ls[chromosome] = []
		chr_id2pos_ls[chromosome].append((position, -1, snp_id))	#-1 is an indicator of snp versus alignment
		chr_id2snp_pos_ls[chromosome].append((position, -1, snp_id))
	
	#sort all pos_ls
	for chr_id in chr_id2pos_ls:	#chr_id2pos_ls.keys() encompass chr_id2alignment_pos_ls.keys() and chr_id2snp_pos_ls.keys()
		chr_id2pos_ls[chr_id].sort()
		if chr_id in chr_id2alignment_pos_ls:
			chr_id2alignment_pos_ls[chr_id].sort()
		if chr_id in chr_id2snp_pos_ls:
			chr_id2snp_pos_ls[chr_id].sort()
	sys.stderr.write("Done.\n")
	return chr_id2pos_ls, chr_id2alignment_pos_ls, chr_id2snp_pos_ls

def get_chr_id2size(curs, chromosome_table='at.chromosome'):
	"""
	2007-10-12
	"""
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
	"""
	2007-10-12
	"""
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

def get_alignment_type_group_ls(curs, alignment_type2alignment_table = 'at.alignment_type2alignment', alignments_cnt_range_ls=[(1,10),(10,50),(100,150),(150,200),(1300,1600)]):
	"""
	2007-10-12
	"""
	sys.stderr.write("Getting alignment_type_group_ls ...")
	alignment_type_group_ls = []
	alignment_type_cnt_group_ls = []
	for cnt_range in alignments_cnt_range_ls:
		curs.execute("select alignment_type, count(alignment_id) as cnt from %s group by alignment_type  order by cnt, alignment_type"%alignment_type2alignment_table)
		rows = curs.fetchall()
		alignment_type_group = []
		alignment_type_cnt_group = []
		for row in rows:
			alignment_type, cnt = row
			if cnt>=cnt_range[0] and cnt<=cnt_range[1]:
				alignment_type_group.append(alignment_type)
				alignment_type_cnt_group.append(cnt)
		alignment_type_group_ls.append(alignment_type_group)
		alignment_type_cnt_group_ls.append(alignment_type_cnt_group)
	sys.stderr.write("Done.\n")
	return alignment_type_group_ls, alignment_type_cnt_group_ls

def is_snp_within_alignment_pos_ls(snp_pos, alignment_pos_ls):
	"""
	2007-10-14
		alignment_pos_ls is already sorted
	"""
	position = snp_pos[0]
	return_value = False
	for alignment_pos in alignment_pos_ls:
		if position>=alignment_pos[0] and position<=alignment_pos[1]:	#found one
			return_value = True
			break
		if position<alignment_pos[0]:	#skip the rest
			break
	return return_value


"""
#some tests for is_snp_within_alignment_pos_ls()
alignment_pos_ls= [ [1,10,1],[100,1000,2]]
is_snp_within_alignment_pos_ls([101,-1], alignment_pos_ls)
"""

def get_alignment_type2no_of_aligns_and_149_snps(curs, alignment_id2pos, snp_id2pos, alignment_type2alignment_table='at.alignment_type2alignment'):
	"""
	2007-10-14
	"""
	sys.stderr.write("Getting alignment_type2no_of_aligns_and_149_snps ...\n")
	alignment_type2no_of_aligns_and_149_snps = {}
	curs.execute("select distinct alignment_type from %s"%alignment_type2alignment_table)
	rows = curs.fetchall()
	for row in rows:
		alignment_type = row[0]
		chr_id2pos_ls, chr_id2alignment_pos_ls, chr_id2snp_pos_ls = get_chr_id2pos_ls(curs, alignment_type, alignment_id2pos, snp_id2pos, alignment_type2alignment_table)
		no_of_overlappings = 0
		for chr_id, snp_pos_ls in chr_id2snp_pos_ls.iteritems():
			if chr_id in chr_id2alignment_pos_ls:
				for snp_pos in snp_pos_ls:
					if is_snp_within_alignment_pos_ls(snp_pos, chr_id2alignment_pos_ls[chr_id]):
						no_of_overlappings += 1
		alignment_type2no_of_aligns_and_149_snps[alignment_type] = [sum(map(len,chr_id2alignment_pos_ls.values())), no_of_overlappings]
	sys.stderr.write("Done.\n")
	return alignment_type2no_of_aligns_and_149_snps

def LatexSummaryAlignmentType(curs, alignments_cnt_range_ls, alignment_type_group_ls, alignment_type_cnt_group_ls, fig_fname_ls, alignment_type2no_of_aligns_and_149_snps, output_fname, accession_table='at.accession', ecotype2accession_table='at.ecotype2accession', stock_schema='stock20071008', accession2alignment_type_table='at.accession2alignment_type'):
	"""
	2007-10-12
	2007-10-15
		beautify the format, but still use the tabular. change the page size in 2010SequenceReport.tex to fit the table in. rather than use supertabular or longtable. cuz the table is too wide.
	"""
	sys.stderr.write("Output alignment summary in latex ...\n")
	from pymodule.latex import escape_characters
	of = open(output_fname, 'w')
	for i in range(len(alignment_type_group_ls)):
		sys.stderr.write("\t alignment group %s "%i)
		of.write('\\begin{table}\n')
		alignment_type_group = alignment_type_group_ls[i]
		cnt_range = alignments_cnt_range_ls[i]
		alignment_type_cnt_group = alignment_type_cnt_group_ls[i]
		flabel = 'fatg%s'%i
		caption = "The number of alignments for these strains is from %s to %s falling into %s. %s distinct alignment combinations. alignment-type is id for one kind of distinct alignment combination. \\#fragments is the number of fragments sequenced in 2010. \\#149snps is the number of the 149 snps falling into that strain's fragments. Figure~\\ref{%s} is an example chromosome chart of this table."%(alignment_type_cnt_group[0], alignment_type_cnt_group[-1], cnt_range, len(alignment_type_group), flabel)
		of.write('\\caption{%s}\n'%caption)
		label = 'tatg%s'%i
		of.write('\\label{%s}\n'%label)
		of.write('\\begin{tabular}{|rrrr|rrrrrrrr|rrr|}\n')
		of.write('\\hline\n')
		#get all accessions within this alignment_type_group
		curs.execute("select aa.accession_id, a.name, a.origin, a.number, aa.alignment_type from %s aa, %s a where a.id=aa.accession_id and aa.alignment_type in (%s) order by name, accession_id "%(accession2alignment_type_table, accession_table, repr(map(int ,alignment_type_group))[1:-1]))
		rows = curs.fetchall()
		of.write('\\multicolumn{15}{|c|}{%s strains} \\\\\n'%(len(rows)) )
		of.write('\\hline\n')
		of.write('\\multicolumn{4}{|c|}{accession table} & \\multicolumn{8}{|c|}{ecotype table} & \\multicolumn{3}{|c|}{fragments info} \\\\\n')
		of.write('\\hline\n')
		of.write('id & name & origin & number & id & name & nativename & stockparent & lat & lon & site & country & alignment-type & \\#fragments & \\#149snps\\\\\n')
		of.write('\\hline\n')
		for row in rows:	#table entry starts here. one by one
			accession_id, name, origin, number, alignment_type = row
			curs.execute("select e.id, e.name, e.nativename, e.stockparent, e.latitude, e.longitude, s.name, c.abbr from %s.ecotype e, %s.address a, %s.site s, %s.country c, %s ea where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and e.id=ea.ecotype_id and ea.accession_id=%s "%(stock_schema, stock_schema, stock_schema, stock_schema, ecotype2accession_table, accession_id))
			ecotype_rows = curs.fetchall()
			if ecotype_rows:
				ecotype_id, ecotype_name, nativename, stockparent, latitude, longitude, site, country = ecotype_rows[0]
			else:
				ecotype_id=ecotype_name= nativename= stockparent= latitude= longitude= site= country = ''
			no_of_alns, no_of_149snps = alignment_type2no_of_aligns_and_149_snps[alignment_type]
			one_table_entry = '%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s  \\\\\n'%(accession_id, name, origin, number, ecotype_id, ecotype_name, nativename, stockparent, latitude, longitude, site, country, alignment_type, no_of_alns, no_of_149snps)
			one_table_entry = escape_characters(one_table_entry)
			of.write(one_table_entry)
		of.write('\\hline\n')
		of.write('\\end{tabular}\n')
		of.write('\\end{table}\n')
		of.write('\n')
		of.write('\\begin{figure}\n')
		of.write('\\includegraphics[width=1\\textwidth, height=1\\textheight]{figures/%s}\n'%fig_fname_ls[i])
		of.write('\\caption{Red ticks are locations of 149 snps. Black blocks are locations of fragments.}\\label{%s}\n'%flabel)
		of.write('\\end{figure}\n')
		of.write('\n')
		sys.stderr.write("\n")
	del of
	sys.stderr.write("Done.\n")


"""
#2007-10-12
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
hostname='localhost'
dbname='stock20071008'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()
alignment_type2alignment_table='at.alignment_type2alignment'

alignment_id2pos = get_alignment_id2pos(curs, alignment_table='at.alignment')
snp_id2pos = get_snp_id2pos(curs, snp_table='%s.snps'%dbname)
chr_id2size = get_chr_id2size(curs, chromosome_table='at.chromosome')


for alignment_type in range(1,106):
	chr_id2pos_ls, chr_id2alignment_pos_ls, chr_id2snp_pos_ls = get_chr_id2pos_ls(curs, alignment_type, alignment_id2pos, snp_id2pos, alignment_type2alignment_table)
	output_fname = '/tmp/chr_2010_alignment_type%s_149snps.txt'%alignment_type
	output_chr_id2pos_ls(chr_id2pos_ls, chr_id2size, output_fname)

alignment_type2no_of_aligns_and_149_snps = get_alignment_type2no_of_aligns_and_149_snps(curs, alignment_id2pos, snp_id2pos, alignment_type2alignment_table)

alignments_cnt_range_ls=[(1,10),(10,50),(100,150),(150,200),(1300,1600)]
alignment_type_group_ls, alignment_type_cnt_group_ls = get_alignment_type_group_ls(curs, alignment_type2alignment_table, alignments_cnt_range_ls)

accession_table='at.accession'
ecotype2accession_table='at.ecotype2accession'
stock_schema='stock20071008'
accession2alignment_type_table='at.accession2alignment_type'

fig_fname_ls = []	#this corresponds to alignments_cnt_range_ls
for alignment_type in [105, 97, 101, 98, 92]:
	fig_fname = 'chr_2010_alignment_type%s_149snps.eps'%alignment_type
	fig_fname_ls.append(fig_fname)

output_fname = 'alignment_type_tables_figures.tex'
LatexSummaryAlignmentType(curs, alignments_cnt_range_ls, alignment_type_group_ls, alignment_type_cnt_group_ls, fig_fname_ls, alignment_type2no_of_aligns_and_149_snps, output_fname, accession_table, ecotype2accession_table, stock_schema, accession2alignment_type_table)
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
