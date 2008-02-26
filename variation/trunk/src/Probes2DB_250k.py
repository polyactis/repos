#!/usr/bin/env python
"""
Usage: Probes2DB_250k.py [OPTIONS] -i input_fname

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock20071008(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-i ...,	input file, snp.RData outputted by write.table() in R
	-p ...,	probes table, 'probes_250k'(default)
	-s ...,	snps table, specify if the input is derived from snp.RData
	-n	both the probes snps tables are new. to be created(IGNORE)
	-c	commit the database submission
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	Probes2DB_250k.py -i /tmp/snp -c
	
	Probes2DB_250k.py -i /tmp/snp -d stock_250k -c -p probes -s snps
	
	Probes2DB_250k.py -i /tmp/atsnptile1 -d stock_250k -c -p probes

Description:
	1st usage (with snps_table specified, -s xxx): read the R's write.table() output of snp.RData
	(250k SNPs, made by Xu Zhang of Justin Borevitz lab)
	and dump them into probes and snps table.
	
	2nd usage (with snps_table unspecified, -s ''): dump the R's write.table() output of atsnptile1.RData
	into table probes
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import sys, getopt, csv

class Probes2DB_250k:
	"""
	2007-12-10
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', \
		input_fname='', probes_table='probes_250k', snps_table='', \
		new_table=0, commit=0, debug=0, report=0):
		"""
		2007-12-10
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.probes_table = probes_table
		self.snps_table = snps_table
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def create_probes_table(self, curs, probes_table):
		sys.stderr.write("Creating %s ..."%probes_table)
		curs.execute("create table %s(\
			id	integer auto_increment primary key,\
			snpid	integer,\
			seq	varchar(200),\
			chromosome	integer,\
			position	integer,\
			allele	varchar(1),\
			strand	varchar(10),\
			xpos	integer,\
			ypos	integer)"%probes_table)
		sys.stderr.write("Done.\n")
	
	def create_snps_table(self, curs, snps_table):
		sys.stderr.write("Creating %s ..."%snps_table)
		curs.execute("create table %s(\
			id	integer auto_increment primary key,\
			snpid	varchar(200),\
			chromosome	integer,\
			position	integer,\
			allele1	varchar(1),\
			allele2	varchar(2))"%snps_table)
		sys.stderr.write("Done.\n")
		
	def get_snps_and_probes(self, input_fname):
		"""
		2007-12-10
			row[1:] => row[1:9] more specific
		"""
		sys.stderr.write("Getting snps and probes ...")
		reader = csv.reader(open(input_fname), delimiter='\t')
		snps_ls = []
		probes_ls = []
		chr_pos2snpid = {}
		snpid2allele_ls = {}
		reader.next()	#toss out the 1st header row
		snpid = 0
		for row in reader:
			snpacc, seq, chr, position, allele, strand, xpos, ypos = row[1:9]	#row[0] is data frame row number. So write.table() of R outputs data row with one more column than its header.
			chr = int(chr)
			position = int(position)
			xpos = int(xpos)
			ypos = int(ypos)
			chr_pos_key = (chr, position)
			if chr_pos_key not in chr_pos2snpid:
				snpid += 1
				chr_pos2snpid[chr_pos_key] = snpid
				snps_ls.append([snpid, snpacc, chr, position])
			snpid = chr_pos2snpid[chr_pos_key]
			if snpid not in snpid2allele_ls:
				snpid2allele_ls[snpid] = []
			snpid2allele_ls[snpid].append(allele)	#each allele has two copies (antisense, sense)
			probes_ls.append([snpid, seq, chr, position, allele, strand, xpos, ypos])
		del reader
		sys.stderr.write("Done.\n")
		return snps_ls, probes_ls, snpid2allele_ls
	
	def read_atsnptile1(self, input_fname):
		"""
		2008-02-25
		"""
		sys.stderr.write("Read in probes from %s..."%os.path.basename(input_fname))
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()	#toss out the 1st header row
		probes_ls = []
		for row in reader:
			xpos, ypos, chr, position, direction, seq, gene, RNA, tu, flank, expressedClones, totalClones,\
			strand, multiTranscript, LerDel, LerCopy, LerSNPdelL, LerSNPdelR, LerSNPpos, promoter, utr5,\
			utr3, intron, intergenic, downstream, cda = row[1:]	#row[0] is data frame row number. So write.table() of R outputs data row with one more column than its header.
			chr = int(chr)
			position = int(position)
			xpos = int(xpos)
			ypos = int(ypos)
			try:
				expressedClones = int(expressedClones)
			except:
				expressedClones = int(round(float(expressedClones)))	#something like 8.001, 8.002, 0.001, 0.003 appear in it. deem it as typo.
			totalClones = int(totalClones)
			vars_need_NULL_transform = [LerCopy, LerSNPdelL, LerSNPdelR, LerSNPpos]	#need to transform 'NA' to 'NULL' in order for db insertion
			for i in range(len(vars_need_NULL_transform)):
				if vars_need_NULL_transform[i] =='NA':
					vars_need_NULL_transform[i] = 'NULL'
				else:
					vars_need_NULL_transform[i] = int(vars_need_NULL_transform[i])
			LerCopy, LerSNPdelL, LerSNPdelR, LerSNPpos = vars_need_NULL_transform
			probes_ls.append([seq, chr, position, strand, xpos, ypos, direction, gene, RNA, tu, flank,\
							expressedClones, totalClones, multiTranscript, LerDel, LerCopy, LerSNPdelL,\
							LerSNPdelR, LerSNPpos, promoter, utr5, utr3, intron, intergenic, downstream, cda])
		del reader
		sys.stderr.write("Done.\n")
		return probes_ls	
	
	def submit_snps_ls(self, curs, snps_ls, snpid2allele_ls, snps_table):
		sys.stderr.write("Submitting snps_ls...")
		for snpid, snpacc, chr, position in snps_ls:
			allele_ls = [snpid2allele_ls[snpid][0], snpid2allele_ls[snpid][2]]	#take the 1st and 3rd position
			curs.execute("insert into %s(id, snpid, chromosome, position, allele1, allele2) values(%s, '%s', %s, %s, '%s', '%s')"%(snps_table, snpid, snpacc, chr, position, allele_ls[0], allele_ls[1]))
		sys.stderr.write("Done.\n")
	
	def submit_probes_ls(self, curs, probes_ls, probes_table):
		"""
		2008-02-18
			snpid in probes table is renamed to snps_id
		"""
		sys.stderr.write("Submitting probes_ls ...")
		for snpid, seq, chr, position, allele, strand, xpos, ypos in probes_ls:
			curs.execute("insert into %s(snps_id, seq, chromosome, position, allele, strand, xpos, ypos) values (%s, '%s', %s, %s, '%s', '%s', %s, %s)"%(probes_table, snpid, seq, chr, position, allele, strand, xpos, ypos) )
		sys.stderr.write("Done.\n")
	
	def submit_atsnptile1(self, curs, probes_ls, probes_table):
		"""
		2008-02-25
			submit atsnptile1
		"""
		sys.stderr.write("Submitting probes_ls of atsnptile1 ...")
		for row in probes_ls:
			seq, chr, position, strand, xpos, ypos, direction, gene, RNA, tu, flank,\
				expressedClones, totalClones, multiTranscript, LerDel, LerCopy, LerSNPdelL,\
				LerSNPdelR, LerSNPpos, promoter, utr5, utr3, intron, intergenic, downstream, cda = row
			curs.execute("insert into %s(seq, chromosome, position, strand, xpos, ypos, direction, gene, RNA, tu, flank,\
				expressedClones, totalClones, multiTranscript, LerDel, LerCopy, LerSNPdelL,\
				LerSNPdelR, LerSNPpos, promoter, utr5, utr3, intron, intergenic, downstream, cda) values \
				('%s', %s, %s, '%s', %s, %s, '%s', '%s', '%s', '%s', '%s',\
				%s, %s, '%s', '%s', %s, %s, \
				%s, %s, %s, %s, %s, %s, %s, %s, %s )"%\
				(probes_table, seq, chr, position, strand, xpos, ypos, direction, gene, RNA, tu, flank,\
				expressedClones, totalClones, multiTranscript, LerDel, LerCopy, LerSNPdelL,\
				LerSNPdelR, LerSNPpos, promoter, utr5, utr3, intron, intergenic, downstream, cda) )
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2008-02-18
			tables are created according to mysql.sql
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		if self.snps_table:
			snps_ls, probes_ls, snpid2allele_ls = self.get_snps_and_probes(self.input_fname)
		else:
			probes_ls = self.read_atsnptile1(self.input_fname)
		if self.commit:
			#if self.new_table:
			#	self.create_probes_table(curs, self.probes_table)
			#	self.create_snps_table(curs, self.snps_table)
			if self.snps_table:
				self.submit_snps_ls(curs, snps_ls, snpid2allele_ls, self.snps_table)
				self.submit_probes_ls(curs, probes_ls, self.probes_table)
			else:
				self.submit_atsnptile1(curs, probes_ls, self.probes_table)
			curs.execute("commit")

if __name__ == '__main__':
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:p:s:ncbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock20071008'
	schema = 'dbsnp'
	input_fname = ''
	probes_table = 'probes_250k'
	snps_table = ''
	new_table = 0
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
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-p",):
			probes_table = arg
		elif opt in ("-s",):
			snps_table = arg
		elif opt in ("-n",):
			new_table = 1
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if hostname and dbname and schema and input_fname and probes_table:
		instance = Probes2DB_250k(hostname, dbname, schema, input_fname,\
			probes_table, snps_table, new_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
