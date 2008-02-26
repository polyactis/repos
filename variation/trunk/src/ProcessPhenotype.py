"""
2008-02-13
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, numpy
import warnings, traceback

class ProcessPhenotype:
	"""
	2008-02-13
		class to process phenotype data from at.phenotype and at.experiment to get flowering time.
		Flowering time is "time of first flower open" - "date counted as germination"
		
	Argument list:
		-z ..., --hostname=...	the hostname, localhost(default)
		-d ..., --dbname=...	the database name, stock20071008(default)
		-k ..., --schema=...	which schema in the database, dbsnp(default)
		-o ..., output_fname (if wanna output)
		-a ..., raw_phenotype_table=argument1, at.phenotype(default)
		-e ..., experiment_table=argument2, at.experiment_table(default)
		-b,	toggle debug
		-r, toggle report
	Examples:
		main.py -y 2 -o /tmp/phenotype.tsv
	"""
	def __init__(self, hostname, dbname, schema, output_fname=None, raw_phenotype_table='at.phenotype', experiment_table='at.experiment', debug=0, report=0):
		"""
		2008-02-20
		"""
		if not hostname or not dbname:
			print self.__doc__
			sys.exit(2)
		if not raw_phenotype_table:
			raw_phenotype_table='at.phenotype'
		if not experiment_table:
			experiment_table='at.experiment'
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.output_fname = output_fname
		self.raw_phenotype_table = raw_phenotype_table
		self.experiment_table = experiment_table
		self.debug = int(debug)
		self.report = int(report)
	
	def get_experiment_id2data(self, curs, experiment_table='at.experiment'):
		"""
		2008-02-20
			
		"""
		sys.stderr.write("Getting experiment_id2data ...\n")
		curs.execute("select id, description, day_counted_as_germination from %s"%experiment_table)
		experiment_id2data = {}
		rows = curs.fetchall()
		for row in rows:
			experiment_id, description, day_counted_as_germination = row
			day_counted_as_germination = str(day_counted_as_germination)	#transform the datetime.date type into string
			time_germination = time.strptime(day_counted_as_germination, '%Y-%m-%d')
			experiment_id2data[experiment_id] = [description, time_germination]
		sys.stderr.write("Done.\n")
		return experiment_id2data
	
	def get_experiment_id2accession_id2reading_ls(self, curs, raw_phenotype_table):
		"""
		2008-02-20
		2008-02-23
			'long ago' now mapped to 'NA'
			the reading_error_dict applied before testing whether it's NA or not
		"""
		sys.stderr.write("Getting experiment_id2accession_id2reading_ls ...\n")
		curs.execute("select experiment, accession, replicate, reading from %s where measure='first flower open'"%raw_phenotype_table)
		rows = curs.fetchall()
		experiment_id2accession_id2reading_ls = {}
		reading_error_dict = {'1/14/2003+C41':'1/14/03',
							'110/20/03':'11/20/03',
							'1/12003':'1/1/03',
							'10//24/03':'10/24/03',
							'long ago':'NA'}
		for row in rows:
			experiment_id, accession_id, replicate, reading = row
			if reading in reading_error_dict:	#correct the typo
				reading = reading_error_dict[reading]
			if reading!='NA' and reading!='N':
				if experiment_id not in experiment_id2accession_id2reading_ls:
					experiment_id2accession_id2reading_ls[experiment_id] = {}
				if accession_id not in experiment_id2accession_id2reading_ls[experiment_id]:
					experiment_id2accession_id2reading_ls[experiment_id][accession_id] = []
				try:
					if experiment_id==3 and reading[-2:]=='99':	#wrong year input
						reading = reading[:-2]+'03'
					if experiment_id==4 and reading[-2:]=='00':	#wrong year input
						reading = reading[:-2]+'02'
					if reading!='DNF':
						time_first_flower_open = time.strptime(reading, '%m/%d/%y')	#transform it into a time tuple
					else:
						time_first_flower_open = 'DNF'
				except:
					print "reading:",reading
					traceback.print_exc()
					print sys.exc_info()
					sys.exit(2)
				experiment_id2accession_id2reading_ls[experiment_id][accession_id].append(time_first_flower_open)
		sys.stderr.write("Done.\n")
		return experiment_id2accession_id2reading_ls
	
	def get_experiment_id2accession_id2FT(self, experiment_id2data, experiment_id2accession_id2reading_ls):
		"""
		2008-02-20
		"""
		sys.stderr.write("Getting experiment_id2accession_id2FT ...\n")
		experiment_id2accession_id2FT = {}
		for experiment_id, accession_id2reading_ls in experiment_id2accession_id2reading_ls.iteritems():
			description, time_germination = experiment_id2data[experiment_id]
			experiment_id2accession_id2FT[experiment_id] = {}
			for accession_id, reading_ls in accession_id2reading_ls.iteritems():
				FT_ls = []
				for reading in reading_ls:
					if reading!='DNF':
						FT = (time.mktime(reading)-time.mktime(time_germination))/(3600.0*24)	#counted as days
					else:
						FT = 200
					if FT<0:
						warnings.warn("Warning: experiment_id=%s, accession_id=%s, FT is negative: %s, reading is %s, time_germination is %s.\n"%\
										(experiment_id, accession_id, FT, reading, time_germination))
					FT_ls.append(FT)
				if len(FT_ls)>0:
					avg_FT = numpy.average(FT_ls)
					if len(FT_ls)>1:
						std_FT = numpy.std(FT_ls)
					else:
						std_FT = 'NA'
					experiment_id2accession_id2FT[experiment_id][accession_id] = [avg_FT, std_FT, len(FT_ls)]
		sys.stderr.write("Done.\n")
		return experiment_id2accession_id2FT
	
	def output_experiment_id2accession_id2FT(self, experiment_id2accession_id2FT, experiment_id2data, output_fname):
		"""
		2008-02-20
		2008-02-25
			stdev and sample size become standalone columns
		"""
		sys.stderr.write("outputting experiment_id2accession_id2FT ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		expt_id_ls = experiment_id2data.keys()
		expt_id_ls.sort()
		accession_id_ls = experiment_id2accession_id2FT[expt_id_ls[0]].keys()
		accession_id_ls.sort()
		header = ['accession id']
		for expt_id in expt_id_ls:
			header.append('%s %s'%(expt_id, experiment_id2data[expt_id][0]))
			header.append('%s %s (stdev)'%(expt_id, experiment_id2data[expt_id][0]))
			header.append('%s %s (sample size)'%(expt_id, experiment_id2data[expt_id][0]))
		writer.writerow(header)
		matrix = []
		for accession_id in accession_id_ls:
			data_row = [accession_id]
			for expt_id in expt_id_ls:
				if accession_id in experiment_id2accession_id2FT[expt_id]:
					avg_FT, std_FT, sample_size = experiment_id2accession_id2FT[expt_id][accession_id]
					data_row += [avg_FT, std_FT, sample_size]
				else:
					data_row += ['NA']*3
			matrix.append(data_row)
		for data_row in matrix:
			writer.writerow(data_row)
		del writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		experiment_id2data = self.get_experiment_id2data(curs, self.experiment_table)
		experiment_id2accession_id2reading_ls = self.get_experiment_id2accession_id2reading_ls(curs, self.raw_phenotype_table)
		experiment_id2accession_id2FT = self.get_experiment_id2accession_id2FT(experiment_id2data, experiment_id2accession_id2reading_ls)
		if self.output_fname:
			self.output_experiment_id2accession_id2FT(experiment_id2accession_id2FT, experiment_id2data, self.output_fname)