#!/usr/bin/env python

import os,sys

class DisplaySNPMatrix(object):
	"""
	2007-03-05
		show an image visualizing SNP data
	2007-06-05
		set aspect='auto' in imshow(), the default (pylab.image.rcParams['image.aspect'])='equal', which is bad
	"""
	def display_snp_matrix(input_fname, output_fname=None, need_sort=0, need_savefig=0, xlabel='', ylabel=''):
		import csv, Numeric, pylab
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		data_matrix = []
		for row in reader:
			data_row = row[2:]
			data_row = map(int, data_row)
			data_matrix.append(data_row)
		del reader
		data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
		if need_sort:
			data_matrix.sort()
		data_matrix = Numeric.array(data_matrix)
		
		pylab.clf()
		pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
		pylab.colorbar()
		if xlabel:
			pylab.xticks([data_matrix.shape[1]/2], [xlabel])
		if ylabel:
			pylab.yticks([data_matrix.shape[0]/2], [ylabel])
		if need_savefig:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
	
	def make_snp_matrix_legend(value_ls, label_ls, output_fname=None):
		"""
		2007-10-25
			to pair with display_snp_matrix()
		"""
		import numpy, pylab
		label_ls_copy = label_ls[:]
		data_matrix = numpy.zeros([len(value_ls), 1], numpy.int)
		for i in range(len(value_ls)):
			data_matrix[i,0] = value_ls[i]
		label_ls_copy.reverse()	#pylab put the label on starting from the bottom of the picture
		pylab.clf()
		pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')
		pylab.yticks(range(len(value_ls)), label_ls_copy, fontsize=60, verticalalignment='bottom', horizontalalignment='right')
		pylab.xticks([],[])
		if output_fname:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
	
	"""
	2007-10-25
	value_ls = range(3)
	label_ls = ['NA', 'A', 'C']
	make_snp_matrix_legend(value_ls, label_ls)
	
	"""


class DB149(object):
	
	"""
	2007-04-09
	"""
	def draw_SNP_gap_histogram(curs, snp_locus_table, output_fname, need_savefig=0):
		SNP_gap_ls = []
		prev_chromosome = None
		prev_position = None
		curs.execute("select chromosome, position from %s order by chromosome, position"%snp_locus_table)
		rows = curs.fetchall()
		for row in rows:
			chromosome, position = row
			if prev_chromosome==None or prev_position==None:
				prev_chromosome = chromosome
				prev_position = position
			elif chromosome==prev_chromosome:
				SNP_gap_ls.append(position-prev_position)
			prev_chromosome = chromosome
			prev_position = position
		import pylab
		pylab.clf()
		pylab.hist(SNP_gap_ls, 20)
		pylab.title("hist of SNP gap")
		if need_savefig:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
		return SNP_gap_ls
	
	@classmethod
	def blast_snp_segment(cls, curs, snp_locus_table, output_fname, database_fname, flanking_seq_length=12, \
						max_no_of_hits_to_be_outputted=3, blast_bin_path=os.path.expanduser('~/bin/blast/bin/blastall'), \
						annot_assembly_table='sequence.annot_assembly', \
						raw_sequence_table='sequence.raw_sequence', tmp_blast_infname='/tmp/blast_input'):
		"""
		2007-04-30
		"""
		import sys
		from Bio.Blast import NCBIXML, NCBIStandalone
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
		from codense.common import get_sequence_segment
		curs.execute("select s.id, s.acc, s.chromosome, s.position, a.gi from %s s, %s a\
			where a.tax_id=s.tax_id and a.chromosome=s.chromosome"%(snp_locus_table, annot_assembly_table))
		rows = curs.fetchall()
		counter = 0
		outf = open(output_fname, 'w')
		for row in rows:
			snp_locus_id, acc, chromosome, position, genomic_gi = row
			genomic_start = position-flanking_seq_length
			genomic_stop = position + flanking_seq_length
			seq = get_sequence_segment(curs, genomic_gi, genomic_start, genomic_stop, annot_assembly_table, raw_sequence_table)
			if seq:
				for candidate_snp in ['A', 'C', 'G', 'T']:
					seq = seq[:flanking_seq_length] + candidate_snp + seq[flanking_seq_length+1:]	#str doesn't support item assignment
					inf = open(tmp_blast_infname, 'w')
					inf.write('>%s(snp_locus_id=%s) candidate_snp=%s\n'%(acc, snp_locus_id, candidate_snp))
					inf.write('%s\n'%seq)
					del inf
					result_handle, error_info = NCBIStandalone.blastall(blast_bin_path, "blastn", database_fname, tmp_blast_infname, align_view=7)	#align_view=7 toggles the xml output
					blast_parser = NCBIXML.BlastParser()	#the parser has to be re-instanced every time otherwise the blast records parsed in the previous round stay in the parser
					blast_records = blast_parser.parse(result_handle)
					del blast_parser
					no_of_hits = min(max_no_of_hits_to_be_outputted, len(blast_records.alignments))
					outf.write("%s (id=%s) candidate_snp=%s\n"%(acc, snp_locus_id, candidate_snp))
					outf.write("query sequence = %s\n"%seq)
					for i in range(no_of_hits):
						outf.write("\t%s\n"%(blast_records.alignments[i].title))
						for hsp in blast_records.alignments[i].hsps:
							outf.write("\t\tquery_start=%s, sbjct_start=%s, frame=%s, identities=%s/%s E=%s\n"%(hsp.query_start, hsp.sbjct_start, hsp.frame, hsp.identities, 2*flanking_seq_length+1, hsp.expect))
							outf.write("\t\t%s\n"%hsp.query)
							outf.write("\t\t%s\n"%hsp.match)
							outf.write("\t\t%s\n"%hsp.sbjct)
							outf.write("\n")
						outf.write("\n")
					outf.write("\n")
			sys.stderr.write("%s%s"%('\x08'*10, counter))
			counter += 1
		del outf
	
	"""
	my_blast_db = os.path.expanduser("~/bin/blast/db/Arabidopsis_thaliana.main_genome.fasta")
	DB149.blast_snp_segment(curs, 'snp_locus', './blast_149snps_vs_thaliana_len_25.txt', my_blast_db)
	
	
	my_blast_db = os.path.expanduser("~/bin/blast/db/Arabidopsis_lyrata.main_genome.scaffolds.fasta")
	DB149.blast_snp_segment(curs, 'snp_locus', './blast_149snps_vs_lyrata_len_51.txt', my_blast_db, flanking_seq_length=25)
	"""
	
	"""
	2007-04-30
	"""
	def fill_snp_locus_table_with_25mer_thaliana_call(curs, snp_locus_table, annot_assembly_table='sequence.annot_assembly', \
		raw_sequence_table='sequence.raw_sequence', need_commit=0):
		import sys
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
		from codense.common import get_sequence_segment
		curs.execute("select s.id, s.acc, s.chromosome, s.position, a.gi from %s s, %s a\
			where a.tax_id=s.tax_id and a.chromosome=s.chromosome"%(snp_locus_table, annot_assembly_table))
		rows = curs.fetchall()
		counter = 0
		for row in rows:
			snp_locus_id, acc, chromosome, position, genomic_gi = row
			genomic_start = position - 12
			genomic_stop = position + 12
			seq = get_sequence_segment(curs, genomic_gi, genomic_start, genomic_stop, annot_assembly_table, raw_sequence_table)
			if seq:
				thaliana_call = seq[12]
				curs.execute("update %s set thaliana_call='%s', flanking_25mer='%s' where id=%s"%(snp_locus_table, thaliana_call, seq, snp_locus_id))
			sys.stderr.write("%s%s"%('\x08'*10, counter))
			counter += 1
		if need_commit:
			curs.execute("end")
	
	"""
	fill_snp_locus_table_with_25mer_thaliana_call(curs, 'dbsnp.snp_locus', need_commit=1)
	"""
	
	
	"""
	2007-03-05
		add the position info given by Yan Li from Borevitz Lab
	"""
	def fill_snp_locus_table_with_position_info(curs, input_fname, output_table, need_commit):
		import csv
		reader = csv.reader(open(input_fname))
		snp_locus_acc_list = reader.next()[1:]
		snp_locus_position_list = reader.next()[1:]
		for i in range(len(snp_locus_acc_list)):
			snp_locus_acc = snp_locus_acc_list[i]
			snp_locus_acc = snp_locus_acc.replace('.', ' ')
			snp_locus_position = snp_locus_position_list[i]
			chromosome, position = snp_locus_position.split('-')
			curs.execute("update %s set chromosome='%s', position=%s where acc='%s'"%(output_table, chromosome, position, snp_locus_acc))
		if need_commit:
			curs.execute("end")
	
	
	"""
	fill_snp_locus_table_with_position_info(curs, './script/variation/data/snp_position.csv', 'dbsnp.snp_locus', need_commit=1)
	"""
	
	"""
	2007-09-21
		check duplicate calls in the database
	"""
	def check_inconsistent_duplicate_calls(curs, ecotype_table, calls_table, debug=0):
		"""
		2007-09-22
			buggy, don't use it. use check_inconsistent_duplicate_calls2()
		"""
		sys.stderr.write("Checking inconsistent duplicate calls ...\n")
		from common import nt2number
		no_of_strain_snp_pairs = 0
		no_of_inconsistent = 0
		inconsistency = 0
		old_call_number = 0
		old_strain_snp_pair = None
		curs.execute("select e.nativename, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
		rows = curs.fetchall()
		counter = 0
		from sets import Set
		strain_snp_pair_set = Set()
		inconsistent_dup_strain_snp_pair_set = Set()
		if debug:
			import pdb
			pdb.set_trace()
		for row in rows:
			nativename, snpid, call, callhet = row
			nativename = nativename.upper()	#bug here. same name appears in >1 forms differing in cases
			call = call.upper()
			if callhet:
				callhet.upper()
				call = call+callhet
			call_number = nt2number[call]
			strain_snp_pair = (nativename, snpid)
			if strain_snp_pair != old_strain_snp_pair:
				no_of_strain_snp_pairs += 1
				no_of_inconsistent += inconsistency
				if inconsistency == 1:
					inconsistent_dup_strain_snp_pair_set.add(old_strain_snp_pair)
				old_strain_snp_pair = strain_snp_pair
				inconsistency = 0
			else:
				if inconsistency == 0 and old_call_number!=0 and call_number!=0 and old_call_number!=call_number:	#if inconsistency==1, don't change it. comparison only between same pairs
					inconsistency = 1
			old_call_number = call_number
			counter += 1
			strain_snp_pair_set.add(strain_snp_pair)
		#don't miss out the last one
		no_of_inconsistent += inconsistency
		old_strain_snp_pair = strain_snp_pair
		sys.stderr.write('%s%s'%('\x08'*20, counter))
		sys.stderr.write("\nDone.\n")
		print "%s strain snp pairs"%(len(strain_snp_pair_set))
		print 'inconsistent ratio: %s/%s=%s'%(no_of_inconsistent, no_of_strain_snp_pairs, float(no_of_inconsistent)/no_of_strain_snp_pairs)
		return strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set
	
	def check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=1, debug=0):
		"""
		use_ecotypeid_as_strainname controls whether it's nativename or ecotypeid in strain_snp_pair
		"""
		sys.stderr.write("Checking inconsistent duplicate calls ...\n")
		if strainname_type==1:
			print 'distinct on (nativename,snpid) pair:'
		elif strainname_type == 2:
			print 'distinct on (ecotypeid,snpid) pair:'
		elif strainname_type == 3:
			print 'distinct on (nativename, stockparent ,snpid) pair:'
		else:
			sys.stderr.write("unsupported strainname_type %s\n"%strainname_type)
			return None, None
		from common import nt2number
		no_of_strain_snp_pairs = 0
		no_of_inconsistent = 0
		old_strain_snp_pair = None
		duplicated_times = 0
		if strainname_type == 1:
			curs.execute("select e.nativename, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
		elif strainname_type == 2:
			curs.execute("select e.id, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
		elif strainname_type == 3:
			curs.execute("select e.nativename, e.stockparent, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, stockparent, snpid"%(ecotype_table, calls_table))
		rows = curs.fetchall()
		counter = 0
		from sets import Set
		strain_snp_pair_set = Set()
		inconsistent_dup_strain_snp_pair_set = Set()
		duplicate_call_type2counter = {}
		if debug:
			import pdb
			pdb.set_trace()
		for row in rows:
			if strainname_type==3:
				nativename, stockparent, snpid, call, callhet = row
			elif strainname_type==1:
				nativename, snpid, call, callhet = row
			elif strainname_type==2:
				ecotypeid, snpid, call, callhet = row
			if strainname_type!=2:
				nativename = nativename.upper()	#bug here. same name appears in >1 forms differing in cases
			call = call.upper()
			if callhet:
				callhet.upper()
				call = call+callhet
			call_number = nt2number[call]
			if strainname_type==1:
				strain_snp_pair = (nativename, snpid)
			elif strainname_type==2:
				strain_snp_pair = (ecotypeid, snpid)
			elif strainname_type==3:
				strain_snp_pair = (nativename, stockparent, snpid)
			if old_strain_snp_pair==None:	#first time
				old_strain_snp_pair = strain_snp_pair
				duplicated_times += 1
			elif strain_snp_pair != old_strain_snp_pair:
				if duplicated_times>1:
					no_of_strain_snp_pairs += 1
					if len(duplicate_call_type2counter)>1:
						no_of_inconsistent += 1
						inconsistent_dup_strain_snp_pair_set.add(old_strain_snp_pair)
				old_strain_snp_pair = strain_snp_pair
				duplicated_times = 1	#the current position is counted as 1
				duplicate_call_type2counter = {}
			else:
				duplicated_times += 1
			if call_number!=0:	#have to be put last
				if call_number not in duplicate_call_type2counter:
					duplicate_call_type2counter[call_number] = 0
				duplicate_call_type2counter[call_number] += 1
			counter += 1
			strain_snp_pair_set.add(strain_snp_pair)
		#don't miss out the last one
		if duplicated_times>1:
			no_of_strain_snp_pairs += 1
			if len(duplicate_call_type2counter)>1:
				no_of_inconsistent += 1
		old_strain_snp_pair = strain_snp_pair
		print "%s strain snp pairs"%counter
		print "%s distinct strain snp pairs"%(len(strain_snp_pair_set))
		print 'inconsistent ratio: %s/%s=%s'%(no_of_inconsistent, no_of_strain_snp_pairs, float(no_of_inconsistent)/no_of_strain_snp_pairs)
		sys.stderr.write("Done.\n")
		return strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set
	
	
	def get_nativename_snp_distinct_set(curs, ecotype_table, calls_table):
		"""
		2007-09-22 find out the bug caused by the case-insensitive problem of sql.
		"""
		curs.execute("select distinct e.nativename, c.snpid from %s e, %s c where e.id=c.ecotypeid" %(ecotype_table, calls_table))
		from sets import Set
		strain_snp_pair_set = Set()
		rows = curs.fetchall()
		for row in rows:
			nativename, snpid = row
			strain_snp_pair = (nativename, snpid)
			strain_snp_pair_set.add(strain_snp_pair)
		return strain_snp_pair_set
	"""
	ecotype_table = 'ecotype'
	calls_table = 'calls'
	strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set = check_inconsistent_duplicate_calls(curs, ecotype_table, calls_table)
	strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
	strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
	strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)
	
	nativename_snp_distinct_set = get_nativename_snp_distinct_set(curs, ecotype_table, calls_table)
	
	strain_snp_pair_set - nativename_snp_distinct_set gives 1937 different pairs. turns out that sql's "distinct e.nativename, c.snpid" ignores the cases , like Kz-9 and KZ-9 are same. however python doesn't do that.
	"""

class Data149_PopulationStructure(object):
	#2007-06-17 function to calculate the great circle distance of a sphere given latitude and longitude
	# http://en.wikipedia.org/wiki/Great-circle_distance
	def cal_great_circle_distance(lat1, lon1, lat2, lon2, earth_radius=6372.795):
		import math
		lat1_rad = lat1*math.pi/180
		lon1_rad = lon1*math.pi/180
		lat2_rad = lat2*math.pi/180
		lon2_rad = lon2*math.pi/180
		long_diff = abs(lon1_rad-lon2_rad)
		sin_lat1 = math.sin(lat1_rad)
		cos_lat1 = math.cos(lat1_rad)
		sin_lat2 = math.sin(lat2_rad)
		cos_lat2 = math.cos(lat2_rad)
		spheric_angular_diff = math.atan2(math.sqrt(math.pow(cos_lat2*math.sin(long_diff),2) + math.pow(cos_lat1*sin_lat2-sin_lat1*cos_lat2*math.cos(long_diff),2)), sin_lat1*sin_lat2+cos_lat1*cos_lat2*math.cos(long_diff))
		return earth_radius*spheric_angular_diff
	
	#2007-06-17 draw location of all strains onto a map
	def test_draw_strain_location(pic_area=[-180,-90,180,90]):
		import pylab
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
	            resolution='l',projection='mill')
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,1,1])
		m.drawmeridians(pylab.arange(-180,180,60), labels=[1,1,1,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		
		import MySQLdb
		db = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
		c = db.cursor()
		c.execute("""select latitude, longitude from ecotype where latitude is not null and longitude is not null""")
		lat_lon_ls = c.fetchall()
		euc_coord1_ls = []
		euc_coord2_ls = []
		for lat, lon in lat_lon_ls:
			euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
		m.scatter(euc_coord1_ls, euc_coord2_ls, 25, marker='o', zorder=10)
		pylab.title("worldwide distribution of %s strains"%(len(lat_lon_ls)))
		pylab.show()
	
	#2007-06-17 calculate pairwise distance among strains by calling cal_great_circle_distance()
	def cal_pairwise_distance_of_strains():
		import MySQLdb, pylab
		db = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
		c = db.cursor()
		c.execute("""select latitude, longitude from ecotype where latitude is not null and longitude is not null""")
		lat_lon_ls = c.fetchall()
		distance_ls = []
		no_of_sites = len(lat_lon_ls)
		for i in range(no_of_sites):
			for j in range(i+1, no_of_sites):
				distance_ls.append(cal_great_circle_distance(lat_lon_ls[i][0], lat_lon_ls[i][1], lat_lon_ls[j][0], lat_lon_ls[j][1]))
		pylab.clf()
		pylab.hist(distance_ls, 20)
		pylab.title("histogram of distances between %s strains(km)"%no_of_sites)
		pylab.show()
	
	"""
	#2007-06-18, calll cal_great_circle_distance()
	import MySQLdb
	conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
	cursor = conn.cursor()
	cursor.execute("select latitude, longitude from ecotype where latitude is not null and longitude is not null")
	lat_lon_ls = cursor.fetchall()
	"""
	
	
	"""
	2007-09-13
	following functions to investigate inter/intra-population identity
	"""
	def construct_identity_pair_ls(input_fname):
		"""
		2007-09-13
			input_fname is a data frame file with snps coded in integers. 1st row is the ecotype id in table ecotype
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		no_of_strains = len(strain_acc_list)
		no_of_snps = len(header)-2
		identity_pair_ls = []
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				no_of_same_cols = 0
				for k in range(no_of_snps):
					if data_matrix[i][k] == data_matrix[j][k] or data_matrix[i][k]==0 or data_matrix[j][k]==0:
						no_of_same_cols += 1
				if no_of_same_cols == no_of_snps:
					identity_pair_ls.append([strain_acc_list[i], strain_acc_list[j]])
		return identity_pair_ls
	
	
	def construct_graph_out_of_identity_pair(identity_pair_ls):
		"""
		2007-09-18
			construct a graph based on identity pairs in identity_pair_ls
			identity_pair_ls's entries will be converted to integers
		"""
		sys.stderr.write("Constructing graph out of identity_pair_ls ...")
		import networkx as nx
		g = nx.XGraph()
		for identity_pair in identity_pair_ls:
			identity_pair = map(int, identity_pair)
			g.add_edge(identity_pair[0], identity_pair[1], 1)
		sys.stderr.write("Done.\n")
		return g
	
	def get_singleton_strain_id_ls(g, input_fname):
		"""
		2007-10-01
			node in g is in ecotype id format
			so strain_acc (1st column) in input_fname should be integer (ecotype id)
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		sys.stderr.write("Getting singleton_strain_id_ls ...")
		from sets import Set
		graph_node_set = Set(g.nodes())
		singleton_strain_id_ls = []
		for strain_acc in strain_acc_list:
			strain_acc = int(strain_acc)
			if strain_acc not in graph_node_set:
				singleton_strain_id_ls.append(strain_acc)
		sys.stderr.write("Done.\n")
		return singleton_strain_id_ls
	
	def construct_site_graph_out_of_strain_graph(g, ecotypeid2pos, node_weight_by_no_of_nodes=0):
		"""
		2007-09-20
			#reduce the drawing by omitting identical edges (linking same GPS positions), construct a graph whose node is (lat,lon) and then plot this graph.
		2007-10-01
			add node_weight_by_no_of_nodes
		"""
		sys.stderr.write("Constructing site graph out of strain graph...")
		site2pos = {}
		site2weight = {}
		import networkx as nx
		site_g = nx.XGraph()
		for e in g.edges():
			pos1 = (round(ecotypeid2pos[e[0]][0],2), round(ecotypeid2pos[e[0]][1],2))
			pos2 = (round(ecotypeid2pos[e[1]][0],2), round(ecotypeid2pos[e[1]][1],2))
			if not node_weight_by_no_of_nodes:	#2007-10-01
				if pos1 not in site2weight:
					site2weight[pos1] = 0
				site2weight[pos1] += 1
				if pos2 not in site2weight:
					site2weight[pos2] = 0
				site2weight[pos2] += 1
			if pos1==pos2:
				site_g.add_node(pos1)
			else:
				if not site_g.has_edge(pos1, pos2):
					site_g.add_edge(pos1, pos2, 0)
				site_g.adj[pos1][pos2] += 1
				site_g.adj[pos2][pos1] += 1
		if node_weight_by_no_of_nodes:	#2007-10-01
			for n in g.nodes():
				pos = (round(ecotypeid2pos[n][0],2), round(ecotypeid2pos[n][1],2))
				if pos not in site2weight:
					site2weight[pos] = 0
				site2weight[pos] += 1
		for n in site_g:
			site2pos[n] = n
		sys.stderr.write("Done.\n")
		return site_g, site2weight, site2pos
	
	def collapseStrainsWithSamePos(ecotypeid_ls, ecotypeid2pos):
		"""
		2007-10-09
		"""
		sys.stderr.write("Collapsing Strains with same GPS ...")
		site2pos = {}
		site2weight = {}
		for ecotypeid in ecotypeid_ls:
			pos = (round(ecotypeid2pos[ecotypeid][0],2), round(ecotypeid2pos[ecotypeid][1],2))
			if pos not in site2weight:
				site2weight[pos] = 0
				site2pos[pos] = pos
			site2weight[pos] += 1
		sys.stderr.write("Done.\n")
		return site2weight, site2pos
	
	
	def draw_strains_on_map(ecotypeid_ls, ecotypeid2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None, label_ls=None, need_collapse_strains_with_same_pos=0):
		"""
		2007-10-09
		2007-11-08
			add label_ls, need_collapse_strains_with_same_pos
		"""
		import os, sys
		sys.stderr.write("Drawing strains on a map ...\n")
		import pylab, math
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		
		
		sys.stderr.write("\tDrawing nodes ...")
		euc_coord1_ls = []
		euc_coord2_ls = []
		diameter_ls = []
		ax=pylab.gca()
		if need_collapse_strains_with_same_pos:
			site2weight, site2pos = collapseStrainsWithSamePos(ecotypeid_ls, ecotypeid2pos)
			for n in site2pos:
				lat, lon = site2pos[n]
				euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
				euc_coord1_ls.append(euc_coord1)
				euc_coord2_ls.append(euc_coord2)
				diameter_ls.append(8*math.sqrt(site2weight[n]))
		else:
			for i in range(len(ecotypeid_ls)):
				ecotypeid = ecotypeid_ls[i]
				lat, lon = ecotypeid2pos[ecotypeid]
				euc_coord1, euc_coord2 = m(lon, lat)
				euc_coord1_ls.append(euc_coord1)
				euc_coord2_ls.append(euc_coord2)
				if label_ls:
					ax.text(euc_coord1, euc_coord2, str(label_ls[i]), size=2, alpha=0.5, horizontalalignment='center', verticalalignment='center', zorder=12)
				diameter_ls.append(3)
		import numpy
		#diameter_ls = numpy.array(diameter_ls)
		m.scatter(euc_coord1_ls, euc_coord2_ls, diameter_ls, marker='o', color='r', alpha=0.4, zorder=10, faceted=False)
		sys.stderr.write("Done.\n")
		
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		
		pylab.title(pic_title)
		if output_fname_prefix:
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=600)
			#pylab.savefig('%s.svg'%output_fname_prefix, dpi=600)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=600)
		sys.stderr.write("Done.\n")
	
	def display_matrix_of_component(input_fname, ecotypeid_ls, ecotypeid2pos, output_fname=None, need_sort=0, need_savefig=0):
		"""
		2007-09-20
			display the data from that component
			
		"""
		import csv, Numeric, pylab
		cc_ecotypeid_pos = []
		for ecotypeid in ecotypeid_ls:
			cc_ecotypeid_pos.append(ecotypeid2pos[ecotypeid])
		import numpy
		argsort_index = numpy.argsort(cc_ecotypeid_pos, 0)	#watch it's two dimensional
		ecotypeid2row_index = {}
		ytick_label_ls = []
		cc_size = len(ecotypeid_ls)
		for i in range(cc_size):
			ecotypeid_ls_index = argsort_index[i][1]	#sort based on longitude
			ecotypeid = ecotypeid_ls[ecotypeid_ls_index]
			ecotypeid2row_index[ecotypeid] = i
			ytick_label_ls.append('%s (%.2f, %.2f)'%(ecotypeid, ecotypeid2pos[ecotypeid][0], ecotypeid2pos[ecotypeid][1]))
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		data_matrix = [0]*cc_size
		for row in reader:
			ecotypeid = int(row[0])
			if ecotypeid in ecotypeid2row_index:
				data_row = row[2:]
				data_row = map(int, data_row)
				data_matrix[ecotypeid2row_index[ecotypeid]] = data_row
		del reader
		data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
		if need_sort:
			data_matrix.sort()
		data_matrix = Numeric.array(data_matrix)
		
		pylab.clf()
		pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
		pylab.colorbar()
		pylab.yticks(range(cc_size), ytick_label_ls)
		if need_savefig:
			pylab.savefig('%s.eps'%output_fname, dpi=300)
			pylab.savefig('%s.svg'%output_fname, dpi=300)
			pylab.savefig('%s.png'%output_fname, dpi=300)
		pylab.show()
	
	def get_ecotypeid2pos(curs, ecotype_table, with_data_affiliated=0):
		"""
		2007-09-18
		2007-10-09
			add with_data_affiliated
		"""
		sys.stderr.write("Getting ecotypeid2pos from %s ..."%ecotype_table)
		ecotypeid2pos = {}
		if with_data_affiliated:
			curs.execute("select distinct e.id, e.latitude, e.longitude from %s e, calls c where e.latitude is not null and e.longitude is not null and e.id=c.ecotypeid"%ecotype_table)
		else:
			curs.execute("select id, latitude, longitude from %s where latitude is not null and longitude is not null"%ecotype_table)
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, latitude, longitude = row
			ecotypeid2pos[ecotypeid] = (latitude, longitude)
		sys.stderr.write("Done.\n")
		return ecotypeid2pos
	
	def get_popid2pos_size_and_ecotypeid2pop_id(curs, popid2ecotypeid_table, selected=0):
		"""
		2007-09-13
		2007-09-24
			doesn't have to be selected
		"""
		sys.stderr.write("Getting popid2pos_size & ecotypeid2pop_id ...")
		popid2pos_size = {}
		ecotypeid2pop_id = {}
		curs.execute("select popid, latitude, longitude, ecotypeid from %s where selected=%s"%(popid2ecotypeid_table, selected))
		rows = curs.fetchall()
		for row in rows:
			popid, latitude, longitude, ecotypeid = row
			if popid not in popid2pos_size:
				popid2pos_size[popid] = [None, 0]
			if popid2pos_size[popid][0]==None:
				popid2pos_size[popid][0] = (latitude, longitude)
			popid2pos_size[popid][1] += 1
			ecotypeid2pop_id[ecotypeid] = popid
		sys.stderr.write("Done.\n")
		return popid2pos_size, ecotypeid2pop_id
	
	def construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls):
		"""
		2007-09-18
		"""
		import os, sys
		popid2pos_size, ecotypeid2pop_id = get_popid2pos_size_and_ecotypeid2pop_id(curs, popid2ecotypeid_table)
		
		no_of_inter_pop_identities = 0
		no_of_valid_identities = 0	#valid means both strains belong to some population
		import networkx as nx
		g = nx.XGraph()
		popid2intra_pop_connections = {}
		for e in identity_pair_ls:
			e = map(int, e)
			popid1 = ecotypeid2pop_id.get(e[0])
			popid2 = ecotypeid2pop_id.get(e[1])
			if popid1 and popid2:	#both are classified into populations
				no_of_valid_identities += 1
				if popid1 not in popid2intra_pop_connections:
					popid2intra_pop_connections[popid1] = 0
					g.add_node(popid1)
				if popid2 not in popid2intra_pop_connections:
					popid2intra_pop_connections[popid2] = 0
					g.add_node(popid2)
				if popid1!=popid2:
					if not g.has_edge(popid1, popid2):
						g.add_edge(popid1, popid2, 0)
					g.adj[popid1][popid2] += 1	#increase the edge weight
					g.adj[popid2][popid1] += 1
					no_of_inter_pop_identities += 1
				else:
					popid2intra_pop_connections[popid1] += 1
		print '%s/%s=%s inter population identities'%(no_of_inter_pop_identities, no_of_valid_identities, float(no_of_inter_pop_identities)/no_of_valid_identities)
		return g, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id
	
	def get_pic_area(pos_ls):
		"""
		2007-10-02
			get a good pic_area
		"""
		min_lat = pos_ls[0][0]
		max_lat = pos_ls[0][0]
		min_lon = pos_ls[0][1]
		max_lon = pos_ls[0][1]
		for i in range(1, len(pos_ls)):
			pos = pos_ls[i]
			if pos[0]<min_lat:
				min_lat = pos[0]
			if pos[0]>max_lat:
				max_lat = pos[0]
			if pos[1]<min_lon:
				min_lon = pos[1]
			if pos[1]>max_lon:
				max_lon = pos[1]
		pic_area = [-180, -90, 180, 90]	#this is the boundary
		if min_lon-5>pic_area[0]:
			pic_area[0] = min_lon-5
		if min_lat-5>pic_area[1]:
			pic_area[1] = min_lat-5
		if max_lon+5<pic_area[2]:
			pic_area[2] = max_lon+5
		if max_lat+5<pic_area[3]:
			pic_area[3] = max_lat+5
		return pic_area
	
	def draw_graph_on_map(g, node2weight, node2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None, need_draw_edge=0):
		"""
		2007-09-13
			identity_pair_ls is a list of pairs of strains (ecotype id as in table ecotype)
		2007-10-08
			correct a bug in 4*diameter_ls, diameter_ls has to be converted to array first.
			sqrt the node weight, 8 times the original weight
		"""
		import os, sys
		sys.stderr.write("Drawing graph on a map ...\n")
		import pylab, math
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		
		sys.stderr.write("\tDrawing nodes ...")
		euc_coord1_ls = []
		euc_coord2_ls = []
		diameter_ls = []
		for n in g.nodes():
			lat, lon = node2pos[n]
			euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
			diameter_ls.append(math.sqrt(node2weight[n]))
		import numpy
		diameter_ls = numpy.array(diameter_ls)
		m.scatter(euc_coord1_ls, euc_coord2_ls, 8*diameter_ls, marker='o', color='r', alpha=0.4, zorder=12, faceted=False)
		sys.stderr.write("Done.\n")
		
		if need_draw_edge:
			sys.stderr.write("\tDrawing edges ...")
			ax=pylab.gca()
			for popid1, popid2, no_of_connections in g.edges():
				lat1, lon1 = node2pos[popid1]
				lat2, lon2 = node2pos[popid2]
				x1, y1 = m(lon1, lat1)
				x2, y2 = m(lon2, lat2)
				ax.plot([x1,x2],[y1,y2], 'g', linewidth=math.log(no_of_connections+1)/2, alpha=0.2, zorder=10)
			sys.stderr.write("Done.\n")
		
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		
		pylab.title(pic_title)
		if output_fname_prefix:
			pylab.savefig('%s_site_network.eps'%output_fname_prefix, dpi=600)
			pylab.savefig('%s_site_network.svg'%output_fname_prefix, dpi=600)
			pylab.savefig('%s_site_network.png'%output_fname_prefix, dpi=600)
		sys.stderr.write("Done.\n")
	
	def check_cross_ocean_components(g, cc, ecotypeid2pos, longi_watershed=-30):
		from sets  import Set
		cross_ocean_cc2lon_pair_set = {}
		for i in range(len(cc)):
			gs = g.subgraph(cc[i])
			for e in gs.edges():
				l1 = ecotypeid2pos[e[0]][1]
				l2 = ecotypeid2pos[e[1]][1]
				if (l1<longi_watershed and l2>longi_watershed) or (l1>longi_watershed and l2<longi_watershed):
					if i not in cross_ocean_cc2lon_pair_set:
						cross_ocean_cc2lon_pair_set[i] = Set()
					cross_ocean_cc2lon_pair_set[i].add((l1, l2))
		return cross_ocean_cc2lon_pair_set
	
	def get_clique_ls_indexed_by_cc_number(g):
		"""
		2007-10-01 each component is partitioned into cliques
			the resulting clique_ls is put in the order of connected components
			like clique_ls_indexed_by_cc_number[0] is clique_ls from connected component 0
		"""
		sys.stderr.write("Getting clique_ls_indexed_by_cc_number ...\n")
		import networkx as nx
		from pymodule.PartitionGraphIntoCliques import PartitionGraphIntoCliques
		PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(algorithm_type=2)
		g_cc_ls = nx.connected_components(g)
		clique_ls_indexed_by_cc_number = []
		for i in range(len(g_cc_ls)):
			sys.stderr.write("%s\t%s"%('\x08'*20, i))
			cc_g = g.subgraph(g_cc_ls[i])
			PartitionGraphIntoCliques_ins.partition(cc_g.copy())
			clique_ls_indexed_by_cc_number.append(PartitionGraphIntoCliques_ins.clique_ls)
		sys.stderr.write("Done.\n")
		return g_cc_ls, clique_ls_indexed_by_cc_number
	
	def get_haplotype_size2number(clique_ls_indexed_by_cc_number):
		"""
		2007-10-02
		"""
		sys.stderr.write("Getting haplotype_size2number ...")
		haplotype_size2number = {}
		for i in range(len(clique_ls_indexed_by_cc_number)):
			for clique in clique_ls_indexed_by_cc_number[i]:
				haplotype_size = len(clique)
				if haplotype_size not in haplotype_size2number:
					haplotype_size2number[haplotype_size] = 0
				haplotype_size2number[haplotype_size] += 1
		sys.stderr.write("Done.\n")
		return haplotype_size2number
	
	def get_ordered_clique_ls(clique_ls_indexed_by_cc_number):
		"""
		2007-10-02
			return cliques with size in descending order.
		"""
		sys.stderr.write("Getting ordered_clique_ls ...")
		clique_ls = []
		ordered_clique_ls = []
		import numpy
		for i in range(len(clique_ls_indexed_by_cc_number)):
			for clique in clique_ls_indexed_by_cc_number[i]:
				clique_ls.append(clique)
		clique_ls_argsort = list(numpy.argsort(map(len, clique_ls)))
		clique_ls_argsort.reverse()	#from big to small
		for index in clique_ls_argsort:
			ordered_clique_ls.append(clique_ls[index])
		sys.stderr.write("Done.\n")
		return ordered_clique_ls
	
	def get_strains_within_longitude_span(ecotypeid_ls, ecotypeid2pos, longitude_span=[-12, 50]):
		"""
		2007-10-04
			check the map roughly
			[-12,50] is europe
			[-130,-60] is north america
			[60,90] is central asia
			[130, 150] is japan
			
		"""
		no_of_strains = 0
		strain_ls = []
		for ecotypeid in ecotypeid_ls:
			pos = ecotypeid2pos[ecotypeid]
			if pos[1]>=longitude_span[0] and pos[1]<=longitude_span[1]:
				no_of_strains += 1
				strain_ls.append(ecotypeid)
		print "#strains with %s is %s"%(longitude_span, no_of_strains)
		return strain_ls
	
	def draw_histogram_of_haplotype_size(haplotype_size2number, output_fname_prefix=None):
		"""
		2007-10-02
		"""
		import pylab
		pylab.clf()
		haplotype_size_ls = []
		count_ls = []
		for haplotype_size, count in haplotype_size2number.iteritems():
			haplotype_size_ls.append(haplotype_size)
			count_ls.append(count)
		pylab.title("haplotype_size bar chart")
		pylab.xlabel("haplotype_size")
		pylab.ylabel("count")
		pylab.bar(haplotype_size_ls, count_ls)
		if output_fname_prefix:
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.show()
	
	def calGenoDistAndGeoDist(data_matrix_fname, ecotypeid2pos, longitude_span=[-180, 180], dist_range=[0.0,1.0]):
		"""
		2008-01-30
			add dist_range to restrict the output pairs within that genotype distance
		2008-01-29
			return strain_pair_ls as well
			add a general cutoff
		2007-10-03
			plot to see relationship between genotype distance and geographic distance
			-cal_great_circle_distance()
		2007-10-04
			add longitude_span
		"""
		sys.stderr.write("Calculating geno, geo distance ...\n")
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		strain_pair_ls = []
		geno_dist_ls = []
		geo_dist_ls = []
		no_of_NA_pairs = 0
		no_of_strains = len(data_matrix)
		no_of_snps = len(data_matrix[0])
		no_of_pairs = 0
		no_of_included_pairs = 0
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				pos1 = ecotypeid2pos[int(strain_acc_list[i])]
				pos2 = ecotypeid2pos[int(strain_acc_list[j])]
				if pos1[1]>=longitude_span[0] and pos1[1]<=longitude_span[1] and pos2[1]>=longitude_span[0] and pos2[1]<=longitude_span[1]:
					no_of_pairs += 1
					no_of_valid_pairs = 0.0
					no_of_matching_pairs = 0.0
					for k in range(no_of_snps):
						if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
							no_of_valid_pairs += 1
							if data_matrix[i][k] == data_matrix[j][k]:
								no_of_matching_pairs += 1
					if no_of_valid_pairs!=0:
						distance = 1 - no_of_matching_pairs/no_of_valid_pairs
						if distance>=dist_range[0] and distance<=dist_range[1]:
							geno_dist_ls.append(distance)
							geo_dist = cal_great_circle_distance(pos1[0], pos1[1], pos2[0], pos2[1])
							geo_dist_ls.append(geo_dist)
							strain_pair_ls.append([int(strain_acc_list[i]), int(strain_acc_list[j])])
							no_of_included_pairs += 1
					else:
						no_of_NA_pairs += 1
		print "out of %s pairs, %s are NA, %s pairs' distance is within range(%s)."%(no_of_pairs, no_of_NA_pairs, no_of_included_pairs, dist_range)
		return strain_pair_ls, geno_dist_ls, geo_dist_ls
	
	def calGenoDistAndGeoDistBetweenTwoAreas(data_matrix_fname, ecotypeid2pos, longitude_span1=[-180, 180], longitude_span2=[-180, 180]):
		"""
		2007-10-05
			modified from calGenoDistAndGeoDist with two longitude spans
		"""
		sys.stderr.write("Calculating geno, geo distance ...\n")
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		geno_dist_ls = []
		geo_dist_ls = []
		no_of_NA_pairs = 0
		no_of_strains = len(data_matrix)
		no_of_snps = len(data_matrix[0])
		no_of_pairs = 0
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				pos1 = ecotypeid2pos[int(strain_acc_list[i])]
				pos2 = ecotypeid2pos[int(strain_acc_list[j])]
				if (pos1[1]>=longitude_span1[0] and pos1[1]<=longitude_span1[1] and pos2[1]>=longitude_span2[0] and pos2[1]<=longitude_span2[1]) or (pos2[1]>=longitude_span1[0] and pos2[1]<=longitude_span1[1] and pos1[1]>=longitude_span2[0] and pos1[1]<=longitude_span2[1]):
					no_of_pairs += 1
					no_of_valid_pairs = 0.0
					no_of_matching_pairs = 0.0
					for k in range(no_of_snps):
						if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
							no_of_valid_pairs += 1
							if data_matrix[i][k] == data_matrix[j][k]:
								no_of_matching_pairs += 1
					if no_of_valid_pairs!=0:
						distance = 1 - no_of_matching_pairs/no_of_valid_pairs
						geno_dist_ls.append(distance)
						geo_dist = cal_great_circle_distance(pos1[0], pos1[1], pos2[0], pos2[1])
						geo_dist_ls.append(geo_dist)
					else:
						no_of_NA_pairs += 1
		print "out of %s pairs, %s are NA"%(no_of_pairs, no_of_NA_pairs)
		return geno_dist_ls, geo_dist_ls
	
	def plotHist(ls, title=None, output_fname_prefix=None, no_of_bins=40):
		import pylab
		pylab.clf()
		pylab.hist(ls, no_of_bins)
		if title:
			pylab.title(title)
		if output_fname_prefix:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.show()
	
	def plotXY(xls, yls, title=None, xlabel=None, ylabel=None, output_fname_prefix=None):
		import pylab
		pylab.clf()
		
		pylab.title(title)
		pylab.xlabel(xlabel)
		pylab.ylabel(ylabel)
		pylab.plot(xls, yls, '.')
		if output_fname_prefix:
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.show()
	
	def outputGraphInCliquerFormat(graph, output_fname):
		"""
		2007-10-05
			to test if cliquer suffers the same problem
		"""
		of = open(output_fname, 'w')
		of.write("p edge %s %s\n"%(graph.number_of_nodes(), graph.number_of_edges()))
		node_name2index = {}
		for n in graph.nodes():
			node_name2index[n] = len(node_name2index)+1
		for e in graph.edges():
			of.write("e %s %s\n"%(node_name2index[e[0]], node_name2index[e[1]]))
		del of
	
	def get_ecotypeid_ls_given_popid(curs, popid2ecotypeid_table, popid):
		ecotypeid_ls = []
		curs.execute("select ecotypeid from %s where popid=%s and selected=1"%(popid2ecotypeid_table, popid))
		rows = curs.fetchall()
		for row in rows:
			ecotypeid_ls.append(row[0])
		return ecotypeid_ls
	
	"""
	2008-01-25
		following functions to investigate how identity spreads against distance in different regions.
	"""
	def find_span_index_out_of_span_ls(pos, longitude_span_ls):
		"""
		2008-01-25
			find the index of the longitude_span in which pos drops in
		"""
		span_index = -1
		for i in range(len(longitude_span_ls)):
			if pos[1]>=longitude_span_ls[i][0] and pos[1]<=longitude_span_ls[i][1]:
				span_index = i
		return span_index
	
	def get_span_pair_count_ls(input_fname, ecotypeid2pos, span_ls):
		"""
		2008-01-30
			calculate the number of possible pairs for a specific span type (either single span or two crossing-spans)
			for normalization purpose in group_identity_pair_according_to_region_distance()
		"""
		sys.stderr.write("Getting span_pair_count_ls ...")
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		strain_acc_list = map(int, strain_acc_list)
		#1st check how many strains in each span
		span_count_ls = [0]*len(span_ls)
		for ecotypeid in strain_acc_list:
			pos = ecotypeid2pos[ecotypeid]
			for i in range(len(span_ls)):
				if pos[1]>=span_ls[i][0] and pos[1]<=span_ls[i][1]:
					span_count_ls[i] += 1
		#2nd calculate the number of pairs in each span
		span_pair_count_ls = []
		for i in range(len(span_count_ls)):
			span_pair_count_ls.append(span_count_ls[i]*span_count_ls[i]/2.0)
		#3rd calculate the number of pairs crossing two spans
		#check all combinations of span_ls
		for i in range(len(span_ls)):
			for j in range(i+1, len(span_ls)):
				span_pair_count_ls.append(span_count_ls[i]*span_count_ls[j])
		sys.stderr.write("Done.\n")
		return span_pair_count_ls
	
	def group_identity_pair_according_to_region_distance(strain_pair_ls, geno_dist_ls, geo_dist_ls, ecotypeid2pos, span_ls, span_label_ls, geo_distance_window_size=100, max_identity_geno_dist=0.0, span_pair_count_ls=[]):
		"""
		2008-01-30
			add span_pair_count_ls to normalize the counts
		2008-01-29
			strain_pair_ls, geno_dist_ls, geo_dist_ls replace identity_pair_ls
			use max_identity_geno_dist to determine whether to call it identity or not
		2008-01-25
			to group the identity_pair_ls into different regions and distance windows.
			geo_distance_window_size is in unit km
		"""
		sys.stderr.write("Grouping identity pairs according to region and geographic distance ...\n")
		import math
		list_of_window_no2count = []
		#1. get the row labels (label for the span type of each identity_pair)
		#2. link the span_index_pair to the index of list_of_window_no2count
		#3. initiate list_of_window_no2count
		row_label_ls =[] #region designation
		span_index_pair2list_index = {}
		for i in range(len(span_ls)):
			row_label_ls.append(span_label_ls[i])
			span_index_pair = (i,i)
			span_index_pair2list_index[span_index_pair] = i
			list_of_window_no2count.append({})
		for i in range(len(span_ls)):
			for j in range(i+1, len(span_ls)):
				row_label_ls.append('%s - %s'%(span_label_ls[i], span_label_ls[j]))
				span_index_pair2list_index[(i,j)] = len(row_label_ls)-1 #both order have the same index
				span_index_pair2list_index[(j,i)] = len(row_label_ls)-1
				list_of_window_no2count.append({})
	
		#record all the counts
		max_window_no = -1 #used to determine the final maxtrix shape
		for i in range(len(strain_pair_ls)):
			strain_pair = strain_pair_ls[i]
			geno_dist = geno_dist_ls[i]
			geo_dist = geo_dist_ls[i]
			if geno_dist<=max_identity_geno_dist:
				pos1 = ecotypeid2pos[strain_pair[0]]
				pos2 = ecotypeid2pos[strain_pair[1]]
				span_index1 = find_span_index_out_of_span_ls(pos1, span_ls)
				span_index2 = find_span_index_out_of_span_ls(pos2, span_ls)
				if span_index1>=0 and span_index2>=0:	#both data are in span_ls
					window_no = math.floor(geo_dist/geo_distance_window_size) #take the floor as the bin number
					if window_no>max_window_no:
						max_window_no = window_no
					list_index = span_index_pair2list_index[(span_index1, span_index2)]
					window_no2count = list_of_window_no2count[list_index]
					if window_no not in window_no2count:
						window_no2count[window_no] = 0
					window_no2count[window_no] += 1
		#transform list_of_window_no2count to matrix
		import numpy
		if max_window_no>=0: #-1 means nothing recorded
			matrix_of_counts_by_region_x_geo_dist = numpy.zeros([len(list_of_window_no2count), max_window_no+1], numpy.float) #+1 because window_no starts from 0
			for i in range(len(list_of_window_no2count)):
				for j in range(max_window_no+1):
					if j in list_of_window_no2count[i]:
						if span_pair_count_ls:
							matrix_of_counts_by_region_x_geo_dist[i][j] = list_of_window_no2count[i][j]/float(span_pair_count_ls[i])	#2008-01-30 normalize
						else:
							matrix_of_counts_by_region_x_geo_dist[i][j] = list_of_window_no2count[i][j]
		else:
			matrix_of_counts_by_region_x_geo_dist = None
		#generate the labels for the columns of matrix_of_counts_by_region_x_geo_dist
		col_label_ls = [] #geo distance window designation
		col_index_ls = [] #the X-axis value
		for i in range(max_window_no+1):
			col_index = (i+0.5)*geo_distance_window_size #right in the middle of the window/bin
			col_index_ls.append(col_index)
			col_label_ls.append('%skm'%int(col_index)) #km is unit
		sys.stderr.write("Done.\n")
		return matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls
	
	def draw_table_as_bar_chart(matrix, row_label_ls, col_label_ls, col_index_ls=None, output_fname_prefix=None):
		"""
		2008-01-30
			bars from different rows no longer sit on one another. get them on the same ground and one after another.
		2008-01-25
			originated from http://matplotlib.sourceforge.net/screenshots/table_demo.py
		"""
		sys.stderr.write("Drawing Table as bar char ...")
		import matplotlib
		import pylab
		from pymodule.colours import get_colours
		pylab.clf()
		pylab.axes([0.1, 0.3, 0.8, 0.6])   # leave room below the axes for the table
		# Get some pastel shades for the colours
		colours = get_colours(len(row_label_ls))
		colours = list(colours)
		#colours.reverse()
		rows = len(matrix) #matrix.shape[0] will also do, but in case matrix is a 2D list
		ind = pylab.arange(len(col_label_ls))  # the x locations for the groups
		cellText = []
		width = 0.2    # the width of the bars
		#width = (col_index_ls[1]-col_index_ls[0])/2     # the width of the bars is half the window size
		yoff = pylab.array([0.0] * len(col_label_ls)) # the bottom values for stacked bar chart
		bar_ins_ls = []
		step = 0.32
		for row in xrange(rows):
		    bar_ins = pylab.bar(ind+step*row, matrix[row], width, bottom=yoff, color=colours[row])
		    bar_ins_ls.append(bar_ins)
		    #yoff = yoff + matrix[row]
		    cellText.append(matrix[row])
		# Add a table at the bottom of the axes
		colours.reverse()
		cellText.reverse()
		row_label_ls_reverse = row_label_ls[:]
		row_label_ls_reverse.reverse()
		the_table = pylab.table(cellText=cellText, rowLabels=row_label_ls_reverse, rowColours=colours, colLabels=col_label_ls)
		pylab.ylabel("Counts")
		pylab.xlabel("Geographic Distance")
		#vals = arange(0, 2500, 500)
		#pylab.yticks(vals*1000, ['%d' % val for val in vals])
		#pylab.xticks([])
		pylab.title('Counts of identity pairs by region X geo distance')
		
		pylab.legend([bar_ins[0] for bar_ins in bar_ins_ls], row_label_ls, shadow=True)
		if output_fname_prefix:
			pylab.savefig('%s_identity_bar_char_by_region_x_geo_dist.svg'%output_fname_prefix, dpi=600)
			pylab.savefig('%s_identity_bar_char_by_region_x_geo_dist.png'%output_fname_prefix, dpi=600)
		pylab.show()
		sys.stderr.write("Done.\n")
	
	"""
	from misc import *
	
	input_fname = 'script/variation/stock20071008/data_d110_c0_5.tsv'
	identity_pair_ls = construct_identity_pair_ls(input_fname)
	#store it in a file as it takes a long time to generate
	import cPickle
	of = open('%s_identity_pair_ls'%os.path.splitext(input_fname)[0], 'w')
	cPickle.dump(identity_pair_ls, of)
	del of
	
	import cPickle
	inf = open('%s_identity_pair_ls'%os.path.splitext(input_fname)[0], 'r')
	identity_pair_ls = cPickle.load(inf)
	del inf
	
	#codes to group strains into populations
	popid2ecotypeid_table = 'popid2ecotypeid_25'
	pop_iden_g, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls)
	
	popid2pos = {}
	for popid, pos_size in popid2pos_size.iteritems():
		popid2pos[popid] = pos_size[0]
	
	import math
	popid2weight = {}
	for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
		popid2weight[popid] = 8*math.log(intra_pop_connections+1)
	
	draw_graph_on_map(pop_iden_g, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map1')
	
	##graph of strains
	strain_iden_g = construct_graph_out_of_identity_pair(identity_pair_ls)
	ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
	ecotypeid2weight = {}
	for n in strain_iden_g:
		ecotypeid2weight[n] = 1
	
	import networkx as nx
	strain_iden_g_cc = nx.connected_components(strain_iden_g)
	
	#2007-10-09 draw all strains on map
	output_fname_prefix = '%s'%os.path.splitext(input_fname)[0]
	strain_map_fname_prefix = '%s_strain_map'%output_fname_prefix
	draw_strains_on_map(ecotypeid2pos.keys(), ecotypeid2pos, 'map of strains',  get_pic_area(ecotypeid2pos.values()), strain_map_fname_prefix)
	
	ecotypeid2pos_with_data = get_ecotypeid2pos(curs, 'ecotype', 1)
	draw_strains_on_map(ecotypeid2pos_with_data.keys(), ecotypeid2pos_with_data, 'map of strains with snp data',  get_pic_area(ecotypeid2pos_with_data.values()), '%s_with_data'%strain_map_fname_prefix)
	
	#draw cc of strain_iden_g
	for i in range(len(strain_iden_g_cc)):
		gs = strain_iden_g.subgraph(strain_iden_g_cc[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
		for site, weight in site2weight.iteritems():
			site2weight[site] =  8*math.log(weight+1)
		
		display_matrix_of_component(input_fname, strain_iden_g_cc[i], ecotypeid2pos, output_fname='ecotype_identity_map_cc%s'%i, need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, pic_area=[-130,34,35,65], output_fname_prefix='ecotype_identity_cc%s'%i)
	
	
	#2007-10-02 use clique to represent unique haplotype
	strain_iden_g_cc, clique_ls_indexed_by_cc_number = get_clique_ls_indexed_by_cc_number(strain_iden_g)
	haplotype_size2number = get_haplotype_size2number(clique_ls_indexed_by_cc_number)
	singleton_strain_id_ls = get_singleton_strain_id_ls(strain_iden_g, input_fname)
	haplotype_size2number[1] += len(singleton_strain_id_ls)	#need to make up the singletons lost in graph
	hist_output_fname_prefix = '%s_haplotype_size_bar'%os.path.splitext(input_fname)[0]
	draw_histogram_of_haplotype_size(haplotype_size2number, hist_output_fname_prefix)
	ordered_clique_ls = get_ordered_clique_ls(clique_ls_indexed_by_cc_number)
	for i in range(len(ordered_clique_ls)):
		gs = strain_iden_g.subgraph(ordered_clique_ls[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos, node_weight_by_no_of_nodes=1)
		display_matrix_of_component(input_fname, ordered_clique_ls[i], ecotypeid2pos, output_fname='haplotype_%s_size%s'%(i,len(ordered_clique_ls[i])), need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'haplotype %s size %s'%(i,len(ordered_clique_ls[i])), pic_area=get_pic_area(site2pos.values()), output_fname_prefix='haplotype_%s_size%s_map'%(i,len(ordered_clique_ls[i])))
	
	europe_lon_span = [-12,50]
	norame_lon_span = [-130,-60]
	centralasia_lon_span = [60,90]
	japan_lon_span = [130, 150]
	norame_strain_ls = get_strains_within_longitude_span(ordered_clique_ls[0], ecotypeid2pos, norame_lon_span)
	japan_strain_ls = get_strains_within_longitude_span(ordered_clique_ls[0], ecotypeid2pos, japan_lon_span)
	
	
	cross_atlantic_cc2lon_pair_set = check_cross_ocean_components(strain_iden_g, ordered_clique_ls, ecotypeid2pos, longi_watershed=-30)
	cross_eurasia_cc2lon_pair_set = check_cross_ocean_components(strain_iden_g, ordered_clique_ls, ecotypeid2pos, longi_watershed=70)
	
	
	#2007-10-03 to see whether geographic distance correlates with genotype distance
	import numpy
	strain_pair_ls, geno_dist_ls, geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos)
	
	#store them into cPickle
	import cPickle
	of = open('%s_strain_pair_geno_dist_geo_dist_ls'%os.path.splitext(input_fname)[0], 'w')
	cPickle.dump([strain_pair_ls, geno_dist_ls, geo_dist_ls], of)
	del of
	
	geo_dist_ls = numpy.array(geo_dist_ls)
	geno_dist_ls = numpy.array(geno_dist_ls)
	
	output_fname_prefix = '%s'%os.path.splitext(input_fname)[0]
	geo_output_fname_prefix = '%s_geo_distance_hist'%output_fname_prefix
	plotHist(geo_dist_ls, title="Histogram of non-NA geographic distances", output_fname_prefix=geo_output_fname_prefix)
	
	import random
	index_ls = range(len(geo_dist_ls))
	index_selected_ls = random.sample(index_ls, 10000)
	geo_dist_selected_ls = geo_dist_ls[index_selected_ls]
	geno_dist_selected_ls = geno_dist_ls[index_selected_ls]
	
	geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist'%output_fname_prefix
	plotXY(geo_dist_selected_ls, geno_dist_selected_ls, title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=geno_vs_geo_output_fname_prefix)
	plot_LD(geo_dist_selected_ls, geno_dist_selected_ls,  title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=geno_vs_geo_output_fname_prefix)
	
	#2007-10-04 divide data into continents to see whether geographic distance correlates with genotype distance
	europe_lon_span = [-12,50]
	norame_lon_span = [-130,-60]
	centralasia_lon_span = [60,90]
	japan_lon_span = [130, 150]
	def sample_geno_geo_correlation(geno_dist_ls, geo_dist_ls, output_fname_prefix):
		import numpy
		geo_dist_ls = numpy.array(geo_dist_ls)
		geno_dist_ls = numpy.array(geno_dist_ls)
		import random
		index_ls = range(len(geo_dist_ls))
		index_selected_ls = random.sample(index_ls, 10000)
		geo_dist_selected_ls = geo_dist_ls[index_selected_ls]
		geno_dist_selected_ls = geno_dist_ls[index_selected_ls]
		plot_LD(geo_dist_selected_ls, geno_dist_selected_ls,  title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=output_fname_prefix)
	
	eur_strain_pair_ls, eur_geno_dist_ls, eur_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, europe_lon_span)
	eur_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_eur'%output_fname_prefix
	sample_geno_geo_correlation(eur_geno_dist_ls, eur_geo_dist_ls, eur_geno_vs_geo_output_fname_prefix)
	
	noramer_strain_pair_ls, noramer_geno_dist_ls, noramer_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, norame_lon_span)
	noramer_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_noramer'%output_fname_prefix
	sample_geno_geo_correlation(noramer_geno_dist_ls, noramer_geo_dist_ls, noramer_geno_vs_geo_output_fname_prefix)
	
	eur_noramer_strain_pair_ls, eur_noramer_geno_dist_ls, eur_noramer_geo_dist_ls =  calGenoDistAndGeoDistBetweenTwoAreas(input_fname, ecotypeid2pos, europe_lon_span, norame_lon_span)
	eur_noramer_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_eur_noramer'%output_fname_prefix
	sample_geno_geo_correlation(eur_noramer_geno_dist_ls, eur_noramer_geo_dist_ls, eur_noramer_geno_vs_geo_output_fname_prefix)
	
	#2007-10-01 parition cc into cliques
	from PartitionGraphIntoCliques import PartitionGraphIntoCliques
	PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(0)
	for for i in range(len(strain_iden_g_cc)):
		g_cc = strain_iden_g.subgraph(strain_iden_g_cc[i])
		PartitionGraphIntoCliques_ins.partition(g_cc.copy())
		map(len, PartitionGraphIntoCliques_ins.clique_ls)
		for j in range(len(PartitionGraphIntoCliques_ins.clique_ls)):
			gs = g_cc.subgraph(PartitionGraphIntoCliques_ins.clique_ls[j])
			site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
			for site, weight in site2weight.iteritems():
				site2weight[site] =  8*math.log(weight+1)
			
			display_matrix_of_component(input_fname, PartitionGraphIntoCliques_ins.clique_ls[j], ecotypeid2pos, output_fname='ecotype_identity_map_cc%sc%s'%(i,j), need_sort=0, need_savefig=1)
			draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s clique%s'%(i,j), pic_area=[-130,34,120,65], output_fname_prefix='ecotype_identity_cc%sc%s'%(i,j))
	
	#2008-01-25 to investigate how identity spreads against distance in different regions.
	europe_lon_span = [-12,50]
	norame_lon_span = [-130,-60]
	centralasia_lon_span = [60,90]
	japan_lon_span = [130, 150]
	
	span_ls = [europe_lon_span, norame_lon_span]
	span_label_ls = ['europe', 'nor america']
	
	dist_range = [0.0, 0.4]
	strain_pair_ls, geno_dist_ls, geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, dist_range=dist_range)
	
	def float2str(input):
		input = repr(input)
		input = input.replace('.', '_')
		return input
	
	import cPickle
	output_fname = '%s_strain_pair_geno_dist_%s_%s_geo_dist_ls'%(os.path.splitext(input_fname)[0], dist_range[0], dist_range[1])
	output_fname = output_fname.replace('.', '_')
	of = open(output_fname, 'w')
	cPickle.dump([strain_pair_ls, geno_dist_ls, geo_dist_ls], of)
	del of
	
	matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls = group_identity_pair_according_to_region_distance(strain_pair_ls, geno_dist_ls, geo_dist_ls, ecotypeid2pos, span_ls, span_label_ls, geo_distance_window_size=100, max_identity_geno_dist=0.0)
	draw_table_as_bar_chart(matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls)
	"""
	
	"""
	2007-09-24
		magnus's idea to reduce fake identities
	"""
	def get_data_matrix_sorted_by_NA_perc(input_fname):
		"""
		2007-09-24
			input_fname is a data frame file with snps coded in integers. 1st row is the ecotype id in table ecotype
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		no_of_strains = len(strain_acc_list)
		no_of_snps = len(header)-2
		data_matrix_sorted_by_NA_perc = []
		for i in range(no_of_strains):
			no_of_NAs = 0
			for j in range(no_of_snps):
				if data_matrix[i][j] == 0:
					no_of_NAs += 1
			data_matrix_sorted_by_NA_perc.append([float(no_of_NAs)/no_of_snps, data_matrix[i], strain_acc_list[i]])
		data_matrix_sorted_by_NA_perc.sort()
		return data_matrix_sorted_by_NA_perc
	
	def get_strain_id2index_from_data_matrix_sorted_by_NA_perc(data_matrix_sorted_by_NA_perc):
		sys.stderr.write("Getting strain_id2index from data_matrix_sorted_by_NA_perc ...")
		strain_id2index = {}
		no_of_strains = len(data_matrix_sorted_by_NA_perc)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		for i in range(no_of_strains):
			ecotypeid = data_matrix_sorted_by_NA_perc[i][2]
			strain_id2index[ecotypeid] = i
		sys.stderr.write("Done.\n")
		return strain_id2index
	
	
	def get_unique_haplotype_set(data_matrix_sorted_by_NA_perc):
		sys.stderr.write("Getting unique haplotype set ...")
		no_of_strains = len(data_matrix_sorted_by_NA_perc)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		from sets import Set
		unique_haplotype_set = Set()
		for i in range(no_of_strains):
			no_of_prev_identity_strains = 0
			for j in range(i):
				no_of_same_cols = 0
				for k in range(no_of_snps):
					call1 = data_matrix_sorted_by_NA_perc[i][1][k]
					call2 = data_matrix_sorted_by_NA_perc[j][1][k]
					if call1 == call2 or call1==0 or call2==0:
						no_of_same_cols += 1
				if no_of_same_cols == no_of_snps:
					no_of_prev_identity_strains += 1
			if no_of_prev_identity_strains==0:
				unique_haplotype_set.add(data_matrix_sorted_by_NA_perc[i][2])
		sys.stderr.write("Done.\n")
		return unique_haplotype_set
	
	
	def assign_remainder_strains_to_unique_haplotype_set(unique_haplotype_set, data_matrix_sorted_by_NA_perc, strain_id2index):
		sys.stderr.write("Assigning remainder strains to unique_haplotype_set ...")
		identity_pair_ls = []
		no_of_strains = len(data_matrix_sorted_by_NA_perc)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		for i in range(no_of_strains):
			ecotypeid = data_matrix_sorted_by_NA_perc[i][2]
			if ecotypeid not in unique_haplotype_set:	#skip strains in unique_haplotype_set
				prev_identity_strains = []
				for seedid in unique_haplotype_set:
					seed_data_row = data_matrix_sorted_by_NA_perc[strain_id2index[seedid]][1]
					no_of_same_cols = 0
					for k in range(no_of_snps):
						call1 = data_matrix_sorted_by_NA_perc[i][1][k]
						call2 = seed_data_row[k]
						if call1 == call2 or call1==0 or call2==0:
							no_of_same_cols += 1
					if no_of_same_cols == no_of_snps:
						prev_identity_strains.append(seedid)
				if len(prev_identity_strains)==1:	#uniquely assigned
					identity_pair_ls.append([ecotypeid, prev_identity_strains[0]])
		sys.stderr.write("Done.\n")
		return identity_pair_ls
	
	def construct_identity_pair_ls_given_strain_id_list(strain_id_list, data_matrix_sorted_by_NA_perc, strain_id2index):
		sys.stderr.write("Constructing identity_pair_ls given strain_id_list ...")
		identity_pair_ls = []
		no_of_strains = len(strain_id_list)
		no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
		for i in range(no_of_strains):
			ecotypeid1 = repr(strain_id_list[i])	#convert int to char
			index_of_ecotypeid1 = strain_id2index[ecotypeid1]
			for j in range(i+1, no_of_strains):
				ecotypeid2 = repr(strain_id_list[j])
				index_of_ecotypeid2 = strain_id2index[ecotypeid2]
				no_of_same_cols = 0
				for k in range(no_of_snps):
					call1 = data_matrix_sorted_by_NA_perc[index_of_ecotypeid1][1][k]
					call2 = data_matrix_sorted_by_NA_perc[index_of_ecotypeid2][1][k]
					if call1 == call2 or call1==0 or call2==0:
						no_of_same_cols += 1
				if no_of_same_cols == no_of_snps:
					identity_pair_ls.append([ecotypeid1, ecotypeid2])
		sys.stderr.write("Done.\n")
		return identity_pair_ls
	
	
	"""
	input_fname = 'script/variation/stock20070919/data_d110_c0_5.tsv'
	data_matrix_sorted_by_NA_perc = get_data_matrix_sorted_by_NA_perc(input_fname)
	
	unique_haplotype_set = get_unique_haplotype_set(data_matrix_sorted_by_NA_perc)
	
	strain_id2index = get_strain_id2index_from_data_matrix_sorted_by_NA_perc(data_matrix_sorted_by_NA_perc)
	
	identity_pair_ls2 = assign_remainder_strains_to_unique_haplotype_set(unique_haplotype_set, data_matrix_sorted_by_NA_perc, strain_id2index)
	
	import cPickle
	of = open('%s_identity_pair_ls2'%os.path.splitext(input_fname)[0], 'w')
	cPickle.dump(identity_pair_ls2, of)
	del of
	popid2ecotypeid_table = 'popid2ecotypeid_25'
	pop_iden_g2, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls2)
	
	popid2pos = {}
	for popid, pos_size in popid2pos_size.iteritems():
		popid2pos[popid] = pos_size[0]
	
	import math
	popid2weight = {}
	for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
		popid2weight[popid] = 8*math.log(intra_pop_connections+1)
	
	draw_graph_on_map(pop_iden_g2, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map2')
	
	
	strain_iden_g2 = construct_graph_out_of_identity_pair(identity_pair_ls2)
	
	ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
	ecotypeid2weight = {}
	for n in strain_iden_g:
		ecotypeid2weight[n] = 1
	
	import networkx as nx
	strain_iden_g_cc2 = nx.connected_components(strain_iden_g2)
	
	i=11
	
	
	for i in range(len(strain_iden_g_cc2)):
		gs = strain_iden_g2.subgraph(strain_iden_g_cc2[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
		for site, weight in site2weight.iteritems():
			site2weight[site] =  8*math.log(weight+1)
		
		display_matrix_of_component(input_fname, strain_iden_g_cc2[i], ecotypeid2pos, output_fname='ecotype_identity_tcc%s'%i, need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, output_fname_prefix='ecotype_identity_tcc%s'%i)
	
	
	#2007-09-25 complete the identity graph
	identity_pair_ls2_complete = construct_identity_pair_ls_given_strain_id_list(strain_iden_g2.nodes(), data_matrix_sorted_by_NA_perc, strain_id2index)
	
	pop_iden_g2_complete, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls2_complete)
	
	popid2pos = {}
	for popid, pos_size in popid2pos_size.iteritems():
		popid2pos[popid] = pos_size[0]
	
	import math
	popid2weight = {}
	for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
		popid2weight[popid] = 8*math.log(intra_pop_connections+1)
	
	draw_graph_on_map(pop_iden_g2_complete, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map2_complete')
	
	strain_iden_g2_complete = construct_graph_out_of_identity_pair(identity_pair_ls2_complete)
	strain_iden_g_cc2_complete = nx.connected_components(strain_iden_g2_complete)
	
	for i in range(len(strain_iden_g_cc2_complete)):
		gs = strain_iden_g2_complete.subgraph(strain_iden_g_cc2_complete[i])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
		for site, weight in site2weight.iteritems():
			site2weight[site] =  8*math.log(weight+1)
		
		display_matrix_of_component(input_fname, strain_iden_g_cc2_complete[i], ecotypeid2pos, output_fname='ecotype_identity_cc2_complete%s'%i, need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, output_fname_prefix='ecotype_identity_cc2_complete%s'%i)
	
	
	cross_ocean_cc_set2_atlantic = check_cross_ocean_components(strain_iden_g2_complete, strain_iden_g_cc2_complete, ecotypeid2pos)
	
	cross_ocean_cc_set2_eurasia = check_cross_ocean_components(strain_iden_g2_complete, strain_iden_g_cc2_complete, ecotypeid2pos, 60)
	"""



class LD:
	"""
	2007-09-26
		calculate LD (r^2)
	2008-01-15
		put everything into class
	"""
	def __init__(self):
		pass
	
	def group_ordered_snps_into_chr_snp_2layer_ls(self, curs, snp_acc_list, snp_locus_table='snps'):
		"""
		2007-09-26
			assume snps are already in chromosome, position order, just need to find out
			where to stop the chromosome
		"""
		sys.stderr.write("Grouping ordered snps into chr_snp_2layer_ls ...")
		old_chromosome = -1
		chr_snp_2layer_ls = []
		snp_position_ls = []	#no chromosome, just position.
		chr_ls = []
		for i in range(len(snp_acc_list)):
			curs.execute("select chromosome, position from %s where snpid='%s'"%(snp_locus_table, snp_acc_list[i]))
			rows = curs.fetchall()
			chromosome, position = rows[0]
			snp_position_ls.append(position)
			if old_chromosome == -1:	#1st encounter
				old_chromosome = chromosome
			elif chromosome!=old_chromosome:
				chr_snp_2layer_ls.append(chr_ls)
				chr_ls = []
				old_chromosome = chromosome
			chr_ls.append(i)
		chr_snp_2layer_ls.append(chr_ls)
		sys.stderr.write("Done.\n")
		return chr_snp_2layer_ls, snp_position_ls
	
	def fill_in_snp_allele2index(self, diploid_allele, allele2index):
		from common import number2nt, nt2number
		if diploid_allele>4:
			nt = number2nt[diploid_allele]
			allele1 = nt2number[nt[0]]
			allele2 = nt2number[nt[1]]
		else:
			allele1 = allele2 = diploid_allele
		if allele1 not in allele2index:
			allele2index[allele1] = len(allele2index)
		if allele2 not in allele2index:
			allele2index[allele2] = len(allele2index)
		return allele1, allele2
	
	def calculate_LD(self, input_fname, curs, snp_locus_table='snps', debug=0, check_bit_ls=[1,0,0,0,0,0,0,0,0,0], chr_ls=[]):
		"""
		exclude pairs with one or two NAs
		exclude pairs both of who are heterozygous calls (can't figure out the phase)
		(only one of pairs is heterozygous is all right)
		2008-01-23 add chr_ls
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		if debug:
			import pdb
			pdb.set_trace()
		import random
		no_of_strains = len(strain_acc_list)
		snp_acc_list = header[2:]
		chr_snp_2layer_ls, snp_position_ls = self.group_ordered_snps_into_chr_snp_2layer_ls(curs, snp_acc_list, snp_locus_table)
		no_of_chrs = len(chr_snp_2layer_ls)
		if len(chr_ls)==0:	#2008-01-23 if there's nothing in chr_ls
			chr_ls = range(no_of_chrs)
		r_square_ls = []
		distance_ls = []
		allele_freq_ls = []
		snp_pair_ls = []
		D_ls = []
		D_prime_ls = []
		import Numeric
		counter = 0
		period = len(check_bit_ls)
		for c in chr_ls:
			no_of_snps = len(chr_snp_2layer_ls[c])
			for i in range(no_of_snps):
				if len(check_bit_ls)>1:	#2008-01-23 no shuffling if there's only 1 bit
					random.shuffle(check_bit_ls)
				for j in range(i+1, no_of_snps):
					snp1_index = chr_snp_2layer_ls[c][i]
					snp2_index = chr_snp_2layer_ls[c][j]
					distance = abs(snp_position_ls[snp1_index]-snp_position_ls[snp2_index])
					if distance>10000:	#2008-01-15 skip SNPs two far
						continue
					counter += 1
					if check_bit_ls[counter%period]==0:	#2008-01-15 skip some portion of them randomly
						continue
					counter_matrix = Numeric.zeros([2,2])
					snp1_allele2index = {}
					snp2_allele2index = {}
					for k in range(no_of_strains):
						snp1_allele = data_matrix[k][snp1_index]
						snp2_allele = data_matrix[k][snp2_index]
						if snp1_allele!=0 and snp2_allele!=0 and not (snp1_allele>4 and snp2_allele>4):
							snp1_allele1, snp1_allele2 = self.fill_in_snp_allele2index(snp1_allele, snp1_allele2index)
							snp2_allele1, snp2_allele2 = self.fill_in_snp_allele2index(snp2_allele, snp2_allele2index)
							counter_matrix[snp1_allele2index[snp1_allele1],snp2_allele2index[snp2_allele1]] += 1
							counter_matrix[snp1_allele2index[snp1_allele2],snp2_allele2index[snp2_allele2]] += 1
					PA = sum(counter_matrix[0,:])
					Pa = sum(counter_matrix[1,:])
					PB = sum(counter_matrix[:,0])
					Pb = sum(counter_matrix[:,1])
					total_num = float(PA+Pa)
					try:
						PA = PA/total_num
						Pa = Pa/total_num
						PB = PB/total_num
						Pb = Pb/total_num
						PAB = counter_matrix[0,0]/total_num
						D = PAB-PA*PB
						PAPB = PA*PB
						PAPb = PA*Pb
						PaPB = Pa*PB
						PaPb = Pa*Pb
						Dmin = max(-PAPB, -PaPb)
						Dmax = min(PAPb, PaPB)
						if D<0:
							D_prime = D/Dmin
						else:
							D_prime = D/Dmax
						r2 = D*D/(PA*Pa*PB*Pb)
					except:	#2008-01-23 exceptions.ZeroDivisionError, Dmin or Dmax could be 0 if one of(-PAPB, -PaPb)  is >0 or <0
						sys.stderr.write('Unknown except, ignore: %s\n'%repr(sys.exc_info()[0]))
						continue
					allele_freq = (min(PA, Pa),min(PB, Pb))
					D_ls.append(D)
					D_prime_ls.append(D_prime)
					r_square_ls.append(r2)
					distance_ls.append(distance)
					allele_freq_ls.append(allele_freq)
					snp_pair_ls.append((snp_acc_list[i], snp_acc_list[j]))
		return D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls
	
	
	def plot_LD(self, x_ls, y_ls, title, xlabel, ylabel, max_dist=0, output_fname_prefix=None):
		import pylab
		import numpy
		import scipy.interpolate
		#x_ls = numpy.array(x_ls)
		#y_ls = numpy.array(y_ls)
		x_argsort_ls = numpy.argsort(x_ls)
		new_x_ls = []
		new_y_ls = []
		for i in range(len(x_argsort_ls)):
			if max_dist>0:
				if x_ls[x_argsort_ls[i]]<=max_dist:
					new_x_ls.append(x_ls[x_argsort_ls[i]])
					new_y_ls.append(y_ls[x_argsort_ls[i]])
			else:
				new_x_ls.append(x_ls[x_argsort_ls[i]])
				new_y_ls.append(y_ls[x_argsort_ls[i]])
		
		sp = scipy.interpolate.UnivariateSpline(new_x_ls,new_y_ls)
		step = (new_x_ls[-1]-new_x_ls[0])/100
		n_x_ls = numpy.arange(new_x_ls[0], new_x_ls[-1], step)
		n_y_ls = map(sp, n_x_ls)
		pylab.clf()
		pylab.title(title)
		pylab.xlabel(xlabel)
		pylab.ylabel(ylabel)
		pylab.plot(new_x_ls, new_y_ls, '.')
		pylab.plot(n_x_ls, n_y_ls)
		if output_fname_prefix:
			if max_dist:
				output_fname_prefix = '%s_%s'%(output_fname_prefix, max_dist)
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=450)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=450)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=450)
	
	def get_snp_acc2LD_data(self, LD_ls, distance_ls, snp_pair_ls):
		"""
		2008-01-24
			the previus model doesn't give good results.
			curve and real data are quite far apart. it's easily disrupted by points with low LD. axis distance=0 (the separating line) sometimes sits in the middle, like 5kb.
			so try Alex Platt's suggestion:
				r^2 = ax^b => log(r^2) = log(a) + b*log(x)

		2008-01-23
			process data for linear model fitting
			r^2 = 1/(a+bx) (page 483 of Hartl2007) => 1/r^2 = a+bx
		"""
		sys.stderr.write("Getting snp_acc2LD_data ...")
		snp_acc2LD_data = {}
		import math
		for i in range(len(snp_pair_ls)):
			snp_pair = snp_pair_ls[i]
			for snp_acc in snp_pair:
				if snp_acc not in snp_acc2LD_data:
					snp_acc2LD_data[snp_acc] = [[], []]
				if LD_ls[i]>0 and distance_ls[i]>0:
					snp_acc2LD_data[snp_acc][0].append(math.log(LD_ls[i]))
					snp_acc2LD_data[snp_acc][1].append([1, math.log(distance_ls[i])])
		sys.stderr.write("Done.\n")
		return snp_acc2LD_data
	
	def LD_linear_fitting(self, snp_acc2LD_data):
		"""
		2008-01-23
			actual linear fitting
		"""
		sys.stderr.write("Linear fitting for LD ...")
		from annot.bin.module_cc.linear_model import linear_model
		linear_model_ins = linear_model()
		snp_acc2LD_linear_fitting = {}
		for snp_acc, LD_data in snp_acc2LD_data.iteritems():
			LD_ls, X_ls = LD_data
			if len(LD_ls)>5:
				linear_model_ins.prepare_data(LD_ls, X_ls)
				linear_model_ins.run()
				coeff_list = linear_model_ins.coefficients()
				chisq_tuple = linear_model_ins.chisq_return()
				snp_acc2LD_linear_fitting[snp_acc] = [coeff_list, chisq_tuple]
				linear_model_ins.cleanup()
		del linear_model_ins
		sys.stderr.write("Done.\n")
		return snp_acc2LD_linear_fitting
	
	def check_one_snp_acc_LD_decay(self, snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc):
		"""
		2008-01-24
			follow the model change in get_snp_acc2LD_data()
		2008-01-23
			unfold the data processed in snp_acc2LD_data
			to see whether the theoretical curve (blue line) fits real data (red dots).
		"""
		import math
		a = snp_acc2LD_linear_fitting[snp_acc][0][0]
		b = snp_acc2LD_linear_fitting[snp_acc][0][1]
		theoretical_curve_func = lambda x: math.exp(a)*math.pow(x,b)
		x_ls = range(500,10000,100)
		y_ls = map(theoretical_curve_func, x_ls)
		import pylab
		pylab.clf()
		pylab.plot(x_ls, y_ls)
		
		inverse_func = lambda x: math.exp(x)
		r2_ls = map(inverse_func, snp_acc2LD_data[snp_acc][0])
		distance_ls = [math.exp(row[1]) for row in snp_acc2LD_data[snp_acc][1]]
		pylab.plot(distance_ls, r2_ls, 'r.')
		pylab.show()
	
	"""
from variation.src.misc import LD
input_fname = os.path.expanduser('~/script/variation/stock20070919/data_d110_c0_5.tsv')
snp_locus_table = 'snps'
LD_ins = LD()
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table)

title = 'LD decay'
xlabel = 'Distance'
ylabel = r'$r^2$'
output_fname_prefix = '%s_LD_r2'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = 'D'
output_fname_prefix = '%s_LD_D'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = "D'"
output_fname_prefix = '%s_LD_D_prime'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)


#2008-01-15 250k data
input_fname = os.path.expanduser('~/script/variation/genotyping/250ksnp/data/data_250k.tsv')
snp_locus_table = 'snps_250k'
LD_ins = LD()
#2008-01-23 chromosome 1 and no random skipping
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table, check_bit_ls=[1], chr_ls=[1])
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table, check_bit_ls=[1]+range(100), chr_ls=[2])

pickle_fname = os.path.expanduser('~/.pickle/LD_ins')
of = open(pickle_fname, 'w')
cPickle.dump([D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls], of)
del of

#2008-01-23 load the computed data by cPickle
import cPickle
inf = open(pickle_fname, 'r')
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = cPickle.load(inf)
del inf

snp_acc2LD_data = LD_ins.get_snp_acc2LD_data(r_square_ls, distance_ls, snp_pair_ls)
snp_acc2LD_linear_fitting = LD_ins.LD_linear_fitting(snp_acc2LD_data)

def check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc):
	import math
	a = snp_acc2LD_linear_fitting[snp_acc][0][0]
	b = snp_acc2LD_linear_fitting[snp_acc][0][1]
	theoretical_curve_func = lambda x: math.exp(a)*math.pow(x,b)
	x_ls = range(500,10000,100)
	y_ls = map(theoretical_curve_func, x_ls)
	import pylab
	pylab.clf()
	pylab.plot(x_ls, y_ls)
	
	inverse_func = lambda x: math.exp(x)
	r2_ls = map(inverse_func, snp_acc2LD_data[snp_acc][0])
	distance_ls = [math.exp(row[1]) for row in snp_acc2LD_data[snp_acc][1]]
	pylab.plot(distance_ls, r2_ls, 'r.')
	pylab.show()


LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting,'1_10219_T_A')
LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting,'1_3102_A_G')

for snp_acc in snp_acc2LD_linear_fitting:
	LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc)
	print snp_acc2LD_linear_fitting[snp_acc]
	raw_input(":")
	"""
	
	def plotLDHist(input_fname, output_fname_prefix, which_LD_statistic=1, min_r2=0.95, discard_perc=0.999, no_of_bins=50, \
				   pairType=0, normed=False, debug=False):
		"""
		2009-4-5
			add argument discard_perc which controls the percentage of pairs discarded randomly to reduce computational load
			add argument no_of_bins
			pairType=0, all
			pairType=1, intra-chromosomal
			pairType=2, inter-chromosomal
		02/20/09
			histogram of LD under cutoff by min_r2
		"""
		import csv, random
		from pymodule import getColName2IndexFromHeader, PassingData
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		xylim_data = PassingData(xlim = [-1,-1], ylim=[-1,-1])
		LD_ls = []
		counter = 0
		real_counter = 0
		for row in reader:
			counter += 1
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = tuple(map(int, snp2))
			if pairType==1 and snp1[0]!=snp2[0]:	#intra-chromosomal, but the two snps are from different chromosomes
				continue
			if pairType==2 and snp1[0]==snp2[0]:	#inter-chromosomal, but the two snps are from the same chromosome
				continue
			if discard_perc==0:
				u = 1
			else:	#sample a uniform unless discard_perc is not 0
				u = random.random()
			if u<discard_perc:
				continue
			
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			#if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
			LD_stat = float(row[col_name2index[LD_statistic.get_name(which_LD_statistic)]])
			LD_stat = abs(LD_stat)
			if LD_stat>=min_r2:
				LD_ls.append(LD_stat)
				real_counter+=1
			if debug and counter==50000:
				break
			if real_counter%5000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader
		
		import pylab
		pylab.clf()
		pylab.title('Histogram of %s'%LD_statistic.get_name(which_LD_statistic))
		pylab.ylabel('frequency')
		pylab.hist(LD_ls, no_of_bins, alpha=0.5, normed=normed)
		pylab.savefig('%s_%s_min_r2_%s_discard_perc_%s_hist.png'%\
					  (output_fname_prefix, LD_statistic.get_name(which_LD_statistic), min_r2, discard_perc), dpi=300)
	
	"""
input_fname = os.path.expanduser('~/panfs/250k/dataset/call_method_29_LD_D_prime_m0.8.tsv')
output_fname_prefix = os.path.expanduser('%s'%os.path.splitext(input_fname)[0])
plotLDHist(input_fname, output_fname_prefix, min_r2=0.95)

input_fname = os.path.expanduser('~/panfs/250k/dataset/call_method_29_LD_D_prime_m0_p.997.tsv')
output_fname_prefix = os.path.expanduser('%s_intra_chr'%os.path.splitext(input_fname)[0])
plotLDHist(input_fname, output_fname_prefix, which_LD_statistic=1, min_r2=0, discard_perc=0.994, no_of_bins=200, pairType=1, normed=True)
output_fname_prefix = os.path.expanduser('%s_inter_chr'%os.path.splitext(input_fname)[0])
plotLDHist(input_fname, output_fname_prefix, which_LD_statistic=1, min_r2=0, discard_perc=0.994, no_of_bins=200, pairType=2, normed=True)

	"""
	
	def filterLD(input_fname, output_fname, which_LD_statistic=1, min_r2=0.95, debug=False):
		"""
		2009-3-6
			select entries from input_fname that are above min_r2 into output_fname
		"""
		import csv
		from pymodule import getColName2IndexFromHeader
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter= 0
		LD_ls = []
		counter = 0
		real_counter = 0
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = tuple(map(int, snp2))
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			#if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
			LD_stat = float(row[col_name2index[LD_statistic.get_name(which_LD_statistic)]])
			LD_stat = abs(LD_stat)
			if LD_stat>=min_r2:
				#LD_ls.append(LD_stat)
				writer.writerow(row)
				real_counter+=1
			counter += 1
			if debug and counter==50000:
				break
			if real_counter%2000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader, writer
		
	"""
input_fname = os.path.expanduser('~/panfs/250k/call_method_29_LD_D_prime_m0.8.tsv')
min_r2=0.85
output_fname = '%s_min_r2_%s.tsv'%(os.path.expanduser('%s'%os.path.splitext(input_fname)[0]), min_r2)
filterLD(input_fname, output_fname, min_r2=min_r2)
	
	"""
	
	def getSNPPairDistWithLDAboveMin(input_fname, which_LD_statistic=1, min_r2=0.95, debug=False):
		"""
		2009-3-6
			
		"""
		import csv
		from pymodule import getColName2IndexFromHeader
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter= 0
		LD_ls = []
		counter = 0
		real_counter = 0
		inter_chr_pair_ls = []
		snp_pair_dist_ls = []
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = tuple(map(int, snp2))
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			#if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
			LD_stat = float(row[col_name2index[LD_statistic.get_name(which_LD_statistic)]])
			LD_stat = abs(LD_stat)
			if LD_stat>=min_r2:
				#LD_ls.append(LD_stat)
				real_counter+=1
				if snp1[0]!=snp2[0]:
					inter_chr_pair_ls.append((snp1, snp2))
				else:
					snp_pair_dist = abs(snp1[1]-snp2[1])
					snp_pair_dist_ls.append(snp_pair_dist)
			
			counter += 1
			if debug and counter==50000:
				break
			if real_counter%2000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader
		return inter_chr_pair_ls, snp_pair_dist_ls
	
	"""
input_fname = os.path.expanduser('~/panfs/250k/call_method_29_LD_D_prime_m0.8_min_r2_0.85.tsv')
min_r2 = 0.95
inter_chr_pair_ls, snp_pair_dist_ls = getSNPPairDistWithLDAboveMin(input_fname, min_r2=min_r2)
	"""
	class Counter(object):
		"""
		2009-4-8
			a class for calSNPPairStats()
		"""
		def __init__(self):
			self.distance_ls = []
			self.lowest_freq_ls = []
			self.no_of_intrachr_pairs = 0
			self.no_of_interchr_pairs = 0

	def calSNPPairStats(input_fname, output_fname_prefix, discard_perc=0.999, no_of_bins=50, \
					normed=False, debug=False):
		"""
		2009-4-8
			given the data which records the frequency of each allelic combination from two SNPs
				count how many SNP pairs with 4 combos, how many with 3, how many with 2.
					and among each category, how many are inter-chromosomal, how many are intra.
						for the intra pairs, draw histogram of distance
					also among each category, draw histogram of the frequency of the least-frequent combo
		"""
		import csv, random
		reader = csv.reader(open(input_fname), delimiter='\t')
		counter_dc_idx = 9
		from pymodule import getColName2IndexFromHeader, PassingData
		from variation.src.DrawSNPRegion import DrawSNPRegion, LD_statistic
		reader = csv.reader(open(input_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter= 0
		
		no_of_combos2counter = {}
		for row in reader:
			counter_dict = eval(row[counter_dc_idx])
			no_of_combos = len(counter_dict)
			if no_of_combos not in no_of_combos2counter:
				no_of_combos2counter[no_of_combos] = Counter()
			
			counter += 1
			snp1 = row[col_name2index['snp1_id']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2_id']].split('_')
			snp2 = tuple(map(int, snp2))
			
			if discard_perc==0:
				u = 1
			else:	#sample a uniform unless discard_perc is not 0
				u = random.random()
			
			if snp1[0]!=snp2[0]:	#intra-chromosomal, but the two snps are from different chromosomes
				no_of_combos2counter[no_of_combos].no_of_interchr_pairs += 1
				
			if snp1[0]==snp2[0]:	#inter-chromosomal, but the two snps are from the same chromosome
				no_of_combos2counter[no_of_combos].no_of_intrachr_pairs += 1
				if u>=discard_perc:
					snp_pair_dist = abs(snp1[1]-snp2[1])
					no_of_combos2counter[no_of_combos].distance_ls.append(snp_pair_dist)
			
			if u>=discard_perc:
				lowest_no_between_two_snps = min(counter_dict.values())
				no_of_total_accessions = sum(counter_dict.values())
				lowest_freq = lowest_no_between_two_snps/float(no_of_total_accessions)
				no_of_combos2counter[no_of_combos].lowest_freq_ls.append(lowest_freq)
				real_counter+=1
			if debug and counter==50000:
				break
			if real_counter%5000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
		del reader
		sys.stderr.write("%s%s\t\t%s\n"%('\x08'*40, real_counter, counter))
		
		import pylab
		for no_of_combos, counter_obj in no_of_combos2counter.iteritems():
			no_of_total_pairs = counter_obj.no_of_interchr_pairs + counter_obj.no_of_intrachr_pairs
			perc = counter_obj.no_of_interchr_pairs/float(no_of_total_pairs)
			print '%s/%s=%.5f inter-chromosomal SNP pairs for pairs with %s combos'%(counter_obj.no_of_interchr_pairs, no_of_total_pairs, perc, no_of_combos)
			pylab.clf()
			pylab.title('Histogram of Intra-chromosal Dist For SNPPairs with %s combos'%no_of_combos)
			pylab.ylabel('frequency')
			pylab.hist(counter_obj.distance_ls, no_of_bins, alpha=0.5, normed=normed)
			pylab.savefig('%s_%s_combos_discard_perc_%s_dist_hist.png'%\
					  (output_fname_prefix, no_of_combos, discard_perc), dpi=300)
			pylab.clf()
			pylab.title('Histogram of Frequency of Least-freq Combo For SNPPairs with %s combos'%no_of_combos)
			pylab.ylabel('frequency')
			pylab.hist(counter_obj.lowest_freq_ls, no_of_bins, alpha=0.5, normed=normed)
			pylab.savefig('%s_%s_combos_discard_perc_%s_lowest_freq_hist.png'%\
					  (output_fname_prefix, no_of_combos, discard_perc), dpi=300)
"""
input_fname = os.path.expanduser('~/panfs/250k/InterSNPCount/SNPpair_7_FT_22C.tsv')
output_fname_prefix = os.path.expanduser('%s'%os.path.splitext(input_fname)[0])
calSNPPairStats(input_fname, output_fname_prefix, discard_perc=0.4, no_of_bins=200, normed=True, debug=True)

"""
			

class Trio(object):
	"""
	2007-03-18
		the input_fname is output of MpiTrioAncestryInference.py
		for each strain, it outputs no_of_pairs, no_of_distinct_parents, avg_no_of_jumps
	"""
	def TrioAncestrySummary(input_fname, output_fname):
		import sys, os, csv
		from sets import Set
		sys.stderr.write("Reading trio ancestry information ...")
		reader = csv.reader(open(input_fname), delimiter='\t')
		strain_index2data = {}	#value is list of [no_of_pairs, Set_of_parents, no_of_jump_list]
		for row in reader:
			i,j,k = row[:3]
			no_of_jumps = int(row[-1])
			i = int(i)
			j = int(j)
			k = int(k)
			if k not in strain_index2data:
				strain_index2data[k] = [0, Set(), []]
			strain_index2data[k][0] += 1
			strain_index2data[k][1].add(i)
			strain_index2data[k][1].add(j)
			strain_index2data[k][2].append(no_of_jumps)
		del reader
		sys.stderr.write('done\n')
		
		sys.stderr.write("Outputting summary data ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for strain_index, data in strain_index2data.iteritems():
			no_of_distinct_parents = len(data[1])
			avg_no_of_jumps = sum(data[2])/float(data[0])
			writer.writerow([strain_index, data[0], no_of_distinct_parents, avg_no_of_jumps])
		del writer
		sys.stderr.write('done\n')
	
	"""
	TrioAncestrySummary('./script/variation/data/justin_data_filtered.trio_ances', './script/variation/data/justin_data_filtered.trio_ances.summary')
	"""
	
	"""
	2007-03-18
	"""
	def draw_histogram_of_data_from_specified_column(input_fname, column=-1, no_of_bins=20):
		import pylab, csv
		data_ls = []
		reader = csv.reader(open(input_fname), delimiter='\t')
		for row in reader:
			data_ls.append(float(row[column]))
		pylab.clf()
		pylab.hist(data_ls, no_of_bins)
		pylab.show()
	
	"""
	draw_histogram_of_data_from_specified_column('./script/variation/data/justin_data_filtered.trio_ances.summary', 1)
	"""
	
	"""
	2007-03-19
	"""
	def get_one_ancestry_line(trio_ances_fname, triplet):
		import csv, sys
		sys.stderr.write("Getting one ancestry line...")
		reader = csv.reader(open(trio_ances_fname), delimiter='\t')
		ancestry_ls = None
		for row in reader:
			this_triplet = row[:3]
			this_triplet = map(int, this_triplet)
			if this_triplet[0] == triplet[0] and this_triplet[1]==triplet[1] and this_triplet[2] == triplet[2]:
				ancestry_ls = row[3:-2]	#-2 to discard the last chromosome separator and no_of_jumps
				break
		sys.stderr.write("Done.\n")
		del reader
		return ancestry_ls
	
	def DrawTrioAncestry(data_matrix_fname, trio_ances_fname, triplet=[]):
		"""
		2007-04-17 double the layers to accomodate the heterozygotes
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		if not triplet:
			print 'Error: triplet is not specified. '
			return
		ancestry_ls = get_one_ancestry_line(trio_ances_fname, triplet)
		import Numeric
		sub_data_matrix = Numeric.zeros([12, len(ancestry_ls)])
		from common import nt2number, number2nt
		no_of_chrs = 0
		for i in range(len(ancestry_ls)):
			if ancestry_ls[i]=='||':
				no_of_chrs += 1
			else:
				#first fill up the last 6 encoding lines
				sub_data_matrix[6,i] = 1
				sub_data_matrix[7,i] = 1
				sub_data_matrix[8,i] = 2
				sub_data_matrix[9,i] = 2
				ancestry_bit = int(ancestry_ls[i])
				if ancestry_bit == 2:
					sub_data_matrix[10,i] = 1
					sub_data_matrix[11,i] = 2
				else:
					sub_data_matrix[10,i] = sub_data_matrix[11,i] = ancestry_bit+1	#the original coding is 0 or 1, but 0 is reserved for NA
				for j in range(len(triplet)):	#second deal with the 1st 3 lines with real data
					SNP_call = data_matrix[triplet[j]][i-no_of_chrs]
					if SNP_call>4:
						SNP_call_1 = nt2number[number2nt[SNP_call][0]]
						sub_data_matrix[2*j,i] = SNP_call_1
						SNP_call_2 = nt2number[number2nt[SNP_call][1]]
						sub_data_matrix[2*j+1,i] = SNP_call_2
					else:
						sub_data_matrix[2*j,i] = sub_data_matrix[2*j+1,i] = SNP_call
		import pylab
		pylab.clf()
		pylab.imshow(sub_data_matrix, interpolation='nearest')
		pylab.colorbar()
		ytick_label_ls = []
		for i in range(2):	#repeat the labels
			for j in range(len(triplet)):
				ytick_label_ls.append(strain_acc_list[triplet[j]])	#2007-04-17 double it
				ytick_label_ls.append(strain_acc_list[triplet[j]])
		ytick_loc_ls = range(12)
		ytick_loc_ls.reverse()
		pylab.yticks(ytick_loc_ls, ytick_label_ls)
		pylab.show()
	
	"""
	2007-10-16
	"""
	def get_snp_index2pos(snp_acc_list, curs, snps_table='stock20071008.snps'):
		sys.stderr.write("Getting snp_index2pos ...")
		snp_index2pos = {}
		for i in range(len(snp_acc_list)):
			curs.execute("select chromosome, position from %s where snpid='%s'"%(snps_table, snp_acc_list[i]))
			rows = curs.fetchall()
			chromosome, position = rows[0]
			snp_index2pos[i] = (chromosome, position)
		sys.stderr.write("Done.\n")
		return snp_index2pos
	
	def get_shared_block_ls(row1, row2, snp_index2pos):
		shared_block_ls = []
		shared_block = []
		for i in range(len(row1)):
			if row1[i] == row2[i]:	#same allele
				if i==0 or snp_index2pos[i-1][0]==snp_index2pos[i][0]:	#1st snp or previous snp and this one are on same chromosome
					shared_block.append(i)
				elif shared_block:	#same allele, but previous snp is on different chromosome
					shared_block_ls.append(shared_block)
					shared_block = []	#clear it up
					shared_block.append(i)
			elif shared_block:	#different allele but shared_block is not empty
				shared_block_ls.append(shared_block)	#store the shared_block
				shared_block = []	#clear it up
		if row1[-1]==row2[-1]:	#last snp is same. then the shared_block is not put in yet.
			shared_block_ls.append(shared_block)
		return shared_block_ls
	
	def get_all_pairs_with_max_shared_block_length_ls(strain_acc_list, snp_index2pos, data_matrix):
		sys.stderr.write("Getting all_pairs_with_max_shared_block_length_ls ...")
		all_pairs_with_max_shared_block_length_ls = []
		no_of_strains = len(strain_acc_list)
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				shared_block_ls = get_shared_block_ls(data_matrix[i], data_matrix[j], snp_index2pos)
				shared_block_length_ls = map(len, shared_block_ls)
				shared_block_length_ls.sort()
				all_pairs_with_max_shared_block_length_ls.append([shared_block_length_ls[-1], (i,j)])	#take the maximum shared block length
		all_pairs_with_max_shared_block_length_ls.sort()
		sys.stderr.write("Done.\n")
		return all_pairs_with_max_shared_block_length_ls
	
	def DrawTwoRowAndSharedBlock(data_matrix, strain_index_pair, strain_acc_list, snp_index2pos):
		"""
		2007-10-16
			used to check get_shared_block_ls() is correct
		"""
		import Numeric
		from common import nt2number, number2nt
		data_row1 = []
		data_row2 = []
		shared_block_row = []
		no_of_chrs = 0
		row1 = data_matrix[strain_index_pair[0]]
		row2 = data_matrix[strain_index_pair[1]]
		for i in range(len(row1)):
			#for shared_block_row: -1 is chromosome separator, 1 is same, 0 is different.
			if i!=0 and snp_index2pos[i-1][0]!=snp_index2pos[i][0]:	#a different chromosome
				shared_block_row.append(-1)
				data_row1.append(-1)	#-1 is chromosome separator
				data_row2.append(-1)
			if row1[i]==row2[i]:
				shared_block_row.append(1)
			else:
				shared_block_row.append(0)
			data_row1.append(row1[i])
			data_row2.append(row2[i])
		sub_data_matrix = Numeric.array([data_row1, data_row2, shared_block_row])
		import pylab
		pylab.clf()
		pylab.imshow(sub_data_matrix, interpolation='nearest')
		pylab.colorbar()
		ytick_label_ls = []
		for strain_index in strain_index_pair:
			ytick_label_ls.append(strain_acc_list[strain_index])
		ytick_label_ls.append('block')
		pylab.yticks([2,1,0], ytick_label_ls)	#it's reversed.
		pylab.show()
	
	
	def DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size):
		"""
		2008-02-01
			add comments
		2007-10-16
		"""
		import pylab
		pylab.clf()
		#draw the chromosome separator as green circles
		for chr_id, cumu_size in chr_id2cumu_size.iteritems():
			pylab.plot([cumu_size], [1], 'go')
		#draw the snp first as red sticks
		for snp_index,pos in snp_index2pos.iteritems():
			chr_id, chr_pos = pos
			cumu_chr_pos = chr_id2cumu_size[chr_id-1]+chr_pos
			pylab.plot([cumu_chr_pos], [1], 'r|')
		#draw the blocks as lines crossing those red sticks
		for shared_block in shared_block_ls:
			if len(shared_block)>1:
				starting_snp_index = shared_block[0]
				ending_snp_index = shared_block[-1]
				
				chr_id = snp_index2pos[starting_snp_index][0]
				starting_chr_pos = snp_index2pos[starting_snp_index][1]
				ending_chr_pos = snp_index2pos[ending_snp_index][1]
				#add the cumulative chromosome size
				cumu_starting_chr_pos = chr_id2cumu_size[chr_id-1]+starting_chr_pos
				cumu_ending_chr_pos = chr_id2cumu_size[chr_id-1]+ending_chr_pos
				#draw the block
				pylab.plot([cumu_starting_chr_pos, cumu_ending_chr_pos], [1,1], c='b')
		#chromosome separator on the x axis
		xtick_label_ls = []
		xtick_loc_ls = []
		chr_id_ls = chr_id2cumu_size.keys()
		chr_id_ls.sort()
		for i in range(1,len(chr_id_ls)):
			chr_id = chr_id_ls[i]
			xtick_loc_ls.append(chr_id2cumu_size[chr_id])
			xtick_label_ls.append('%s:%s'%(chr_id, chr_id2size[chr_id]))
		pylab.xticks(xtick_loc_ls, xtick_label_ls)	#it's reversed.
		pylab.show()

	"""
	#2007-10-16 check shared block between pairwise strains
	
	input_fname = 'script/variation/stock20071008/data_d110_c0_5_d001.tsv'
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	snp_acc_list = header[2:]
	snp_index2pos = get_snp_index2pos(snp_acc_list, curs, snps_table='stock20071008.snps')
	all_pairs_with_max_shared_block_length_ls = get_all_pairs_with_max_shared_block_length_ls(strain_acc_list, snp_index2pos, data_matrix)
	
	#2007-10-16 test to draw the diagrams
	from common import get_chr_id2size, get_chr_id2cumu_size
	chr_id2size = get_chr_id2size(curs)
	data_ls = get_chr_id2cumu_size(chr_id2size)
	chr_id2cumu_size, chr_gap = data_ls[:2]
	shared_block_ls = get_shared_block_ls(data_matrix[1], data_matrix[2], snp_index2pos)
	DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size)
	
	"""
	"""
	2007-03-21
		The stat is d(i,j)-|d(i,k)-d(j,k)| for k is a child of i and j
		the bigger this value is, the more divergent between two parents and child is in equal distance to two parents
	"""
	def cal_trio_stat(trio_ances_fname, distance_matrix, output_fname, pic_fname_prefix, need_savefig=0, report=0):
		import csv
		import Numeric
		no_of_strains = distance_matrix.shape[0]
		triplet2trio_stat = {}
		trio_stat_ls = []
		reader = csv.reader(open(trio_ances_fname), delimiter='\t')
		counter = 0
		for row in reader:
			triplet = row[:3]
			triplet = map(int, triplet)
			i,j,k = triplet
			trio_stat = distance_matrix[i,j] - abs(distance_matrix[i,k] - distance_matrix[j,k])
			triplet2trio_stat[tuple(triplet)] = trio_stat
			trio_stat_ls.append(trio_stat)
			counter += 1
			if report and counter%3000 == 0:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
		if report:
			sys.stderr.write("%s%s\n"%('\x08'*20, counter))
		del reader
		sys.stderr.write("outputting trio_stat...")
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		for triplet, trio_stat in triplet2trio_stat.iteritems():
			writer.writerow(list(triplet)+[trio_stat])
		del writer
		sys.stderr.write("Done\n")
		import pylab
		pylab.clf()
		pylab.hist(trio_stat_ls, 20)
		pylab.title("hist of d(i,j)-|d(i,k)-d(j,k)|")
		if need_savefig:
			pylab.savefig('%s_trio_stat_his.eps'%pic_fname_prefix, dpi=300)
			pylab.savefig('%s_trio_stat_hist.svg'%pic_fname_prefix, dpi=300)
			pylab.savefig('%s_trio_stat_hist.png'%pic_fname_prefix, dpi=300)
		pylab.show()
		return triplet2trio_stat
	
	"""
	triplet2trio_stat = cal_trio_stat('./script/variation/data/justin_data_filtered.trio_ances', distance_matrix, './script/variation/data/justin_data_filtered.trio_ances.trio_stat', './script/variation/data/justin_data_filtered.trio_ances', need_savefig=1, report=1)
	"""

class Data2010(object):
	"""
	2007-11-08
		get 192 strains for suzi and draw them on the map
	"""
	def get_ecotypeid_ls_nativename_ls_of_192_strains(curs):
		ecotypeid_ls = []
		nativename_ls = []
		curs.execute("select b.ecotypeid, e.nativename from ecotype e, batch_ecotype b where b.ecotypeid=e.id and b.batchid=2")
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, nativename = row
			ecotypeid_ls.append(ecotypeid)
			nativename_ls.append(nativename)
		return ecotypeid_ls, nativename_ls
	
	
	"""
	ecotypeid_ls, nativename_ls = get_ecotypeid_ls_nativename_ls_of_192_strains(curs)
	ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
	draw_strains_on_map(ecotypeid_ls, ecotypeid2pos, pic_title='192 strains',  pic_area=[-130,10,140,70], output_fname_prefix='/tmp/suzi_192', label_ls=nativename_ls, need_collapse_strains_with_same_pos=0)
	
	"""
	
	@classmethod
	def getFRIAlignment(cls, output_fname, alignment_id=1843):
		"""
		2008-07-31
			output alignment in a format for DrawSNPMatrix.py to draw SNP matrix plot
		"""
		from variation.src.AtDB import AtDB, Sequence, Alignment
		db = AtDB(hostname='localhost')
		rows = Sequence.query.filter_by(alignment=alignment_id).order_by(Sequence.accession).all()
		"""
		#2008-07-31 output in fasta format
		outf = open(output_fname, 'w')
		is_target_alignment_outputted = 0
		for row in rows:
			if not is_target_alignment_outputted:
				outf.write('>ref\n')
				outf.write('%s\n'%row.alignment_obj.target)
				is_target_alignment_outputted = 1
			outf.write('>%s %s\n'%(row.accession, row.accession_obj.name))
			outf.write('%s\n'%(row.bases))
		del outf
		"""
		#2008-08-01 output in a yh SNP matrix format for DrawSNPMatrix.py
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		is_target_alignment_outputted = 0
		for row in rows:
			if not is_target_alignment_outputted:
				header_row = ['name', 1]
				one_row = ['ref', 1]
				for i in range(len(row.alignment_obj.target)):
					base = row.alignment_obj.target[i]
					one_row.append(base)
					header_row.append(i+1)
				writer.writerow(header_row)
				writer.writerow(one_row)
				is_target_alignment_outputted = 1
			one_row = ['%s %s'%(row.accession, row.accession_obj.name), 1]
			for base in row.bases:
				one_row.append(base)
			writer.writerow(one_row)
		del writer
	
	"""
	getFRIAlignment('/tmp/alignment_1843.fasta')
	"""
	
class GWA(object):
	
	
	"""
	2009-4-29
		check the distance between top genes and its most significant nearby SNP
		top genes are partitioned into two categories, candidate and non-candidate 
	"""
	def cmpGeneSNPDistance(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, output_fname_prefix, \
						min_score=3.5, no_of_top_lines=None):
		
		sys.stderr.write("Comparing %s %s %s ...\n"%(call_method_id, phenotype_method_id, analysis_method_id))
		from GeneListRankTest import GeneListRankTest 
		candidate_gene_list = GeneListRankTest.dealWithCandidateGeneList(list_type_id)
		candidate_gene_set = set(candidate_gene_list)
		
		rm = Stock_250kDB.ResultsByGene.results_method
		rbg = Stock_250kDB.ResultsByGene.query.filter(rm.has(call_method_id=call_method_id)).\
				filter(rm.has(phenotype_method_id=phenotype_method_id)).filter(rm.has(analysis_method_id=analysis_method_id)).\
				filter_by(min_distance=40000).filter_by(get_closest=0).first()
		
		result_fname = getattr(rbg, 'filename', None)
		if result_fname is None or not os.path.isfile(result_fname):
			sys.stderr.write("%s doesn't exist.\n"%result_fname)
			return None
		import csv
		from pymodule import getColName2IndexFromHeader
		reader = csv.reader(open(result_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		candidate_gene_snp_dist_ls = []
		non_candidate_gene_snp_dist_ls = []
		counter = 0
		no_of_candidate_genes = 0
		no_of_touches = 0
		for row in reader:
			gene_id = int(row[col_name2index['gene_id']])
			score = float(row[col_name2index['score']])
			disp_pos = int(row[col_name2index['disp_pos']])
			snps_id = int(row[col_name2index['snps_id']])
			snps_context = Stock_250kDB.SnpsContext.query.filter_by(snps_id=snps_id).filter_by(gene_id=gene_id).first()
			if snps_context.disp_pos_comment=='touch':	#touch, distance shall be zero
				disp_pos = 0
				no_of_touches += 1
			if no_of_top_lines is None and score<min_score:
				break
			elif counter>no_of_top_lines:
				break
			if gene_id in candidate_gene_set:
				candidate_gene_snp_dist_ls.append(abs(disp_pos))
				no_of_candidate_genes += 1
			else:
				non_candidate_gene_snp_dist_ls.append(abs(disp_pos))
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, no_of_touches))
		sys.stderr.write("%s\t%s\t%s\n"%('\x08'*40, counter, no_of_touches))
		
		import pylab
		pylab.clf()
		n1 = pylab.hist(candidate_gene_snp_dist_ls, min(max(10, no_of_candidate_genes/10),20), alpha=0.4, normed=1)
		n2 = pylab.hist(non_candidate_gene_snp_dist_ls, 20, alpha=0.4, normed=1, facecolor='r')
		pylab.title("Histogram of SNP-gene dist Top %s %s"%(no_of_top_lines, rbg.short_name))
		pylab.legend([n1[2][0], n2[2][0]], ['%s candidate genes'%no_of_candidate_genes, '%s non-candidate genes'%(counter-no_of_candidate_genes)])
		pylab.xlabel("gene-SNP distance")
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
		
		
	"""
call_method_id = 32
phenotype_method_id = 1
analysis_method_id = 7
list_type_id = 145
no_of_top_lines = 10000
output_fname_prefix = '/tmp/SNP_gene_dist_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_Top%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, no_of_top_lines)
cmpGeneSNPDistance(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, output_fname_prefix, no_of_top_lines=no_of_top_lines)
for phenotype_method_id in range(1,8):
	for analysis_method_id in [1,7]:
		output_fname_prefix = '/tmp/SNP_gene_dist_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_Top%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, no_of_top_lines)
		cmpGeneSNPDistance(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, output_fname_prefix, no_of_top_lines=no_of_top_lines)		
	"""
	
	"""
	2009-4-29
		check the number of SNPs per gene and see whether it's different between candidate and non-candidate 
			top genes are partitioned into two categories, candidate and non-candidate 
	"""
	def checkNoOfSNPsPerGene(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, snps_context_wrapper, output_fname_prefix, \
						min_score=3.5, no_of_top_lines=None):
		sys.stderr.write("Checking %s %s %s ...\n"%(call_method_id, phenotype_method_id, analysis_method_id))
		from GeneListRankTest import GeneListRankTest 
		candidate_gene_list = GeneListRankTest.dealWithCandidateGeneList(list_type_id)
		candidate_gene_set = set(candidate_gene_list)
		
		rm = Stock_250kDB.ResultsMethod.query.filter_by(call_method_id=call_method_id).\
			filter_by(phenotype_method_id=phenotype_method_id).filter_by(analysis_method_id=analysis_method_id).\
			filter_by(results_method_type_id=1).first()
		
		result_fname = getattr(rm, 'filename', None)
		if result_fname is None or not os.path.isfile(result_fname):
			sys.stderr.write("%s doesn't exist.\n"%result_fname)
			return None
		import csv
		from pymodule import PassingData
		genome_wide_result = GeneListRankTest.getResultMethodContent(rm, min_MAF=0.0)
		
		candidate_gene_snp_dist_ls = []
		non_candidate_gene_snp_dist_ls = []
		counter = 0
		no_of_candidate_genes = 0
		no_of_non_candidate_genes = 0
		if genome_wide_result is None:
			return None
		candidate_gene_id2no_of_snps = {}
		non_candidate_gene_id2no_of_snps = {}
		counter = 0
		for data_obj in genome_wide_result.data_obj_ls:
			score = data_obj.value
			snps_context_matrix = snps_context_wrapper.returnGeneLs(data_obj.chromosome, data_obj.position)
			if score<min_score:
				continue
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					dc = candidate_gene_id2no_of_snps
					no_of_candidate_genes += 1
				else:
					dc = non_candidate_gene_id2no_of_snps
					no_of_non_candidate_genes += 1	
				if gene_id not in dc:
					dc[gene_id] = 0
				dc[gene_id] += 1
			counter += 1
			if counter%5000==0:
				sys.stderr.write("%s\t%s\t%s\t%s"%('\x08'*40, counter, no_of_candidate_genes, no_of_non_candidate_genes))
		sys.stderr.write("%s\t%s\t%s\t%s\n"%('\x08'*40, counter, no_of_candidate_genes, no_of_non_candidate_genes))
		import pylab
		pylab.clf()
		
		#recalculate two counters below because previous calculations are not unique gene
		no_of_candidate_genes = len(candidate_gene_id2no_of_snps)
		no_of_non_candidate_genes = len(non_candidate_gene_id2no_of_snps)
		n1 = pylab.hist(candidate_gene_id2no_of_snps.values(), min(max(10, no_of_candidate_genes/20),20), alpha=0.4, normed=1)
		n2 = pylab.hist(non_candidate_gene_id2no_of_snps.values(), 20, alpha=0.4, normed=1, facecolor='r')
		pylab.title("Histogram of #SNPs per gene (score>=%s) %s"%(min_score, rm.short_name))
		pylab.legend([n1[2][0], n2[2][0]], ['%s candidate genes'%no_of_candidate_genes, '%s non-candidate genes'%(no_of_non_candidate_genes)])
		pylab.xlabel("Number of SNPs per gene")
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	
	"""
call_method_id = 32
phenotype_method_id = 1
analysis_method_id = 7
list_type_id = 145
min_score = 3
output_fname_prefix = '/tmp/no_of_SNPs_per_gene_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_min_score_%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, min_score)

from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
min_distance = 40000
get_closest = 0
snps_context_picklef = '/Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m40000'
snps_context_wrapper = GeneListRankTest.dealWithSnpsContextWrapper(snps_context_picklef, min_distance, get_closest)

checkNoOfSNPsPerGene(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, snps_context_wrapper, output_fname_prefix, min_score=min_score)
for phenotype_method_id in range(1,8):
	for analysis_method_id in [1,7]:
		if analysis_method_id==1:
			min_score = 5
		else:
			min_score = 3
		output_fname_prefix = '/tmp/no_of_SNPs_per_gene_call_method_id_%s_phenotype_%s_analysis_%s_list_type_%s_min_score_%s'%\
			(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, min_score)
		checkNoOfSNPsPerGene(call_method_id, phenotype_method_id, analysis_method_id, list_type_id, snps_context_wrapper, output_fname_prefix, min_score=min_score)
	"""
	
	def simulatePvalue(curs, snps_table, output_fname):
		"""
		2008-01-30
			simulate uniform p-values for snps selected from db
		"""
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		curs.execute("select chromosome, position from %s order by chromosome, position"%snps_table)
		rows = curs.fetchall()
		import random
		for row in rows:
			chromosome, position = row
			writer.writerow([chromosome, position, random.expovariate(1)])
		del writer
	
	"""
	snps_table = 'snps_250k'
	output_fname = '/tmp/simulate.pvalue'
	simulatePvalue(curs, snps_table, output_fname)
	"""
	
	def minusLogPvalue(input_fname, output_fname):
		"""
		2008-02-14
			take log of pvalue in input_fname
		"""
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		import math
		for row in reader:
			chromosome, position, pvalue = row
			pvalue = float(pvalue)
			if pvalue>0:
				pvalue = -math.log(float(pvalue))
				writer.writerow([chromosome, position, pvalue])
		del writer, reader
	
	def exponentialMinusLogPvalue(input_fname, output_fname):
		"""
		2008-02-14
			take log of pvalue in input_fname
		"""
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()	#skip the header
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		import math
		for row in reader:
			chromosome, position, pvalue = row
			pvalue = float(pvalue)
			if pvalue>0:
				pvalue = math.exp(-float(pvalue))
				writer.writerow([chromosome, position, pvalue])
		del writer, reader
	
	@classmethod
	def estimateFDRofGWA(cls, results_id, output_fname_prefix, top_p_value_cutoff=1, lambda_gap=0.1):
		"""
		2008-11-04
			try estimating FDR on our genome association results
		"""
		from Stock_250kDB import Stock_250kDB, ResultsMethod
		rm = ResultsMethod.get(results_id)
		from GeneListRankTest import GeneListRankTest
		from pymodule import PassingData
		pd = PassingData(do_log10_transformation=False)
		gwr = GeneListRankTest.getResultMethodContent(rm, min_MAF=0, pdata=pd)
		pvalue_ls = [data_obj.value for data_obj in gwr.data_obj_ls]
		pvalue_ls.sort()
		from transfac.src.AnalyzeTRANSFACHits import AnalyzeTRANSFACHits
		AnalyzeTRANSFACHits_ins = AnalyzeTRANSFACHits()
		AnalyzeTRANSFACHits_ins.remove_top_p_values(pvalue_ls, top_p_value_cutoff)
		figure_fname = '%s_p_value_hist.png'%output_fname_prefix
		AnalyzeTRANSFACHits_ins.draw_pvalue_histogram(pvalue_ls, figure_fname)
		figure_fname = '%s_pi0Tolambda.png'%output_fname_prefix
		lambda_list, pi0_list = AnalyzeTRANSFACHits_ins.calculate_pi0_list(pvalue_ls, figure_fname, \
																		top_p_value_cutoff=top_p_value_cutoff, lambda_gap=lambda_gap)
		
		estimated_pi0 = AnalyzeTRANSFACHits_ins.estimate_pi0(lambda_list, pi0_list)
		
		AnalyzeTRANSFACHits_ins.cal_q_value_list(pvalue_ls, estimated_pi0, top_p_value_cutoff, output_fname_prefix)
		
	"""
	results_id = 2318
	output_fname_prefix = '/tmp/2318_FDR'
	GWA.estimateFDRofGWA(results_id, output_fname_prefix)
	"""
	
	"""
	2008-01-04 upgrade it to output a specified beta or its pvalue
	2008-12-18 take a genome-wide-result file, replace score with beta1 and output into a new file
	"""
	def outputGWABetaPvalue(input_fname, output_fname, pos_index=1, need_beta=True, min_value_cutoff=None, do_log10_transformation=False):
		from pymodule import getGenomeWideResultFromFile
		genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff, do_log10_transformation)
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['chromosome', 'position', 'score', 'MAF', 'MAC', 'genotype_var_perc', 'beta0']
		writer.writerow(header)
		for data_obj in genome_wide_result.data_obj_ls:
			value = None
			if data_obj.comment:
				beta_and_pvalue_ls = data_obj.comment.split(',')	#a list of betas are all separated by ','
				if len(beta_and_pvalue_ls)>pos_index:
					beta, pvalue = beta_and_pvalue_ls[pos_index].split(':')	#beta and its pvalue are separted by ':'
					if need_beta:
						value = abs(float(beta)*10)	#increase it 10 fold to match pvalue, also take absolute value
					else:
						value = abs(float(pvalue))
			if value is not None:
				row = [data_obj.chromosome, data_obj.position, value, data_obj.maf, data_obj.mac, data_obj.genotype_var_perc] + data_obj.extra_col_ls
				writer.writerow(row)
		del writer
	
	"""
	input_fname = '/Network/Data/250k/db/results/type_1/2898_results.tsv'	#LM_with_PC12 on LD
	output_fname = '/tmp/LM_with_PC12_on_LD_beta1.tsv'
	outputGWABetaPvalue(input_fname, output_fname)
	
	input_fname = '/Network/Data/250k/db/results/type_1/2736_results.tsv'	#Emma on LD
	output_fname = '/tmp/Emma_on_LD_beta1.tsv'
	outputGWABetaPvalue(input_fname, output_fname)
	
	#2008-01-04 check the gene-environ interaction pvalue
	input_fname = '/Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6_pheno_2_LD+V_1_LD.tsv'
	output_fname = '/Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6_pvalue_by_int_pheno_2_LD+V_1_LD.tsv'
	outputGWABetaPvalue(input_fname, output_fname, pos_index=3, need_beta=False)
	
	pheno_pair_ls = ['2_LD+V_1_LD', '1_LD_3_SD', '1_LD_4_SD+V', '2_LD+V_3_SD', '2_LD+V_4_SD+V', '4_SD+V_3_SD', '5_FT_10C_6_FT_16C', '5_FT_10C_7_FT_22C', '6_FT_16C_7_FT_22C', '39_JIC0W_42_JIC8W']
	common_input_fname = '/Network/Data/250k/tmp-yh/association_results/lm/call_method_17_y6'
	for pheno_pair in pheno_pair_ls:
		input_fname = '%s_pheno_%s.tsv'%(common_input_fname, pheno_pair)
		output_fname = '%s_pvalue_by_int_pheno_%s.tsv'%(common_input_fname, pheno_pair)
		outputGWABetaPvalue(input_fname, output_fname, pos_index=3, need_beta=False)
	"""
	
	"""
	2008-1-11 draw nicer genome wide plots
	"""
	@classmethod
	def drawGWANicer(cls, db, genome_wide_result, output_fname_prefix, min_value=2.5, need_svg=False):
		chr2xy_ls = {}
		for data_obj in genome_wide_result.data_obj_ls:
			if data_obj.value>=min_value:
				chr = data_obj.chromosome
				if chr not in chr2xy_ls:
					chr2xy_ls[chr] = [[],[]]
				chr2xy_ls[chr][0].append(data_obj.position)
				chr2xy_ls[chr][1].append(data_obj.value)
		
		from variation.src.common import get_chr_id2size, get_chr_id2cumu_size
		chr_id_int2size = get_chr_id2size(db.metadata.bind)
		chr_id2cumu_size, chr_gap, chr_id_ls = get_chr_id2cumu_size(chr_id_int2size, chr_gap=0)
		
		import pylab
		pylab.clf()
		fig = pylab.figure(figsize=(10,2))
		#ax = pylab.axes()
		ax = fig.gca()
		import numpy
		chr_ls = chr2xy_ls.keys()
		chr_ls.sort()
		for chr in chr_ls:
			xy_ls = chr2xy_ls[chr]
			x_ls = numpy.array(xy_ls[0])
			x_ls += chr_id2cumu_size[chr]-chr_id_int2size[chr]
			if xy_ls:
				ax.plot(x_ls, xy_ls[1], '.', markeredgewidth=0, markersize=5, alpha=0.8)
		
		#separate each chromosome
		#for chr in chr_ls[:-1]:
		#	print chr
		#	ax.axvline(chr_id2cumu_size[chr], linestyle='--', color='k', linewidth=0.8)
		
		
		#draw the bonferroni line
		bonferroni_value = -math.log10(0.01/len(genome_wide_result.data_obj_ls))
		ax.axhline(bonferroni_value, linestyle='--', color='k', linewidth=0.8)
		
		#ax.set_ylabel("-log(P-value)")
		#ax.set_xlabel('Chromosomal Position')
		#ax.set_xlim([0, chr_id2cumu_size[chr_ls[-1]]])
		y_lim = ax.get_ylim()
		ax.set_ylim([0, y_lim[1]])
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		if need_svg:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
	
	
	"""
	input_fname = '/Network/Data/250k/db/results/type_1/3018_results.tsv'	#KW on LD, call method 22
	from pymodule import getGenomeWideResultFromFile
	genome_wide_result = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
	output_fname_prefix = '/tmp/call_22_KW_on_LD'
	GWA.drawGWANicer(db_250k, genome_wide_result, output_fname_prefix)
	
	input_fname = '/Network/Data/250k/db/results/type_1/3025_results.tsv'	#Emma on LD, call method 22
	from pymodule import getGenomeWideResultFromFile
	genome_wide_result2 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
	output_fname_prefix = '/tmp/call_22_Emma_on_LD'
	GWA.drawGWANicer(db_250k, genome_wide_result2, output_fname_prefix, min_value=1)
	
	input_fname = '/Network/Data/250k/db/results/type_1/3024_results.tsv'	#KW on FT_22C, call method 22
	from pymodule import getGenomeWideResultFromFile
	genome_wide_result3 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
	output_fname_prefix = '/tmp/call_22_KW_on_7_FT_22C'
	GWA.drawGWANicer(db_250k, genome_wide_result3, output_fname_prefix)
	
	input_fname = '/Network/Data/250k/db/results/type_1/3031_results.tsv'	#Emma on FT_22C, call method 22
	from pymodule import getGenomeWideResultFromFile
	genome_wide_result4 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
	output_fname_prefix = '/tmp/call_22_Emma_on_7_FT_22C'
	GWA.drawGWANicer(db_250k, genome_wide_result4, output_fname_prefix, min_value=1)
	
	input_fname_ls = ['/Network/Data/250k/db/results/type_1/729_results.tsv',\
					'/Network/Data/250k/db/results/type_1/3554_results.tsv',\
					'/Network/Data/250k/db/results/type_1/3881_results.tsv']
	output_fname_prefix_ls = [os.path.expanduser('~/doc/20090930EckerLabVisit/figures/call_7_KW_on_Germ10'),\
							os.path.expanduser('~/doc/20090930EckerLabVisit/figures/call_32_KW_on_Germ10'),\
							os.path.expanduser('~/doc/20090930EckerLabVisit/figures/call_43_KW_on_Germ10')]
	from pymodule import getGenomeWideResultFromFile
	for i in range(len(input_fname_ls)):
		input_fname = input_fname_ls[i]
		output_fname_prefix = output_fname_prefix_ls[i]
		genome_wide_result2 = getGenomeWideResultFromFile(input_fname, min_value_cutoff=None, do_log10_transformation=True)
		GWA.drawGWANicer(db_250k, genome_wide_result2, output_fname_prefix, min_value=1)
	
	"""
	
	"""
	2008-06-27 read in pvalues from a file
	"""
	def plot_maf_vs_pvalue(maf_vector, input_fname, do_log10_transformation=True):
		from GenomeBrowser import GenomeBrowser
		genome_wide_result = GenomeBrowser.getGenomeWideResultFromFile(input_fname, do_log10_transformation=do_log10_transformation)
		pvalue_ls = [genome_wide_result.data_obj_ls[i].value for i in range(len(genome_wide_result.data_obj_ls))]
		import pylab
		pylab.clf()
		pylab.plot(maf_vector, pvalue_ls, '.')
		pylab.show()
	
	"""
pvalue_fname = '/Network/Data/250k/db/results/type_1/394_results.tsv'
plot_maf_vs_pvalue(maf_vector, pvalue_fname)
for i in range(389, 758):
	pvalue_fname = '/Network/Data/250k/db/results/type_1/%s_results.tsv'%i
	plot_maf_vs_pvalue(maf_vector, pvalue_fname)
	"""
	
	
	@classmethod
	def plotCmpEnrichmentOfTwoAnalysisMethods(cls, db, output_fname_prefix):
		"""
		2009-4-17
			plot the enrichment ratios in each cell of the 2X2 partition between two analysis methods across all phenotypes
		"""
		import matplotlib
		import pylab, numpy
		pylab.clf()
		#Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods
		rows = db.metadata.bind.execute("select c.*, r.phenotype_method_id, r.call_method_id from %s c, %s r, %s p where c.results_id1=r.id and \
								p.id=r.phenotype_method_id and p.biology_category_id=1 order by r.phenotype_method_id, c.r1_min_score"%\
								(Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods.table.name, Stock_250kDB.ResultsMethod.table.name,\
								Stock_250kDB.PhenotypeMethod.table.name))
		axe_height = 0.2
		ax1 = pylab.axes([0.1, 0.1, 0.8, axe_height], frameon=False)
		axe_gap = 0.02
		ax2 = pylab.axes([0.1, 0.1+axe_height+axe_gap, 0.8, axe_height], frameon=False, sharey=ax1)
		ax3 = pylab.axes([0.1, 0.1+axe_height*2+axe_gap*2, 0.8, axe_height], frameon=False, sharey=ax1)
		ax4 = pylab.axes([0.1, 0.1+axe_height*3+axe_gap*3, 0.8, axe_height], frameon=False, sharey=ax1)
		non_significant_enrich_ratio_ls = []
		emma_enrich_ratio_ls = []
		kw_enrich_ratio_ls = []
		double_enrich_ratio_ls = []
		xlabel_ls = []
		for row in rows:
			if row.enrichment_ratio is None:
				enrichment_ratio = 0
			else:
				enrichment_ratio = row.enrichment_ratio
			if row.r1_max_score is None and row.r2_max_score is None:
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
				double_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is None and row.r2_max_score is not None:
				emma_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is not None and row.r2_max_score is None:
				kw_enrich_ratio_ls.append(enrichment_ratio)
			else:
				non_significant_enrich_ratio_ls.append(enrichment_ratio)
		N = len(xlabel_ls)
		ind = numpy.arange(N)    # the x locations for the groups
		width = 0.35       # the width of the bars: can also be len(x) sequenc
		p1 = ax1.bar(ind, double_enrich_ratio_ls, width)
		ax1.set_ylabel('double')
		p2 = ax2.bar(ind, emma_enrich_ratio_ls, width)
		ax2.set_ylabel('Emma')
		p3 = ax3.bar(ind, kw_enrich_ratio_ls, width)
		ax3.set_ylabel('KW')
		p4 = ax4.bar(ind, non_significant_enrich_ratio_ls, width)
		
		ax1.set_xticks(ind+width/2.,  xlabel_ls)
		ax1.set_ylim([0,5])
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)

	"""
output_fname_prefix = '/tmp/enrich_ratio'
plotCmpEnrichmentOfTwoAnalysisMethods(db_250k, output_fname_prefix)
	"""
	
	@classmethod
	def generateCommonSQLClause(cls, biology_category_id=1, threshold_type_id=12, type_id=56):
		"""
		2009-5-2
			used by plot3DCmpEnrichmentOfTwoAnalysisMethods() + cmpRankAtThresholdOfTwoAnalysisMethods()
		"""
		common_sql_sentence = "from %s c, %s r, %s p where c.results_id1=r.id and c.type_id=%s and \
							p.id=r.phenotype_method_id and p.biology_category_id=%s and c.threshold_type_id=%s order by \
							r.phenotype_method_id, c.results_id1, c.results_id2, c.r1_min_score"%\
							(Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods.table.name, Stock_250kDB.ResultsMethod.table.name,\
							Stock_250kDB.PhenotypeMethod.table.name, type_id, biology_category_id, threshold_type_id)
		return common_sql_sentence
	
	@classmethod
	def plot3DCmpEnrichmentOfTwoAnalysisMethods(cls, db, biology_category_id=1, threshold_type_id=12, type_id=56, \
											min_non_candidate_sample_size=20):
		"""
		2009-4-17
			3D-version plot of plotCmpEnrichmentOfTwoAnalysisMethods 
		"""
		import matplotlib
		import pylab, numpy
		pylab.clf()
		
		#Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods
		common_sql_sentence = cls.generateCommonSQLClause(biology_category_id, threshold_type_id, type_id)
		
		rows = db.metadata.bind.execute("select distinct r.phenotype_method_id %s"%common_sql_sentence)
		no_of_phenotypes = rows.rowcount*2
		
		rows = db.metadata.bind.execute("select c.*, r.phenotype_method_id, r.call_method_id %s"%common_sql_sentence)

		x, y = numpy.mgrid[0:2*no_of_phenotypes:1, 0:4:1]	#added a gap of 1 column between two phenotypes. one phenotype occupies two rows & two columns.
		
		#remove the gap in x & y
		needed_index_ls = []
		for i in range(0, no_of_phenotypes):
			needed_index_ls.append(2*i)
			#needed_index_ls.append(3*i+1)
			#y[3*i+1][1]=2
		x = x[needed_index_ls]
		y = y[needed_index_ls]
		enrichment_matrix = numpy.ones(x.shape, numpy.float)
		
		non_significant_enrich_ratio_ls = []
		emma_enrich_ratio_ls = []
		kw_enrich_ratio_ls = []
		double_enrich_ratio_ls = []
		xlabel_ls = []
		i = 0	#
		old_results_id_tuple = None
		for row in rows:
			results_id_tuple = (row.results_id1, row.results_id2)
			if old_results_id_tuple == None:
				old_results_id_tuple = results_id_tuple
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
			elif old_results_id_tuple!=results_id_tuple:
				i += 1
				old_results_id_tuple = results_id_tuple
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
			
			if row.enrichment_ratio is None or row.non_candidate_sample_size<=min_non_candidate_sample_size:
				enrichment_ratio = 1
			else:
				enrichment_ratio = row.enrichment_ratio
				
			if row.r1_max_score is None and row.r2_max_score is None:
				enrichment_matrix[i, 3] = enrichment_ratio
			elif row.r1_max_score is None and row.r2_max_score is not None:
				enrichment_matrix[i, 2] = enrichment_ratio
				#emma_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is not None and row.r2_max_score is None:
				enrichment_matrix[i, 1] = enrichment_ratio
				#kw_enrich_ratio_ls.append(enrichment_ratio)
			else:
				enrichment_matrix[i, 0] = enrichment_ratio
				#non_significant_enrich_ratio_ls.append(enrichment_ratio)
		
		from enthought.mayavi import mlab
		mlab.clf()
		bar = mlab.barchart(x, y , enrichment_matrix, lateral_scale=0.8, opacity=1.)
		#mlab.ylabel("KW")
		#mlab.xlabel("Emma")
		#mlab.zlabel("Enrichment Ratio")
		from pymodule.DrawMatrix import get_font 
		font = get_font()
		
		for i in range(len(xlabel_ls)):
			label = xlabel_ls[i]
			char_width, char_height = font.getsize(label)	#W is the the biggest(widest)
			
			mlab.text(2*i, 0, label, z=0, width=char_width/1500.)	#min(0.0075*len(label), 0.04))
		
		s = numpy.ones(x.shape, numpy.int)
		surf = mlab.surf(x, y, s, opacity=0.6, extent=[-1, 2*no_of_phenotypes, -1, 4, 1.0, 1.0])
		
	
	"""
	biology_category_id = 1
	threshold_type_id = 1
	type_id = 57
	min_non_candidate_sample_size = 10
	output_fname_prefix = '/tmp/enrich_ratio_biology_category_id_%s_threshold_type_id_%s_type_id_%s'%(biology_category_id, threshold_type_id, type_id)
	#GWA.cmpRankAtThresholdOfTwoAnalysisMethods(db_250k, output_fname_prefix, biology_category_id, threshold_type_id, type_id)
	GWA.plot3DCmpEnrichmentOfTwoAnalysisMethods(db_250k, biology_category_id, threshold_type_id, type_id, min_non_candidate_sample_size)
	"""
	
	
	@classmethod
	def cmpRankAtThresholdOfTwoAnalysisMethods(cls, db, output_fname_prefix, biology_category_id=1, threshold_type_id=12, type_id=56):
		"""
		2009-5-1
		GWA based on two analysis methods are selected at different cutoffs to do enrichment.
			this function looks at the pvalue rank at that cutoff (=number of SNPs above that cutoff) in different phenotypes.
			to see how different two analysis methods differ in terms of that.
			
			an output table of all values allows further invesitgation of whether the number of SNPs above that cutoff
			has any correlation with enrichment ratio or other correlations.
		"""
		sys.stderr.write("Comparing ...")
		import matplotlib
		import pylab, numpy
		pylab.clf()
		
		#Stock_250kDB.CmpEnrichmentOfTwoAnalysisMethods
		common_sql_sentence = cls.generateCommonSQLClause(biology_category_id, threshold_type_id, type_id)
		
		rows = db.metadata.bind.execute("select distinct r.phenotype_method_id %s"%common_sql_sentence)
		no_of_phenotypes = rows.rowcount
		
		rows = db.metadata.bind.execute("select c.*, r.phenotype_method_id, r.call_method_id %s"%common_sql_sentence)
		
		x, y = numpy.mgrid[0:3*no_of_phenotypes:1, 0:2:1]	#added a gap of 1 column between two phenotypes. one phenotype occupies two rows & two columns.
		
		non_significant_enrich_ratio_ls = []
		emma_enrich_ratio_ls = []
		kw_enrich_ratio_ls = []
		double_enrich_ratio_ls = []
		
		non_significant_sample_size_ls = []
		emma_sample_size_ls = []
		kw_sample_size_ls = []
		double_sample_size_ls = []
		xlabel_ls = []
		i = 0	#
		old_results_id_tuple = None
		r1_threshold = None
		r2_threshold = None
		for row in rows:
			results_id_tuple = (row.results_id1, row.results_id2)
			if old_results_id_tuple == None:
				old_results_id_tuple = results_id_tuple
			elif old_results_id_tuple!=results_id_tuple:
				i += 1
				old_results_id_tuple = results_id_tuple
			if row.candidate_sample_size is None:
				candidate_sample_size = 0
			else:
				candidate_sample_size = row.candidate_sample_size
				
			if row.non_candidate_sample_size is None:
				non_candidate_sample_size = 0
			else:
				non_candidate_sample_size = row.non_candidate_sample_size
			if row.enrichment_ratio is None:
				enrichment_ratio = 1
			else:
				enrichment_ratio = row.enrichment_ratio
			
			#if enrichment_ratio==0:
			#	enrichment_ratio==0.01
			sample_size = candidate_sample_size+non_candidate_sample_size
			if row.r1_max_score is None and row.r2_max_score is None:
				pm = Stock_250kDB.PhenotypeMethod.get(row.phenotype_method_id)
				xlabel_ls.append(pm.short_name)
				double_sample_size_ls.append(sample_size)
				double_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is None and row.r2_max_score is not None:
				emma_sample_size_ls.append(sample_size)
				emma_enrich_ratio_ls.append(enrichment_ratio)
			elif row.r1_max_score is not None and row.r2_max_score is None:
				kw_sample_size_ls.append(sample_size)
				kw_enrich_ratio_ls.append(enrichment_ratio)
			else:
				r1_threshold = row.r1_max_score
				r2_threshold = row.r2_max_score
				non_significant_sample_size_ls.append(sample_size)
				non_significant_enrich_ratio_ls.append(enrichment_ratio)
				#enrichment_matrix[2*i, 0] = enrichment_ratio
		import numpy
		double_sample_size_ls = numpy.array(double_sample_size_ls, numpy.float)
		emma_sample_size_ls = numpy.array(emma_sample_size_ls, numpy.float)
		kw_sample_size_ls = numpy.array(kw_sample_size_ls, numpy.float)
		real_emma_sample_size_ls = emma_sample_size_ls + double_sample_size_ls
		real_kw_sample_size_ls = kw_sample_size_ls + double_sample_size_ls
		real_emma_div_kw_sample_size_ls = real_emma_sample_size_ls/real_kw_sample_size_ls
		emma_div_kw_sample_size_ls = emma_sample_size_ls/kw_sample_size_ls
		emma_div_kw_enrichment_ratio_ls = numpy.array(emma_enrich_ratio_ls)/numpy.array(kw_enrich_ratio_ls)
		
		#output a table for xy plot investigation
		import csv
		writer = csv.writer(open('%s.tsv'%output_fname_prefix, 'w'), delimiter='\t')
		header = ['phenotype_name', 'phenotype_name', \
				'double_sample_size', 'emma_only_sample_size', 'kw_only_sample_size', \
				'emma_div_kw_sample_size', \
				'real_emma_sample_size', 'real_kw_sample_size', 'real_emma_div_kw_sample_size', \
				'double_enrich_ratio', 'emma_enrich_ratio', 'kw_enrich_ratio', 'emma_div_kw_enrichment_ratio']
		writer.writerow(header)
		for i in range(len(xlabel_ls)):
			row = [xlabel_ls[i], xlabel_ls[i], double_sample_size_ls[i], emma_sample_size_ls[i], kw_sample_size_ls[i], emma_div_kw_sample_size_ls[i],\
				real_emma_sample_size_ls[i], real_kw_sample_size_ls[i], real_emma_div_kw_sample_size_ls[i], 
				double_enrich_ratio_ls[i],\
				emma_enrich_ratio_ls[i], kw_enrich_ratio_ls[i], emma_div_kw_enrichment_ratio_ls[i]]
			writer.writerow(row)
		del writer
		import pylab
		pylab.clf()
		n1 = pylab.hist(kw_sample_size_ls, 10, alpha=0.4)
		n2 = pylab.hist(emma_sample_size_ls, 10, alpha=0.4, facecolor='r')
		pylab.legend([n1[2][0], n2[2][0]], ['#SNPs above %s in KW'%r2_threshold, '#SNPs above %s in Emma'%r1_threshold])
		pylab.xlabel("#SNPs above a threshold")
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
biology_category_id = 1
output_fname_prefix = '/tmp/NoOfSNPsAboveThreshold_Flowering'
GWA.cmpRankAtThresholdOfTwoAnalysisMethods(db_250k, output_fname_prefix, biology_category_id=biology_category_id)
	"""
	
	def drawRandomEffectAndResidual(association_output, chromosome=None, position=None):
		"""
		2009-3-24
			association_output is special output by Association.py with y-xb trailling in each row
			get (u+\epsilon) = y-xb, in order to check its distribution
		"""
		import csv, pylab
		reader = csv.reader(open(association_output), delimiter='\t')
		for row in reader:
			chr = int(row[0])
			pos = int(row[1])
			if chromosome is not None and position is not None and (chr!=chromosome or pos!=position):
				continue
			random_effect_residual_ls = row[8:]
			random_effect_residual_ls = map(float, random_effect_residual_ls)
			pylab.clf()
			pylab.title('SNP %s %s mean=%s'%(row[0], row[1], pylab.mean(random_effect_residual_ls)))
			pylab.hist(random_effect_residual_ls, 20)
			pylab.show()
			to_continue = raw_input("Conitnue? (Y/n)")
			if to_continue=='n' or to_continue=='N':
				break
			_chromosome = raw_input("SNP chromosome: (%s)"%(chromosome))
			if _chromosome:
				chromosome = int(_chromosome)
			_position = raw_input("SNP position: (%s)"%(position))
			if _position:
				position = int(_position)
		del reader

	"""
association_output = os.path.expanduser('~/panfs/250k/association_results/call_method_29_y3_pheno_41_JIC4W.tsv')
drawRandomEffectAndResidual(association_output)
	"""
	
	@classmethod
	def contrastPvalueFromTwoGWA(cls, input_fname1, input_fname2, output_fname_prefix, list_type_id=132):
		"""
		2009-7-8
			contrast pvalue from two methods. each dot is a SNP. x-axis is one method. y-axis is the other.
			output the figure to output_fname_prefix
			
			- for input_fname1,  do_log10_transformation=False
			- for input_fname2,  do_log10_transformation=False
		"""
		
		from pymodule import PassingData
		from pymodule.SNP import getGenomeWideResultFromFile
		import Stock_250kDB
		from PlotCmpTwoAnalysisMethods import PlotCmpTwoAnalysisMethods
		from GeneListRankTest import GeneListRankTest, SnpsContextWrapper
		
		param_data = PassingData(need_the_value=1, construct_chr_pos2index=True)
		gwar1 = getGenomeWideResultFromFile(input_fname1, do_log10_transformation=False, pdata=param_data)
		gwar2 = getGenomeWideResultFromFile(input_fname2, do_log10_transformation=False, pdata=param_data)
		
		genome_wide_result_ls = [gwar1, gwar2]
		
		#snps_context_wrapper is not to be used,
		snps_context_picklef = '/Network/Data/250k/tmp-yh/snps_context/snps_context_g0_m5000'
		snps_context_wrapper = PlotCmpTwoAnalysisMethods.dealWithSnpsContextWrapper(snps_context_picklef, min_distance=5000, get_closest=0)
		candidate_gene_set = PlotCmpTwoAnalysisMethods.dealWithCandidateGeneList(list_type_id, return_set=True)
		pvalue_matching_data = PlotCmpTwoAnalysisMethods.matchPvaluesFromTwoResults(genome_wide_result_ls, snps_context_wrapper, \
																				candidate_gene_set)
		import pylab
		pylab.clf()
		"""
		pylab.subplots_adjust(left=0.08, right=0.92,bottom = 0.05)
		#calculate the number of rows needed according to how many score_rank_data, always two-column
		no_of_rows = 1
		no_of_cols = 1
		rm1, rm2 = rm_ls
		ax_scatter_pvalue = pylab.subplot(no_of_rows, no_of_cols, 1, frameon=False)
		"""
		ax = pylab.axes([0.1, 0.1, 0.8,0.8], frameon=False)
		rm1 = Stock_250kDB.ResultsMethod.get(3847)
		rm2 = Stock_250kDB.ResultsMethod.get(3847)
		rm_ls = [rm1, rm2]
		PlotCmpTwoAnalysisMethods.plot_scatter(ax, pvalue_matching_data, rm_ls, data_type=2)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		
	"""
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper	#otherwise,  snps_context_wrapper = cPickle.load(picklef) won't work
input_fname1 = '/Network/Data/250k/db/results/type_1/3847_results.tsv' #RF_Emwa1_32
input_fname2 = '/Network/Data/250k/tmp-atarone/Adnane new analyses/One node 20K trees/RF_nT20K_10_Emwa1.imp'
output_fname_prefix = '/tmp/RF_10k_vs_20k'
GWA.contrastPvalueFromTwoGWA(input_fname1, input_fname2, output_fname_prefix, list_type_id=132)
	"""
	
	@classmethod
	def cofactorLM(cls, genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, cofactors=[],\
					cofactor_phenotype_id_ls=[], report=1, run_type=1,):
		"""
		2010-1-11
			add argument cofactor_phenotype_id_ls to have the functionality of treating some phenotypes as cofactors
		2009-8-26
			one phenotype at a time:
				1. create a new SNP matrix which includes accessions whose phenotypes (this phenotype + cofactor_phenotype_id_ls) are non-NA
				2. one SNP at a time
					1. add cofactor SNP matrix if cofactors are present
					2. add cofactor phenotype matrix if cofactor_phenotype_id_ls exist
					3. run association
			
			parameter run_type is passed to Association.linear_model()
				run_type 1: pure linear model by python
				run_type 2: EMMA
				run_type 3: pure linear model by R (Field test shows run_type 3 is same as 1.)
			
			start_snp and end_snp are on the same chromosome
		"""
		sys.stderr.write("Running association (pure linear model or EMMA) with cofactor ... \n")
		from Association import Association
		
		start_chr, start_pos = start_snp.split('_')[:2]
		start_chr = int(start_chr)
		start_pos = int(start_pos)
		
		stop_chr, stop_pos = stop_snp.split('_')[:2]
		stop_chr = int(stop_chr)
		stop_pos = int(stop_pos)
		
		eigen_vector_fname = ''
		
		test_type = 3 #Emma
		initData = Association.readInData(phenotype_fname, genotype_fname, eigen_vector_fname, phenotype_method_id_ls, test_type=test_type)
		
		which_phenotype_index_ls = initData.which_phenotype_ls
		environment_matrix = None
		data_matrix = initData.snpData.data_matrix
		min_data_point = 3
		
		# 2010-1-11 create a cofactor_phenotype_index_ls
		from PlotGroupOfSNPs import PlotGroupOfSNPs
		cofactor_phenotype_index_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(initData.phenData, set(cofactor_phenotype_id_ls))
		
		#create an index list of cofactor SNPs
		cofactors_indices = []
		cofactors_set = set(cofactors)
		for i in range(len(initData.snpData.col_id_ls)):
			col_id = initData.snpData.col_id_ls[i]
			if col_id in cofactors_set:
				cofactors_indices.append(i)
		
		import numpy, rpy
		rpy.r.source(os.path.expanduser('~/script/variation/src/gwa/emma/R/emma.R'))
		for which_phenotype in which_phenotype_index_ls:
			phenotype_name = initData.phenData.col_id_ls[which_phenotype]
			phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
			output_fname='%s_pheno_%s.tsv'%(os.path.splitext(output_fname_prefix)[0], phenotype_name)	#make up a new name corresponding to this phenotype
			
			#create non NA phenotype
			phenotype_ls = initData.phenData.data_matrix[:, which_phenotype]
			non_phenotype_NA_row_index_ls = []
			non_NA_phenotype_ls = []
			for i in range(len(phenotype_ls)):
				# 2010-1-11 make sure no NA in the cofactor phenotype matrix
				this_row_has_NA_phenotype = False
				if cofactor_phenotype_index_ls:
					for phenotype_index in cofactor_phenotype_index_ls:
						if numpy.isnan(initData.phenData.data_matrix[i, phenotype_index]):
							this_row_has_NA_phenotype = True
							break
				
				if numpy.isnan(phenotype_ls[i]):
					this_row_has_NA_phenotype = True
				if not this_row_has_NA_phenotype:
					non_phenotype_NA_row_index_ls.append(i)
					non_NA_phenotype_ls.append(phenotype_ls[i])
			
			non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
			new_data_matrix = data_matrix[non_phenotype_NA_row_index_ls,:]
			if run_type==2:
				kinship_matrix = Association.get_kinship_matrix(new_data_matrix)
				eig_L = rpy.r.emma_eigen_L(None, kinship_matrix)	#to avoid repeating the computation of eig_L inside emma.REMLE
			else:
				kinship_matrix = None
				eig_L = None

			no_of_rows, no_of_cols = new_data_matrix.shape
			results = []
			counter = 0
			real_counter = 0
			
			#create the cofactor matrix based on new SNP data matrix
			if len(cofactor_phenotype_index_ls)>0:
				cofactor_phenotype_matrix = new_data_matrix[:, cofactor_phenotype_index_ls]
			else:
				cofactor_phenotype_matrix = None
			
			#create the cofactor matrix based on new SNP data matrix
			if len(cofactors_indices)>0:
				cofactor_data_matrix = new_data_matrix[:, cofactors_indices]
			else:
				cofactor_data_matrix = None
			
			#do association
			for j in range(no_of_cols):
				col_id = initData.snpData.col_id_ls[j]
				if col_id not in cofactors_set:	#same SNP appearing twice would cause singular design matrix
					chr, pos = col_id.split('_')[:2]
					chr = int(chr)
					pos = int(pos)
					if chr>=start_chr and chr<=stop_chr and pos>=start_pos and pos<=stop_pos:
						genotype_matrix = new_data_matrix[:,j]
						if cofactor_data_matrix is not None:
							genotype_matrix = genotype_matrix.reshape([len(genotype_matrix), 1])
							genotype_matrix = numpy.hstack((genotype_matrix, cofactor_data_matrix))
						if cofactor_phenotype_matrix is not None:
							if len(genotype_matrix.shape)<2:
								genotype_matrix = genotype_matrix.reshape([len(genotype_matrix), 1])
							genotype_matrix = numpy.hstack((genotype_matrix, cofactor_phenotype_matrix))
						pdata = Association.linear_model(genotype_matrix, non_NA_phenotype_ls, min_data_point, snp_index=j, \
												kinship_matrix=kinship_matrix, eig_L=eig_L, run_type=run_type)
						if pdata is not None:
							results.append(pdata)
							real_counter += 1
				
				counter += 1
				if report and counter%2000==0:
					sys.stderr.write("%s\t%s\t%s"%('\x08'*40, counter, real_counter))
			if report:
				sys.stderr.write("%s\t%s\t%s\n"%('\x08'*40, counter, real_counter))
			
			#output
			Association.output_lm_results(results, initData.snpData.col_id_ls, output_fname)
		
		sys.stderr.write("Done.\n")
	
	"""
	genotype_fname = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_method_id_ls = [285]
	output_fname_prefix = '/tmp/cofactorLM.tsv'
	start_snp = '1_2005921'
	stop_snp = '1_20054408'
	cofactors = ['1_3832974']
	cofactor_phenotype_id_ls = [77]
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
		cofactors=cofactors, cofactor_phenotype_id_ls=cofactor_phenotype_id_ls)
	
	genotype_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_fname = '/tmp/phenotype_43_FLC.tsv'
	phenotype_method_id_ls = [43]
	start_snp = '4_1'
	stop_snp = '4_1750000'
					
	run_type = 2
	cofactors = ['4_268809']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_268809', '4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0], cofactors[1])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	
	run_type = 1
	cofactors = ['4_268809']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	cofactors = ['4_268809', '4_269962']
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s_Cofactor_%s_%s.tsv'%(run_type, start_snp, stop_snp, cofactors[0], cofactors[1])
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=cofactors, run_type=run_type)
	
	run_type =3
	output_fname_prefix = '/tmp/Runtype_%s_SNP_%s_%s.tsv'%(run_type, start_snp, stop_snp)
	GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
					cofactors=[], run_type=run_type)
	
	genotype_fname = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_method_id_ls = [285]
	output_fname_prefix = '/tmp/cofactorLM.tsv'
	start_snp = '1_2005921'
	stop_snp = '1_20054408'
	cofactors = ['1_3832974']
	
	genotype_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv'
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype.tsv'
	phenotype_method_id_ls = [285]
	start_snp = '1_1'
	stop_snp = '5_60000000'
	cofactor_phenotype_id_ls = [77]
	for run_type in [1,2]:
		output_fname_prefix = os.path.expanduser('~/Runtype_%s_Cofactor_Phenotype_%s.tsv'%(run_type, cofactor_phenotype_id_ls[0]))
		GWA.cofactorLM(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, start_snp, stop_snp, \
			cofactors=[], cofactor_phenotype_id_ls=cofactor_phenotype_id_ls, run_type=run_type)
	
	"""
	
	@classmethod
	def phenotypePCA(cls, db, phenotype_fname, output_fname_prefix):
		"""
		2009-9-1
			Run PCA on the phenotype matrix either ecotype-wise or phenotype-wise.
			
			output the data in a matrix fashion that the web MotionChartAppMCPanel app would recognize 
		"""
		from pymodule import read_data, SNPData
		import numpy, csv
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=2, matrix_data_type=float)
		data_matrix_phen = numpy.array(data_matrix_phen)
		phenData = SNPData(header=header_phen, strain_acc_list=strain_acc_list_phen, category_list=category_list_phen, \
						data_matrix=data_matrix_phen)
		
		
		min_data_point = 4
		sys.stderr.write("Removing ecotypes with too few phenotype points (<%s) ..."%min_data_point)
		import pca_module
		rows_to_be_kept = []
		row_labels = []
		for i in range(len(phenData.row_id_ls)):
			row_id = phenData.row_id_ls[i]
			
			no_of_valid_points = sum(numpy.isfinite(phenData.data_matrix[i,:]))
			if no_of_valid_points>=min_data_point:
				rows_to_be_kept.append(row_id)
				row_labels.append('%s %s'%(row_id[1], row_id[0]))
		filteredPhenData = SNPData.keepRowsByRowID(phenData, rows_to_be_kept)
			#order in rows_to_be_kept is kept because it's in the same order as phenData.row_id_ls
		phenData = filteredPhenData
		sys.stderr.write("Done.\n")
		
		phenData_trans = SNPData(row_id_ls=phenData.col_id_ls, col_id_ls=phenData.row_id_ls, \
								data_matrix=numpy.transpose(phenData.data_matrix))
		
		sys.stderr.write("Carrying out phenotype-wise PCA ...")
		# phenotype-wise PCA
		import Stock_250kDB 
		phenotypePCA_fname = '%s_phenotype.tsv'%output_fname_prefix
		phenotypePCA_writer = csv.writer(open(phenotypePCA_fname, 'w'), delimiter='\t')
		
		import pca_module
		from pymodule.PCA import PCA
		#T, P, explained_var = pca_module.PCA_svd(phenData_trans.data_matrix, standardize=True)
		T, P, explained_var = PCA.eig(phenData_trans.data_matrix, normalize=False)	#normalize=True causes missing value in the covariance matrix
		# get the category information for each phenotype
		header = ['phenotype_label', 'category|string', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
		phenotypePCA_writer.writerow(header)
		for i in range(len(phenData_trans.row_id_ls)):
			row_id = phenData_trans.row_id_ls[i]
			row_id_tuple = row_id.split('_')
			phenotype_method_id=int(row_id_tuple[0])
			pm = Stock_250kDB.PhenotypeMethod.get(phenotype_method_id)
			if pm.biology_category:
				category = pm.biology_category.short_name
			else:
				category = 'other'
			data_row = [row_id, category] + list(T[i,0:6])
			phenotypePCA_writer.writerow(data_row)
		del phenotypePCA_writer
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Carrying out ecotype-wise PCA ... \n")
		# ecotype-wise PCA
		
		# run PCA
		#T, P, explained_var = pca_module.PCA_svd(phenData.data_matrix, standardize=True)	#SVD doesn't converge
		T, P, explained_var = PCA.eig(phenData.data_matrix, normalize=False)	#normalize=True gives rise to missing value in the covariance matrix
		from common import getEcotypeInfo
		ecotype_info = getEcotypeInfo(db)
		
		# output
		ecotypePCA_fname = '%s_ecotype.tsv'%output_fname_prefix
		ecotypePCA_writer = csv.writer(open(ecotypePCA_fname, 'w'), delimiter='\t')
		header = ['ecotype_label', 'region|string', 'country|string', 'lat', 'lon', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']
		ecotypePCA_writer.writerow(header)
		for i in range(len(row_labels)):
			row_label = row_labels[i]
			ecotype_id = int(rows_to_be_kept[i][0])
			ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
			if ecotype_obj is not None:
				region = ecotype_obj.region
				country = ecotype_obj.country
				lat = ecotype_obj.latitude
				lon = ecotype_obj.longitude
			else:
				region, country, lat, long = None, None, None, None
			data_row = [row_label, region, country, lat, lon] + list(T[i,0:6])
			ecotypePCA_writer.writerow(data_row)
		del ecotypePCA_writer
		sys.stderr.write("Done.\n")
	
	"""
	phenotype_fname = os.path.expanduser('~/mnt/panfs/250k/phenotype20090902.tsv')
	output_fname_prefix = '/tmp/phenotypePCA'
	GWA.phenotypePCA(db_250k, phenotype_fname, output_fname_prefix)
	
	phenotype_fname = os.path.expanduser('~/mnt/panfs/250k/phenotype20090902_cypress_col_1_9.tsv')	#~/mnt/panfs/250k/phenotype20090902_cypress.tsv')
	output_fname_prefix = '/tmp/phenotypePCA_cypress_col_1_9'
	GWA.phenotypePCA(db_250k, phenotype_fname, output_fname_prefix)
	
	phenotype_fname = os.path.expanduser('~/mnt/panfs/250k/phenotype20090902.tsv')
	output_fname_prefix = '/tmp/phenotypePCA'
	GWA.phenotypePCA(db_250k, phenotype_fname, output_fname_prefix)
	"""
	
	@classmethod
	def plotEnrichmentRatioVsCutoff(cls, db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[1], list_type_id=149, type_id=2, stat_on_x_axis_type=1,\
								call_method_id=32, min_distance=20000):
		"""
		2009-10-22
			Fetch enrichment data from candidate_gene_top_snp_test_rm, plot the enrichment ratio against the cutoff or no of top SNPs.
				for one analysis method, and one/multiple phenotypes
			
			stat_on_x_axis_type=1: min_score
			stat_on_x_axis_type=2: no_of_top_snps
			
			
			The alternative way to get data matrix thru raw sql:
			
			select r.phenotype_method_id, p.short_name, r.analysis_method_id, newt.* from
				(select c.results_id, c.list_type_id, c.pvalue, c.candidate_sample_size, c.candidate_gw_size, 
				c.non_candidate_sample_size, c.non_candidate_gw_size, c.no_of_top_snps, c.max_score, c.min_score 
				from candidate_gene_top_snp_test_rm c where c.type_id=2 and c.results_id in (select id from 
				results_method where call_method_id =32 and phenotype_method_id >=9 and phenotype_method_id <=13 
				and analysis_method_id in (1,6)) and c.list_type_id=149 and c.min_distance=20000) as newt , 
				results_method r, phenotype_method p where r.id=newt.results_id and r.phenotype_method_id=p.id 
				order by phenotype_method_id, analysis_method_id, no_of_top_snps;
		"""
		sys.stderr.write("Plotting enrichment ratio vs cutoff ...")
		from Stock_250kDB import CandidateGeneTopSNPTestRM, ResultsMethod
		import math
		rows = ResultsMethod.query.filter(ResultsMethod.analysis_method_id.in_(analysis_method_id_ls)).\
			filter_by(call_method_id=call_method_id).\
			filter(ResultsMethod.phenotype_method_id.in_(phenotype_method_id_ls))
		results_id_ls = []
		for row in rows:
			results_id_ls.append(row.id)
		
		rows = CandidateGeneTopSNPTestRM.query.filter_by(type_id=type_id).filter_by(list_type_id=list_type_id).\
			filter_by(min_distance=min_distance).\
			filter(CandidateGeneTopSNPTestRM.results_id.in_(results_id_ls)).order_by(CandidateGeneTopSNPTestRM.results_id).\
			order_by(CandidateGeneTopSNPTestRM.no_of_top_snps)
		phenotype_id2name = {}
		phenotype_id2xy_ls = {}
		for row in rows:
			phenotype = row.result.phenotype_method
			phenotype_id = phenotype.id
			if phenotype_id not in phenotype_id2name:
				phenotype_id2name[phenotype_id] = phenotype.short_name
				phenotype_id2xy_ls[phenotype_id] = [[], []]
			candidate_ratio = row.candidate_sample_size/float(row.candidate_gw_size)
			non_candidate_ratio = row.non_candidate_sample_size/float(row.non_candidate_gw_size)
			if non_candidate_ratio >0:
				enrichment_ratio = candidate_ratio/non_candidate_ratio
				if stat_on_x_axis_type==1:
					stat_on_x_axis = row.min_score
					analysis_method = row.result.analysis_method
					if analysis_method.smaller_score_more_significant==1:
						stat_on_x_axis = math.pow(10, -stat_on_x_axis)	# the score is -log(pvalue), need to convert it back
				elif stat_on_x_axis_type==2:
					stat_on_x_axis = row.no_of_top_snps
				phenotype_id2xy_ls[phenotype_id][0].append(stat_on_x_axis)
				phenotype_id2xy_ls[phenotype_id][1].append(enrichment_ratio)
		
		import pylab
		pylab.clf()
		phenotype_id_ls = phenotype_id2xy_ls.keys()
		phenotype_id_ls.sort()
		legend_ls = []
		for phenotype_id in phenotype_id_ls:
			x_ls, y_ls = phenotype_id2xy_ls[phenotype_id]
			phenotype_name = phenotype_id2name[phenotype_id]
			legend_ls.append(phenotype_name)
			pylab.semilogx(x_ls, y_ls)
			
			#pylab.loglog(x_ls, y_ls, basex=10)
		pylab.ylabel("Enrichment Ratio")
		pylab.legend(legend_ls)
		if output_fname_prefix:
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	output_fname_prefix = '/tmp/diseaseResistance_Wilcoxon'
	GWA.plotEnrichmentRatioVsCutoff(db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[1], stat_on_x_axis_type=1)
	
	output_fname_prefix = '/tmp/diseaseResistance_Wilcoxon_no_of_top_SNPs_on_X'
	GWA.plotEnrichmentRatioVsCutoff(db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[1], stat_on_x_axis_type=2)
	
	output_fname_prefix = '/tmp/diseaseResistance_RF'
	GWA.plotEnrichmentRatioVsCutoff(db_250k, output_fname_prefix, phenotype_method_id_ls=[9,10,11,12,13], \
								analysis_method_id_ls=[6], stat_on_x_axis_type=2)
	
	"""
	
	@classmethod
	def predictByRpart(cls, data_matrix, phenotype_ls, col_id_ls=None, output_fname=None):
		"""
		2009-11-18
			call rpart (=randomForest) to predict phenotype_ls based on data_matrix
		"""
				
		import rpy
		from rpy import r
		r.library('rpart')
		rpart_cp = 0.01
		rpy.set_default_mode(rpy.NO_CONVERSION)
			
		pre_data_frame_dict = {}
		for i in range(len(col_id_ls)):
			col_id = col_id_ls[i]
			pre_data_frame_dict[col_id] = rpy.r.as_factor(data_matrix[:,i])	# all SNPs are treated as factor
		
		pre_data_frame_dict["phenotype"] = phenotype_ls
		data_frame = r.as_data_frame(pre_data_frame_dict)
		fit = r.rpart(r("phenotype~."), data=data_frame, method="anova", control=r.rpart_control(cp=rpart_cp))
			#,\
			#	parms=r.list(prior=prior_prob, loss=r.matrix(loss_matrix) ) )
		r.postscript('%s_rsq.ps'%output_fname)
		a=r.rsq_rpart(fit)
		r.dev_off()
		r.postscript('%s_tree.ps'%output_fname)
		a=r.plot(fit)
		r.text(fit, use_n=rpy.r.TRUE)
		r.dev_off()
		r.postscript('%s_cp.ps'%output_fname)
		a=r.plotcp(fit)
		r.dev_off()
	
	@classmethod
	def predictBySVM(cls, data_matrix, phenotype_ls, col_id_ls=None, output_fname=None):
		"""
		2009-11-18
			use SVM (support vector machine) to predict phenotype_ls based on data_matrix
		"""
		from svm import svm_problem, svm_parameter, svm_model, cross_validation, LINEAR, POLY, RBF
		import numpy
		problem = svm_problem(phenotype_ls, data_matrix)
		size = len(phenotype_ls)
		
		total_correct = 0.
		kernels2kname = {LINEAR:'linear', POLY:'polynomial', RBF:'rbf'}
		kernels2kname = {POLY:'polynomial'}	# 2009- 11-18 only polynomial
		import sys
		mute_device = open('/dev/null', 'w')
		
		param = svm_parameter(C = 10)	#,nr_weight = 2,weight_label = [1,0],weight = [10,1])
		C_ls = [100, 10, 1, 0.1, 0.01]	# tried 1000, 10, 1, 0.1 ,0.01
		#C_ls = [100]
		#gamma_ls = [100, 10, 1, 0.1, 0.0001, 0]
		gamma_ls = [0]	# test shows that gamma affects rbf quite a bit
		nr_fold = 10	# 10 fold cross validation
		for k, kname in kernels2kname.iteritems():
			for C in C_ls:
				for gamma in gamma_ls:
					param.gamma = gamma
					param.C = C
					param.kernel_type = k;
					sys.stderr = mute_device
					sys.stdout = mute_device
					target = cross_validation(problem, param, nr_fold)
					#model = svm_model(problem,param)
					sys.stderr = sys.__stderr__
					sys.stdout = sys.__stdout__
					errors = 0.
					for i in range(size):
							#prediction = model.predict(data_matrix[i])
							prediction = target[i]
							#probability = model.predict_probability
							errors += (prediction-phenotype_ls[i])*(prediction-phenotype_ls[i])
							"""
							if (phenotype_ls[i] != prediction):
									errors = errors + 1
							"""
					errors = errors/size
					print "##########################################"
					print " kernel %s, C %s, gamma %s: mse = %s" % (kname, C, gamma, errors)
					print "##########################################"
	
	@classmethod
	def predictPhenotypeBasedOnGenotype(cls, genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, cofactors=[],\
									report=1, run_type=1):
		"""
		2009-11-18
			predict phenotype based on genotype matrix using machine-learning approaches (rpart, svm)
		"""
		sys.stderr.write("Predicting phenotype based on genotype ... \n")
		from Association import Association
		eigen_vector_fname = ''
		
		test_type = 3 #Emma
		initData = Association.readInData(phenotype_fname, genotype_fname, eigen_vector_fname, phenotype_method_id_ls, test_type=test_type)
		
		which_phenotype_index_ls = initData.which_phenotype_ls
		which_phenotype_index_ls.sort()
		environment_matrix = None
		data_matrix = initData.snpData.data_matrix
		min_data_point = 3
		
		#create an index list of cofactor SNPs
		cofactors_indices = []
		cofactors_set = set(cofactors)
		for i in range(len(initData.snpData.col_id_ls)):
			col_id = initData.snpData.col_id_ls[i]
			if col_id in cofactors_set:
				cofactors_indices.append(i)
		import numpy
		for which_phenotype in which_phenotype_index_ls:
			phenotype_name = initData.phenData.col_id_ls[which_phenotype]
			phenotype_name = phenotype_name.replace('/', '_')	#'/' will be recognized as directory in output_fname
			print phenotype_name
			output_fname='%s_pheno_%s'%(os.path.splitext(output_fname_prefix)[0], phenotype_name)	#make up a new name corresponding to this phenotype
			
			#create non NA phenotype
			phenotype_ls = initData.phenData.data_matrix[:, which_phenotype]
			non_phenotype_NA_row_index_ls = []
			non_NA_phenotype_ls = []
			for i in range(len(phenotype_ls)):
				if not numpy.isnan(phenotype_ls[i]):
					non_phenotype_NA_row_index_ls.append(i)
					non_NA_phenotype_ls.append(phenotype_ls[i])
			
			non_NA_phenotype_ar = numpy.array(non_NA_phenotype_ls)
			new_data_matrix = data_matrix[non_phenotype_NA_row_index_ls,:]
			if run_type==1:
				#new_data_matrix = new_data_matrix.tolist()
				cls.predictByRpart(new_data_matrix, non_NA_phenotype_ls, col_id_ls=initData.snpData.col_id_ls, output_fname=output_fname)			
			elif run_type==2:
				cls.predictBySVM(new_data_matrix, non_NA_phenotype_ls, col_id_ls=initData.snpData.col_id_ls, output_fname=output_fname)
			else:
				sys.stderr.write("Run type %s not supported.\n"%run_type)
		sys.stderr.write("Done.\n")
	"""
	genotype_fname = "/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv"
	genotype_fname = "/Network/Data/250k/db/dataset/call_method_32.tsv"
	phenotype_fname = '/Network/Data/250k/tmp-yh/phenotype20091113.tsv'
	phenotype_method_id_ls = [1,2]
	phenotype_method_id_ls = range(1,8)
	output_fname_prefix = '/tmp/predictPhenotype'
	run_type =2 
	GWA.predictPhenotypeBasedOnGenotype(genotype_fname, phenotype_fname, phenotype_method_id_ls, output_fname_prefix, run_type=run_type)
	"""
	
	
class FileFormatExchange(object):
	"""
	2008-03-18
		check chiamo output
	"""
	def convertChiamoOutput(chiamo_infname, chiamo_outfname, ref_250k_infname, output_fname, posterior_min=0.95):
		from variation.src.common import nt2number
		import csv
		import numpy
		
		reader = csv.reader(open(ref_250k_infname), delimiter='\t')
		SNP_acc_ls = reader.next()[2:]
		del reader
		
		SNP_id2allele_ls = {}
		SNP_id2info = {}
		for SNP_acc in SNP_acc_ls:
			chr, pos, alleleA, alleleB = SNP_acc.split('_')
			SNP_id = '%s_%s'%(chr, pos)
			SNP_id2allele_ls[SNP_id] = [alleleA, alleleB]
			SNP_id2info[SNP_id] = [SNP_acc, len(SNP_id2allele_ls)-1]	#[name, index]
		
		reader = csv.reader(open(chiamo_infname), delimiter='\t')
		strain_id_dup_ls = reader.next()[5:]
		strain_id_ls = []
		for i in range(0, len(strain_id_dup_ls), 2):
			strain_id_ls.append(strain_id_dup_ls[i][:-2])
		del reader
		
		sys.stderr.write("Reading chiamo output ...\n")
		reader = csv.reader(open(chiamo_outfname), delimiter=' ')
		snp_acc_ls = []
		data_matrix = numpy.zeros([len(strain_id_ls), len(SNP_id2allele_ls)], numpy.integer)
		counter = 0
		for row in reader:
			SNP_id = row[0]
			if SNP_id not in SNP_id2allele_ls:
				continue
			allele_ls = SNP_id2allele_ls[SNP_id]
			snp_acc, col_index = SNP_id2info[SNP_id]
			snp_acc_ls.append(snp_acc)
			for i in range(5, len(row)-1, 3):
				posterior_ls = [float(row[i]), float(row[i+1]), float(row[i+2])]	#posterior for 3 classes
				max_posterior_index = numpy.argmax(posterior_ls)
				if posterior_ls[max_posterior_index]>=posterior_min:
					if max_posterior_index==2:	#heterozygous
						call = '%s%s'%(allele_ls[0], allele_ls[1])
					else:
						call = allele_ls[max_posterior_index]
					row_index = (i-5)/3
					try:
						data_matrix[row_index][col_index] = nt2number[call]
					except:
						print SNP_id, snp_acc, col_index, allele_ls, row_index, strain_id_ls[row_index], call
						return
			counter += 1
			if counter%2000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
					
		del reader
		sys.stderr.write("Done.\n")
		
		from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header = ['ecotypeid', 'ecotypeid'] + snp_acc_ls
		FilterStrainSNPMatrix_instance.write_data_matrix(data_matrix, output_fname, header, strain_id_ls, [1]*len(strain_id_ls))
		sys.stderr.write("Finished.\n")
	
	"""
	chiamo_infname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo_14SNPs.in')
	chiamo_outfname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo.out_0_mcmc')
	ref_250k_infname = os.path.expanduser('~/script/variation/genotyping/250ksnp/data/data_250k.tsv')
	output_fname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo_out.tsv')
	
	convertChiamoOutput(chiamo_infname, chiamo_outfname, ref_250k_infname, output_fname, posterior_min=0.95)
	"""
	
	

	@classmethod
	def removeRowsBasedOnSNPAllele(input_fname, output_fname, SNP_label, allele='-'):
		"""
		2008-08-04 investigate whether throwing off some rows help to increase significance
		"""
		from pymodule import read_data, write_data_matrix, nt2number
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		snp_acc_ls = header[2:]
		snp_index = -1
		for i in range(len(snp_acc_ls)):
			if snp_acc_ls[i]==SNP_label:
				snp_index = i
		
		allele = nt2number[allele]
		rows_to_be_tossed_out = set()
		for i in range(len(data_matrix)):
			if data_matrix[i][snp_index]==allele:
				rows_to_be_tossed_out.add(i)
		
		write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=rows_to_be_tossed_out)
	
	"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_2.tsv')
SNP_label = '4_268809_0'
FileFormatExchange.removeRowsBasedOnSNPAllele(input_fname, output_fname, SNP_label, allele='-')
	"""
	
	@classmethod
	def removeSNPsWithMoreThan2Alleles(cls, input_fname, output_fname):
		"""
		2008-08-05
			NPUTE can't work with SNPs with >2 alleles
		"""
		from pymodule import SNPData
		snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1)
		newSNPData = snpData.removeSNPsWithMoreThan2Alleles(snpData)
		newSNPData.tofile(output_fname)
	
	"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_only_2_alleles.tsv')
FileFormatExchange.removeSNPsWithMoreThan2Alleles(input_fname, output_fname)
	"""
	
	
	@classmethod
	def turnNPUTEOutputIntoYuFormat(cls, input_fname, output_fname):
		"""
		2008-08-05
			NPUTE output format is SNPXStrain by and large.
				1st and 2nd column are same as input's 1st row. 1st row is input's 1st column. 2nd row is input's 2nd column.
			
		"""
		from pymodule import SNPData
		snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1, double_header=1, ignore_2nd_column=1)
		snpData.col_id_ls = snpData.header[0][2:]	#take the first header
		
		from pymodule.SNP import transposeSNPData
		newSNPData = transposeSNPData(snpData)
		newSNPData.strain_acc_list = newSNPData.row_id_ls
		newSNPData.category_list = snpData.header[1][2:]
		newSNPData.tofile(output_fname)
	"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.NPUTE.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.NPUTE.yh.tsv')
FileFormatExchange.turnNPUTEOutputIntoYuFormat(input_fname, output_fname)

input_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.impute.tsv')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.impute.Yu.format.tsv')
FileFormatExchange.turnNPUTEOutputIntoYuFormat(input_fname, output_fname)
	"""
	

	@classmethod
	def outputSNPmatrixGivenRegion(cls, snpData, output_fname, chr, start_pos, stop_pos):
		"""
		2008-12-15 output a chromosome region of the SNP matrix
		"""
		sys.stderr.write("Outputting a selected region of SNP matrix ...")
		from pymodule import read_data, write_data_matrix, nt2number
		#header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		snp_acc_ls = snpData.col_id_ls
		cols_to_be_tossed_out = set()
		for i in range(len(snp_acc_ls)):
			chr_pos = snp_acc_ls[i]
			chr_pos = chr_pos.split('_')
			chr_pos = map(int, chr_pos)
			pos = chr_pos[1]
			if not (chr_pos[0]==chr and pos>=start_pos and pos<=stop_pos):
				cols_to_be_tossed_out.add(i)
		
		write_data_matrix(snpData.data_matrix, output_fname, snpData.header, snpData.strain_acc_list, snpData.category_list, cols_to_be_tossed_out=cols_to_be_tossed_out)
	
	"""
input_fname= '/Network/Data/250k/tmp-yh/call_method_17.tsv'
from pymodule import SNPData
snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1)
chr=4
start_pos=100000
stop_pos=700000
output_fname = '/Network/Data/250k/tmp-yh/250k_data/call_method_17_chr%s_%s_%s.tsv'%(chr, start_pos, stop_pos)
FileFormatExchange.outputSNPmatrixGivenRegion(snpData, output_fname, chr, start_pos, stop_pos)
	"""

	def reduce250kOnlyToCertainSNPs(cls, snp_id_ls):
		#2008-10-07 form a smaller 250k test dataset for PlotGroupOfSNPs.py, snp_id_ls is the top 200 snps from (KW,LD)
		
		from pymodule import SNPData, write_data_matrix
		input_fname = '/Network/Data/250k/tmp-yh/call_method_17.tsv'
		snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1)
		output_fname = '/Network/Data/250k/tmp-yh/call_method_17_test.tsv'
		from sets import Set
		good_snp_id_set = Set(snp_id_ls)
		col_index_to_be_tossed_set = Set()
		for snp_id, col_index in snpData.col_id2col_index.iteritems():
			if snp_id not in good_snp_id_set:
				col_index_to_be_tossed_set.add(col_index)
		write_data_matrix(snpData.data_matrix, output_fname, snpData.header, snpData.strain_acc_list, snpData.category_list, 
						rows_to_be_tossed_out=None, cols_to_be_tossed_out=col_index_to_be_tossed_set)
	"""
FileFormatExchange.reduce250kOnlyToCertainSNPs = classmethod(reduce250kOnlyToCertainSNPs)
	"""
	
	"""
	2007-03-28
	"""
	@classmethod
	def strip_2010_strain_info(cls, input_fname, output_fname):
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		#skip the 1st two sentences
		reader.next()
		reader.next()
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for row in reader:
			for i in range(len(row)):
				row[i] = row[i].strip()
			writer.writerow(row)
		del reader, writer
	
	"""
FileFormatExchange.strip_2010_strain_info('./script/variation/data/2010/2010_strain_info.csv', './script/variation/data/2010/2010_strain_info_stripped.csv')
	"""
	
	"""
	2007-03-15
	"""
	def reformat_data_for_chris(curs, input_fname, output_fname, snp_locus_table):
		from common import number2nt
		
		sys.stderr.write("Getting snp_acc2pos ...")
		snp_acc2pos = {}
		curs.execute("select acc, chromosome, position from %s"%snp_locus_table)
		rows = curs.fetchall()
		for row in rows:
			snp_acc, chromosome, position = row
			snp_acc2pos[snp_acc] = [chromosome, position]
		sys.stderr.write("Done.\n")
		
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = reader.next()
		new_header = [header[0]] + header[2:]
		writer.writerow(new_header)
		
		chr_ls = ['']
		pos_ls = ['']
		for snp_acc in new_header[1:]:
			chr_ls.append(snp_acc2pos[snp_acc][0])
			pos_ls.append(snp_acc2pos[snp_acc][1])
		writer.writerow(chr_ls)
		writer.writerow(pos_ls)
		
		for row in reader:
			new_row = [row[0]]
			for call in row[2:]:
				nt = number2nt[int(call)]
				if nt=='NA':	#chris wants 'N' instead of 'NA'
					nt = 'N'
				new_row.append(nt)
			writer.writerow(new_row)
		del writer
	
	"""
	hostname='zhoudb'
	dbname='graphdb'
	schema = 'dbsnp'
	conn, curs = db_connect(hostname, dbname, schema)
	reformat_data_for_chris(curs, '../data/justin_data_filtered.csv', '../data/justin_data_filtered_for_chris.csv', 'snp_locus')
	"""


	@classmethod
	def removeRowsBasedOnPhenotypeValue(input_fname, phenotype_fname, phenotype_method_id, output_fname, min_pheno_value=None, max_pheno_value=None):
		"""
		2009-4-13 investigate whether throwing off some rows help to increase significance.
			similar purpose to removeRowsBasedOnSNPAllele()
			throw away accessions whose phenotypes are not in the phenotype range specified
		"""
		from pymodule import read_data, write_data_matrix, SNPData
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
						data_matrix=data_matrix)
		header_phen, strain_acc_list_phen, category_list_phen, data_matrix_phen = read_data(phenotype_fname, turn_into_integer=0)
		from variation.src.Association import Association
		data_matrix_phen = Association.get_phenotype_matrix_in_data_matrix_order(strain_acc_list, strain_acc_list_phen, data_matrix_phen)
		phenData = SNPData(header=header_phen, strain_acc_list=snpData.strain_acc_list, data_matrix=data_matrix_phen)
		
		from variation.src.PlotGroupOfSNPs import PlotGroupOfSNPs
		which_phenotype_ls = PlotGroupOfSNPs.findOutWhichPhenotypeColumn(phenData, set([phenotype_method_id]))
		which_phenotype = which_phenotype_ls[0]
		
		import numpy
		rows_to_be_tossed_out = set()
		for i in range(len(phenData.data_matrix)):
			phenotype_value = phenData.data_matrix[i][which_phenotype]
			if numpy.isnan(phenotype_value):
				rows_to_be_tossed_out.add(i)
				continue
			
			if min_pheno_value is not None and phenotype_value<min_pheno_value:
				rows_to_be_tossed_out.add(i)
				continue
			if max_pheno_value is not None and phenotype_value>max_pheno_value:
				rows_to_be_tossed_out.add(i)
				continue
		
		write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=rows_to_be_tossed_out)
		
	"""
	input_fname='/Network/Data/250k/tmp-yh/250k_data/call_method_17_test.tsv'
	input_fname='/Network/Data/250k/tmp-yh/250k_data/call_method_29.tsv'
	phenotype_fname='/Network/Data/250k/tmp-yh/phenotype_no_transform.tsv'
	phenotype_method_id=44
	min_pheno_value = 0.45
	output_fname = '/tmp/call_method_29_pheno_id_%s_above_%s.tsv'%(phenotype_method_id, min_pheno_value)
	FileFormatExchange.removeRowsBasedOnPhenotypeValue(input_fname, phenotype_fname, phenotype_method_id, output_fname, min_pheno_value)
	"""
	
	@classmethod
	def reverseComplementFastaInput(cls, input_fname, output_fname):
		"""
		2009-3-25
			read in sequences from input_fname (in fasta format), reverse complement each sequence and output them 
		"""
		inf = open(input_fname)
		of = open(output_fname, 'w')
		from Bio import SeqIO
		for seq_record in SeqIO.parse(inf, "fasta") :
			seq_rc = seq_record.seq.reverse_complement()
			of.write('>%s\n'%seq_record.id)
			of.write('%s\n'%seq_rc.tostring())
	
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC.seq')
	output_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC_rc.seq')
	AnalyzeSNPData.reverseComplementFastaInput(input_fname, output_fname)
	"""
	
	
	"""
	2009-3-20
		output fasta format in DrawSNPMatrix.py input format.
		also partition the input into blocks each with number of columns <= maxNumberColsPerBlock.
		too many columns render some parts of the image generated by DrawSNPMatrix.py unviewable.
	"""
	@classmethod
	def fasta2DrawSNPMatrixFormat(cls ,input_fname, output_fname, maxNumberColsPerBlock=2000):
		import csv, math
		blockIndex2writer = {}
		inf = open(input_fname)
		from Bio import SeqIO
		isHeaderOutputted = False
		header_row = []
		for seq_record in SeqIO.parse(inf, "fasta") :
			one_row = []
			no_of_blocks = int(math.ceil(len(seq_record.seq)/float(maxNumberColsPerBlock)))
			for i in range(len(seq_record.seq)):
				nt = seq_record.seq[i]
				if not isHeaderOutputted:
					header_row.append(i+1)
				one_row.append(nt)
			
			#output header
			if not isHeaderOutputted:
				for i in range(no_of_blocks):
					start_index = i*maxNumberColsPerBlock
					stop_index = (i+1)*maxNumberColsPerBlock
					row_for_this_block = header_row[start_index:stop_index]
					if i not in blockIndex2writer:
						output_fname_ls = os.path.splitext(output_fname)
						block_output_fname = '%s_%s_%s%s'%(output_fname_ls[0], start_index, stop_index, output_fname_ls[1])
						writer = csv.writer(open(block_output_fname, 'w'), delimiter='\t')
						blockIndex2writer[i] = writer
					blockIndex2writer[i].writerow(['name',1]+row_for_this_block)
				isHeaderOutputted = True
			
			#output data
			for i in range(no_of_blocks):
				start_index = i*maxNumberColsPerBlock
				stop_index = (i+1)*maxNumberColsPerBlock
				row_for_this_block = one_row[start_index:stop_index]
				blockIndex2writer[i].writerow([seq_record.id, 1] + row_for_this_block)
		
		del blockIndex2writer, inf
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLCAlign.fasta')
	output_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLCAlign.tsv')
	AnalyzeSNPData.fasta2DrawSNPMatrixFormat(input_fname, output_fname)
	"""
	
class Data250k(object):
	def __init__(self):
		"""
		2009-11-14
			keep the old name alive
		"""
		choose192OutOf250k = self.retainSNPDataFromCertainEcotypes
	
	"""
	2008-08-16
		check if the array_ids and ecotype_ids in call files generated by bjarni match the ones in db
	"""
	@classmethod
	def checkBjarniFile(cls, input_fname, curs):
		import csv
		reader = csv.reader(open(input_fname))
		array_id_ls = reader.next()[2:]
		array_id_ls = map(int, array_id_ls)
		print "%s arrays"%(len(array_id_ls))
		ecotype_id_ls = reader.next()[2:]
		ecotype_id_ls = map(int, ecotype_id_ls)
		print "%s ecotypes"%(len(ecotype_id_ls))
		for i in range(len(array_id_ls)):
			array_id = array_id_ls[i]
			ecotype_id = ecotype_id_ls[i]
			curs.execute("select tg_ecotypeid from stock.ecotypeid2tg_ecotypeid where ecotypeid=%s"%ecotype_id)
			rows = curs.fetchall()
			if rows:
				bjarni_ecotype_id = rows[0][0]
			else:
				sys.stderr.write( "ecotype_id %s has no tg_ecotypeid in ecotypeid2tg_ecotypeid.\n"%(ecotype_id))
				bjarni_ecotype_id = ecotype_id
			curs.execute("select maternal_ecotype_id from stock_250k.array_info where id=%s"%array_id)
			rows = curs.fetchall()
			maternal_ecotype_id = rows[0][0]
			if bjarni_ecotype_id!=maternal_ecotype_id:
				print i, array_id, bjarni_ecotype_id, maternal_ecotype_id
		del reader
	
	"""
	input_fname = '/Network/Data/250k/dataFreeze_080608/250K_f8_080608.csv'
	checkBjarniFile(input_fname, curs)
	"""
	
	@classmethod
	def recoverArrayID2ndCol(cls, input_fname, data_with_array_id_fname, output_fname):
		"""
		2008-01-04 recover array id (2nd column) of a StrainXSNP data file based on its old version
			the 1st column (strain id) of input_fname has to be unique,
			otherwise a random one among duplicates would get a possibly wrong array id assigned.
		"""
		from pymodule import read_data, SNPData
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, data_matrix=data_matrix)
		#
		header2, strain_acc_list2, category_list2, data_matrix2 = read_data(data_with_array_id_fname)
		#put the array_id from category_list2 into corresponding position in category_list
		for i in range(len(strain_acc_list2)):
			ecotype_id = strain_acc_list2[i]
			array_id = category_list2[i]
			row_index = snpData.row_id2row_index.get(ecotype_id)	#ecotype_id might not be in input_fname
			if row_index is not None:
				category_list[row_index] = array_id
		
		newSNPData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list,\
							data_matrix=data_matrix)
		newSNPData.tofile(output_fname)
	
	"""
	input_fname = os.path.expanduser('~/panfs/NPUTE_data/input/250k_l3_y.85_uniq_ecotype_20080919_3_FRI_del_no_array_id.tsv')
	data_with_array_id_fname = os.path.expanduser('~/panfs/NPUTE_data/input/250k_l3_y.85_uniq_ecotype_20080919.tsv')
	output_fname = os.path.expanduser('~/panfs/NPUTE_data/input/250k_l3_y.85_uniq_ecotype_20080919_3_FRI_del.tsv')
	Data250k.recoverArrayID2ndCol(input_fname, data_with_array_id_fname, output_fname)
	"""
	
	@classmethod
	def chooseFirst96(cls, input_fname, db, output_fname):
		"""
		2009-2-17
			pick 1st 96 2010 accessions out of 250k data
		"""
		from pymodule import SNPData
		snp_data = SNPData(input_fname=input_fname, turn_into_array=1, ignore_2nd_column=1)
		rows = db.metadata.bind.execute("select * from at.accession2tg_ecotypeid where accession_id<=96")
		rows_to_be_preserved = set()
		for row in rows:
			row_id = '%s'%row.ecotype_id
			if row_id in snp_data.row_id2row_index:
				rows_to_be_preserved.add(snp_data.row_id2row_index[row_id])
		rows_to_be_tossed_out = set(range(len(snp_data.row_id_ls))) - rows_to_be_preserved
		snp_data.tofile(output_fname, rows_to_be_tossed_out=rows_to_be_tossed_out)
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/call_method_29.tsv')
	output_fname = os.path.expanduser('~/panfs/250k/call_method_29_only96.tsv')
	Data250k.chooseFirst96(input_fname, db_250k, output_fname)
	"""
	@classmethod
	def filterGWAToRetainSegregatingSNPs(cls, db_250k, call_method_id, phenotype_method_id, analysis_method_id, snpData, \
										 ecotype_id1, ecotype_id2, output_fname):
		"""
		2009-3-12
			read in a GWA result, and its affiliated genotype matrix
			and retain SNPs in GWA result that are segregating between ecotypes
		"""
		gwar = db_250k.getGWA(call_method_id, phenotype_method_id, analysis_method_id)
		row_index1 = snpData.row_id2row_index[str(ecotype_id1)]
		row_index2 = snpData.row_id2row_index[str(ecotype_id2)]
		
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['chromosome', 'position', 'pvalue'])
		for data_obj in gwar.data_obj_ls:
			snp_id = '%s_%s'%(data_obj.chromosome, data_obj.position)
			col_index = snpData.col_id2col_index[snp_id]
			allele1 = snpData.data_matrix[row_index1][col_index]
			allele2 = snpData.data_matrix[row_index2][col_index]
			if allele1!=allele2:
				writer.writerow([data_obj.chromosome, data_obj.position, data_obj.value])
		del writer

	"""
	call_method_id = 29
	phenotype_method_id=1
	analysis_method_id=7
	ecotype_id1=6917	'Fab-2'
	ecotype_id2=6904	'Br-0'
	output_fname = '/tmp/gwa_call_%s_phenotype_%s_analysis_%s_segregating_in_%s_%s.tsv'%(call_method_id, phenotype_method_id, analysis_method_id, ecotype_id1, ecotype_id2)
	snpData = db_250k.getSNPMatrix(call_method_id)
	Data250k.filterGWAToRetainSegregatingSNPs(db_250k, call_method_id, phenotype_method_id, analysis_method_id, snpData, ecotype_id1, ecotype_id2, output_fname)
	"""
	
	@classmethod
	def retainSNPDataFromCertainEcotypes(cls, input_fname, ecotype_id_fname, output_fname):
		"""
		2009-11-14
			function renamed from choose192OutOf250k()
			
			retain 250k data (input_fname) that whose ecotype IDs are within ecotype_id_fname (one-line-header)
			output the filtered data into he new output_fname.
		2009-3-27
		"""
		from pymodule import SNPData
		import csv
		sys.stderr.write("Getting ecotype IDs from %s ..."%ecotype_id_fname)
		reader = csv.reader(open(ecotype_id_fname), delimiter=',')
		reader.next()	#toss the header
		ecotype_set = set()
		for row in reader:
			ecotype_set.add(row[0])
		del reader
		sys.stderr.write("%s ecotypes.\n"%len(ecotype_set))
		
		snp_data = SNPData(input_fname=input_fname, turn_into_array=1)	#2nd column array id is also needed
		
		rows_to_be_tossed_out = set()
		for row_id in snp_data.row_id_ls:
			ecotype_id = row_id[0]	#1st column is ecotype_id, 2nd is array id
			if ecotype_id not in ecotype_set:
				rows_to_be_tossed_out.add(snp_data.row_id2row_index[row_id])
		snp_data.tofile(output_fname, rows_to_be_tossed_out=rows_to_be_tossed_out)


	"""
	input_fname = '/Network/Data/250k/db/dataset/call_method_29.tsv'
	fname_with_192_ids = '/tmp/192_accessions_031009.csv' # email from bjarni. coma delimited, 1st column is ecotype id. actually 199 accessions.
	output_fname = '/tmp/call_method_29_with_192_accession.tsv'
	Data250k.choose192OutOf250k(input_fname, fname_with_192_ids, output_fname)
	"""
	
	@classmethod
	def filterCallMethod1ToHaveSameArraysAsCallMethod2(cls, call_method_id1, call_method_id2, output_fname):
		"""
		2009-5-19
			 generate a new snp dataset based on call_method_id1 but without the arrays that are not in call_method_id2
		"""
		import Stock_250kDB
		call_method1 = Stock_250kDB.CallMethod.get(call_method_id1)
		call_method2 = Stock_250kDB.CallMethod.get(call_method_id2)
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=call_method1.filename, turn_into_array=1)
		
		sys.stderr.write("Selecting overlapping arrays ...")
		array_id_in_call_method2_set = set()
		for call_info in call_method2.call_info_ls:
			array_id_in_call_method2_set.add(call_info.array_id)
		
		#print array_id_in_call_method2_set
		
		row_id_wanted_ls = []
		for row_id in snpData1.row_id_ls:
			array_id = int(row_id[1])
			if array_id in array_id_in_call_method2_set:
				row_id_wanted_ls.append(row_id)
		#print row_id_wanted_ls
		sys.stderr.write("Done.\n")
		
		newSNPData = SNPData.keepRowsByRowID(snpData1, row_id_wanted_ls)
		newSNPData.tofile(output_fname)
		del newSNPData
		del snpData1
	
	"""
call_method_id1 = 34
call_method_id2 = 33
output_fname = '/tmp/call_method_%s_arrays_from_%s.tsv'%(call_method_id1, call_method_id2)
Data250k.filterCallMethod1ToHaveSameArraysAsCallMethod2(call_method_id1, call_method_id2, output_fname)
	"""
	
class BooleanInteraction(object):
	
	"""
	2008-08-05
		test boolean relationship between SNPs
	"""
	
	def returnTop2Allele(snp_allele2count):
		"""
		2008-08-06 remove redundant argument snp_allele_ls
		2008-08-05 in the descending order of count for each allele, assign index ascending
		"""
		snp_allele_count_ls = []
		snp_allele_ls = snp_allele2count.keys()
		for snp_allele in snp_allele_ls:
			snp_allele_count_ls.append(snp_allele2count[snp_allele])
		import numpy
		argsort_ls = numpy.argsort(snp_allele_count_ls)
		new_snp_allele2index = {}
		for i in [-1, -2]:
			snp_index = argsort_ls[i]	#-1 is index for biggest, -2 is next biggest
			new_snp_allele2index[snp_allele_ls[snp_index]] = -i-1
		return new_snp_allele2index
	
	@classmethod
	def booleanMergeSNPs(cls, input_fname, output_fname, SNP_label1, SNP_label2, operator_type=1):	#1 is and, 2 is or
		"""
		2008-08-05
			alleles not in the top 2 are taken as NA. major allele is coded as 0. minor allele is coded as 1.
		"""
		from pymodule import read_data, write_data_matrix, nt2number
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
		
		snp_acc_ls = header[2:]
		snp_index1 = -1
		snp_index2 = -1
		for i in range(len(snp_acc_ls)):
			if snp_acc_ls[i]==SNP_label1:
				snp_index1 = i
			if snp_acc_ls[i]==SNP_label2:
				snp_index2 = i
		
		snp_allele2count1 = {}
		snp_allele2count2 = {}
		no_of_rows = len(data_matrix)
		
		for i in range(no_of_rows):
			snp1_allele = data_matrix[i][snp_index1]
			snp2_allele = data_matrix[i][snp_index2]
			if snp1_allele!=0:
				if snp1_allele not in snp_allele2count1:
					snp_allele2count1[snp1_allele] = 0
				snp_allele2count1[snp1_allele] += 1
			if snp2_allele!=0:
				if snp2_allele not in snp_allele2count2:
					snp_allele2count2[snp2_allele] = 0
				snp_allele2count2[snp2_allele] += 1
		print snp_allele2count1
		print snp_allele2count2
		snp_allele2index1 = returnTop2Allele(snp_allele2count1)
		
		snp_allele2index2 = returnTop2Allele(snp_allele2count2)
		print snp_allele2index1
		print snp_allele2index2
		
		no_of_cols = 1
		new_data_matrix = data_matrix	#replace the 1st SNP's data with the new boolean result
		for i in range(no_of_rows):
			snp1_allele = data_matrix[i][snp_index1]
			snp2_allele = data_matrix[i][snp_index2]
			if snp1_allele in snp_allele2index1:
				snp_code1 = snp_allele2index1[snp1_allele]
			else:
				snp_code1 = 0
			
			if snp2_allele in snp_allele2index2:
				snp_code2 = snp_allele2index2[snp2_allele]
			else:
				snp_code2 = 0
			if operator_type==1:
				if snp1_allele in snp_allele2index1 and snp2_allele in snp_allele2index2:
					new_data_matrix[i][snp_index1] = (snp_code1 and snp_code2) + 1
			elif operator_type==2:
				if snp1_allele in snp_allele2index1 or snp2_allele in snp_allele2index2:
					new_data_matrix[i][snp_index1] = (snp_code1 or snp_code2) + 1
		write_data_matrix(new_data_matrix, output_fname, header, strain_acc_list, category_list)
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
	output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2.tsv')
	SNP_label1 = '4_268809_0'
	SNP_label2 = '4_269962_8'
	booleanMergeSNPs(input_fname, output_fname, SNP_label1, SNP_label2, operator_type=2)
	input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2.tsv')
	output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2_or_3.tsv')
	SNP_label3 = '4_270712_0'
	booleanMergeSNPs(input_fname, output_fname, SNP_label1, SNP_label3, operator_type=2)
	"""
	
	
	@classmethod
	def handleMostSignificantOperatorPerSNPPair(cls, input_fname, function_handler, param_obj):
		"""
		2009-2-19
			a general function which calls function_handler to do a bit processing for the most significant operator of each SNP pair.
				(results from the same SNP pair are always cluttered together.)
		"""
		import csv
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()
		prev_snp2_id = None
		prev_snp1_id = None
		min_pvalue = None
		row_with_min_pvalue = None
		i = 0
		#input_fname_basename = os.path.splitext(os.path.basename(input_fname))[0]	#function calls this general function should supply this
		#param_obj.input_fname_basename = input_fname_basename
		counter = 0
		real_counter = 0
		for row in reader:
			snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
			counter += 1
			pvalue = float(pvalue)
			if prev_snp2_id is None:	#first snp2
				prev_snp1_id = snp1_id
				prev_snp2_id = snp2_id
				min_pvalue = pvalue
				row_with_min_pvalue = row
			elif snp2_id == prev_snp2_id and snp1_id==prev_snp1_id:	#same snp2
				i += 1
				if pvalue<min_pvalue:
					min_pvalue=pvalue
					row_with_min_pvalue = row
			else:	#new pairs
				i = 0
				prev_snp1_id = row_with_min_pvalue[0]
				function_handler(row_with_min_pvalue, param_obj)
				real_counter += 1
				prev_snp1_id = snp1_id
				prev_snp2_id = snp2_id
				min_pvalue = pvalue
				row_with_min_pvalue = row
				if real_counter%50000==0:
					sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
	
	"""
	2009-1-25
		input_fname is output of MpiIntraGeneSNPPairAsso.py or MpiInterGeneSNPPairAsso.py or MpiInterSNPPairAsso.py
		
		separate one snp-pairwise association results into single-SNP GWA files according to the snp1_id.
	"""
	@classmethod
	def outputRow(cls, row, param_obj):
		snp_id2writer = getattr(param_obj, "snp_id2writer")
		output_dir = getattr(param_obj, "output_dir")
		input_fname_basename = getattr(param_obj, "input_fname_basename")
		header = getattr(param_obj, "header")
		
		snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
		if snp1_id not in snp_id2writer:
			output_fname = os.path.join(output_dir, '%s_vs_%s.tsv'%(input_fname_basename, snp1_id))
			writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
			writer.writerow(header)
			snp_id2writer[snp1_id] = writer
		writer = snp_id2writer[snp1_id]
		chr, pos = snp2_id.split('_')
		count1 = int(count1)
		count2 = int(count2)
		mac = min(count1, count2)
		maf = float(mac)/(count1+count2)
		row = [chr, pos, pvalue, maf, mac, bool_type]
		writer.writerow(row)
	
	@classmethod
	def outputBooleanPairIntoGWAFormat(cls, input_fname, output_dir, pos_index=1, need_beta=True, min_value_cutoff=None, do_log10_transformation=False):
		"""
		2009-1-25
			call outputRow()
		"""
		import csv
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		snp_id2writer = {}	#each SNP1_ID would have a separate file for results from all SNPs that tested with that SNP
		
		input_fname_basename = os.path.splitext(os.path.basename(input_fname))[0]
		from pymodule import PassingData
		param_obj = PassingData(output_dir=output_dir, input_fname_basename=input_fname_basename, snp_id2writer=snp_id2writer,\
							header=header)
		cls.handleMostSignificantOperatorPerSNPPair(input_fname, cls.outputRow, param_obj)
		del snp_id2writer
	
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/boolean_gene_vs_snp_DOG1/SNPpair_41_JIC4W.tsv')
	output_dir = os.path.expanduser('~/250k/boolean_gene_vs_snp_DOG1_in_gwr/')
	outputBooleanPairIntoGWAFormat(input_fname, output_dir)
	"""
	
	@classmethod
	def getPvaluePerOperator(cls, input_fname, no_of_lines_to_skip=0):
		"""
		2009-2-19
			take boolean KW output as input_fname
			return pvalue list for each operator, no_of_lines_to_skip is the average number of lines to skip after each line to control intake.
		"""
		import csv, sys, traceback, random
		operator2pvalue_ls = {}
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()
		prev_snp2_id = None
		prev_snp1_id = None
		min_pvalue = None
		row_with_min_pvalue = None
		
		random_no_pool = []
		if no_of_lines_to_skip>0:
			pop_pool = range(0,no_of_lines_to_skip*2)
			for i in range(100):
				random_no_pool.append(random.sample(pop_pool, 1)[0])
		
		counter = 0
		real_counter = 0
		for row in reader:
			snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
			bool_type = int(bool_type)
			pvalue = float(pvalue)
			counter += 1
			real_counter += 1
			if bool_type not in operator2pvalue_ls:
				operator2pvalue_ls[bool_type] = []
			operator2pvalue_ls[bool_type].append(pvalue)
			if real_counter%50000==0:
				sys.stderr.write("%s%s\t\t%s"%('\x08'*40, real_counter, counter))
			if no_of_lines_to_skip>0:
				random_no = random_no_pool[counter%len(random_no_pool)]
				for i in range(random_no):
					try:
						reader.next()
						counter += 1
					except:
						traceback.print_exc()
						sys.stderr.write('%s.\n'%repr(sys.exc_info()))
						break
		del reader
		return operator2pvalue_ls
	
	@classmethod
	def stuffPvalueIntoDict(cls, row, param_obj):
		"""
		2009-2-19
			called by getPvaluePerOperatorOnlyTopInSNPPair()
		"""
		operator2pvalue_ls = getattr(param_obj, "operator2pvalue_ls")
		snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
		bool_type = int(bool_type)
		if bool_type not in operator2pvalue_ls:
			operator2pvalue_ls[bool_type] = []
		pvalue = float(pvalue)
		operator2pvalue_ls[bool_type].append(pvalue)
	
	@classmethod
	def getPvaluePerOperatorOnlyTopInSNPPair(cls, input_fname, no_of_lines_to_skip=0):
		"""
		2009-2-19
			take boolean KW output as input_fname
			return pvalue list for each operator, no_of_lines_to_skip is the average number of lines to skip after each line to control intake.
			
			call stuffPvalueIntoDict()
		"""
		import csv, sys, traceback, random
		operator2pvalue_ls = {}
		input_fname_basename = os.path.splitext(os.path.basename(input_fname))[0]
		from pymodule import PassingData
		param_obj = PassingData(operator2pvalue_ls=operator2pvalue_ls, input_fname_basename=input_fname_basename)
		cls.handleMostSignificantOperatorPerSNPPair(input_fname, cls.stuffPvalueIntoDict, param_obj)
		return operator2pvalue_ls
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/IntraGeneSNPPair_BooleanKW_m10k_call29/SNPpair_7_FT_22C.tsv')
	output_fname_prefix = os.path.expanduser('%s_top_bool_type_pvalue_hist'%os.path.splitext(input_fname)[0])
	operator_top2pvalue_ls = BooleanInteraction.getPvaluePerOperatorOnlyTopInSNPPair(input_fname, no_of_lines_to_skip=0)
	BooleanInteraction.drawPvalueHist(operator_top2pvalue_ls, output_fname_prefix, minus_log_pvalue=True)
	"""
	
	@classmethod
	def drawPvalueHist(cls, operator2pvalue_ls, output_fname_prefix, minus_log_pvalue=False):
		"""
		2009-2-19
			draw histogram distribution for each each operator's pvalue list
		"""
		import pylab
		pylab.clf()
		from variation.src.PlotGenePairAssoResult import PlotGenePairAssoResult
		#bool_type2marker_and_name
		import random
		hist_handler_ls = []
		legend_ls = []
		bool_type_ls = operator2pvalue_ls.keys()
		bool_type_ls.sort()
		color_ls = ['r', 'g', 'c', 'b', 'y', 'k']
		log_func = lambda x: -math.log(x)
		for i in range(len(bool_type_ls)):
			pylab.clf()
			bool_type = bool_type_ls[i]
			pvalue_ls = operator2pvalue_ls[bool_type]
			
			pvalue_ls = random.sample(pvalue_ls, 10000)
			
			if minus_log_pvalue:
				xlabel = '-logPvalue'
				pvalue_ls = map(log_func, pvalue_ls)
			else:
				xlabel = 'pvalue'
			n1 = pylab.hist(pvalue_ls, 50, alpha=0.3, normed=1, facecolor=color_ls[i])
			hist_handler_ls.append(n1[2][0])
			marker, bool_type_name = PlotGenePairAssoResult.bool_type2marker_and_name[bool_type][:2]
			legend_ls.append(bool_type_name)
			#pylab.legend(hist_handler_ls, legend_ls)
			pylab.ylabel('frequency')
			
			pylab.title(bool_type_name)
			pylab.savefig('%s_%s_%s.png'%(output_fname_prefix, bool_type, bool_type_name), dpi=300)
	
	"""
	input_fname = os.path.expanduser('~/panfs/250k/IntraGeneSNPPair_BooleanKW_m10k_call29/SNPpair_7_FT_22C.tsv')
	output_fname_prefix = os.path.expanduser('%s_bool_type_pvalue_hist'%os.path.splitext(input_fname)[0])
	operator2pvalue_ls = BooleanInteraction.getPvaluePerOperator(input_fname, no_of_lines_to_skip=0)
	BooleanInteraction.drawPvalueHist(operator2pvalue_ls, output_fname_prefix, minus_log_pvalue=True)
	"""
	
	@classmethod
	def countBoolTypeTopForEachSNPPair(cls, row, param_obj):
		"""
		2009-2-19
		"""
		import math, sys, traceback
		operator_int_pvalue2count = getattr(param_obj, "operator_int_pvalue2count")
		int_pvalue2count = getattr(param_obj, "int_pvalue2count")
		snp1_id, gene1_id, snp2_id, gene2_id, bool_type, pvalue, count1, count2 = row[:8]
		bool_type = int(bool_type)
		try:
			minus_log_pvalue = -math.log(float(pvalue))
		except:
			traceback.print_exc()
			sys.stderr.write('%s. minus_log_pvalue=50.\n'%repr(sys.exc_info()))
			minus_log_pvalue = 50
		int_pvalue = int(minus_log_pvalue)
		if int_pvalue not in int_pvalue2count:
			int_pvalue2count[int_pvalue] = 0
		int_pvalue2count[int_pvalue] += 1
		
		operator_int_pvalue = (bool_type, int_pvalue)
		if operator_int_pvalue not in operator_int_pvalue2count:
			operator_int_pvalue2count[operator_int_pvalue] = 0
		operator_int_pvalue2count[operator_int_pvalue] += 1
	
	@classmethod
	def plot_operator_int_pvalue2count(cls, operator_int_pvalue2count, int_pvalue2count, output_fname_prefix, \
									bool_type_ls=[1,2,4,6,7], int_pvalue_step=2):
		"""
		2009-2-19
		"""
		#from variation.src.PlotGenePairAssoResult import PlotGenePairAssoResult
		#PlotGenePairAssoResult.bool_type2marker_and_name[bool_type]
		int_pvalue_ls = int_pvalue2count.keys()
		int_pvalue_ls.sort()
		
		#group int_pvalue into groups according to int_pvalue_step
		row_id_ls = []
		for i in range(len(int_pvalue_ls)):
			int_pvalue = int_pvalue_ls[i]
			row_id = int_pvalue/int_pvalue_step
			if len(row_id_ls)==0:
				row_id_ls.append(row_id)
			else:
				if row_id!=row_id_ls[-1]:	#only add when it's new
					row_id_ls.append(row_id)
		
		col_id_ls = bool_type_ls
		import numpy
		#col_id_ls = numpy.array(col_id_ls)
		no_of_rows = len(row_id_ls)
		no_of_cols = len(col_id_ls)
		data_matrix = numpy.zeros([no_of_rows, no_of_cols], numpy.int)
		from pymodule import SNPData
		countData = SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=data_matrix)
		
		for int_pvalue in int_pvalue_ls:
			row_id = int_pvalue/int_pvalue_step
			row_index = countData.row_id2row_index[row_id]
			for bool_type in col_id_ls:
				operator_int_pvalue = (bool_type, int_pvalue)
				if operator_int_pvalue in operator_int_pvalue2count:
					count = operator_int_pvalue2count[operator_int_pvalue]
				else:
					count = 0
				col_index = countData.col_id2col_index[bool_type]
				countData.data_matrix[row_index][col_index] += count
		
		perc_matrix = numpy.zeros([no_of_rows, no_of_cols], numpy.float)
		percData = SNPData(row_id_ls=row_id_ls, col_id_ls=col_id_ls, data_matrix=perc_matrix)
		row_sum_array = numpy.sum(countData.data_matrix, 1, numpy.float)
		for j in range(no_of_cols):
			percData.data_matrix[:,j] = numpy.sum(countData.data_matrix[:,range(j+1)], 1)/row_sum_array	#cumulative percentage
		
		import pylab
		pylab.clf()
		pylab.grid(True, alpha=0.3)
		from variation.src.PlotGenePairAssoResult import PlotGenePairAssoResult
		color_ls = ['r', 'g', 'c', 'b', 'y', 'k']
		plot_ls = []
		legend_ls = []
		row_id_ls = numpy.array(row_id_ls)
		for j in range(no_of_cols):
			bool_type = percData.col_id_ls[j]
			marker, bool_type_name = PlotGenePairAssoResult.bool_type2marker_and_name[bool_type][:2]
			x_ls = percData.data_matrix[:,j]
			y_ls = row_id_ls*int_pvalue_step
			a = pylab.plot(x_ls, y_ls, '.-', color=color_ls[j])
			plot_ls.append(a)
			legend_ls.append(bool_type_name)
		
		pylab.legend(plot_ls, legend_ls)
		pylab.ylabel('-logPvalue')
		pylab.xlabel('percentage')
		pylab.title('Percentage of each bool type at certain log pvalue range')
		pylab.savefig('%s_perc_bool_by_log_pvalue.png'%(output_fname_prefix), dpi=300)
			
	
	"""
operator_int_pvalue2count = {}
int_pvalue2count = {}
from pymodule import PassingData
param_obj = PassingData(operator_int_pvalue2count=operator_int_pvalue2count, int_pvalue2count=int_pvalue2count)
BooleanInteraction.handleMostSignificantOperatorPerSNPPair(input_fname, BooleanInteraction.countBoolTypeTopForEachSNPPair, param_obj)
BooleanInteraction.plot_operator_int_pvalue2count(operator_int_pvalue2count, int_pvalue2count, output_fname_prefix)
	"""


class FRIDeletion(object):
	"""
		in processing FRI deletion data from Shindo2005. check plone doc, /research/variation/log-2008-07.
	2008-08-05
	"""
	@classmethod
	def outputShindo2005(cls, input_fname, output_fname, which_type_of_id_to_output=1):
		"""
		2008-08-06
			in output. if which_type_of_id_to_output==1, output in ecotypeid; else output in accession_id
		"""
		from variation.src.AtDB import AtDB, Sequence, Alignment
		db = AtDB(hostname='localhost')
		from pymodule import read_data, write_data_matrix
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=0)
		
		no_of_rows = len(strain_acc_list)
		import numpy
		new_data_matrix = numpy.ones([no_of_rows, 2], numpy.int8)
		new_strain_acc_list = []
		new_category_list = []
		for i in range(len(strain_acc_list)):
			strain_acc = strain_acc_list[i].upper()
			#2008-08-06 for magnus's data.csv	no swedish letters.
			rows = db.metadata.bind.execute("select a2e.ecotype_id, a2e.nativename, a2e.accession_id from accession2tg_ecotypeid a2e, magnus_192_vs_accession m2a where upper(m2a.linename)='%s' and m2a.accession_id=a2e.accession_id"%\
									(strain_acc))
			
			#2008-08-06 for stuff extracted from Supplemental Figure 4A, Figure 4B. this query doesn't help to recognize those swedish letters
			#rows = db.metadata.bind.execute("select a2e.ecotype_id, a2e.nativename, a2e.accession_id from accession2tg_ecotypeid a2e where upper(a2e.accession_name)='%s'"%\
			#						(strain_acc))
			try:
				row = rows.fetchone()
				if which_type_of_id_to_output==1:
					new_strain_acc_list.append(row.ecotype_id)
				else:
					new_strain_acc_list.append(row.accession_id)
				new_category_list.append(row.nativename)
				deletion_code = int(category_list[i])
				deletion_code_index = deletion_code-2
				if deletion_code_index>=0:
					new_data_matrix[i][deletion_code_index] = 2	#allele 2
			except:
				print i, strain_acc
				new_strain_acc_list.append(strain_acc)
				new_category_list.append(strain_acc)
			
		new_header = header[:2] + ['4_268809_0', '4_269962_8']
		write_data_matrix(new_data_matrix, output_fname, new_header, new_strain_acc_list, new_category_list)
	
	"""
	input_fname = os.path.expanduser('~/script/variation/doc/FRI/Shindo2005_data.csv')
	output_fname = os.path.expanduser('~/script/variation/doc/FRI/Shindo2005_data_SNP.tsv')
	outputShindo2005(input_fname, output_fname)
	"""


class Data149_Haplotype(object):
	"""
	2008-09-10
		Alex sent me some files of 149SNP haplotype consensus seq and country info to get pairwise picture. convert it to my format.
	"""
	
	def get_haplotype_group_name2country_ls(file_country):
		import csv, os, sys
		sys.stderr.write("Reading in country info of haplotypes ...")
		reader = csv.reader(open(file_country), delimiter='\t')
		reader.next()
		haplotype_group_name2country_ls = {}
		for row in reader:
			g_name = row[0][:-1]
			if g_name in haplotype_group_name2country_ls:
				sys.stderr.write("Error: %s already existed in haplotype_group_name2country_ls.\n"%g_name)
			haplotype_group_name2country_ls[g_name] = row[1:]
		del reader
		sys.stderr.write("Done.\n")
		return haplotype_group_name2country_ls
		
	def processAlexHaplotypeFile(file_country, file2, output_fname):
		haplotype_group_name2country_ls = get_haplotype_group_name2country_ls(file_country)
		import csv, os, sys
		sys.stderr.write("Reading in consensus of haplotypes ...")
		from pymodule import nt2number
		reader = csv.reader(open(file2), delimiter='\t')
		reader.next()
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		i = 0
		data_row = []
		for row in reader:
			g_name = row[0]
			if g_name in haplotype_group_name2country_ls:	#2008-09-12 has to be in the country file
				country_str = '_'.join(haplotype_group_name2country_ls[g_name])
				label = '%s_%s_%s'%(g_name, country_str, row[-1])
				data_row.append(label)
				data_row.append(label)
				for j in range(1, len(row)-2):
					for k in range(len(row[j])):
						data_row.append(nt2number[row[j][k]])
				if i==0:	#need to write a header
					writer.writerow(['','']+range(1, len(data_row)-1))
				writer.writerow(data_row)
				data_row = []
				i += 1
		sys.stderr.write("Done.\n")
		del reader, writer
	
	
	def getMatrixOutofCrossMatchResult(curs, cross_match_outfile, file_country, output_fname, max_mismatch_rate=1):
		haplotype_group_name2country_ls = get_haplotype_group_name2country_ls(file_country)
		import csv, os, sys
		longitude_g_name_ls = []
		for g_name, country_ls in haplotype_group_name2country_ls.iteritems():
			if len(country_ls)==1:
				curs.execute("select latitude, longitude, abbr from stock.country where abbr='%s'"%(country_ls[0]))
				rows = curs.fetchall()
				longitude = rows[0][1]
				longitude_g_name_ls.append((longitude, country_ls[0], g_name))
			else:
				longitude_g_name_ls.append((-360, '', g_name))
		
		haplotype_group_name2index = {}
		longitude_g_name_ls.sort()
		prev_country = None
		no_of_countries = 0
		label_ls = []
		for i in range(len(longitude_g_name_ls)):
			country = longitude_g_name_ls[i][1]
			if prev_country==None:
				prev_country = country
			elif country!=prev_country:	#insert a separator
				no_of_countries += 1
				g_name = '-%s'%no_of_countries
				haplotype_group_name2index[g_name] = len(haplotype_group_name2index)
				label_ls.append('')
				prev_country = country
			g_name = longitude_g_name_ls[i][-1]
			haplotype_group_name2index[g_name] = len(haplotype_group_name2index)
			label_ls.append(g_name)	#change it to a fuller one later
		
		sys.stderr.write("Getting data matrix ... \n")
		import numpy
		data_matrix = numpy.zeros([len(haplotype_group_name2index), len(haplotype_group_name2index)], numpy.float)
		data_matrix[:,:] = -1	#mark everything as NA
		reader = csv.reader(open(cross_match_outfile), delimiter='\t')
		#figure out which variable is in which column
		header = reader.next()
		col_name2index = {}
		for i in range(len(header)):
			column_name = header[i]
			col_name2index[column_name] = i
		
		for row in reader:
			source_id = row[col_name2index['source_id']]
			source_g_name = source_id.split('_')[0]
			target_id = row[col_name2index['target_id']]
			target_g_name = target_id.split('_')[0]
			mismatch_rate = float(row[col_name2index['mismatch_rate']])
			no_of_mismatches = int(row[col_name2index['no_of_mismatches']])
			no_of_non_NA_pairs = int(row[col_name2index['no_of_non_NA_pairs']])
			row_index = haplotype_group_name2index[source_g_name]
			col_index = haplotype_group_name2index[target_g_name]
			label_ls[row_index] = source_id
			if no_of_non_NA_pairs>=20 and mismatch_rate<=max_mismatch_rate:
				data_matrix[row_index][col_index] = data_matrix[col_index][row_index] = mismatch_rate
		
		for g_name, boundary_index in haplotype_group_name2index.iteritems():
			if g_name[0]=='-':
				data_matrix[boundary_index,:] = -3
				data_matrix[:,boundary_index] = -3
		sys.stderr.write("Done.\n")
		from pymodule import write_data_matrix
		write_data_matrix(data_matrix, output_fname, ['','']+label_ls, label_ls, label_ls)
			
	"""
file_country = os.path.expanduser('~/script/variation/data/149SNP/haps.countries.tsv')
file2 = os.path.expanduser('~/script/variation/data/149SNP/consensus_seqs.HG(0.005.S).tsv')
output_fname = os.path.expanduser('~/script/variation/data/149SNP/haps_in_149.tsv')
processAlexHaplotypeFile(file_country, file2, output_fname)

file_country = os.path.expanduser('~/panfs/149CrossMatch/haps.countries.tsv')
cross_match_outfile = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_cross_match.tsv')
output_fname = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_cross_match_matrix.tsv')
max_mismatch_rate=1
getMatrixOutofCrossMatchResult(curs, cross_match_outfile, file_country, output_fname, max_mismatch_rate=1)

#2008-09-12
file_country = os.path.expanduser('~/script/variation/data/149SNP/haps.countries.x11.tsv')
output_fname = os.path.expanduser('~/script/variation/data/149SNP/haps_in_149_no_4_bad_plates.tsv')
processAlexHaplotypeFile(file_country, file2, output_fname)

file_country = os.path.expanduser('~/banyan_home/script/variation/data/149SNP/haps.countries.x11.tsv')
cross_match_outfile = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_no_4_bad_plates_cross_match.tsv')
output_fname = os.path.expanduser('~/panfs/149CrossMatch/alex_hap_no_4_bad_plates_cross_match_matrix_a0.3.tsv')
getMatrixOutofCrossMatchResult(curs, cross_match_outfile, file_country, output_fname, max_mismatch_rate=0.3)
	"""
	
	"""
	2008-09-11
		found out sequenom group 149SNPs into 4 blocks (38, 37, 37, 37) and then genotype each block on their plate separately
		the SNPs are not in chromosome,position order. It's in id of table snps order.
		now partition the data file into 4 blocks accordingly.
	"""
	def partition149SNPDataInto4Blocks(input_fname, db_149, output_fname_prefix):
		from pymodule import read_data, SNPData, write_data_matrix
		header, strain_acc_list, category_list, data_matrix = read_data(input_fname)
		snpData = SNPData(header=header, strain_acc_list=strain_acc_list, category_list=category_list, data_matrix=data_matrix, turn_into_array=1)
		
		block1_snp_id_ls = []
		block2_snp_id_ls = []
		block3_snp_id_ls = []
		block4_snp_id_ls = []
		sys.stderr.write("Grouping SNPs ...")
		rows = db_149.metadata.bind.execute("select * from snps order by id")
		for row in rows:
			snp_id = '%s_%s'%(row.chromosome, row.position)
			if row.id<=38:
				block1_snp_id_ls.append(snp_id)
			elif row.id>38 and row.id<=75:
				block2_snp_id_ls.append(snp_id)
			elif row.id>75 and row.id<=112:
				block3_snp_id_ls.append(snp_id)
			else:
				block4_snp_id_ls.append(snp_id)
		sys.stderr.write("Done.\n")
		sys.stderr.write("Splitting matrix ... \n")
		no_of_rows = snpData.data_matrix.shape[0]
		import numpy
		block_snp_id_ls_ls = [block1_snp_id_ls, block2_snp_id_ls, block3_snp_id_ls, block4_snp_id_ls]
		for i in range(len(block_snp_id_ls_ls)):
			output_fname = '%s_%s.tsv'%(output_fname_prefix, i)
			no_of_snps_in_this_block = len(block_snp_id_ls_ls[i])
			d_matrix = numpy.zeros([no_of_rows, no_of_snps_in_this_block], numpy.int8)
			for j in range(no_of_snps_in_this_block):
				snp_id = block_snp_id_ls_ls[i][j]
				col_index = snpData.col_id2col_index[snp_id]
				d_matrix[:,j] = snpData.data_matrix[:,col_index]
			header = ['', ''] + block_snp_id_ls_ls[i]
			write_data_matrix(d_matrix, output_fname, header, strain_acc_list, category_list)
		sys.stderr.write("Done.\n")
	"""
input_fname = os.path.expanduser('~/panfs/149CrossMatch/stock_149SNP_y0000110101.tsv')
output_fname_prefix = os.path.expanduser('~/panfs/149CrossMatch/149SNPSequenomBlock')
partition149SNPDataInto4Blocks(input_fname, db_149, output_fname_prefix)
	"""
	
	def getDistanceMatrixOutofCrossMatchFile(input_fname, max_no_of_strains, strainid2index, strainid_ls, min_no_of_non_NA_pairs=10):
		import os, sys,csv
		sys.stderr.write("Getting distance matrix from %s ... "%input_fname)
		reader = csv.reader(open(input_fname), delimiter='\t')
		header = reader.next()
		col_name2index = {}
		for i in range(len(header)):
			column_name = header[i]
			col_name2index[column_name] = i
		import numpy
		data_matrix = numpy.zeros([max_no_of_strains, max_no_of_strains], numpy.float)
		data_matrix[:] = -1	#set default to NA
		i = 0
		for row in reader:
			strainid = int(row[col_name2index['strainid']])
			target_id = int(row[col_name2index['target_id']])
			mismatch_rate = float(row[col_name2index['mismatch_rate']])
			no_of_mismatches = int(row[col_name2index['no_of_mismatches']])
			no_of_non_NA_pairs = int(row[col_name2index['no_of_non_NA_pairs']])
			if no_of_non_NA_pairs>=min_no_of_non_NA_pairs:
				if strainid not in strainid2index:
					strainid2index[strainid] = len(strainid2index)
					strainid_ls.append(strainid)
				if target_id not in strainid2index:
					strainid2index[target_id] = len(strainid2index)
					strainid_ls.append(target_id)
				row_index = strainid2index[strainid]
				col_index = strainid2index[target_id]
				if row_index<max_no_of_strains and col_index<max_no_of_strains:
					data_matrix[row_index][col_index] = mismatch_rate
					data_matrix[col_index][row_index] = mismatch_rate
				else:
					sys.stderr.write("no of strains in this data exceeds max_no_of_strains %s.\n"%(max_no_of_strains))
		del reader
		sys.stderr.write("%s strains. Done.\n"%(len(strainid_ls)))
		return data_matrix
	
	def findStrainShowMsmatchRateVariationIn4Blocks(block_fname_ls, output_fname, max_no_of_strains=7000, min_no_of_non_NA_pairs=10, min_var=0.03):
		import os, sys,csv
		from pymodule import PassingData
		no_of_blocks = len(block_fname_ls)
		strainid2index = {}
		strainid_ls = []
		data_matrix_ls = []
		for i in range(no_of_blocks):
			data_matrix = getDistanceMatrixOutofCrossMatchFile(block_fname_ls[i], max_no_of_strains, strainid2index, strainid_ls, min_no_of_non_NA_pairs)
			data_matrix_ls.append(data_matrix)
		
		sys.stderr.write("Looking for pairs having mismatch_rate var >= %s ..."%(min_var))
		import rpy
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		no_of_strains = len(strainid_ls)
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				mismatch_ls = []
				mismatch_non_NA_ls = []
				for data_matrix in data_matrix_ls:
					if data_matrix[i][j]!=-1:
						mismatch_non_NA_ls.append(data_matrix[i][j])
					mismatch_ls.append(data_matrix[i][j])
				mismatch_var = rpy.r.var(mismatch_non_NA_ls)
				if mismatch_var>=min_var:
					writer.writerow([strainid_ls[i], strainid_ls[j], mismatch_var]+mismatch_ls)
		del writer
	
	"""
block_fname_ls = []
for i in range(4):
	block_fname_ls.append(os.path.expanduser('~/panfs/149CrossMatch/149SNPSequenomBlock_%s_cross_match.tsv'%(i)))

output_fname = os.path.expanduser('~/panfs/149CrossMatch/149SNPSequenomBlock_high_var.tsv')
findStrainShowMsmatchRateVariationIn4Blocks(block_fname_ls, output_fname, max_no_of_strains=7000, min_no_of_non_NA_pairs=10, min_var=0.03)
	"""
	
class CmpDifferentData(object):
	"""
	2009-2-17 like QC
	"""
	
	"""
	2008-10-11
		window average of SNP mismatch rates.
		input_fname is output of TwoSNPData.output_col_id2NA_mismatch_rate_InGWRFormat() (invoked in Qc.py)
	"""
	@classmethod
	def average_SNP_mismatch_rate_within_window(cls, input_fname, output_fname, window_size=100000):
		"""
		"""
		from pymodule import figureOutDelimiter
		import csv
		delimiter = figureOutDelimiter(input_fname)
		reader = csv.reader(open(input_fname), delimiter=delimiter)
		
		chr_pos2mismatch_rate = {}
		for row in reader:
			chr, pos, mismatch_rate = row
			chr = int(chr)
			pos = int(pos)
			mismatch_rate = float(mismatch_rate)
			chr_pos2mismatch_rate[(chr, pos)] = mismatch_rate
		
		prev_chr_pos = None
		start_pos = None
		no_of_mismatches = 0.
		no_of_snps = 0
		chr_pos_ls = chr_pos2mismatch_rate.keys()
		chr_pos_ls.sort()
		chr_start_stop2mismatch_rate_ls = {}
		for chr_pos in chr_pos_ls:
			chr, pos = chr_pos
			if prev_chr_pos==None:
				prev_chr_pos = chr_pos
			if start_pos == None:
				start_pos = pos
			
			if chr==prev_chr_pos[0]:
				if pos-start_pos<=window_size:
					no_of_mismatches += chr_pos2mismatch_rate[chr_pos]
					no_of_snps += 1
				else:	#out of window
					chr_start_stop2mismatch_rate_ls[(prev_chr_pos[0], start_pos, prev_chr_pos[1])] = [no_of_mismatches/no_of_snps, no_of_snps]
					
					#reset
					start_pos = pos
					no_of_mismatches = chr_pos2mismatch_rate[chr_pos]
					no_of_snps = 1
			elif chr!=prev_chr_pos[0]:
				chr_start_stop2mismatch_rate_ls[(prev_chr_pos[0], start_pos, prev_chr_pos[1])] = [no_of_mismatches/no_of_snps, no_of_snps]
				
				
				#reset
				start_pos = pos
				no_of_mismatches = chr_pos2mismatch_rate[chr_pos]
				no_of_snps = 1
			
			prev_chr_pos=chr_pos
		
		writer = csv.writer(open(output_fname, 'w'), delimiter=delimiter)
		chr_start_stop_ls = chr_start_stop2mismatch_rate_ls.keys()
		chr_start_stop_ls.sort()
		for chr_start_stop in chr_start_stop_ls:
			mismatch_rate_ls = chr_start_stop2mismatch_rate_ls[chr_start_stop]
			chr, start_pos, stop_pos = chr_start_stop
			writer.writerow([chr, start_pos, mismatch_rate_ls[0], stop_pos, mismatch_rate_ls[1]])
		del writer
	
	"""
	input_fname = '/tmp/chromosmal_SNP_mismatch_139'
	output_fname = '%s_window_100k'%(input_fname)
	average_SNP_mismatch_rate_within_window(input_fname, output_fname, window_size=100000)
	"""


class DB250k(object):
	@classmethod
	def dumpFTGene2File(cls, db, output_fname, list_type_id=28):
		"""
		2008-11-25 dump all flowering time genes from table ft_gene into file. (for MpiInterGeneSNPPairAsso.py)
		"""
		#rows = db.metadata.bind.execute("select distinct f.gene_id, g.gene_symbol from genome.gene g, ft_gene f where f.gene_id=g.gene_id")
		rows = db.metadata.bind.execute("select distinct l.gene_id, g.gene_symbol from genome.gene g, candidate_gene_list l where l.gene_id=g.gene_id and l.list_type_id=%s order by gene_symbol"%(list_type_id))
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for row in rows:
			writer.writerow([row.gene_id, row.gene_symbol])
		del writer
	
	"""
output_fname = '/tmp/ft_gene.tsv'
dumpFTGene2File(db_250k, output_fname)
	"""
	
	@classmethod
	def outputSNPCandidateGeneAssociation(cls, snps_context_picklef, list_type_id, output_fname):
		"""
		2008-11-13
		output SNP-gene association (from a file containing the pickled snps_context_wrapper) into a file
		for magnus
		"""
		import cPickle, csv
		from GeneListRankTest import GeneListRankTest
		candidate_gene_set = GeneListRankTest.dealWithCandidateGeneList(list_type_id, return_set=True)
		from Stock_250kDB import SnpsContext
		
		picklef = open(snps_context_picklef)
		snps_context_wrapper = cPickle.load(picklef)
		del picklef
		
		chr_pos_ls = snps_context_wrapper.chrpos2snps_id.keys()
		chr_pos_ls.sort()
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['gene_id', 'snps_id', 'chromosome', 'position', 'disp_pos', 'left_or_right', 'disp_pos_comment']
		writer.writerow(header)
		for chr, pos in chr_pos_ls:
			snps_context_matrix = snps_context_wrapper.returnGeneLs(chr, pos)
			assign_snp_candidate_gene = 0
			assign_snp_non_candidate_gene = 0
			for snps_context in snps_context_matrix:
				snps_id, disp_pos, gene_id = snps_context
				if gene_id in candidate_gene_set:
					snps_context = SnpsContext.query.filter_by(snps_id=snps_id).filter_by(gene_id=gene_id).first()
					left_or_right = getattr(snps_context, 'left_or_right', '')
					disp_pos_comment = getattr(snps_context, 'disp_pos_comment', '')
					writer.writerow([gene_id, snps_id, chr, pos, disp_pos, left_or_right, disp_pos_comment])
		del writer
	
	"""
snps_context_picklef = './mnt2/panfs/250k/snps_context_g0_m500'
list_type_id = 28
output_fname = '/tmp/snps_context_g0_m500_list28.tsv'
outputSNPCandidateGeneAssociation(snps_context_picklef, list_type_id, output_fname)
	"""
	
	
	"""
	2008-04-11
		this is a one-time 250k pipeline fix to avoid memory-hefty intensity re-output.
		
		form 2 lists of intensity matrix filenames based on the new and old array_info_table
		output them as two-column (old_fname, new_fname) into a output_fname.
		shell/file_batch_move.py reads the output_fname and handle name changing.
	
	2008-04-11
		turns out to be useless because the header in each intensity matrix file has array_id embedded.
		have to output each array into intensity matrix.
	"""
	@classmethod
	def output_intensity_fname(cls, curs, new_array_info_table, old_array_info_table, output_fname):
		"""
		"""
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		curs.execute("select n.id, o.id from %s n, %s o where n.original_filename=o.filename"%(new_array_info_table, old_array_info_table))
		rows = curs.fetchall()
		for row in rows:
			new_array_id, old_array_id = row
			new_fname = '%s_array_intensity.tsv'%new_array_id
			old_fname = '%s_array_intensity.tsv'%old_array_id
			writer.writerow([old_fname, new_fname])
		del writer
	
	"""
	new_array_info_table='stock_250k.array_info'
	old_array_info_table='stock_250k.array_info_2008_04_11'
	output_fname = '/tmp/intensity_fname.rename.tsv'
	output_intensity_fname(curs, new_array_info_table, old_array_info_table, output_fname)
	"""
	
	
	"""
	2008-05-31
		output results from stock_250k.results, totally db based, not reading results.filename
	"""
	@classmethod
	def outputResults(cls, db, results_method_id, output_fname):
		import csv
		import sqlalchemy as sql
		
		conn = db.connection	#establish the connection before referring db.tables (in case it hasn't been setup)
		sys.stderr.write("Getting marker_id2position ... ")
		marker_id2position = {}
		snps_table = db.tables['snps'].alias()
		results = conn.execute(sql.select([snps_table.c.id, snps_table.c.chromosome, snps_table.c.position, snps_table.c.end_position]))
		for row in results:
			marker_id2position[row.id] = (row.chromosome, row.position, row.end_position)
		del results
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Outputting results_method_id=%s ... "%results_method_id)
		results_table = db.tables['results'].alias()
		results = conn.execute(sql.select([results_table.c.snps_id, results_table.c.score], results_table.c.results_method_id==results_method_id,\
										  order_by=[results_table.c.snps_id]))
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['chromosome', 'position', 'score'])
		for row in results:
			chromosome, position, end_position = marker_id2position[row.snps_id]
			writer.writerow([chromosome, position, row.score])
		del writer
		sys.stderr.write("Done.\n")
	
	"""
	from variation.src.db import Stock_250kDatabase
	db = Stock_250kDatabase(username='nordborglab',
					   password='papaya', hostname='papaya.usc.edu', database='stock_250k')
	conn = db.connection	#establish the connection before referring db.tables (it needs to be setup)
	outputResults(db, 5, '/tmp/5_results.tsv')
	outputResults(db, 6, '/tmp/6_results.tsv')
	"""
	
	"""
	2008-07-21
	"""
	@classmethod
	def drawScoreHistogram(cls, curs, results_method_id, list_type_id=1, do_log10_transformation=True, min_or_max_func='min'):
		from Stock_250kDB import ResultsMethod, GeneList, ResultsByGene, CandidateGeneRankSumTestResult
		rm = ResultsMethod.get(results_method_id)
		db_rows = GeneList.query.filter_by(list_type_id=list_type_id)
		from sets import Set
		gene_set = Set()
		for gl in db_rows:
			gene_set.add(gl.gene_id)
		
		import math
		score_ls1 = []
		score_ls2 = []
		#db_rows = ResultsByGene.query.filter_by(results_method_id=results_method_id)
		curs.execute("select r.gene_id, %s(r.score) as score from results_by_gene r  where r.results_method_id=%s group by r.gene_id"%(min_or_max_func, results_method_id))
		db_rows = curs.fetchall()
		for row in db_rows:
			gene_id, score = row
			if do_log10_transformation:
				score = -math.log(score)
			if gene_id in gene_set:
				score_ls1.append(score)
			else:
				score_ls2.append(score)
		import pylab
		pylab.clf()
		pylab.title('results_method_id=%s, (%s on %s) by list type: %s.'%(rm.id, rm.analysis_method.short_name, rm.phenotype_method.short_name, gl.list_type.short_name))
		n1 = pylab.hist(score_ls1, 100, alpha=0.4, normed=1)
		n2 = pylab.hist(score_ls2, 100, alpha=0.4, normed=1, facecolor='r')
		pylab.legend([n1[2][0], n2[2][0]], ['candidate gene list', 'non-candidate gene list'])
		pylab.show()
		return score_ls1, score_ls2
	
	"""
	score_ls1, score_ls2 = drawScoreHistogram(curs, 23, 1)
	"""
	
	@classmethod
	def plotArrayMismatchVsMedianIntensity(cls, db, array_id_ls, take_log=False):
		"""
		2009-3-11		
			plot array mismatch rate vs median intensity of all probes
			
			array_id_ls is comma or dash-separated list of array ids
		"""
		import math
		if array_id_ls:
			from pymodule import getListOutOfStr
			array_id_ls = getListOutOfStr(array_id_ls, data_type=str)
		else:
			array_id_ls = range(500)
			array_id_ls = map(str, array_id_ls)
		rows = db.metadata.bind.execute("select mismatch_rate, median_intensity from stock_250k.view_qc where array_id in (%s)"%(','.join(array_id_ls)))
		y_ls = []
		x_ls = []
		for row in rows:
			if take_log:
				median_intensity = math.log10(row.median_intensity)
			else:
				median_intensity = row.median_intensity
			x_ls.append(median_intensity)
			y_ls.append(row.mismatch_rate)
		
		import pylab
		pylab.plot(x_ls, y_ls, '.', alpha=0.5)
		pylab.title('mismatch rate vs median_intensity of arrays')
		pylab.xlabel('median_intensity')
		pylab.ylabel('mismatch rate')
		pylab.show()
	
	"""
plotArrayMismatchVsMedianIntensity(db_250k, '616-1045')
	"""
	
	@classmethod
	def saveArrayQCIntoTableFile(cls, db, array_id_ls, output_fname, take_log=False):
		"""
		2009-3-13		
			output array QC information in a matrix format which pymodule/DataMatrixGuiXYProbe.py can recognize and do clickable scatterplots 
		"""
		import math
		if array_id_ls:
			from pymodule import getListOutOfStr
			array_id_ls = getListOutOfStr(array_id_ls, data_type=str)
		else:
			array_id_ls = range(500)
			array_id_ls = map(str, array_id_ls)
		rows = db.metadata.bind.execute("select * from stock_250k.view_qc where array_id in (%s)"%(','.join(array_id_ls)))
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['nativename', 'original_array_name', 'ecotype_id', 'array_id', 'mismatch_rate', 'median_intensity', \
					'QC_NA_rate', 'no_of_mismatches', 'no_of_non_NA_pairs', 'qc_method_id', 'call_info_id', 'array_original_filename']
		writer.writerow(header)
		for row in rows:
			if take_log:
				median_intensity = math.log10(row.median_intensity)
			else:
				median_intensity = row.median_intensity
			original_array_name = os.path.basename(row.array_original_filename)
			output_row = []
			for col_name in header:
				if col_name=='original_array_name':
					output_row.append(original_array_name)
				elif col_name=='median_intensity':
					output_row.append(median_intensity)
				else:
					output_row.append(getattr(row, col_name))
			writer.writerow(output_row)
		del writer
		
	"""
output_fname = '/tmp/CHLA_2009_01_QC.tsv'
saveArrayQCIntoTableFile(db_250k, '616-1045', output_fname)
	"""

	"""
	2009-4-13 function to check all 250k calls to see which is bad , which is good.
		for bad ones, further check cross_match results to see if mis-labelling happens
		
		check view_qc
		check qc_cross_match_table to see if any cross-labeling
			no_of_non_NA_pairs>=40
			mismatch_rate<2%
		
		quality meaning:
			0: bad
			1: good
			2: mis-labelled
			
	"""
	@classmethod
	def reportBadArrays(cls, db, call_method_id, qc_method_id, output_fname, max_mismatch_rate=0.1, min_no_of_non_NA_pairs=40, \
					max_mislabel_mismatch_rate=0.02,\
					view_qc_table='view_qc', qc_cross_match_table='qc_cross_match'):
		from variation.src.common import get_ecotypeid2nativename
		ecotypeid2nativename = get_ecotypeid2nativename(db.metadata.bind, ecotype_table='stock.ecotype')
		sys.stderr.write("Reporting 250k arrays ... \n")
		rows = db.metadata.bind.execute("select * from %s where qc_method_id=%s and call_method_id=%s order by nativename"%\
									(view_qc_table, qc_method_id, call_method_id))
		
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header_row = ['ecotype_id', 'nativename', 'array_original_filename', 'array_created', 'quality', 'true ecotype_id', 'true nativename']
		writer.writerow(header_row)
		results = []
		counter = 0
		no_of_bad_ones = 0
		no_of_mis_labelled_ones = 0
		
		good_ecotype_id_set = set()
		bad_ecotype_id_set = set()
		true_ecotype_id_set = set()
		
		for row in rows:
			quality = 0
			counter += 1
			true_ecotype_id_ls = []
			true_nativename_ls = []
			if row.no_of_non_NA_pairs>=min_no_of_non_NA_pairs and row.mismatch_rate<=max_mismatch_rate:
				quality = 1
			else:
				no_of_bad_ones += 1
				cross_match_rows = db.metadata.bind.execute("select * from %s where call_info_id=%s and no_of_non_NA_pairs>=%s order by mismatch_rate"%\
										(qc_cross_match_table, row.call_info_id, min_no_of_non_NA_pairs))
				for cross_match_row in cross_match_rows:
					if cross_match_row.mismatch_rate<=max_mislabel_mismatch_rate:
						if cross_match_row.vs_ecotype_id in ecotypeid2nativename:	#make sure ecotypeid is still alive
							true_ecotype_id_ls.append(cross_match_row.vs_ecotype_id)
							true_nativename_ls.append(ecotypeid2nativename[cross_match_row.vs_ecotype_id])
				if len(true_ecotype_id_ls)>0:
					quality = 2
					no_of_mis_labelled_ones += 1
			if quality==1:
				good_ecotype_id_set.add(row.ecotype_id)
			elif quality==0 or quality==2:
				bad_ecotype_id_set.add(row.ecotype_id)
			
			if quality==2:
				true_ecotype_id_set.update(set(true_ecotype_id_ls))
			
			true_ecotype_id_ls = map(str, true_ecotype_id_ls)	#in order to use ','.join()
			output_row = [row.ecotype_id, row.nativename, row.array_original_filename, row.array_created, quality, \
						','.join(true_ecotype_id_ls), ','.join(true_nativename_ls)]
			writer.writerow(output_row)
			if counter%100==0:
				sys.stderr.write("%s\t%s\t%s\t%s"%('\x08'*80, no_of_mis_labelled_ones, no_of_bad_ones, counter))
		del writer
		
		rescued_bad_ecotype_id_set = set()
		for ecotype_id in bad_ecotype_id_set:
			if ecotype_id in good_ecotype_id_set or ecotype_id in true_ecotype_id_set:
				rescued_bad_ecotype_id_set.add(ecotype_id)
		
		# assign quality 2 to entries in true_ecotype_id_set if the ecotype_id is not in good_ecotype_id_set
		true_ecotype_id_quality_ls = []
		for ecotype_id in true_ecotype_id_set:
			if ecotype_id not in good_ecotype_id_set:
				true_ecotype_id_quality_ls.append([ecotype_id, 2])
		# remove the rescued bad ones from bad_ecotype_id_set
		bad_ecotype_id_set.difference_update(rescued_bad_ecotype_id_set)
		
		output_fname = '%s_simple%s'%(os.path.splitext(output_fname)[0], os.path.splitext(output_fname)[1])
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header_row = ['ecotype_id', 'nativename', 'quality']
		writer.writerow(header_row)
		for ecotype_id in good_ecotype_id_set:
			writer.writerow([ecotype_id, ecotypeid2nativename[ecotype_id], 1])
		for ecotype_id, quality in true_ecotype_id_quality_ls:
			writer.writerow([ecotype_id, ecotypeid2nativename[ecotype_id], quality])
		
		for ecotype_id in bad_ecotype_id_set:
			writer.writerow([ecotype_id, ecotypeid2nativename[ecotype_id], 0])
		del writer
		sys.stderr.write("%s mis-labelled out of %s bad ones, out of %s in total. Done.\n"%(no_of_mis_labelled_ones, no_of_bad_ones, counter))
		
	"""
call_method_id = 3
qc_method_id = 9
output_fname = '/tmp/250k_good_bad_cross_label_arrays.tsv'
reportBadArrays(db_250k, call_method_id, qc_method_id, output_fname)
	"""
	
	@classmethod
	def outputProbeIntensityOfOneSNPAcrossArrays(cls, db_250k, intensity_matrix_fname, snpData, snp_id, output_fname):
		"""
		2009-5-19
			intensity_matrix_fname is output of affy/CelQuantileNorm/gtype_cel_to_pq (plone doc: log/microarray)
				format is SNP X array. 1st row is name for the array and allele, with array_id embedded, like 67_raw_data_A.
					one array occupies two columns, allele A and allele B.
					1st column is snp_id (chr_pos like  1_657).
			
		"""
		sys.stderr.write("Outputting probe intensity of one SNP ... ")
		
		array_id2ecotype_id_name_ls = {}
		rows = db_250k.metadata.bind.execute("select * from view_array")
		for row in rows:
			array_id2ecotype_id_name_ls[row.array_id] = [row.maternal_ecotype_id, row.maternal_nativename]
		
		import csv
		reader = csv.reader(open(intensity_matrix_fname), delimiter='\t')
		header = reader.next()
		array_id_allele_label_ls = header[1:]
		nativename_ls = []
		ecotype_id_ls = []
		array_id_ls = []
		
		for i in range(0, len(array_id_allele_label_ls), 2):	# every other column
			label = array_id_allele_label_ls[i]
			array_id = int(label.split('_')[0])
			array_id_ls.append(array_id)
			ecotype_id_nativename_ls = array_id2ecotype_id_name_ls.get(array_id)
			if ecotype_id_nativename_ls is not None:
				ecotype_id, nativename = ecotype_id_nativename_ls
			else:
				ecotype_id = 0
				nativename = 'None'
			ecotype_id_ls.append(ecotype_id)
			nativename_ls.append(nativename)
		
		writer = csv.writer(open(output_fname,'w'), delimiter='\t')
		writer.writerow(['nativename', 'ecotype_id', 'array_id', 'alleleA-intensity', 'alleleB-intensity', 'call'])
		for row in reader:
			if row[0]==snp_id:
				intensity_ls = row[1:]
				for i in range(0, len(intensity_ls), 2):
					alleleA_intensity = intensity_ls[i]
					alleleB_intensity = intensity_ls[i+1]
					array_index = i/2
					ecotype_id = ecotype_id_ls[array_index]
					snp_row_index = snpData.row_id2row_index['%s'%ecotype_id]
					snp_col_index = snpData.col_id2col_index[snp_id]
					snp_allele = snpData.data_matrix[snp_row_index][snp_col_index]
					output_row = [nativename_ls[array_index], ecotype_id_ls[array_index], array_id_ls[array_index], alleleA_intensity, \
								alleleB_intensity, snp_allele]
					writer.writerow(output_row)
				break
		del writer, reader
		sys.stderr.write("Done.\n")
		
	"""
intensity_matrix_fname = os.path.expanduser('~/script/affy/250k_test/call_method_32_arrays.QN.tsv')
snp_id = '2_13265102'
output_fname = '/tmp/%s_intensity_matrix.tsv'%snp_id

from pymodule import SNPData
snp_fname = '/Network/Data/250k/db/dataset/call_method_32.tsv' 
snpData = SNPData(input_fname=snp_fname, turn_into_array=1, ignore_2nd_column=1)
DB250k.outputProbeIntensityOfOneSNPAcrossArrays(db_250k, intensity_matrix_fname, snpData, snp_id, output_fname)
	"""
	
	@classmethod
	def cleanUpTablePhenotype(cls, db, make_replicate_continuous=False):
		"""
		2009-8-25
			add an option (make_replicate_continuous) to make the replicate variable continuous, like 1,2,5,7 => 1,2,3,4.
			The unique key in phenotype table must be removed before applying this option. 
		2009-8-12
			increase replicate from one of the identical entries based on (method_id, ecotype_id, replicate)
		"""
		sys.stderr.write("cleaning up table phenotype ...\n")
		db.session.begin()
		
		sys.stderr.write("Getting data from db ...")
		import Stock_250kDB
		
		method_id_ecotype_id2db_id_ls = {}
		method_id_ecotype_id_replicate2count = {}
		method_id_ecotype_id2max_replicate = {}	#2009-8-25 added to find out (method_id, ecotype_id)s that have the discontinuous replicate numbers
		
		rows = Stock_250kDB.Phenotype.query()
		for row in rows:
			method_id_ecotype_id = (row.method_id, row.ecotype_id)
			method_id_ecotype_id_replicate = (row.method_id, row.ecotype_id, row.replicate)
			if method_id_ecotype_id not in method_id_ecotype_id2db_id_ls:
				method_id_ecotype_id2db_id_ls[method_id_ecotype_id] = []
			method_id_ecotype_id2db_id_ls[method_id_ecotype_id].append(row.id)
			
			if method_id_ecotype_id_replicate not in method_id_ecotype_id_replicate2count:
				method_id_ecotype_id_replicate2count[method_id_ecotype_id_replicate] = 0
			method_id_ecotype_id_replicate2count[method_id_ecotype_id_replicate] += 1
			
			#2009-8-25
			if method_id_ecotype_id not in method_id_ecotype_id2max_replicate:
				method_id_ecotype_id2max_replicate[method_id_ecotype_id] = 0
			if row.replicate>method_id_ecotype_id2max_replicate[method_id_ecotype_id]:
				method_id_ecotype_id2max_replicate[method_id_ecotype_id] = row.replicate
		sys.stderr.write("Done.\n")
		
		method_id_ecotype_id_need_fix_set = set()
		for method_id_ecotype_id_replicate, count in method_id_ecotype_id_replicate2count.iteritems():
			if count>1:
				method_id_ecotype_id = (method_id_ecotype_id_replicate[0], method_id_ecotype_id_replicate[1])
				method_id_ecotype_id_need_fix_set.add(method_id_ecotype_id)
		
		#2009-8-25
		if make_replicate_continuous:
			for method_id_ecotype_id, max_replicate in method_id_ecotype_id2max_replicate.iteritems():
				if max_replicate>len(method_id_ecotype_id2db_id_ls[method_id_ecotype_id]):
					method_id_ecotype_id_need_fix_set.add(method_id_ecotype_id)
		
		sys.stderr.write("%s (method_id, ecotype_id) pairs to be fixed.\n"%len(method_id_ecotype_id_need_fix_set))
		
		#re-assign the replicate value in the order of db ids
		for method_id_ecotype_id in method_id_ecotype_id_need_fix_set:
			db_id_ls = method_id_ecotype_id2db_id_ls.get(method_id_ecotype_id)
			if db_id_ls is None:
				continue
			
			db_id_ls.sort()
			for i in range(len(db_id_ls)):
				db_id = db_id_ls[i]
				db_entry = Stock_250kDB.Phenotype.get(db_id)
				db_entry.replicate = i +1
				db.session.save_or_update(db_entry)
				db.session.flush()
		db.session.commit()
		sys.stderr.write("Done.\n")
	
	"""	
DB250k.cleanUpTablePhenotype(db_250k)

DB250k.cleanUpTablePhenotype(db_250k, make_replicate_continuous=True)
	"""
	
	@classmethod
	def updatePhenotypeAvgEntry(cls, db, db_entry, individual_value_ls):
		"""
		2009-8-13
			
		"""
		import numpy
		db_entry.value = numpy.average(individual_value_ls)
		db_entry.sample_size = len(individual_value_ls)
		if db_entry.sample_size>1:
			db_entry.stddev = numpy.std(individual_value_ls) 
		db.session.save_or_update(db_entry)
		
	
	@classmethod
	def cleanUpTablePhenotypeAvg(cls, db):
		"""
		2009-8-12
			enforce the unique constraint, (method_id, ecotype_id) 
			
		"""
		sys.stderr.write("Cleaning up table phenotype_avg ... \n")
		db.session.begin()
		
		sys.stderr.write("Getting data from db ...")
		import Stock_250kDB
		
		method_id_ecotype_id2db_id_ls = {}
		
		rows = Stock_250kDB.PhenotypeAvg.query()
		for row in rows:
			method_id_ecotype_id = (row.method_id, row.ecotype_id)
			if method_id_ecotype_id not in method_id_ecotype_id2db_id_ls:
				method_id_ecotype_id2db_id_ls[method_id_ecotype_id] = []
			method_id_ecotype_id2db_id_ls[method_id_ecotype_id].append(row.id)
		
		sys.stderr.write("Done.\n")
		import numpy
		no_of_replicate_entries = 0
		no_of_replicate_entries_but_no_individual_data = 0
		no_of_replicate_entries_with_individual_data = 0
		for method_id_ecotype_id, db_id_ls in method_id_ecotype_id2db_id_ls.iteritems():
			if len(db_id_ls)>1:
				no_of_replicate_entries  += 1
				method_id, ecotype_id = method_id_ecotype_id
				rows = Stock_250kDB.Phenotype.query.filter_by(method_id=method_id).filter_by(ecotype_id=ecotype_id)
				db_id_ls.sort()
				if rows.count()>0:
					no_of_replicate_entries_with_individual_data += 1
					for i in range(len(db_id_ls)-1):	#delete all but the last one
						db_id = db_id_ls[i] 
						db_entry = Stock_250kDB.PhenotypeAvg.get(db_id)
						
						db.session.delete(db_entry)
					#update the last PhenotypeAvg entry
					db_entry = Stock_250kDB.PhenotypeAvg.get(db_id_ls[-1])
					individual_value_ls = [row.value for row in rows]
					cls.updatePhenotypeAvgEntry(db, db_entry, individual_value_ls)
				else:
					no_of_replicate_entries_but_no_individual_data += 1
					individual_value_ls = []
					for i in range(len(db_id_ls)-1):	#insert into Stock_250kDB.Phenotype and delete all but the last one from Stock_250kDB.PhenotypeAvg
						db_id = db_id_ls[i] 
						db_entry = Stock_250kDB.PhenotypeAvg.get(db_id)
						individual_value_ls.append(db_entry.value)
						phenotype_entry = Stock_250kDB.Phenotype(method_id=db_entry.method_id, ecotype_id=db_entry.ecotype_id, \
																value=db_entry.value, replicate=i+1)
						db.session.save(phenotype_entry)
						db.session.delete(db_entry)
					#insert the last phenotype_avg into Stock_250kDB.Phenotype and update it in  phenotype_avg
					db_entry = Stock_250kDB.PhenotypeAvg.get(db_id_ls[-1])
					individual_value_ls.append(db_entry.value)
					phenotype_entry = Stock_250kDB.Phenotype(method_id=db_entry.method_id, ecotype_id=db_entry.ecotype_id, \
																value=db_entry.value, replicate=len(db_id_ls))
					db.session.save(phenotype_entry)
					cls.updatePhenotypeAvgEntry(db, db_entry, individual_value_ls)
					#insert data from db_id_ls into Stock_250kDB.Phenotype
		db.session.commit()
		sys.stderr.write("%s total replicate entries. %s has individual data. %s has no individual data. Done.\n"%\
						(no_of_replicate_entries, no_of_replicate_entries_with_individual_data, no_of_replicate_entries_but_no_individual_data))
		
	"""
DB250k.cleanUpTablePhenotypeAvg(db)
	"""
	
	@classmethod
	def updatePhenotypeAvgBasedOnPhenotype(cls, db, phenotype_condition=None):
		"""
		2009-8-25
			After values in table phenotype are updated, corresponding ones in table phenotype_avg shall be updated too.
		"""
		sys.stderr.write("Updating phenotype_avg entries based on the ones in table phenotype ...\n")
		db.session.begin()
		
		method_id_ecotype_id2individual_value_ls = {}
		rows = db.metadata.bind.execute("select method_id, ecotype_id, value from phenotype %s"%phenotype_condition)
		replicate_count = 0
		for row in rows:
			method_id_ecotype_id = (row.method_id, row.ecotype_id)
			if method_id_ecotype_id not in method_id_ecotype_id2individual_value_ls:
				method_id_ecotype_id2individual_value_ls[method_id_ecotype_id] = []
			method_id_ecotype_id2individual_value_ls[method_id_ecotype_id].append(row.value)
			replicate_count += 1
		
		sys.stderr.write("need to update %s phenotype_avg entries.\n"%(len(method_id_ecotype_id2individual_value_ls)))
		
		avg_count = 0
		new_avg_count = 0
		multi_avg_count = 0
		import Stock_250kDB
		for method_id_ecotype_id, individual_value_ls in method_id_ecotype_id2individual_value_ls.iteritems():
			
			method_id, ecotype_id = method_id_ecotype_id
			#rows = db.metadata.bind.execute("select ecotype_id, avg(value) as avg_value, stddev(value) as stddev, count(value) as sample_size, \
			#	method_id from phenotype where \
			#	method_id=%s and ecotype_id=%s group by method_id, ecotype_id"%(method_id, ecotype_id))
			rows = Stock_250kDB.PhenotypeAvg.query.filter_by(method_id=method_id).filter_by(ecotype_id=ecotype_id)
			if rows.count()==1:
				phenotype_avg_entry = rows.one()
			elif rows.count()==0:
				phenotype_avg_entry = Stock_250kDB.PhenotypeAvg(method_id=method_id, ecotype_id=ecotype_id)
				new_avg_count += 1
			else:
				multi_avg_count += 1
				sys.stderr.write("Error: method_id=%s, ecotype_id=%s has multiple phenotype_avg entries.\n"%(method_id, ecotype_id))
				phenotype_avg_entry = rows.first()
				sys.exit(3)
			DB250k.updatePhenotypeAvgEntry(db, phenotype_avg_entry, individual_value_ls)
			#phenotype_avg_entry.value = row.avg_value
			#phenotype_avg_entry.stddev = row.stddev
			#phenotype_avg_entry.sample_size = row.sample_size
			#db.session.save_or_update(phenotype_avg_entry)
			avg_count += 1
		db.session.commit()
		sys.stderr.write("%s(%s new, %s multi) average values from %s replicates were updated. Done.\n"%\
						(avg_count, new_avg_count, multi_avg_count, replicate_count))
	
	"""
DB250k.updatePhenotypeAvgBasedOnPhenotype(db_250k, phenotype_condition='where date_created>"2009-08-23" order by method_id, ecotype_id')

#update all phenotype_avg entries that have inidividual values in table phenotype 
DB250k.updatePhenotypeAvgBasedOnPhenotype(db_250k);
	"""	
	
class CNV(object):
	@classmethod
	def drawCNVAmpHist(cls, snpData, output_dir, max_counter=None, start_pos=None, stop_pos=None, chromosome_chosen=None):
		"""
		2009-2-20
			add arguments start_pos, stop_pos
		2008-12-12 draw histogram of CNV amplitudes probe by probe. input file is the amplitude output of RunGADA.py.
		"""
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		import pylab
		counter = 0
		for col_id in snpData.col_id_ls:
			col_id_split = col_id.split('_')
			col_id_split = map(int, col_id_split)
			chromosome, position = col_id_split
			
			#filter according to chromosome
			if chromosome_chosen is not None and chromosome!=chromosome_chosen:
				continue
			#filter according to position
			if start_pos is not None and stop_pos is not None and (position<start_pos-12 or position>stop_pos+12):
				continue
			
			col_index = snpData.col_id2col_index[col_id]
			sys.stderr.write("%s\t%s"%('\x08'*20, counter))
			output_fname = os.path.join(output_dir, '%s_amp_hist.png'%col_id)
			amp_ls = snpData.data_matrix[:,col_index]
			pylab.clf()
			pylab.hist(amp_ls, 40, alpha=0.6)
			pylab.title(col_id)
			pylab.xlabel('amplitude')
			pylab.ylabel('frequency')
			pylab.xlim([-1,1])
			pylab.savefig(output_fname, dpi=200)
			counter += 1
			if max_counter and counter>max_counter:
				break
	"""
input_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_17_CNV_array_intensity_chr4_line_no_888148_1107622_norm_GADA_out_amp.tsv'
from pymodule import SNPData
snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1, ignore_2nd_column=1, matrix_data_type=float)
output_dir = '/Network/Data/250k/tmp-yh/CNV/amp_hist/'
CNV.drawCNVAmpHist(snpData, output_dir, max_counter=1000)
	"""
	
	@classmethod
	def outputCNVMatrixIntoMemmapFormat(cls, input_fname, output_fname):
		"""
		2009-9-27
			input format is output of CNVNormalize.py & DB_250k2Array.py
			output is binary, probe X sample. test if matlab sees it as a memmapfile.
		"""
		from variation.src.CNVNormalize import CNVNormalize
		data_matrix, probe_id_ls, chr_pos_ls, header = CNVNormalize.get_input(input_fname)
		outf = open(output_fname, 'wb')
		import numpy
		data_matrix = numpy.transpose(data_matrix)	#somehow, without this, matlab will read it as transposed to the original matrix 
		data_matrix.tofile(outf, format='%.8f')	# format seems to not matter.
		"""
		no_of_rows, no_of_cols = data_matrix.shape
		for i in range(no_of_rows):
			for j in range(no_of_cols):
				outf.write(struct.pack('d', data_matrix[i,j])
		"""
		del outf
	"""
	input_fname = os.path.expanduser('~/panfs/250k/CNV/call_method_17_CNV_array_inensity_norm_chr4_n101.tsv')
	output_fname = os.path.expanduser('~/panfs/250k/CNV/call_method_17_CNV_array_intensity_norm_chr4_n101.memmap')
	CNV.outputCNVMatrixIntoMemmapFormat(input_fname, output_fname)
	
	for i in range(2,6):
		input_fname = os.path.expanduser('~/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.tsv'%i)
		output_fname = os.path.expanduser('~/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.memmap'%i)
		CNV.outputCNVMatrixIntoMemmapFormat(input_fname, output_fname)
	"""
	
	@classmethod
	def outputSelectedArraysInGWAFormat(cls, input_fname, output_dir, array_id_ls):
		"""
		2009-10-11
			input_fname is DB_250k2Array.py's output.
			output the intensity from specified arrays, one array one file to inspect genome-wide pattern:
				chromosome	position	amplitude 
		"""
		sys.stderr.write("Outputting selected arrays in GWA format ...\n")
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		import csv
		from pymodule import figureOutDelimiter
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		
		
		array_id_set = set(array_id_ls)
		header = reader.next()
		array_index2writer = {}
		for i in range(1, len(header)-2):
			array_id = int(header[i])
			if array_id in array_id_set:
				output_fname = os.path.join(output_dir, 'array_%s_CNV.tsv'%array_id)
				writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				array_index2writer[i] = writer
		
		counter = 0
		for row in reader:
			chr, pos = row[-2:]
			for i in range(1, len(row)-2):
				if i in array_index2writer:
					new_row = [chr, pos, row[i]]
					array_index2writer[i].writerow(new_row)
			counter += 1
			if counter%10000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		del reader, array_index2writer
		sys.stderr.write("Done.\n")
	
	
	"""
	Col_array_id_ls = [1, 2, 43, 139, 145]	# Col-0 arrays
	Ler_array_id_ls = [3, 4, 41, 150, 151,]
	array_id_ls = Col_array_id_ls + Ler_array_id_ls
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_48_CNV_QNormalize.tsv')
	output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV_QNormalize')
	CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')
	output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV')
	CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_intensity.tsv')
	output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_43_CNV')
	CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	
	
	for i in range(1,6):
		input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.tsv'%i)
		output_dir = os.path.expanduser('~/mnt/banyan/tmp/call_43_CNV_norm_intensity_chr%s'%i)
		CNV.outputSelectedArraysInGWAFormat(input_fname, output_dir, array_id_ls)
	"""
	
	@classmethod
	def outputMedianIntensityOfSelectedArraysInGWAFormat(cls, input_fname, output_fname, array_id_ls):
		"""
		2009-10-27
			similar to outputSelectedArraysInGWAFormat() but take median of probes across all arrays
			input_fname is DB_250k2Array.py's output.
			output the intensity from specified arrays, one array one file to inspect genome-wide pattern:
				chromosome	position	amplitude 
		"""
		sys.stderr.write("Outputting selected arrays in GWA format ...\n")
		
		import csv, numpy
		from pymodule import figureOutDelimiter
		reader = csv.reader(open(input_fname), delimiter=figureOutDelimiter(input_fname))
		
		
		array_id_set = set(array_id_ls)
		header = reader.next()
		array_index2writer = {}
		for i in range(1, len(header)-2):
			array_id = int(header[i])
			if array_id in array_id_set:
				#output_fname = os.path.join(output_dir, 'array_%s_CNV.tsv'%array_id)
				#writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
				array_index2writer[i] = None
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		counter = 0
		for row in reader:
			chr, pos = row[-2:]
			intensity_ls = []
			for i in range(1, len(row)-2):
				if i in array_index2writer:
					intensity_ls.append(float(row[i]))
			median = numpy.median(intensity_ls)
			new_row = [chr, pos, median]
			writer.writerow(new_row)
			counter += 1
			if counter%10000==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
		del reader, array_index2writer, writer
		sys.stderr.write("Done.\n")
		
	"""
	Col_array_id_ls = [1, 2, 43, 139, 145]	# Col-0 arrays
	Ler_array_id_ls = [3, 4, 41, 150, 151,]
	
	input_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_48_CNV_QNormalize.tsv')
	output_fname = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV_QNormalize_Col_median.tsv')
	CNV.outputMedianIntensityOfSelectedArraysInGWAFormat(input_fname, output_fname, Col_array_id_ls)
	output_fname = os.path.expanduser('~/mnt/banyan/tmp/call_48_CNV_QNormalize_Ler_median.tsv')
	CNV.outputMedianIntensityOfSelectedArraysInGWAFormat(input_fname, output_fname, Ler_array_id_ls)
	
	"""
	
	@classmethod
	def compareCNVSegmentsAgainstQCHandler(cls, input_fname_ls, ecotype_id2cnv_qc_call_data, function_handler, param_obj, \
										deletion_cutoff=None, max_boundary_diff=10000, \
										max_diff_perc=0.10, min_no_of_probes=5, count_embedded_segment_as_match=False, \
										min_reciprocal_overlap=0.6, report=True):
		"""
		2010-1-26
			value of the ecotype_id2cnv_qc_call_data dictionary is a RBDict (RBTree dictionary) structure.
		2009-12-8
			add argument min_reciprocal_overlap
		2009-11-4
			a general handler to compare CNV segments from input_fname_ls with cnv_qc_call_data.
			Upon a match between a CNV segment from input_fname_ls and cnv_qc_call_data, function_handler would be called with param_obj as argument.
			
			If deletion_cutoff is None, all segments who have matches in ecotype_id2cnv_qc_call_data would be considered in function_handler().
			If deletion_cutoff is not None (some float), only those segments whose amplitude is below this value would be considered in function_handler().
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		from pymodule.CNV import is_reciprocal_overlap, CNVSegmentBinarySearchTreeKey
		sys.stderr.write("Getting probe amplitude from %s ... \n"%repr(input_fname_ls))
		amp_ls = []
		array_id2array = {}
		counter = 0
		real_counter = 0
		no_of_deletions = 0
		no_of_valid_deletions = 0
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		for line in input_handler:
			if line.find("array_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			ecotype_id_idx = col_name2index.get('ecotype_id', col_name2index.get('array_id'))
			cnv_ecotype_id = int(row[ecotype_id_idx])
			array_id = int(row[col_name2index.get('array_id')])
			#row[ecotype_id_idx] = cnv_ecotype_id
			counter += 1
			if cnv_ecotype_id in ecotype_id2cnv_qc_call_data:	# array is in CNVQCDat
				cnv_qc_call_data = ecotype_id2cnv_qc_call_data.get(cnv_ecotype_id)
				start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				stop_probe = row[col_name2index['end_probe']].split('_')
				stop_probe = map(int, stop_probe)
				no_of_probes = int(row[col_name2index['length']])
				if no_of_probes<min_no_of_probes:
					continue
				amplitude = float(row[col_name2index['amplitude']])
				segment_chromosome = start_probe[0]
				segment_start_pos = start_probe[1]-12
				segment_stop_pos = stop_probe[1]+12
				segment_length = abs(segment_stop_pos-segment_start_pos+1)
				if deletion_cutoff is not None and amplitude>deletion_cutoff:
					continue
				cnv_segment_obj = PassingData(ecotype_id=cnv_ecotype_id, start_probe=start_probe, stop_probe=stop_probe,\
											no_of_probes=no_of_probes, amplitude=amplitude, segment_length=segment_length,\
											segment_chromosome=segment_chromosome, array_id=array_id)
				cnvSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=segment_chromosome, span_ls=[segment_start_pos, segment_stop_pos], min_reciprocal_overlap=min_reciprocal_overlap)
				
				no_of_deletions+=1
				
				targetSegmentValue = cnv_qc_call_data.get(cnvSegmentKey)
				if targetSegmentValue:
					no_of_valid_deletions += 1
					cnv_qc_call = targetSegmentValue
					function_handler(param_obj, cnv_segment_obj, cnv_qc_call)
				else:
					function_handler(param_obj, cnv_segment_obj, cnv_qc_call=None)	# this stage, call the handler if it wants to record the number of deletions.
				
				"""
				for cnv_qc_call in cnv_qc_call_data:
					qc_chromosome, qc_start, qc_stop = cnv_qc_call[:3]
					cnv_qc_call_id = cnv_qc_call[-1]
					valid_match = False
					if qc_chromosome==segment_chromosome:
						boundary_diff1 = abs(segment_start_pos-qc_start)
						boundary_diff2 = abs(segment_stop_pos-qc_stop)
						diff1_perc = boundary_diff1/float(segment_length)
						diff2_perc = boundary_diff2/float(segment_length)
						
						is_overlap = is_reciprocal_overlap([segment_start_pos, segment_stop_pos], [qc_start, qc_stop], \
													min_reciprocal_overlap=min_reciprocal_overlap)
						
						if is_overlap:
							no_of_valid_deletions += 1
							valid_match = True
						
						#if boundary_diff1<=max_boundary_diff and boundary_diff2<=max_boundary_diff and diff1_perc<=max_diff_perc and \
						#diff2_perc<=max_diff_perc:
						#	no_of_valid_deletions += 1
						#	valid_match = True
						#elif count_embedded_segment_as_match and segment_start_pos>=qc_start and segment_stop_pos<=qc_stop:	#the segment doesn't match the criteria but very small and within
						#	no_of_valid_deletions += 1
						#	valid_match = True
						
						if valid_match:
							function_handler(param_obj, cnv_segment_obj, cnv_qc_call, )
					elif qc_chromosome>segment_chromosome:
						break
				"""
			if report and counter%10000==0:
				sys.stderr.write('%s%s\t%s\t%s'%('\x08'*80, counter, no_of_deletions, no_of_valid_deletions))
		setattr(param_obj, "no_of_deletions", no_of_deletions)
		setattr(param_obj, "no_of_valid_deletions", no_of_valid_deletions)
		sys.stderr.write("\n")
	
	@classmethod
	def getCNVQCDataFromDB(cls, data_source_id=1, ecotype_id=None, cnv_type_id=None, \
								min_QC_segment_size=None, min_no_of_probes=None, min_reciprocal_overlap=0.6):
		"""
		2010-1-26
			replace the list structure of cnv_qc_call_data in ecotype_id2cnv_qc_call_data with binary_tree structure
		2009-12-9
			add no_of_probes_covered into returning data
			add cnv_type_id
		2009-11-4
			get CNV QC data from database
		"""
		sys.stderr.write("Getting CNV QC data ... \n")
		import Stock_250kDB
		sql_string = "select a.ecotype_id, c.chromosome, c.start, c.stop, c.size_affected, c.no_of_probes_covered, c.copy_number, c.id from %s c,\
						%s a where c.accession_id=a.id and a.data_source_id=%s order by RAND()"%\
						(Stock_250kDB.CNVQCCalls.table.name, Stock_250kDB.CNVQCAccession.table.name, data_source_id)
						# 2010-1-26 random ordering to optimize the tree
		if cnv_type_id is not None:
			sql_string += " and c.cnv_type_id=%s"%cnv_type_id
		if ecotype_id is not None:
			sql_string += " and a.ecotype_id=%s"%ecotype_id
		if min_no_of_probes is not None:
			sql_string += " and c.no_of_probes_covered>=%s"%min_no_of_probes
		if min_QC_segment_size is not None:
			sql_string += " and c.size_affected>=%s"%min_QC_segment_size
		rows = db_250k.metadata.bind.execute(sql_string)
		count = 0
		ecotype_id2cnv_qc_call_data = {}
		#from pymodule.BinarySearchTree import binary_tree	#2010-1-26 replace the list structure of cnv_qc_call_data \
															# in ecotype_id2cnv_qc_call_data with binary_tree structure
		from pymodule.RBTree import RBDict	# 2010-1-26 RBDict is more efficiency than binary_tree.
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey
		for row in rows:
			if row.ecotype_id not in ecotype_id2cnv_qc_call_data:
				#ecotype_id2cnv_qc_call_data[row.ecotype_id] = binary_tree()
				ecotype_id2cnv_qc_call_data[row.ecotype_id] = RBDict()
			
			segmentKey = CNVSegmentBinarySearchTreeKey(chromosome=row.chromosome, span_ls=[row.start, row.stop], \
													min_reciprocal_overlap=min_reciprocal_overlap)
			ecotype_id2cnv_qc_call_data[row.ecotype_id][segmentKey] = (row.chromosome, row.start, row.stop, row.size_affected, row.no_of_probes_covered, row.copy_number, row.id)
			
			#cnv_qc_call_data = ecotype_id2cnv_qc_call_data[row.ecotype_id]
			#cnv_qc_call_data.append((row.chromosome, row.start, row.stop, row.size_affected, row.no_of_probes_covered, row.copy_number, row.id))
			count += 1
		
		
		import math
		for ecotype_id, tree in ecotype_id2cnv_qc_call_data.iteritems():
			print "\tDepth of Ecotype %s's tree: %d" % (ecotype_id, tree.depth())
			print "\tOptimum Depth: %f (%d) (%f%% depth efficiency)" % (tree.optimumdepth(), math.ceil(tree.optimumdepth()),
															  math.ceil(tree.optimumdepth()) / tree.depth())
			#cnv_qc_call_data.sort()
			#ecotype_id2cnv_qc_call_data[ecotype_id] = cnv_qc_call_data
		
		sys.stderr.write("\t%s cnv qc calls for %s ecotypes. Done.\n"%(count, len(ecotype_id2cnv_qc_call_data)))
		return ecotype_id2cnv_qc_call_data
	
	@classmethod
	def countMatchedDeletionsFunctor(cls, param_obj, cnv_segment_obj=None, cnv_qc_call=None):
		"""
		2009-12-9
			store qc data in param_obj.array_id2qc_data
		2009-11-4
			a functor to be called in 
		"""
		from pymodule import PassingData
		if not hasattr(param_obj, 'no_of_valid_deletions'):
			setattr(param_obj, 'no_of_valid_deletions', 0)
		if not hasattr(param_obj, "array_id2qc_data"):
			param_obj.array_id2qc_data = {}
		if not hasattr(param_obj, "array_id2no_of_probes2qc_data"):	# for FPR per no_of_probes
			param_obj.array_id2no_of_probes2qc_data = {}
		if not hasattr(param_obj, "array_id2qc_no_of_probes2qc_data"):	# for FNR per no_of_probes
			param_obj.array_id2qc_no_of_probes2qc_data = {}
		
		array_id = cnv_segment_obj.array_id
		no_of_probes = cnv_segment_obj.no_of_probes
		if array_id not in param_obj.array_id2qc_data:
			param_obj.array_id2qc_data[array_id] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															no_of_valid_deletions=0,\
															no_of_deletions=0,\
															cnv_qc_call_id_set=set())
			param_obj.array_id2no_of_probes2qc_data[array_id] = {}
			param_obj.array_id2qc_no_of_probes2qc_data[array_id] = {}
		if no_of_probes not in param_obj.array_id2no_of_probes2qc_data[array_id]:
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															no_of_valid_deletions=0,\
															no_of_deletions=0,\
															cnv_qc_call_id_set=set())
		# 2010-1-26 increase the no_of_deletions counter no matter whether there's a corresponding cnv_qc_call or not
		param_obj.array_id2qc_data[array_id].no_of_deletions += 1
		param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].no_of_deletions += 1
		if cnv_qc_call is not None:
			qc_chromosome, qc_start, qc_stop = cnv_qc_call[:3]
			cnv_qc_call_id = cnv_qc_call[-1]
			param_obj.array_id2qc_data[array_id].cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.array_id2qc_data[array_id].no_of_valid_deletions += 1
			
			qc_no_of_probes = cnv_qc_call[4]
			if qc_no_of_probes not in param_obj.array_id2qc_no_of_probes2qc_data[array_id]:
				param_obj.array_id2qc_no_of_probes2qc_data[array_id][qc_no_of_probes] = PassingData(ecotype_id=cnv_segment_obj.ecotype_id, \
															cnv_qc_call_id_set=set())
				param_obj.array_id2qc_no_of_probes2qc_data[array_id][qc_no_of_probes].cnv_qc_call_id_set.add(cnv_qc_call_id)
			
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.array_id2no_of_probes2qc_data[array_id][no_of_probes].no_of_valid_deletions += 1
	
	@classmethod
	def outputFalseNegativeRate(cls, param_obj):
		"""
		2009-12-9
			calculate FNR for each class with same number of probes
		2009-11-4
		"""
		for array_id, qc_data in param_obj.array_id2qc_data.iteritems():
			no_of_QCCalls_matched = len(qc_data.cnv_qc_call_id_set)
			no_of_total_QCCalls = len(param_obj.ecotype_id2cnv_qc_call_data[qc_data.ecotype_id])
			false_negative_rate = (no_of_total_QCCalls-no_of_QCCalls_matched)/float(no_of_total_QCCalls)
			sys.stderr.write("Array %s false negative rate: %s/%s(%s).\n"%(array_id, \
																		no_of_total_QCCalls-no_of_QCCalls_matched,\
																		no_of_total_QCCalls, false_negative_rate))
			
			if getattr(param_obj, 'array_id2qc_no_of_probes2qc_data', None):
				qc_no_of_probes2qc_data = param_obj.array_id2qc_no_of_probes2qc_data[array_id]
				no_of_probes_ls = qc_no_of_probes2qc_data.keys()
				no_of_probes_ls.sort()
				for no_of_probes in no_of_probes_ls:
					qc_data = qc_no_of_probes2qc_data[no_of_probes]
					no_of_QCCalls_matched = len(qc_data.cnv_qc_call_id_set)
					no_of_total_QCCalls = len(param_obj.ecotype_id2qc_no_of_probes2cnv_qc_call_id_set[qc_data.ecotype_id][no_of_probes])
					false_negative_rate = (no_of_total_QCCalls-no_of_QCCalls_matched)/float(no_of_total_QCCalls)
					sys.stderr.write("\t%s\t%s\t%s\t%s\n"%(no_of_probes, \
														no_of_total_QCCalls-no_of_QCCalls_matched,\
														no_of_total_QCCalls, false_negative_rate))
	
	@classmethod
	def outputFalsePositiveRate(cls, param_obj):
		"""
		2009-12-9
			calculate FNR for each class with same number of probes
		2009-11-4
		"""
		for array_id, qc_data in param_obj.array_id2qc_data.iteritems():
			no_of_valid_deletions = qc_data.no_of_valid_deletions
			no_of_deletions = qc_data.no_of_deletions
			no_of_non_valid_deletions = no_of_deletions-no_of_valid_deletions
			false_positive_rate = no_of_non_valid_deletions/float(no_of_deletions)
			sys.stderr.write("Array %s false positive rate: %s/%s(%s).\n"%(array_id, \
							no_of_non_valid_deletions, no_of_deletions, false_positive_rate))
			if getattr(param_obj, 'array_id2no_of_probes2qc_data', None):
				no_of_probes2qc_data = param_obj.array_id2no_of_probes2qc_data[array_id]
				no_of_probes_ls = no_of_probes2qc_data.keys()
				no_of_probes_ls.sort()
				for no_of_probes in no_of_probes_ls:
					qc_data = no_of_probes2qc_data[no_of_probes]
					no_of_valid_deletions = qc_data.no_of_valid_deletions
					no_of_deletions = qc_data.no_of_deletions
					no_of_non_valid_deletions = no_of_deletions-no_of_valid_deletions
					false_positive_rate = no_of_non_valid_deletions/float(no_of_deletions)
					sys.stderr.write("\t%s\t%s\t%s\t%s\n"%(no_of_probes, \
									no_of_non_valid_deletions, no_of_deletions, false_positive_rate))
		
	@classmethod
	def countNoOfCNVDeletionsMatchQC(cls, db_250k, input_fname_ls, ecotype_id=6909, data_source_id=3, cnv_type_id=1, \
								min_QC_segment_size=200, deletion_cutoff=-0.33, max_boundary_diff=10000, \
								max_diff_perc=0.10, min_no_of_probes=5,\
								count_embedded_segment_as_match=True, min_reciprocal_overlap=0.6):
		"""
		2010-1-26
			pass min_reciprocal_overlap to cls.getCNVQCDataFromDB()
		2009-12-9
			calculate FNR for each class with same number of probes
		2009-10-29
			for all CNV deletions, check how many are in QC dataset.
		"""
		ecotype_id2cnv_qc_call_data = cls.getCNVQCDataFromDB(data_source_id, ecotype_id, cnv_type_id, min_QC_segment_size, min_no_of_probes,\
															min_reciprocal_overlap=min_reciprocal_overlap)
		
		from pymodule import PassingData
		param_obj = PassingData(no_of_valid_deletions=0, cnv_qc_call_id_set=set(), array_id2qc_data={})
		cls.compareCNVSegmentsAgainstQCHandler(input_fname_ls, ecotype_id2cnv_qc_call_data, cls.countMatchedDeletionsFunctor, param_obj, \
											deletion_cutoff, max_boundary_diff, max_diff_perc, min_no_of_probes, \
											count_embedded_segment_as_match=count_embedded_segment_as_match, \
											min_reciprocal_overlap=min_reciprocal_overlap, report=False)
		sys.stderr.write("For ecotype_id %s, data_source_id %s, min_QC_segment_size %s, deletion_cutoff: %s, min_no_of_probes: %s, min_reciprocal_overlap: %s.\n"%\
						(ecotype_id, \
						data_source_id, min_QC_segment_size, deletion_cutoff, min_no_of_probes, min_reciprocal_overlap))
		param_obj.ecotype_id2cnv_qc_call_data = ecotype_id2cnv_qc_call_data
		
		param_obj.ecotype_id2qc_no_of_probes2cnv_qc_call_id_set = {}
		for ecotype_id, cnv_qc_call_data in ecotype_id2cnv_qc_call_data.iteritems():
			qc_no_of_probes2cnv_qc_call_id_set = {}
			for cnv_qc_call in  cnv_qc_call_data:
				no_of_probes_covered = cnv_qc_call[4]
				cnv_qc_call_id = cnv_qc_call[-1]
				if no_of_probes_covered not in qc_no_of_probes2cnv_qc_call_id_set:
					qc_no_of_probes2cnv_qc_call_id_set[no_of_probes_covered] = set()
				qc_no_of_probes2cnv_qc_call_id_set[no_of_probes_covered].add(cnv_qc_call_id)
			param_obj.ecotype_id2qc_no_of_probes2cnv_qc_call_id_set[ecotype_id] = qc_no_of_probes2cnv_qc_call_id_set
		
		cls.outputFalsePositiveRate(param_obj)
		cls.outputFalseNegativeRate(param_obj)
		
		
	"""
	# 2009-11-4
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, ecotype_id=8215, data_source_id=3, cnv_type_id=1, \
								min_QC_segment_size=200, deletion_cutoff=-0.33, max_boundary_diff=10000, max_diff_perc=0.10)
	
	# 2009-11-5
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	#input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr1.GADA_A0.5T4M5.tsv')]
	ecotype_id_data_source_id_ls = [(6911, 5), (8215, 3), (9169, 3), (6962, 3)]	# Cvi from Bob, Fei-0 from Schneeberger
	min_QC_segment_size = 5
	min_no_of_probes = 5
	count_embedded_segment_as_match = True
	for ecotype_id, data_source_id in ecotype_id_data_source_id_ls:
		for deletion_cutoff in [-0.33, -0.5]:
			for min_reciprocal_overlap in [0.4, 0.6, 0.8]:
				CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, ecotype_id=ecotype_id, data_source_id=data_source_id, \
										cnv_type_id=1,\
										min_QC_segment_size=min_QC_segment_size, deletion_cutoff=deletion_cutoff, \
										min_no_of_probes=min_no_of_probes,\
										count_embedded_segment_as_match=count_embedded_segment_as_match,\
										min_reciprocal_overlap=min_reciprocal_overlap)
	
	
	# 2009-12-09 use min_reciprocal_overlap
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	#input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr1.GADA_A0.5T4M5.tsv')]
	#input_fname_ls = []
	#for i in range(1,6):
	#	input_fname_ls.append(os.path.expanduser('~/mnt2/panfs//250k/CNV/GADA_output/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s_A2.0T8.0M5.tsv'%i))
	
	ecotype_id_data_source_id_ls = [(6932, 7)]	# two different types of Ler-1 deletions from Ler contigs
	min_QC_segment_size = 5
	min_no_of_probes = 5
	count_embedded_segment_as_match = True
	for ecotype_id, data_source_id in ecotype_id_data_source_id_ls:
		for deletion_cutoff in [-0.33, -0.5]:
			for min_reciprocal_overlap in [0.4, 0.8]:
				CNV.countNoOfCNVDeletionsMatchQC(db_250k, input_fname_ls, ecotype_id=ecotype_id, data_source_id=data_source_id, \
										cnv_type_id=1,\
										min_QC_segment_size=min_QC_segment_size, deletion_cutoff=deletion_cutoff, \
										min_no_of_probes=min_no_of_probes,\
										count_embedded_segment_as_match=count_embedded_segment_as_match,\
										min_reciprocal_overlap=min_reciprocal_overlap)
	"""
	
	@classmethod
	def addAmplitudeFunctor(cls, param_obj, cnv_segment_obj, cnv_qc_call=None):
		"""
		2009-12-9
			adjust argument order and process only if cnv_qc_call is not None
		2009-11-4
		"""
		if not hasattr(param_obj, 'amp_ls'):
			setattr(param_obj, 'amp_ls', [])
		if cnv_qc_call is not None:
			qc_chromosome, qc_start, qc_stop = cnv_qc_call[:3]
			cnv_qc_call_id = cnv_qc_call[-1]
			param_obj.cnv_qc_call_id_set.add(cnv_qc_call_id)
			param_obj.amp_ls.append(cnv_segment_obj.amplitude)
	
	@classmethod
	def drawHistOfAmpOfValidatedDeletions(cls, db_250k, input_fname_ls, output_fname_prefix, data_source_id=1, cnv_type_id=1, \
								min_QC_segment_size=200, min_no_of_probes=5, max_boundary_diff=10000, \
								max_diff_perc=0.10, count_embedded_segment_as_match=True, min_reciprocal_overlap=0.6):
		"""
		2009-11-4
			draw histogram of amplitude of segments who are validated according to certain QC data
		"""
		ecotype_id2cnv_qc_call_data = cls.getCNVQCDataFromDB(data_source_id, cnv_type_id=cnv_type_id, \
															min_QC_segment_size=min_QC_segment_size,\
															min_no_of_probes=min_no_of_probes, \
															min_reciprocal_overlap=min_reciprocal_overlap)
		
		from pymodule import PassingData
		param_obj = PassingData(amp_ls=[], cnv_qc_call_id_set=set())
		cls.compareCNVSegmentsAgainstQCHandler(input_fname_ls, ecotype_id2cnv_qc_call_data, cls.addAmplitudeFunctor, param_obj, \
											deletion_cutoff=None, max_boundary_diff=max_boundary_diff, max_diff_perc=max_diff_perc,\
											min_no_of_probes=min_no_of_probes, \
											count_embedded_segment_as_match=count_embedded_segment_as_match)

		sys.stderr.write("data_source_id %s, min_QC_segment_size %s, min_no_of_probes: %s, max_boundary_diff: %s, max_diff_perc: %s, count_embedded_segment_as_match: %s.\n"%\
						(data_source_id, min_QC_segment_size, min_no_of_probes, max_boundary_diff, max_diff_perc,\
						count_embedded_segment_as_match))
		param_obj.ecotype_id2cnv_qc_call_data = ecotype_id2cnv_qc_call_data
		cls.outputFalsePositiveRate(param_obj)
		cls.outputFalseNegativeRate(param_obj)
				
		sys.stderr.write("Drawing ...")
		import pylab
		pylab.clf()
		pylab.title("Histogram of amplitude of known deletions from source %s"%data_source_id)
		pylab.hist(param_obj.amp_ls, 20)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	# 2009-11-4
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.GADA_A0.5T4M5.tsv'%i))
	#input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr1.GADA_A0.5T4M5.tsv')]
	min_QC_segment_size = 5
	min_no_of_probes = 5
	count_embedded_segment_as_match = True
	data_source_id = 1
	for max_boundary_diff in [10000]:
		for max_diff_perc in [0.20, 0.3]:
			output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/Clark2007aDeletionsAmplitudeHistCall48_mxdiff_%s_mxperc_%s'%(max_boundary_diff, max_diff_perc))	
			CNV.drawHistOfAmpOfValidatedDeletions(db_250k, input_fname_ls, output_fname_prefix, data_source_id=data_source_id, \
									cnv_type_id=1, \
								min_QC_segment_size=min_QC_segment_size, \
								max_boundary_diff=max_boundary_diff, max_diff_perc=max_diff_perc, \
								min_no_of_probes=min_no_of_probes,\
								count_embedded_segment_as_match=count_embedded_segment_as_match)
	
	"""
	
	@classmethod
	def drawDeletedCNVProbeHistogram(cls, db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20):
		"""
		2009-10-12
			for the probes that are deleted according to CNV_qc_fname, look at the histogram of their amplitudes outputted from GADA (input_fname_ls)
		"""
		sys.stderr.write("Drawing deleted CNV probe intensity histogram ... \n")
		import fileinput
		
		from pymodule import SNPData
		CNVQCData = SNPData(input_fname=CNV_qc_fname, turn_into_array=1, ignore_2nd_column=1, data_starting_col=2, \
						turn_into_integer=False)
		new_col_id2col_index = {}
		new_col_id_ls = []
		for col_id in CNVQCData.col_id_ls:
			new_col_id = int(col_id.split('_')[0])	# it's ecotypeid, as in "6899_Bay-0".
			new_col_id_ls.append(new_col_id)
			new_col_id2col_index[new_col_id] = CNVQCData.col_id2col_index[col_id]
		CNVQCData.col_id_ls = new_col_id_ls
		CNVQCData.col_id2col_index = new_col_id2col_index
		
		from Stock_250kDB import Probes, ArrayInfo
		new_row_id2row_index = {}
		new_row_id_ls = []
		for row_id in CNVQCData.row_id_ls:
			probe_id = int(row_id)
			probe = Probes.get(probe_id)
			new_row_id = (probe.chromosome, probe.position)
			new_row_id_ls.append(new_row_id)
			new_row_id2row_index[new_row_id] = CNVQCData.row_id2row_index[row_id]
		CNVQCData.row_id_ls = new_row_id_ls
		CNVQCData.row_id2row_index = new_row_id2row_index
		
		sys.stderr.write("Getting probe amplitude from %s ... \n"%repr(input_fname_ls))
		amp_ls = []
		array_id2array = {}
		counter = 0
		real_counter = 0
		for line in fileinput.input(input_fname_ls):
			if line.find("array_id")==0:
				continue
			line = line.strip()
			row = line.split('\t')
			array_id = int(row[0])
			if array_id not in array_id2array:
				array = ArrayInfo.get(array_id)
				array_id2array[array_id] = array
			array = array_id2array[array_id]
			counter += 1
			if array.maternal_ecotype_id in CNVQCData.col_id2col_index:	# array is in CNVQCData
				start_probe = row[1].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				stop_probe = row[2].split('_')
				stop_probe = map(int, stop_probe)
				amplitude = float(row[4])
				col_index = CNVQCData.col_id2col_index[array.maternal_ecotype_id]
				for i in range(len(CNVQCData.row_id_ls)):
					row_id = CNVQCData.row_id_ls[i]
					chr, pos = row_id
					no_of_bp_deleted = CNVQCData.data_matrix[i][col_index]
					if no_of_bp_deleted!='?':	# not NA
						no_of_bp_deleted = int(no_of_bp_deleted)
						if no_of_bp_deleted>=min_no_of_bp_deleted and chr==start_probe[0] and pos>=start_probe[1] and pos<=stop_probe[1]:
							amp_ls.append(amplitude)
							real_counter += 1
			if counter%10000==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*80, counter, real_counter))
		sys.stderr.write("Done.\n")
		
		import pylab
		pylab.clf()
		pylab.title("Histogram of intensity of probes with >=20 bases deleted")
		pylab.hist(amp_ls, 20)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
	
	""""
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
	
	CNV_qc_fname = '/Network/Data/250K/tmp-yh/CNV/dazhe_CNV_probes_w_del.csv'
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/ProbeFromPerlegenCNVQC_AmplitudeHistCall43')
	CNV.drawDeletedCNVProbeHistogram(db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20)
	
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
		CNV_qc_fname = '/Network/Data/250K/tmp-yh/CNV/dazhe_CNV_probes_w_del.csv'
		output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/ProbeFromPerlegenCNVQC_AmplitudeCall43Chr%sHist'%i)
		CNV.drawDeletedCNVProbeHistogram(db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20)
	
	input_fname_ls = []
	for i in [2,4]:
		input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_48_CNV_QNormalize_chr%s.GADAJRN_A0.5T4M10.tsv'%i)]
		CNV_qc_fname = '/Network/Data/250K/tmp-yh/CNV/dazhe_CNV_probes_w_del.csv'
		output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/ProbeFromPerlegenCNVQC_AmplitudeCall48Chr%sHist'%i)
		CNV.drawDeletedCNVProbeHistogram(db_250k, input_fname_ls, CNV_qc_fname, output_fname_prefix, min_no_of_bp_deleted=20)

	
	"""
	
	
	@classmethod
	def drawDeletedCNVDeletionSizeHist(cls, db_250k, input_fname_ls, output_fname_prefix, max_amplitude=-0.33):
		"""
		2009-10-12
			
		"""
		sys.stderr.write("Drawing deleted CNV probe intensity histogram ... \n")
		import fileinput
		
		from Stock_250kDB import Probes, ArrayInfo
		from pymodule import getColName2IndexFromHeader
		sys.stderr.write("Getting probe amplitude from %s ... \n"%repr(input_fname_ls))
		array_id2array = {}
		length_ls = []
		
		counter = 0
		real_counter = 0
		array_id2length_ls = {}
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		for line in input_handler:
			if line.find("array_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			
			ecotype_id_idx = col_name2index.get('ecotype_id', col_name2index.get('array_id'))
			#cnv_ecotype_id = int(row[ecotype_id_idx])
			array_id = int(row[col_name2index.get('array_id')])
			if array_id not in array_id2array:
				array = ArrayInfo.get(array_id)
				array_id2array[array_id] = array
				array_id2length_ls[array_id] = []
			array = array_id2array[array_id]
			counter += 1			
			amplitude = float(row[col_name2index['amplitude']])
			if amplitude<=max_amplitude:
				start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				stop_probe = row[col_name2index['end_probe']].split('_')
				stop_probe = map(int, stop_probe)
				segment_chromosome = start_probe[0]
				segment_start_pos = start_probe[1]-12
				segment_stop_pos = stop_probe[1]+12
				segment_length = abs(segment_stop_pos-segment_start_pos+1)
				
				#length = stop_probe[1]-start_probe[1]+1
				
				length_ls.append(segment_length)
				array_id2length_ls[array_id].append(segment_length)
				real_counter += 1
			if counter%10000==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*80, counter, real_counter))
		sys.stderr.write("Done.\n")
		
		import pylab
		for array_id, length_ls in array_id2length_ls.iteritems():
			if len(length_ls)>20:
				array = array_id2array[array_id]
				output_fname_prefix_array = '%s_array_%s_ecotype_%s'%(output_fname_prefix, array_id, array.maternal_ecotype_id)
				pylab.clf()
				pylab.title("Histogram of length of deletions with max_amplitude=%s"%max_amplitude)
				pylab.hist(length_ls, 40)
				pylab.savefig('%s.png'%output_fname_prefix_array, dpi=300)
				pylab.savefig('%s.svg'%output_fname_prefix_array, dpi=300)
	
	
	"""
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
	max_amplitude = -0.33
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/call_43_deletion_length_histogram_max_amp_%s'%max_amplitude)
	CNV.drawDeletedCNVDeletionSizeHist(db_250k, input_fname_ls, output_fname_prefix, max_amplitude=max_amplitude)
	
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_43_CNV_norm_intensity_chr%s.GADAJRN_A0.5T4M5.tsv'%i))
	
	max_amplitude = -0.33
	output_fname_prefix = os.path.expanduser('~/script/variation/data/CNV/figures/call_43_deletion_length_histogram_max_amp_%s'%max_amplitude)
	CNV.drawDeletedCNVDeletionSizeHist(db_250k, input_fname_ls, output_fname_prefix, max_amplitude=max_amplitude)
	
	"""
	
	@classmethod
	def fillDeletionInCNVQCProbeCalls(cls, db_250k, data_source_id=1, commit=False):
		"""
		2009-10-29
			construct CNVQCProbeCalls based on CNVQCCalls, only deletions
		"""
		sys.stderr.write("Filling in CNVQCProbeCalls based on CNVQCCalls ...\n")
		session = db_250k.session
		#session.begin()
		
		import Stock_250kDB
		from sqlalchemy import desc, asc
		qc_query = Stock_250kDB.CNVQCCalls.query.filter_by(cnv_type_id=1)	# deletion only
		qc_call_data = []
		for qc_call in qc_query:
			qc_call_data.append((qc_call.chromosome, qc_call.start, qc_call.stop, qc_call.id, qc_call.accession_id))
		qc_call_data.sort()
		sys.stderr.write("%s CNVQCCalls entries.\n"%(len(qc_call_data)))
		
		
		query = db_250k.metadata.bind.execute("select id, chromosome, position, xpos, ypos, allele, strand from %s where direction is not null order by chromosome, position"%(Stock_250kDB.Probes.table.name))
		# 2009-10-29 sqlalchemy takes up too much memory
		# query = Stock_250kDB.Probes.query.filter(Stock_250kDB.Probes.direction!=None).order_by(asc(Stock_250kDB.Probes.chromosome)).\
		#	order_by(asc(Stock_250kDB.Probes.position))
		starting_index_for_next_probe = 0
		no_of_probes = 0
		no_of_encounters = 0
		for probe in query:
			probe_start_pos = probe.position-12
			probe_stop_pos = probe.position+12
			no_of_segments_containing_this_probe = 0
			for i in range(starting_index_for_next_probe, len(qc_call_data)):
				chromosome, start, stop, qc_call_id, accession_id= qc_call_data[i]
				if chromosome==probe.chromosome and start<=probe_stop_pos and stop>=probe_start_pos:	# probe within
					no_of_segments_containing_this_probe += 1
					if no_of_segments_containing_this_probe==1:	# first segment encountering this probe, next probe should at least start from here.
						starting_index_for_next_probe = i
					target_position = probe_start_pos-start
					if start>=probe_start_pos:
						size_affected = 25-(start-probe_start_pos)
					elif stop<=probe_stop_pos:
						size_affected = 25-(probe_stop_pos-stop)
					else:
						size_affected = 25
					db_entry = Stock_250kDB.CNVQCProbeCalls(accession_id=accession_id, probe_id=probe.id, chromosome=probe.chromosome, \
														position=probe.position, \
														size_affected=size_affected, target_position=target_position, \
														cnv_qc_call_id = qc_call_id, cnv_type_id=1)
					session.save(db_entry)
					session.flush()
					no_of_encounters += 1
				elif chromosome<probe.chromosome:	# skip
					continue
				elif chromosome>probe.chromosome or (chromosome==probe.chromosome and start>probe_stop_pos):	# all following segments is beyond this probe. exit
					break
			if no_of_segments_containing_this_probe==0:
				starting_index_for_next_probe = i-1		# this probe is not in any segment. Next probe should start from the (i-1)th segment,
											# which is the last segment ahead of the current probe.
			no_of_probes += 1
			if no_of_probes%5000==0:
				sys.stderr.write("%sNo of probes: %s\tNo of encounters: %s"%('\x08'*80, no_of_probes, no_of_encounters))
		#if commit:
		#	session.commit()
		sys.stderr.write("No of probes: %s\tNo of encounters: %s. Done.\n"%(no_of_probes, no_of_encounters))
	
	"""
	CNV.fillDeletionInCNVQCProbeCalls(db_250k, commit=True)
	"""
	
	@classmethod
	def processColName2IndexForArrayInfo(cls, col_name2index):
		"""
		2009-11-20
		# construct ecotype_id2array_id_ls mapper
		"""
		import sys, os
		sys.stderr.write("Constructing ecotype_id2array_id_ls ...")
		from pymodule import PassingData
		import Stock_250kDB, re
		all_number_pattern = re.compile(r'^\d+$')
		array_id2index = {}
		ecotype_id2array_id_ls = {}
		array_label_ls = []
		for col_name in col_name2index:
			if all_number_pattern.search(col_name):
				array_id = int(col_name)
				array = Stock_250kDB.ArrayInfo.get(array_id)
				array_id2index[array_id] = col_name2index[col_name]
				ecotype_id = array.maternal_ecotype_id
				array_label_ls.append('e-%s %s'%(ecotype_id, array_id))
				if ecotype_id not in ecotype_id2array_id_ls:
					ecotype_id2array_id_ls[ecotype_id] = []
				ecotype_id2array_id_ls[ecotype_id].append(array_id)
		array_index_ls = array_id2index.values()
		array_index_ls.sort()
		min_array_index = min(array_index_ls)
		max_array_index = max(array_index_ls)
		returnData = PassingData(array_id2index=array_id2index, ecotype_id2array_id_ls=ecotype_id2array_id_ls,\
								min_array_index=min_array_index, max_array_index=max_array_index, array_label_ls=array_label_ls)
		sys.stderr.write("Done.\n")
		return returnData
	
	@classmethod
	def calculateProbeSTDDEV(cls, db_250k, input_fname_ls, output_fname, min_no_of_replicates=5):
		"""
		2009-11-19
			calculate stddev of one probe in multiple replicate arrays
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		array_id2array = {}
		counter = 0
		real_counter = 0
		
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		returnData = cls.processColName2IndexForArrayInfo(col_name2index)
		ecotype_id2array_id_ls = returnData.ecotype_id2array_id_ls
		min_array_index = returnData.min_array_index
		max_array_index = returnData.max_array_index
		array_id2index = returnData.array_id2index
		
		# only to calculate stddev of probes of ecotypes who have >=min_no_of_replicates arrays.
		ecotype_id_to_check = []
		for ecotype_id, array_id_ls in ecotype_id2array_id_ls.iteritems():
			if len(array_id_ls)>=min_no_of_replicates:
				ecotype_id_to_check.append(ecotype_id)
		
		sys.stderr.write("Calculating stddev of probe intensity from replicate arrays from %s ... \n"%repr(input_fname_ls))
		import csv
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['probe_id', 'chr_pos', 'ecotype_id', 'sample_size', 'stddev']
		writer.writerow(header)
		import numpy
		from annot.bin.codense.common import dict_map
		for line in input_handler:
			if line.find("probes_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			probe_id_idx = col_name2index.get('probes_id')
			probe_id = int(row[probe_id_idx])
			chr = row[col_name2index.get('chromosome')]
			pos = row[col_name2index.get('position')]
			chr_pos = '%s_%s'%(chr, pos)
			total_probe_intensity_ls = row[min_array_index:max_array_index+1]
			total_probe_intensity_ls = map(float, total_probe_intensity_ls)
			total_probe_intensity_ls = numpy.array(total_probe_intensity_ls)
			for ecotype_id in ecotype_id_to_check:
				array_id_ls = ecotype_id2array_id_ls[ecotype_id]
				array_index_ls = dict_map(array_id2index, array_id_ls)
				probe_intensity_ls = total_probe_intensity_ls[array_index_ls]
				stddev = numpy.std(probe_intensity_ls)
				output_row = [probe_id, chr_pos, ecotype_id, len(probe_intensity_ls), stddev]
				writer.writerow(output_row)
			
			stddev = numpy.std(total_probe_intensity_ls)
			output_row = [probe_id, chr_pos, 'all', len(total_probe_intensity_ls), stddev]
			writer.writerow(output_row)
			
			counter += 1			
			if counter%10000==0:
				sys.stderr.write('%s%s'%('\x08'*80, counter))
		del writer
		sys.stderr.write("Done.\n")
	"""
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_48_CNV_probe_stddev.tsv'
	min_no_of_replicates = 5
	CNV.calculateProbeSTDDEV(db_250k, input_fname_ls, output_fname, min_no_of_replicates=min_no_of_replicates)
	
	
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname = '/Network/Data/250k/tmp-yh/CNV/call_method_48_CNV_QNorm_sub_ref_probe_stddev.tsv'
	min_no_of_replicates = 5
	CNV.calculateProbeSTDDEV(db_250k, input_fname_ls, output_fname, min_no_of_replicates=min_no_of_replicates)
	"""
	
	@classmethod
	def calculateProbeQuartilePerChromosome(cls, db_250k, output_fname_prefix=None):
		"""
		2009-11-25
			a newer and expanded version of this function is moved to DB_250k2Array.py
		2009-11-20
			calculate CNV probe intensity quartile (1st, median, 3rd) for each chromosome and store them in database
		"""
		import sys, os, numpy, math
		from scipy import stats	# for scoreatpercentile/percentileatscore to get quartiles
		import Stock_250kDB
		import rpy
		rpy.r.library('affy')
		from DB_250k2Array import DB_250k2Array
		session = db_250k.session
		
		probes, xy_ls, chr_pos_ls, probes_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, Stock_250kDB.Probes.table.name, \
																		snps=None, run_type=2)
		
		
		sys.stderr.write("Getting probes into each chromosome ...")
		chr2xy_ls = {}
		chr2probe_id_ls = {}
		for i in range(len(xy_ls)):
			chr,pos = chr_pos_ls[i]
			if chr not in chr2xy_ls:
				chr2xy_ls[chr] = []
				chr2probe_id_ls[chr] = []	#initialize with the start_probe_id
			chr2xy_ls[chr].append(xy_ls[i])
			chr2probe_id_ls[chr].append(probes_id_ls[i])
		sys.stderr.write("Done.\n")
		
		
		rows = Stock_250kDB.ArrayInfo.query()
		no_of_arrays = rows.count()
		array_width = None
		counter = 0
		for row in rows:
			array_id = row.id
			filename = row.filename
			ecotype_id = row.maternal_ecotype_id
						
			sys.stderr.write("\t%d/%d: Extracting intensity from %s ... "%(counter+1, no_of_arrays, array_id))
			
			#read array by calling R
			array = rpy.r.read_affybatch(filenames=filename)
			intensity_array = rpy.r.intensity(array)	#return a lengthX1 2-Dimensional array.
			if array_width == None:
				intensity_array_size = len(intensity_array)
				array_width = int(math.sqrt(intensity_array_size))	#assume it's square array
			
			chr2intensity_ls = {}
			for chr, chr_xy_ls in chr2xy_ls.iteritems():
				chr2intensity_ls[chr] = []
				for xpos, ypos in chr_xy_ls:
					intensity_array_index = array_width*(array_width - xpos - 1) + ypos
					intensity = math.log10(intensity_array[intensity_array_index][0])
					chr2intensity_ls[chr].append(intensity)
				array_quartile = Stock_250kDB.ArrayQuartile(array_id=array_id, start_probe_id=chr2probe_id_ls[chr][0],\
														stop_probe_id=chr2probe_id_ls[chr][-1], \
														no_of_probes=len(chr2intensity_ls[chr]))
				array_quartile.minimum = numpy.min(chr2intensity_ls[chr])
				array_quartile.first_decile = stats.scoreatpercentile(chr2intensity_ls[chr], 10)
				array_quartile.lower_quartile = stats.scoreatpercentile(chr2intensity_ls[chr], 25)
				array_quartile.median = stats.scoreatpercentile(chr2intensity_ls[chr], 50)
				array_quartile.upper_quartile = stats.scoreatpercentile(chr2intensity_ls[chr], 75)
				array_quartile.last_decile = stats.scoreatpercentile(chr2intensity_ls[chr], 90)
				array_quartile.maximum = numpy.max(chr2intensity_ls[chr])
				
				session.save_or_update(array_quartile)
				
				# find and store the outliers
				IQR = array_quartile.upper_quartile - array_quartile.lower_quartile
				lower_whisker = array_quartile.lower_quartile - 1.5*IQR
				upper_whisker = array_quartile.upper_quartile + 1.5*IQR
				for i in range(len(chr2intensity_ls[chr])):
					intensity = chr2intensity_ls[chr][i]
					if intensity<lower_whisker or intensity>upper_whisker:
						probe_id = chr2probe_id_ls[chr][i]
						array_quartile_outlier = Stock_250kDB.ArrayQuartileOutlier(probe_id=probe_id, intensity=intensity)
						array_quartile_outlier.array_quartile = array_quartile
						session.save_or_update(array_quartile_outlier)
				session.flush()
			counter += 1
			sys.stderr.write("Done.\n")
		
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		array_id2array = {}
		counter = 0
		real_counter = 0
		
		input_handler = fileinput.input(input_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		returnData = cls.processColName2IndexForArrayInfo(col_name2index)
		min_array_index = returnData.min_array_index
		max_array_index = returnData.max_array_index
		
		output_dir = os.path.split(output_fname_prefix)[0]
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		
		sys.stderr.write("Drawing box plots of CNV arrays for each chromosome from %s ... \n"%repr(input_fname_ls))
		probe_intensity_matrix = []
		for i in range(len(returnData.array_label_ls)):
			probe_intensity_matrix.append([])
		import numpy, pylab
		from annot.bin.codense.common import dict_map
		prev_chromosome = None
		
		for line in input_handler:
			if line.find("probes_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			probe_id_idx = col_name2index.get('probes_id')
			probe_id = int(row[probe_id_idx])
			chr = row[col_name2index.get('chromosome')]
			if prev_chromosome==None:
				prev_chromosome = chr
			#pos = row[col_name2index.get('position')]
			#chr_pos = '%s_%s'%(chr, pos)
			if (chr!=prev_chromosome and len(probe_intensity_matrix[0])>20 and counter>0) or (len(probe_intensity_matrix[0])%90000==0 and counter>0):
				sys.stderr.write("boxplot chr %s at %s ..."%(prev_chromosome, counter))
				for i in range(len(returnData.array_label_ls)):
					pylab.clf()
					pylab.boxplot(probe_intensity_matrix[i])
					pylab.title(returnData.array_label_ls[i])
					pylab.ylim([1,5])
					pylab.savefig('%s_%s_chr%s_%s.png'%(output_fname_prefix, returnData.array_label_ls[i], prev_chromosome, counter), dpi=300)
				sys.stderr.write("Done.\n")
				
				#clean up the matrix
				#del probe_intensity_matrix
				#probe_intensity_matrix = []
				for i in range(len(returnData.array_label_ls)):
					#probe_intensity_matrix.append([])
					probe_intensity_matrix[i] = []
			
			if chr!=prev_chromosome:
				prev_chromosome = chr
				for i in range(len(returnData.array_label_ls)):
					probe_intensity_matrix[i] = []
				
			total_probe_intensity_ls = row[min_array_index:max_array_index+1]
			total_probe_intensity_ls = map(float, total_probe_intensity_ls)
			for i in range(len(total_probe_intensity_ls)):
				probe_intensity_matrix[i].append(total_probe_intensity_ls[i])
			
			counter += 1			
			if counter%10000==0:
				sys.stderr.write('%s%s'%('\x08'*80, counter))
		sys.stderr.write("boxplot chr %s at %s ..."%(chr, counter))
		for i in range(len(returnData.array_label_ls)):
			pylab.clf()
			pylab.boxplot(probe_intensity_matrix[i])
			pylab.title(returnData.array_label_ls[i])
			pylab.ylim([1,5])
			pylab.savefig('%s_%s_chr%s_%s.png'%(output_fname_prefix, returnData.array_label_ls[i], prev_chromosome, counter), dpi=300)
		sys.stderr.write("Done.\n")
		"""
	
	"""
	output_fname_prefix = '/Network/Data/250k/tmp-yh/CNV/CNV_intensity_boxplot'
	CNV.calculateProbeQuartilePerChromosome(db_250k, output_fname_prefix)
	"""
	
	class OneCNV(object):
		def __init__(self, max_boundary_diff=20000, max_diff_perc=0.2, min_reciprocal_overlap=0.6):
			self.chromosome = None
			self.start = None
			self.stop = None
			self.segment_length = None
			self.cnv_ls = []
			self.array_id_ls = []
			self.max_boundary_diff = max_boundary_diff
			self.max_diff_perc = max_diff_perc
			self.min_reciprocal_overlap = min_reciprocal_overlap
		
		def addFirstCNV(self, chromosome, start, stop, array_id=None):
			self.addOneCNV(chromosome, start, stop, array_id)
		
		def addOneCNV(self, chromosome, start, stop, array_id=None):
			"""
			"""
			if self.chromosome is not None and chromosome!=self.chromosome:
				return False
			
			if self.chromosome is None:
				self.chromosome = chromosome
			
			self.cnv_ls.append((start, stop))
			self.array_id_ls.append(array_id)
			self.adjustBoundary()
			return True
		
		def adjustBoundary(self):
			"""
			"""
			import numpy
			if len(self.cnv_ls)>0:
				self.start, self.stop = numpy.mean(self.cnv_ls, 0)
				self.segment_length = abs(self.stop-self.start)
		
		def __len__(self):
			return len(self.cnv_ls)
		
		def addNewCNV(self, chromosome, start, stop, array_id=None):
			"""
			"""
			if self.chromosome is None:
				self.addOneCNV(chromosome, start, stop, array_id)
			elif self.chromosome is not None and chromosome!=self.chromosome:
				return False
			else:
				"""
				boundary_diff1 = abs(start-self.start)
				boundary_diff2 = abs(stop-self.stop)
				diff1_perc = boundary_diff1/float(self.segment_length)
				diff2_perc = boundary_diff2/float(self.segment_length)
				if boundary_diff1<=self.max_boundary_diff and boundary_diff2<=self.max_boundary_diff and \
					diff1_perc<=self.max_diff_perc and diff2_perc<=self.max_diff_perc:
					self.addOneCNV(chromosome, start, stop, array_id)
				
				else:
					return False
				"""
				
				is_overlap = is_reciprocal_overlap([start, stop], [self.start, self.stop], \
													min_reciprocal_overlap=self.min_reciprocal_overlap)
				if is_overlap:
					self.addOneCNV(chromosome, start, stop, array_id)
				else:
					return False
			
	
	@classmethod
	def plotCNVOccurrenceInReplicatesHist(cls, db_250k, cnv_intensity_fname, GADA_output_fname_ls, output_fname_prefix, \
										min_no_of_replicates=5, min_no_of_probes=5,\
										deletion_cutoff=None, max_boundary_diff=10000, \
										max_diff_perc=0.10, min_reciprocal_overlap=0.6):
		"""
		2009-11-19
			calculate stddev of one probe in multiple replicate arrays
		"""
		import fileinput
		from pymodule import getColName2IndexFromHeader, PassingData
		array_id2array = {}
		counter = 0
		real_counter = 0
		
		input_handler = fileinput.input([cnv_intensity_fname])
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		
		returnData = cls.processColName2IndexForArrayInfo(col_name2index)
		ecotype_id2array_id_ls = returnData.ecotype_id2array_id_ls
		min_array_index = returnData.min_array_index
		max_array_index = returnData.max_array_index
		array_id2index = returnData.array_id2index
		
		# only to calculate stddev of probes of ecotypes who have >=min_no_of_replicates arrays.
		ecotype_id_to_check = []
		for ecotype_id, array_id_ls in ecotype_id2array_id_ls.iteritems():
			if len(array_id_ls)>=min_no_of_replicates:
				ecotype_id_to_check.append(ecotype_id)
		ecotype_id_to_check_set = set(ecotype_id_to_check)
		fileinput.close()	# must. otherwise next fileinput won't work.
		
		sys.stderr.write("Getting occurrence from replicate arrays from %s ... \n"%repr(GADA_output_fname_ls))
		counter = 0
		real_counter = 0
		no_of_deletions = 0
		input_handler = fileinput.input(GADA_output_fname_ls)
		header = input_handler.readline().strip().split('\t')
		col_name2index = getColName2IndexFromHeader(header)
		ecotype_id2cnv_ls = {}
		for line in input_handler:
			if line.find("array_id")!=-1:
				continue
			line = line.strip()
			row = line.split('\t')
			ecotype_id_idx = col_name2index.get('ecotype_id', col_name2index.get('array_id'))
			cnv_ecotype_id = int(row[ecotype_id_idx])
			array_id = int(row[col_name2index.get('array_id')])
			#row[ecotype_id_idx] = cnv_ecotype_id
			counter += 1
			if cnv_ecotype_id in ecotype_id_to_check_set:	# array is in CNVQCDat
				
				start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
				start_probe = map(int, start_probe)
				stop_probe = row[col_name2index['end_probe']].split('_')
				stop_probe = map(int, stop_probe)
				no_of_probes = int(row[col_name2index['length']])
				if no_of_probes<min_no_of_probes:
					continue
				amplitude = float(row[col_name2index['amplitude']])
				segment_chromosome = start_probe[0]
				segment_start_pos = start_probe[1]-12
				segment_stop_pos = stop_probe[1]+12
				segment_length = abs(segment_stop_pos-segment_start_pos+1)
				if deletion_cutoff is not None and amplitude>deletion_cutoff:
					continue
				no_of_deletions+=1
				if cnv_ecotype_id not in ecotype_id2cnv_ls:
					ecotype_id2cnv_ls[cnv_ecotype_id] = []
				cnv_ls = ecotype_id2cnv_ls[cnv_ecotype_id]
				addToExistingCNV = False
				for one_cnv in cnv_ls:
					if one_cnv.addNewCNV(segment_chromosome, segment_start_pos, segment_stop_pos, array_id):
						addToExistingCNV = True
				if not addToExistingCNV:
					one_cnv = cls.OneCNV(max_boundary_diff =max_boundary_diff, max_diff_perc=max_diff_perc, \
										min_reciprocal_overlap=min_reciprocal_overlap)
					one_cnv.addFirstCNV(segment_chromosome, segment_start_pos, segment_stop_pos, array_id)
					cnv_ls.append(one_cnv)
				
			if counter%10000==0:
				sys.stderr.write('%s%s\t%s'%('\x08'*80, counter, no_of_deletions))
		
		for ecotype_id, cnv_ls in ecotype_id2cnv_ls.iteritems():
			occurrence_ls = map(len, cnv_ls)
			import pylab
			pylab.clf()
			no_of_replicates = len(ecotype_id2array_id_ls[ecotype_id])
			pylab.title('%s(%s reps), del-c %s, min #probes %s, max dist %s, perc %s, overlap %s'%(ecotype_id, no_of_replicates, deletion_cutoff, min_no_of_probes, \
																max_boundary_diff, max_diff_perc, min_reciprocal_overlap))
			pylab.hist(occurrence_ls, no_of_replicates)
			
			pylab.xlabel('occurrence of one deletion')
			pylab.ylabel('count')
			pylab.savefig('%s_e_%s.png'%(output_fname_prefix, ecotype_id))
		
		
		sys.stderr.write("Done.\n")
	
	"""
	
	cnv_intensity_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')
	max_boundary_diff = 20000
	max_diff_perc = 0.2
	
	GADA_output_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/GADA_output/call_method_48_CNV_intensity_QNorm_sub_ref_chr1_A0.8T8.0M10.tsv')]
	deletion_cutoff = -0.33
	min_reciprocal_overlap = 0.6
	min_no_of_probes = 5
	output_fname_prefix = '/tmp/CNVOccurrence_call_48_chr1_A0.8T8.0M10_delCutoff%s_min_p%s_mdist%s_mperc%s_moverlap%s'%\
						(deletion_cutoff, min_no_of_probes, max_boundary_diff, max_diff_perc, min_reciprocal_overlap)
	CNV.plotCNVOccurrenceInReplicatesHist(db_250k, cnv_intensity_fname, GADA_output_fname_ls, output_fname_prefix, min_no_of_replicates=5, \
										min_no_of_probes=min_no_of_probes, \
										deletion_cutoff=deletion_cutoff, max_boundary_diff=max_boundary_diff, \
										max_diff_perc=max_diff_perc, min_reciprocal_overlap=min_reciprocal_overlap)
	
	
	cnv_intensity_fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')
	aAlpha_ls = [0.2, 0.5, 0.8, 1.0, 2.0]
	T_ls = [2.0,4.0,8.0,12.0,14.0]
	M_ls = [5,10]
	deletion_cutoff_ls = [-0.33, -0.4, -0.5, -0.6]
	
	max_boundary_diff = 20000
	max_diff_perc = 0.2
	
	for aAlpha in aAlpha_ls:
		for T in T_ls:
			for M in M_ls:
				for deletion_cutoff in deletion_cutoff_ls:
					for min_no_of_probes in [5,10,20,40]:
						for min_reciprocal_overlap in [0.4, 0.6, 0.8]:
							GADA_output_fname_ls = []
							for chr in range(1,6):
								fname = os.path.expanduser('~/mnt2/panfs/250k/CNV/GADA_output/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s_A%sT%sM%s.tsv'%(chr, aAlpha, T, M))
								if os.path.isfile(fname):
									GADA_output_fname_ls.append(fname)
							output_fname_prefix = '/Network/Data/250k/tmp-yh/CNV/CNVOccurrenceByOverlap/call_48_A%sT%sM%s_delCutoff%s_min_p%s_mdist%s_mperc%s_moverlap%s'%\
								(aAlpha, T, M, deletion_cutoff, min_no_of_probes, max_boundary_diff, max_diff_perc, min_reciprocal_overlap)
							CNV.plotCNVOccurrenceInReplicatesHist(db_250k, cnv_intensity_fname, GADA_output_fname_ls, \
																output_fname_prefix, \
																min_no_of_replicates=5, \
															min_no_of_probes=min_no_of_probes, deletion_cutoff=deletion_cutoff, \
															max_boundary_diff=max_boundary_diff, \
															max_diff_perc=max_diff_perc, min_reciprocal_overlap=min_reciprocal_overlap)
	
	"""
	
	@classmethod
	def drawHistOfDLRSpread(cls, db_250k, output_fname_prefix, call_method_id=None):
		"""
		2009-11-25
		"""
		import os, sys
		sys.stderr.write("Drawing histogram of DLRSpread for call %s ..."%call_method_id)
		import Stock_250kDB
		if call_method_id is not None:
			rows = Stock_250kDB.CallInfo.query.filter_by(method_id=call_method_id)
		else:
			rows = Stock_250kDB.ArrayQuartile.query.filter(Stock_250kDB.ArrayQuartile.dlrspread!=None)
		dlrspread_ls = []
		for row in rows:
			if call_method_id is not None:
				for array_quartile in row.array.array_quartile_ls: 
					if array_quartile.dlrspread is not None:
						dlrspread_ls.append(array_quartile.dlrspread)
			else:
				 dlrspread_ls.append(row.dlrspread)
		
		import pylab
		pylab.clf()
		if call_method_id is not None:
			title = 'arrays from call %s'%call_method_id
		else:
			title = 'all arrays'
		pylab.title(title)
		pylab.hist(dlrspread_ls, 30)
		pylab.xlabel('DLRSpread')
		pylab.ylabel('Count')
		pylab.xlim([0.1,0.5])
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
	
	"""
	output_fname_prefix = '/tmp/DLRSpreadHistAllArrays'
	CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix)
	
	call_method_id = 48
	output_fname_prefix = '/tmp/DLRSpreadHistCall%sArrays'%call_method_id
	CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix, call_method_id=call_method_id)
	
	output_fname_prefix = '/tmp/DLRSpreadHistAllArrays'
	CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix)
	
	for call_method_id in (43, 32, 29, 48, 49):
		output_fname_prefix = '/tmp/DLRSpreadHistCall%sArrays'%call_method_id
		CNV.drawHistOfDLRSpread(db_250k, output_fname_prefix, call_method_id=call_method_id)
	
	"""
	
	@classmethod
	def drawLerContigSizeHist(cls, fasta_input_fname, output_fname_prefix):
		"""
		2009-12-7
			This function draws a histogram of size of Ler contigs.			
				Ler contigs were downloaded from http://www.arabidopsis.org/browse/Cereon/index.jsp, in fasta format.
				a png file with prefix as output_fname_prefix will be generated.
		"""
		import os, sys
		sys.stderr.write("Drawing histogram of Ler contig size ... ")
		inf = open(fasta_input_fname)
		from Bio import SeqIO
		contig_size_ls = []
		for seq_record in SeqIO.parse(inf, "fasta") :
			contig_size_ls.append(len(seq_record.seq))
			#of.write('>%s\n'%seq_record.id)
		
		import pylab
		pylab.clf()
		title = 'Histogram of Ler Contig Size'
		pylab.title(title)
		pylab.hist(contig_size_ls, 30)
		pylab.xlabel('contig size')
		pylab.ylabel('Count')
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
		sys.stderr.write("Done.\n")
		
	"""
	fasta_input_fname = os.path.expanduser('~/script/variation/data/CNV/Cereon_Ath_Ler.fasta')
	output_fname_prefix = '/tmp/Ler-contig-size-hist'
	CNV.drawLerContigSizeHist(fasta_input_fname, output_fname_prefix)
	"""
	
	@classmethod
	def putLerContigsIntoDB(cls, db, data_source, accession_name, fasta_input_fname):
		"""
		2010-1-27
			put ler contigs sequences and ids into table Stock_250kDB.SequenceFragment
		"""
		import os, sys
		sys.stderr.write("Putting Ler contigs into db ... ")
		session = db.session
		import Stock_250kDB
		from CNVQCConstruct import CNVQCConstruct
		data_source_obj = CNVQCConstruct.getDBObj(session, Stock_250kDB.DataSource, data_source)
		acc_obj = CNVQCConstruct.getCNVQCAccessionObj(session, accession_name, data_source_obj)
		inf = open(fasta_input_fname)
		from Bio import SeqIO
		for seq_record in SeqIO.parse(inf, "fasta"):
			sequence_fragment = Stock_250kDB.SequenceFragment(short_name=seq_record.id, size=len(seq_record.seq),\
															description=seq_record.description,\
															sequence=seq_record.seq.tostring())
			sequence_fragment.accession = acc_obj
			session.save(sequence_fragment)
		session.flush()
		sys.stderr.write("Done.\n")
	
	"""
	fasta_input_fname = os.path.expanduser('~/script/variation/data/CNV/Cereon_Ath_Ler.fasta')
	data_source = 'LerContig'
	accession_name = 'Ler-1'
	CNV.putLerContigsIntoDB(db_250k, data_source, accession_name, fasta_input_fname)
	"""
	
	@classmethod
	def discoverLerDeletionDuplication(cls, db_250k, ler_blast_result_fname, output_fname, deletion_only=True, min_no_of_matches=25):
		"""
		2009-12-7
			ler_blast_result_fname is the output of blasting all CNV probes against Ler contigs http://www.arabidopsis.org/browse/Cereon/index.jsp.
			Two functions:
				1. deletion_only=True. make sure the deletion is covered by the sequencing.
					one naive criteria is if the boundary (two adjacent non-deleted probes) is within the same contig, then yes.
				2. deletion_only=False, detect copy number changes. If two adjacent probes have different number of contigs, 
					then it's a copy number change point.
		"""
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		sys.stderr.write("Reading from %s ... \n"%ler_blast_result_fname)
		counter = 0
		real_counter = 0
		import csv
		reader = csv.reader(open(ler_blast_result_fname), delimiter=figureOutDelimiter(ler_blast_result_fname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		probe_id2contig_id_ls = {}
		for row in reader:
			contig_label = row[col_name2index['Alignment_title']]
			probe_id = int(row[col_name2index['Probe_ID']])
			no_of_matches = int(row[col_name2index['Number_matches']])
			if no_of_matches>=min_no_of_matches:
				contig_id = ' '.join(contig_label.split()[1:])
				if probe_id not in probe_id2contig_id_ls:
					probe_id2contig_id_ls[probe_id] = []
				probe_id2contig_id_ls[probe_id].append(contig_id)
		sys.stderr.write("Done.\n")
		del reader
		
		import Stock_250kDB
		from DB_250k2Array import DB_250k2Array
		probes, xy_ls, chr_pos_ls, total_probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, Stock_250kDB.Probes.table.name, \
																		snps=None, run_type=2)
		
		chr2xy_ls, chr2probe_id_ls = DB_250k2Array.organizeProbesIntoChromosome(xy_ls, chr_pos_ls, total_probe_id_ls)
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header_row = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', 'length', 'copy_number']
		writer.writerow(header_row)
		
		sys.stderr.write("Discovering deletions ...\n")
		counter = 0
		real_counter = 0
		for chr, probe_id_ls in chr2probe_id_ls.iteritems():
			no_of_probes = len(probe_id_ls)
			index_of_prev_probe_within_a_contig = None
			index_of_prev_probe_with_a_different_copy_number = None
			for i in range(no_of_probes):
				probe_id = probe_id_ls[i]
				contig_id_ls = probe_id2contig_id_ls.get(probe_id,[])
				copy_number = len(contig_id_ls)
				if i==0:
					index_of_prev_probe_with_a_different_copy_number = -1	# set before the first probe (index=0)
				else:
					prev_probe_contig_id_ls = probe_id2contig_id_ls.get(probe_id_ls[i-1],[])
					prev_copy_number = len(prev_probe_contig_id_ls)
					if not deletion_only:
						if copy_number != prev_copy_number:	# a change point of copy number
							if index_of_prev_probe_with_a_different_copy_number is not None:
								start_probe_id = probe_id_ls[index_of_prev_probe_with_a_different_copy_number+1]
								start_probe = probes.get_one_probe(start_probe_id)
								start_chr_pos = '%s_%s'%(start_probe.chr, start_probe.pos)
								stop_probe_id = probe_id_ls[i-1]
								stop_probe = probes.get_one_probe(stop_probe_id)
								stop_chr_pos = '%s_%s'%(stop_probe.chr, stop_probe.pos)
								row = [start_probe_id, start_chr_pos, stop_probe_id, stop_chr_pos, \
										i-index_of_prev_probe_with_a_different_copy_number-1, stop_probe.pos-start_probe.pos, prev_copy_number]
								writer.writerow(row)
								real_counter += 1
							index_of_prev_probe_with_a_different_copy_number = i-1
					else:	# look for deleton only. The only difference from above is make sure the deletion is covered by the sequencing.
						#one naive criteria is if the boundary is within the same contig, then yes.
						if prev_copy_number>0 and copy_number==0:	# from non-deletion to deletion
							index_of_prev_probe_within_a_contig = i-1
						elif prev_copy_number==0 and copy_number>0:	# from deletion to non-deletion
							if index_of_prev_probe_within_a_contig is not None:	# found a potential deletion
								current_contig_id_set = set(contig_id_ls)
								prev_non_deleted_probe_id = probe_id_ls[index_of_prev_probe_within_a_contig]
								prev_non_deleted_probe_contig_id_ls = probe_id2contig_id_ls.get(prev_non_deleted_probe_id, [])
								prev_non_deleted_probe_contig_id_set = set(prev_non_deleted_probe_contig_id_ls)
								if len(prev_non_deleted_probe_contig_id_set&current_contig_id_set)>0:	#share at least one contig. deletion confirmed
									deletion_start_probe_id = probe_id_ls[index_of_prev_probe_within_a_contig+1]
									deletion_start_probe = probes.get_one_probe(deletion_start_probe_id)
									deletion_start_chr_pos = '%s_%s'%(deletion_start_probe.chr, deletion_start_probe.pos)
									deletion_stop_probe_id = probe_id_ls[i-1]
									deletion_stop_probe = probes.get_one_probe(deletion_stop_probe_id)
									deletion_stop_chr_pos = '%s_%s'%(deletion_stop_probe.chr, deletion_stop_probe.pos)
									row = [deletion_start_probe_id, deletion_start_chr_pos, deletion_stop_probe_id, deletion_stop_chr_pos, \
										i-index_of_prev_probe_within_a_contig-1, deletion_stop_probe.pos-deletion_start_probe.pos, prev_copy_number]
									writer.writerow(row)
									real_counter += 1
							index_of_prev_probe_within_a_contig = i
						elif prev_copy_number>0 and copy_number>0:	# from non-deletion to non-deletion
							index_of_prev_probe_within_a_contig = i
				
				counter += 1
				if counter%10000==0:
					sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
			
			# don't forget the last segment if it's not in the deletion_only mode.
			if not deletion_only and index_of_prev_probe_with_a_different_copy_number is not None:
				start_probe_id = probe_id_ls[index_of_prev_probe_with_a_different_copy_number+1]
				start_probe = probes.get_one_probe(start_probe_id)
				start_chr_pos = '%s_%s'%(start_probe.chr, start_probe.pos)
				stop_probe_id = probe_id_ls[i]	# watch: not i-1.
				stop_probe = probes.get_one_probe(stop_probe_id)
				stop_chr_pos = '%s_%s'%(stop_probe.chr, stop_probe.pos)
				row = [start_probe_id, start_chr_pos, stop_probe_id, stop_chr_pos, \
						i-index_of_prev_probe_with_a_different_copy_number, stop_probe.pos-start_probe.pos, copy_number]	# watch no -1, and it's copy_number
				writer.writerow(row)
				real_counter += 1
		sys.stderr.write("Done.\n")
		del writer
	"""
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/ler_raw_CNV_QC.csv'
	output_fname = '/tmp/Ler-deletions.tsv'
	CNV.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname)
	
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/ler_raw_CNV_QC.csv'
	output_fname = '/tmp/Ler-copy-number.tsv'
	CNV.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname, deletion_only=False)
	
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/tair9_raw.csv'
	output_fname = '/tmp/Col-copy-number.tsv'
	CNV.discoverLerDeletionDuplication(db_250k, ler_blast_result_fname, output_fname, deletion_only=False)
	
	"""
	
	@classmethod
	def discoverLerContigSpanOverCol(cls, ler_blast_result_fname, output_fname, min_no_of_matches=25,
									max_delta_ratio=0.4, max_length_delta=50000):
		"""
		2010-1-27
			discover the span of Ler Contigs in terms of Col reference genome
			
			The algorithm:
				for each contig, get all the probes that perfectly match some part of it.
					sort the probes in chromosomal order
					set the start-probe to the 1st probe
						for stop-probe in [last probe, ..., 2nd probe]:
							if [start-probe, stop-probe] and [alignment-start, alignment-stop] meet the condition:
								col-span of this ler contig = the chromosomal positions of [start-probe, stop-probe]
			
			The condition is either the length delta is <= max_length_delta, or either of the two delta ratios
				(delta/length1, delta/length2) <= max_delta_ratio.
			
		"""
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		sys.stderr.write("Reading from %s ... \n"%ler_blast_result_fname)
		counter = 0
		real_counter = 0
		import csv
		reader = csv.reader(open(ler_blast_result_fname), delimiter=figureOutDelimiter(ler_blast_result_fname))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		contig_id2probe_chr_pos_ls = {}
		for row in reader:
			contig_label = row[col_name2index['Alignment_title']]
			probe_id = int(row[col_name2index['Probe_ID']])
			no_of_matches = int(row[col_name2index['Number_matches']])
			if no_of_matches>=min_no_of_matches:
				contig_id = contig_label.split()[1]	# 2010-1-28, take "ATL8C9990" in the case of "ATL8C9990 ATL7C121_1"
				if contig_id not in contig_id2probe_chr_pos_ls:
					contig_id2probe_chr_pos_ls[contig_id] = []
				chr = int(row[col_name2index['Chromosome']])
				pos = int(row[col_name2index['Position']])
				alignment_start = int(row[col_name2index['Alignment_start']])
				probe_chr_pos = (chr, pos)
				probe_id = int(probe_id)
				contig_id2probe_chr_pos_ls[contig_id].append((probe_chr_pos, probe_id, alignment_start))
		del reader
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Discovering Ler Contig Span over Col ...\n")		
		import csv, Stock_250kDB
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header_row = ['start_probe_id', 'start_chr_pos', 'stop_probe_id', 'stop_chr_pos', 'no_of_probes', 'length', 'copy_number', \
					'contig_id', 'contig_start', 'contig_stop', 'size_difference']
		writer.writerow(header_row)
		counter = 0
		for contig_id, probe_chr_pos_ls in contig_id2probe_chr_pos_ls.iteritems():
			sys.stderr.write('%sContig %s, %s'%('\x08'*80, counter, contig_id))
			probe_chr_pos_ls.sort()	# sorted according to probe position on Col genome
			n = len(probe_chr_pos_ls)
			span_start_probe_index = 0
			while span_start_probe_index<n:
				old_span_start_probe_index = span_start_probe_index
				stop_index_candidate_ls = range(span_start_probe_index+1, n)
				stop_index_candidate_ls.reverse()	# starting from the end to cover as much as it can
				for i in stop_index_candidate_ls:
					start_probe_chr_pos, start_probe_id, start_alignment_pos = probe_chr_pos_ls[span_start_probe_index]
					stop_probe_chr_pos, stop_probe_id, stop_alignment_pos = probe_chr_pos_ls[i]
					start_chr, start_pos = start_probe_chr_pos
					stop_chr, stop_pos = stop_probe_chr_pos
					if start_chr == stop_chr:
						col_length = abs(stop_pos-start_pos)
						ler_length = abs(stop_alignment_pos-start_alignment_pos)
						length_delta = col_length-ler_length
						if col_length>0:
							col_ratio = abs(length_delta)/float(col_length)
						else:
							continue	# 2010-1-28 ignore if col_length = 0
						if ler_length>0:
							ler_ratio = abs(length_delta)/float(ler_length)
						else:
							continue	# 2010-1-28 ignore if ler_length = 0
						if abs(length_delta)<=max_length_delta or col_ratio<=max_delta_ratio or ler_ratio<=max_delta_ratio:
							no_of_probes = i-span_start_probe_index+1
							row = [start_probe_id, '%s_%s'%(start_chr, start_pos-12), stop_probe_id, '%s_%s'%(stop_chr, stop_pos+12), no_of_probes, stop_pos-start_pos+25, \
								'', contig_id, start_alignment_pos, stop_alignment_pos, length_delta]
							writer.writerow(row)
							span_start_probe_index = i+1	# set the new starting index
							break
				if old_span_start_probe_index == span_start_probe_index:	# no change in the for loop above. or nothing found to meet the min_reciprocal_overlap
					span_start_probe_index += 1	# auto increment it by 1
			counter += 1
		del writer
		sys.stderr.write("Done.\n")
	
	"""
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/ler_raw_CNV_QC.csv'
	max_delta_ratio = 0.4
	max_length_delta = 10000
	for max_length_delta in [1,2,3,4,5,6,7,8,9,10,12,15,20]:
		max_length_delta = max_length_delta*10000
		output_fname = '/tmp/Ler-span-over-Col-mdr%s-mld%s.tsv'%(max_delta_ratio, max_length_delta)
		CNV.discoverLerContigSpanOverCol(ler_blast_result_fname, output_fname, min_no_of_matches=25,
										max_delta_ratio=max_delta_ratio, max_length_delta=max_length_delta)
	"""
	
	@classmethod
	def drawIntensityVsProbeTrueCopyNumber(cls, db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=3, cnv_type_id=None, \
										min_reciprocal_overlap=0.6, report=True, max_copy_number=16):
		"""
		2010-1-26
			value of the ecotype_id2cnv_qc_call_data dictionary is a RBDict (RBTree dictionary) structure.
			remove the frame from the final plot, add a grid instead.
			add argument max_copy_number to limit the plot only to probes whose copy numbers are <= this threshold.
		2009-12-9
			purpose is to check whether probe copy number based on Col/Ler blast results actually mean something.
			input_fname_ls is a list of filenames which contains CNV intensity, raw or normalized.
			
			argument min_reciprocal_overlap is not used.
			cnv_type_id = None means all cnv types.
		"""
		ecotype_id2cnv_qc_call_data = cls.getCNVQCDataFromDB(data_source_id, ecotype_id, cnv_type_id, \
															min_reciprocal_overlap=min_reciprocal_overlap)
		
		import Stock_250kDB
		from DB_250k2Array import DB_250k2Array
		probes, xy_ls, chr_pos_ls, total_probe_id_ls = DB_250k2Array.get_probes(db_250k.metadata.bind, Stock_250kDB.Probes.table.name, \
																		snps=None, run_type=2)		
		
		sys.stderr.write("Establish true copy number of each probe ...")	#2009 need to improve the running time by using CNVSegmentBinarySearchTreeKey from pymodule/CNV.py
		cnv_qc_call_data = ecotype_id2cnv_qc_call_data[ecotype_id]
		probe_id2copy_number = {}
		from pymodule.CNV import CNVSegmentBinarySearchTreeKey
		for i in range(len(total_probe_id_ls)):
			probe_id = total_probe_id_ls[i]
			chr, pos= chr_pos_ls[i]
			probeSegmentKey = CNVSegmentBinarySearchTreeKey(chromosome=chr, span_ls=[pos-12, pos+12], min_reciprocal_overlap=min_reciprocal_overlap)
			cnv_qc_call = cnv_qc_call_data.get(probeSegmentKey)
			if cnv_qc_call:
				copy_number = cnv_qc_call[5]
				probe_id2copy_number[probe_id] = copy_number
			"""
			for cnv_qc_call in cnv_qc_call_data:	# not efficient, shall use the binary search tree
				qc_chromosome, qc_start, qc_stop = cnv_qc_call[:3]
				if chr == qc_chromosome and pos>=qc_start and pos<=qc_stop:
					copy_number = cnv_qc_call[5]
					probe_id2copy_number[probe_id] = copy_number
					break
			"""
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Getting intensity data from %s ..."%repr(input_fname_ls))
		import fileinput, csv
		from pymodule import getColName2IndexFromHeader, PassingData, figureOutDelimiter
		input_handler = fileinput.input(input_fname_ls)
		reader = csv.reader(input_handler, delimiter=figureOutDelimiter(input_fname_ls[0]))
		header = reader.next()
		col_name2index = getColName2IndexFromHeader(header)
		
		no_of_cols = len(header)-3
		
		import Stock_250kDB
		array_id_ls = header[1:-2]
		col_index_to_be_checked = set()
		for array_id in array_id_ls:	# array_id is of type 'str'
			array = Stock_250kDB.ArrayInfo.get(int(array_id))
			if array.maternal_ecotype_id==ecotype_id and array.paternal_ecotype_id==ecotype_id:
				col_index_to_be_checked.add(col_name2index[array_id])
		
		array_id2copy_number2intensity_ls = {}
		probe_id_index = col_name2index['probes_id']
		chr_index  = col_name2index['chromosome']
		pos_index =  col_name2index['position']
		counter = 0
		real_counter = 0
		for row in reader:
			if row[0].find("probes_id")!=-1:	# encounter header in one of the files of input_fname_ls
				continue
			probe_id = int(row[probe_id_index])
			#chr = int(row[chr_index])
			#pos = int(row[pos_index])
			copy_number = probe_id2copy_number.get(probe_id)
			if copy_number is not None:
				for j in col_index_to_be_checked:
					array_id = header[j]
					if array_id not in array_id2copy_number2intensity_ls:
						array_id2copy_number2intensity_ls[array_id] = {}
					if copy_number not in array_id2copy_number2intensity_ls[array_id]:
						array_id2copy_number2intensity_ls[array_id][copy_number] = []
					intensity = float(row[j])
					array_id2copy_number2intensity_ls[array_id][copy_number].append(intensity)
					real_counter += 1
			counter += 1
			if counter%10000==0 and report:
				sys.stderr.write("%s%s\t%s"%('\x08'*80, counter, real_counter))
		
		del reader, input_handler
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Drawing intensity vs probe copy number ...")
		import pylab
		for array_id, copy_number2intensity_ls in array_id2copy_number2intensity_ls.iteritems():
			pylab.clf()
			ax1 = pylab.axes([0.05, 0.05, 0.9, 0.9], frameon=False)	# 2010-1-26 remove the frame
			ax1.grid(True, alpha=0.2)
			ax1.title.set_text("Array %s of Ecotype %s"%(array_id, ecotype_id))
			copy_number_ls = copy_number2intensity_ls.keys()
			copy_number_ls.sort()
			for copy_number in copy_number_ls:
				if copy_number<=max_copy_number:
					intensity_ls = copy_number2intensity_ls[copy_number]
					if len(intensity_ls)>10:	# more than 10 values, draw a boxplot.
						ax1.boxplot(intensity_ls, positions=[copy_number])
					else:
						x_value_ls = [copy_number]*len(intensity_ls)
						ax1.plot(x_value_ls, intensity_ls, '.')
			output_fname = '%s_array_%s.png'%(output_fname_prefix, array_id)
			ax1.set_xlim([-1, max_copy_number+1])
			pylab.savefig(output_fname, dpi=200)
			
		sys.stderr.write("Done.\n")
	
	"""
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	
	input_fname_ls = [os.path.expanduser('~/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = '/tmp/call_48_Col_intensity_vs_true_copy_number'
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None)
	
	max_copy_number = 16
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_intensity_QNorm_sub_ref_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number)
	
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_intensity_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number)
	
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number)
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_QNorm_sub_ref_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number)
	
	"""
	
class AnalyzeSNPData(object):
	@classmethod
	def DrawStrainSNP_NA_PercHist(cls, data_matrix_fname, output_fname, need_savefig=0):
		"""
		2007-03-20
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		strain_NA_perc_ls = []
		for i in range(data_matrix.shape[0]):
			strain_NA_perc_ls.append(sum(data_matrix[i,:]==0)/float(data_matrix.shape[1]))	#0 is NA
		SNP_NA_perc_ls = []
		for i in range(data_matrix.shape[1]):
			SNP_NA_perc_ls.append(sum(data_matrix[:,i]==0)/float(data_matrix.shape[0]))	#0 is NA
		import pylab,os
		base_fname = os.path.basename(data_matrix_fname)
		pylab.clf()
		pylab.hist(strain_NA_perc_ls, 20)
		pylab.title("%s Strain NA perc histogram"%base_fname)
		if need_savefig:
			pylab.savefig('%s_strain_NA_perc.eps'%output_fname, dpi=300)
			pylab.savefig('%s_strain_NA_perc.svg'%output_fname, dpi=300)
			pylab.savefig('%s_strain_NA_perc.png'%output_fname, dpi=300)
		pylab.show()
		
		pylab.clf()
		pylab.hist(SNP_NA_perc_ls, 20)
		pylab.title("%s SNP NA perc histogram"%base_fname)
		if need_savefig:
			pylab.savefig('%s_SNP_NA_perc.eps'%output_fname, dpi=300)
			pylab.savefig('%s_SNP_NA_perc.svg'%output_fname, dpi=300)
			pylab.savefig('%s_SNP_NA_perc.png'%output_fname, dpi=300)
		pylab.show()
	
	"""
	DrawStrainSNP_NA_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)
	"""
	@classmethod
	def DrawStrain_Heterozygotes_PercHist(cls, data_matrix_fname, output_fname, need_savefig=0):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		strain_Hetero_perc_ls = []
		for i in range(data_matrix.shape[0]):
			strain_Hetero_perc_ls.append(sum(data_matrix[i,:]>4)/float(data_matrix.shape[1]))	#bigger than 4 is heterozygotes
		
		import pylab,os
		base_fname = os.path.basename(data_matrix_fname)
		pylab.clf()
		pylab.hist(strain_Hetero_perc_ls, 20)
		pylab.title("%s Strain Heterozygote perc histogram"%base_fname)
		if need_savefig:
			pylab.savefig('%s_strain_hz_perc.eps'%output_fname, dpi=300)
			pylab.savefig('%s_strain_hz_perc.svg'%output_fname, dpi=300)
			pylab.savefig('%s_strain_hz_perc.png'%output_fname, dpi=300)
		pylab.show()
		return strain_Hetero_perc_ls
	
	"""
	strain_Hetero_perc_ls = DrawStrain_Heterozygotes_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)
	"""
	
	"""
	2007-09-24 increase the #bins of histogram to 40
	2007-03-21 draw histogram of pairwise accession genetic distance
	"""
	@classmethod
	def DrawDistanceHistogram(cls, data_matrix_fname, output_fname, need_savefig=0):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
		import Numeric
		distance_ls = []
		no_of_NA_pairs = 0
		no_of_strains = len(data_matrix)
		no_of_snps = len(data_matrix[0])
		distance_matrix = Numeric.zeros([no_of_strains, no_of_strains], Numeric.Float)
		for i in range(no_of_strains):
			for j in range(i+1, no_of_strains):
				no_of_valid_pairs = 0.0
				no_of_matching_pairs = 0.0
				for k in range(no_of_snps):
					if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
						no_of_valid_pairs += 1
						if data_matrix[i][k] == data_matrix[j][k]:
							no_of_matching_pairs += 1
				if no_of_valid_pairs!=0:
					distance = 1 - no_of_matching_pairs/no_of_valid_pairs
					distance_matrix[i,j] = distance_matrix[j,i] = distance
					distance_ls.append(distance)
				else:
					no_of_NA_pairs += 1
		print "out of %s pairs, %s are NA"%((no_of_strains*(no_of_strains-1))/2, no_of_NA_pairs)
		import pylab
		pylab.clf()
		pylab.hist(distance_ls, 40)
		pylab.title("Histogram of non-NA distances")
		if need_savefig:
			pylab.savefig('%s_distance_hist.eps'%output_fname, dpi=300)
			pylab.savefig('%s_distance_hist.svg'%output_fname, dpi=300)
			pylab.savefig('%s_distance_hist.png'%output_fname, dpi=300)
		pylab.show()
		return distance_matrix
	
	"""
	distance_matrix = DrawDistanceHistogram('./script/variation/data/justin_data_filtered.csv', './script/variation/data/justin_data_filtered' , 1)
	"""
	
	"""
	2007-09-17
		check allele frequency
	"""
	@classmethod
	def cal_maf_vector(cls, input_fname):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		from EstimateSelfingGeneration import EstimateSelfingGeneration
		EstimateSelfingGeneration_instance = EstimateSelfingGeneration()
		locus_allele_prob_vector = EstimateSelfingGeneration_instance.cal_locus_allele_prob_vector(data_matrix)
		maf_vector = Numeric.zeros(locus_allele_prob_vector.shape[0], Numeric.Float)	#the minor allele frequency vector
		for i in range(locus_allele_prob_vector.shape[0]):
			maf_vector[i] = min(locus_allele_prob_vector[i])
		import pylab
		pylab.hist(maf_vector, 10)
		return maf_vector
	"""
	input_fname = 'script/variation/stock20070829/data_row_na_col_na_bad_snps.tsv'
	input_fname = '/mnt/hpc-cmb/KW/input/250K_method_5_after_imputation_noRedundant_051908.tsv'
	maf_vector = cal_maf_vector(input_fname)
	import pylab
	pylab.clf()
	pylab.plot(range(len(maf_vector)), maf_vector)
	"""
	
	"""
	2007-09-23
		calculate the percentage of NAs in a data_matrix
	"""
	@classmethod
	def calculate_NA_perc(cls, input_fname):
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		import Numeric
		data_matrix = Numeric.array(data_matrix)
		no_of_NAs = 0
		for i in range(data_matrix.shape[0]):
			for j in range(data_matrix.shape[1]):
				if data_matrix[i][j] == 0:
					no_of_NAs += 1
		no_of_totals = data_matrix.shape[0]*data_matrix.shape[1]
		print '%s/%s=%s'%(no_of_NAs, no_of_totals, float(no_of_NAs)/no_of_totals)
	
	"""
	calculate_NA_perc('.//script/variation/stock20070919/data.tsv')
	"""
	
	
	#2008-2-16
	#check the frequency of all allele-combinations between two SNPs
	@classmethod
	def get_allele_combo2freq(cls, snp_data, snp_id1, snp_id2):
		snp_index1 = snp_data.col_id2col_index[snp_id1]
		snp_index2 = snp_data.col_id2col_index[snp_id2]
		no_of_rows, no_of_cols = snp_data.data_matrix.shape
		allele_combo2freq = {}
		for i in range(no_of_rows):
			allele_combo = (snp_data.data_matrix[i][snp_index1], snp_data.data_matrix[i][snp_index2])
			if allele_combo not in allele_combo2freq:
				allele_combo2freq[allele_combo] = 0
			allele_combo2freq[allele_combo] += 1
		return allele_combo2freq
	"""
	from pymodule import SNPData
	snp_data = SNPData(input_fname='panfs/250k/call_method_29_binary.tsv', turn_into_array=1)
	get_allele_combo2freq(snp_data, '1_18234094', '5_18607474')
	"""
	
	@classmethod
	def getAlignmentMatrixFromFasta(cls, fasta_input_fname, output_fname, chromosome=1, start=1, pickPolymorphicColumns=True):
		"""
		2009-5-28
			add argument pickPolymorphicColumns
		2009-3-26
			generate alignment matrix out of a sequence file in fasta format
				(if ID includes an underscore, the part behind underscore is regarded as real accession id.)
			chromosome & start are the position of the 1st nucleotide in the alignment
		"""
		import os, sys, numpy
		sys.stderr.write("Getting alignment matrix from %s ..."%(fasta_input_fname))
		snp_pos_ls = []
		accession_id_ls = []
		name_ls = []
		data_matrix = []
		inf = open(fasta_input_fname)
		from Bio import SeqIO
		from DiscoverSNPFromAlignment import DiscoverSNPFromAlignment
		from pymodule import dict_map, nt2number, PassingData
		counter = 0
		for seq_record in SeqIO.parse(inf, "fasta"):
			if counter == 0:
				snp_pos_ls = DiscoverSNPFromAlignment.get_snp_pos_ls(seq_record.seq.tostring().upper(), chromosome, start)
			record_id_split = seq_record.id.split('_')
			if len(record_id_split)==2:
				record_id = record_id_split[1]
			else:
				record_id = record_id_split[0]
			accession_id_ls.append(record_id)
			name_ls.append(seq_record.id)
			data_row = dict_map(nt2number, seq_record.seq.tostring().upper())
			data_matrix.append(data_row)
			counter += 1
		data_matrix = numpy.array(data_matrix, numpy.int8)
		passingdata = PassingData(snp_pos_ls=snp_pos_ls, accession_id_ls=accession_id_ls, name_ls=name_ls, data_matrix=data_matrix)
		sys.stderr.write(' %s accessions, %s bases. Done.\n'%(len(accession_id_ls), len(snp_pos_ls)))
		
		if pickPolymorphicColumns:
			DiscoverSNPFromAlignment.pickPolymorphicColumns(passingdata)
		
		header = ['id', 'name']
		for snp_pos in passingdata.snp_pos_ls:
			header.append('%s_%s_%s'%snp_pos)
		passingdata.header = header
		return passingdata
	
	@classmethod
	def outputSNPsOutOfFastaAlignMatrixInYuFormat(cls, fasta_input_fname, output_fname, chromosome=1, start=1):
		"""
		2009-5-28
			call getAlignmentMatrixFromFasta(), then write the returned data into file
		"""
		passingdata = cls.getAlignmentMatrixFromFasta(fasta_input_fname, output_fname, chromosome, start)
		from pymodule import write_data_matrix
		write_data_matrix(passingdata.data_matrix, output_fname, passingdata.header, \
						passingdata.accession_id_ls, passingdata.name_ls)
		return passingdata
	
	
	"""
fasta_input_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC_rc_with_ref.align.fasta')
output_fname = os.path.expanduser('~/script/variation/doc/FLC/CarolineDeanFLC_rc_with_ref_snp_matrix.tsv')
snpData = AnalyzeSNPData.outputSNPsOutOfFastaAlignMatrixInYuFormat(fasta_input_fname, output_fname, chromosome=5, start=3170860)
	"""
	
	@classmethod
	def filterColumnsOfFasta(cls, fasta_input_fname, output_fname, chromosome=1, start=1, pickPolymorphicColumns=False):
		"""
		2009-5-28
			read the fasta_input_fname (in fasta format, if ID includes an underscore, the part behind underscore is regarded as real accession id.)
			purpose of this program is to remove columns full of 'N' or '-', which could cause trouble for alignment with other sequences 
			output in the same fasta format for further data fiddling
		"""
		passingdata = cls.getAlignmentMatrixFromFasta(fasta_input_fname, output_fname, chromosome, start, pickPolymorphicColumns=False)
		sys.stderr.write("Filtering all-NA or all-deletion columns ... ")
		from pymodule import nt2number, number2nt, dict_map, number2single_char_nt
		import numpy
		no_of_accessions = len(passingdata.accession_id_ls)
		all_del_col = [nt2number['-']]*no_of_accessions
		all_NA_col = [nt2number['N']]*no_of_accessions
		non_del_NA_col_indices = []
		all_NA_counter = 0
		all_del_counter = 0
		for j in range(passingdata.data_matrix.shape[1]):
			all_del_truth_vector = passingdata.data_matrix[:,j]==all_del_col
			all_NA_truth_vector = passingdata.data_matrix[:,j]==all_NA_col
			
			if all_del_truth_vector.all():
				all_del_counter += 1
			elif all_NA_truth_vector.all():
				all_NA_counter += 1
			else:
				non_del_NA_col_indices.append(j)
		
		
		outf = open(output_fname, 'w')
		for i in range(no_of_accessions):
			sequence_id = passingdata.name_ls[i]
			outf.write('>%s\n'%sequence_id)
			nt_ls = dict_map(number2single_char_nt, passingdata.data_matrix[i, non_del_NA_col_indices])
			outf.write('%s\n'%''.join(nt_ls))
		del outf
		sys.stderr.write("%s columns are all NA. %s columns are all del. Done.\n"%(all_NA_counter, all_del_counter))
	
	"""
fasta_input_fname = os.path.expanduser('~/script/variation/data/ACD6/ACD6.non_Kz10.sequenced.fasta')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/ACD6.non_Kz10.sequenced.del_NA_trimed.fasta')
AnalyzeSNPData.filterColumnsOfFasta(fasta_input_fname, output_fname, chromosome=4, start=8293290)
	"""
	
	@classmethod
	def removeInsertionSNPsFromSNPMatrix(cls, input_fname, output_fname):
		"""
		2009-5-29
			The most comprehensive way to represent a SNP is 'chr_pos_offset'.
			However, some of my programs can only handle 'chr_pos'.
			This function removes those SNPs embedded in insertion (relative to reference genome).
			
		"""
		sys.stderr.write("Removing SNPs embedded in insertions ")
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=input_fname, turn_into_array=1, ignore_2nd_column=1)
		col_id_to_be_kept_ls = []
		for col_id in snpData1.col_id_ls:
			col_id_tuple = col_id.split('_')
			offset = int(col_id_tuple[2])
			if offset==0:
				col_id_to_be_kept_ls.append(col_id)
		snpData2 = SNPData.keepColsByColID(snpData1, col_id_to_be_kept_ls)
		snpData2.tofile(output_fname)
	
	"""
input_fname = os.path.expanduser('~/script/variation/data/ACD6/ACD6.non_Kz10.sequenced.del_NA_trimed.with_ref.clustalx.snp_matrix.84_1530.tsv')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.no_insert.tsv')
AnalyzeSNPData.removeInsertionSNPsFromSNPMatrix(input_fname, output_fname)
	"""
	
	@classmethod
	def filterSNPMatrixBeforeImputation(cls, input_fname, output_fname):
		"""
		2009-5-29
			1. convert het to NA
			2. remove SNPs with >25% missing calls
			3. remove SNPs with >2 alleles cuz NPUTE doesn't support >2 alleles.
			4. remove SNPs with MAF<4%
		"""
		sys.stderr.write("Removing SNPs in preparation for imputation ...")
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=input_fname, turn_into_array=1, ignore_2nd_column=1)
		snpData2 = SNPData.convertHetero2NA(snpData1)
		snpData3 = SNPData.removeColsByNARate(snpData2, max_NA_rate=0.25)
		snpData4 = SNPData.removeSNPsWithMoreThan2Alleles(snpData3)
		snpData5 = SNPData.removeColsByMAF(snpData4, min_MAF=0.04)
		snpData5.tofile(output_fname)
		sys.stderr.write("Done.\n")
	
	"""
input_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.no_insert.tsv')
output_fname = os.path.expanduser('~/script/variation/data/ACD6/snp_matrix.84_1530.no_insert.no_het.maxNA0.25.noSNPsMoreThan2Alleles.minMAF0.04.tsv')
AnalyzeSNPData.filterSNPMatrixBeforeImputation(input_fname, output_fname)
	"""
	
	@classmethod
	def generateDatasetWithImputedCallsOnly(cls, unImputedFname, imputedFname, output_fname):
		"""
		2009-5-22
			generate a dataset which only contains imputed calls and masks everything else as NA (missing) 
				unImputedFname and imputedFname are of same shape, with same accessions & SNPs.
				The former contains the data before imputation. The latter contains the one after imputation. 
		"""
		sys.stderr.write("Calculating stats of imputed calls in the dataset ... \n")
		from pymodule import SNPData
		snpData1 = SNPData(input_fname=unImputedFname, turn_into_array=1)
		snpData2 = SNPData(input_fname=imputedFname, turn_into_array=1)
		row_id2stats = {}
		
		import numpy, copy
		newSnpData = SNPData(col_id_ls=copy.deepcopy(snpData1.col_id_ls), row_id_ls=snpData1.row_id_ls)
		newSnpData.data_matrix = numpy.zeros(snpData1.data_matrix.shape, numpy.int8)
		row_index = 0
		for row_id, row_index1 in snpData1.row_id2row_index.iteritems():
			row_index2 = snpData2.row_id2row_index[row_id]
			if row_id not in row_id2stats:
				row_id2stats[row_id] = 0
			for col_id, col_index1 in snpData1.col_id2col_index.iteritems():
				col_index2 = snpData2.col_id2col_index[col_id]
				allele_unimputed = snpData1.data_matrix[row_index1][col_index1]
				allele_imputed = snpData2.data_matrix[row_index2][col_index2]
				if (allele_unimputed==0 or allele_unimputed==-2) and allele_imputed!=0 and allele_imputed!=-2:
					#If before imputation it's NA (missing), keep the imputed call.
					row_id2stats[row_id] += 1
					newSnpData.data_matrix[row_index1][col_index1] = allele_imputed
		newSnpData.tofile(output_fname)
		sys.stderr.write("Done.\n")
		return row_id2stats
	
	"""
unImputedFname = '/Network/Data/250k/db/dataset/call_method_35.tsv'
imputedFname = '/Network/Data/250k/db/dataset/call_method_33.tsv'
output_fname = '/tmp/call_method_33_imputed_only.tsv'
row_id2no_of_imputed = AnalyzeSNPData.generateDatasetWithImputedCallsOnly(unImputedFname, imputedFname, output_fname)
	"""
	
	@classmethod
	def cmpTwoSNPDatasets(cls, inputFname1, inputFname2):
		"""
		2009-6-12
			compare two SNP datasets, report:
				#rows deleted/added
				#columns deleted/added
				#SNPs changed (from which to which)
		"""
		sys.stderr.write("Comparing two SNP datasets ... \n")
		from pymodule import SNPData, TwoSNPData, PassingData
		snpData1 = SNPData(input_fname=inputFname1, turn_into_array=1)
		snpData2 = SNPData(input_fname=inputFname2, turn_into_array=1)
		row_id_deleted = []
		row_id_added = []
		col_id_deleted = []
		col_id_added = []
		total_row_id_set = set(snpData1.row_id_ls) | set(snpData2.row_id_ls)
		total_col_id_set = set(snpData1.col_id_ls) | set(snpData2.col_id_ls)
		for row_id in total_row_id_set:
			if row_id not in snpData1.row_id2row_index:
				row_id_added.append(row_id)
			elif row_id not in snpData2.row_id2row_index:
				row_id_deleted.append(row_id)
		for col_id in total_col_id_set:
			if col_id not in snpData1.col_id2col_index:
				col_id_added.append(col_id)
			elif col_id not in snpData2.col_id2col_index:
				col_id_deleted.append(col_id)
		
		twoSNPData = TwoSNPData(SNPData1=snpData1, SNPData2=snpData2)
		diff_data = twoSNPData.get_diff_matrix()
		diff_matrix = diff_data[0]
		import numpy
		sth_index_ls = [1] + range(3, diff_matrix.shape[1])	#"deletion" + 10 calls
		print sth_index_ls
		print diff_matrix
		counter_NA_to_sth = numpy.sum(diff_matrix[[0,2],:][:, sth_index_ls])
		counter_sth_to_NA = numpy.sum(diff_matrix[sth_index_ls,:][:, [0,2]])
		sub_diff_matrix = diff_matrix[sth_index_ls,:][:, sth_index_ls]
		counter_sth_to_sth_diff = numpy.sum(sub_diff_matrix) - numpy.sum(numpy.diagonal(sub_diff_matrix))
		
		print "%s rows deleted."%len(row_id_deleted)
		print "%s rows added."%len(row_id_added)
		print "%s cols deleted."%len(col_id_deleted)
		print "%s cols added."%len(col_id_added)
		print "%s SNPs from NA to sth."%counter_NA_to_sth
		print "%s SNPs from sth to NA."%counter_sth_to_NA
		print "%s SNPs from sth to sth different."%counter_sth_to_sth_diff
		return_data = PassingData(row_id_deleted=row_id_deleted, row_id_added=row_id_added, col_id_deleted=col_id_deleted, col_id_added=col_id_added)
		return_data.counter_NA_to_sth = counter_NA_to_sth
		return_data.counter_sth_to_NA = counter_sth_to_NA
		return_data.counter_sth_to_sth_diff = counter_sth_to_sth_diff
		sys.stderr.write("Done.\n")
		return return_data
	
	"""
	inputFname1 = '/Network/Data/250k/db/dataset/call_method_35.tsv'
	inputFname2 = '/Network/Data/250k/db/dataset/call_method_33.tsv'
	return_data = AnalyzeSNPData.cmpTwoSNPDatasets(inputFname1, inputFname2)
	"""
	
	@classmethod
	def cmpOneRowToTheOther(cls, inputFname, row_id1, row_id2):
		"""
		2009-6-17
			compare SNP data of one accession to the other in the same dataset
		"""
		sys.stderr.write("Comparing one row to the other ... \n")
		from pymodule import SNPData, TwoSNPData, PassingData
		row_id_key_set = set([row_id1, row_id2])
		snpData = SNPData(input_fname=inputFname, turn_into_array=1, row_id_key_set=row_id_key_set)
		twoSNPData = TwoSNPData(SNPData1=snpData, SNPData2=snpData)
		print twoSNPData.cmpOneRow(row_id1, row_id2)
	
	"""
	inputFname = '/Network/Data/250k/db/dataset/call_method_29.tsv'
	row_id1 = ('6910', '62')
	row_id2 = ('8290', '181')
	AnalyzeSNPData.cmpOneRowToTheOther(inputFname, row_id1, row_id2)
		
	inputFname = os.path.expanduser('~/mnt2/panfs/NPUTE_data/input/250k_l3_y.85_20091208.tsv')
	row_id1 = ('7034', '1338')	# Blh-1 from Versailles plate
	#row_id2 = ('8265', '243')	# another Blh-1
	row_id2 = ('7035', '336')	# Blh-2
	AnalyzeSNPData.cmpOneRowToTheOther(inputFname, row_id1, row_id2)
	"""
	
	@classmethod
	def cmpAllDuplicatesOfOneEcotype(cls, inputFname, ecotype_id_ls):
		"""
		2009-12-11
			For each ecotype_id in ecotype_id_ls, compare mismatch-rates between duplicates
		"""
		sys.stderr.write("Comparing one row to the other ... \n")
		from pymodule import SNPData, TwoSNPData, PassingData
		ecotype_id_set = set(ecotype_id_ls)
		def row_id_hash_func(row_id):
			return int(row_id[0])
		snpData = SNPData(input_fname=inputFname, turn_into_array=1, row_id_key_set=ecotype_id_set, row_id_hash_func=row_id_hash_func)
		
		ecotype_id2row_id_to_check_ls = {}
		for row_id in snpData.row_id_ls:
			ecotype_id = int(row_id[0])
			if ecotype_id in ecotype_id_set:
				if ecotype_id not in ecotype_id2row_id_to_check_ls:
					ecotype_id2row_id_to_check_ls[ecotype_id] = []
				ecotype_id2row_id_to_check_ls[ecotype_id].append(row_id)
		twoSNPData = TwoSNPData(SNPData1=snpData, SNPData2=snpData)
		for ecotype_id, row_id_to_check_ls in ecotype_id2row_id_to_check_ls.iteritems():
			if len(row_id_to_check_ls)>1:
				print "ecotype_id: %s"%ecotype_id
				no_of_arrays = len(row_id_to_check_ls)
				for i in range(no_of_arrays):
					for j in range(i+1, no_of_arrays):
						row_id1 = row_id_to_check_ls[i]
						row_id2 = row_id_to_check_ls[j]
						print "row_id1 %s vs row_id2 %s"%(row_id1, row_id2)
						print twoSNPData.cmpOneRow(row_id1, row_id2)
	
	"""
	inputFname = os.path.expanduser('~/mnt2/panfs/NPUTE_data/input/250k_l3_y.85_20091208.tsv')
	ecotype_id_ls = [8297, 7317, 6910, 8274, 6911, 6905, 7034, 6909, 6962, 7373, 7270, 6983, 6899]
	AnalyzeSNPData.cmpAllDuplicatesOfOneEcotype(inputFname, ecotype_id_ls)
	
	"""
	
	@classmethod
	def linkEcotypeIDFromSuziPhenotype(cls, fname_with_ID, fname_with_phenotype, output_fname):
		"""
		2009-7-31
			she gave me two files
				one has phenotype data and accession names but with no ecotype ID
				2nd is a map from accession name to ecotype ID
		"""
		sys.stderr.write("Linking accession names to ecotype ID ... ")
		import csv
		inf_phenotype = csv.reader(open(fname_with_phenotype, 'r'), delimiter='\t')
		accession_name_ls = []
		#skip two lines
		inf_phenotype.next()
		inf_phenotype.next()
		for row in inf_phenotype:
			accession_name_ls.append(row[0])
		del inf_phenotype
		
		inf_with_ID = csv.reader(open(fname_with_ID), delimiter='\t')
		inf_with_ID.next()
		accession_name2ecotype_id = {}
		for row in inf_with_ID:
			ecotype_id = row[0]
			accession_name = row[5]
			accession_name2ecotype_id[accession_name] = ecotype_id
		del inf_with_ID
		print accession_name2ecotype_id
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for accession_name in accession_name_ls:
			ecotype_id = accession_name2ecotype_id.get(accession_name)
			writer.writerow([accession_name, ecotype_id])
		del writer
		
		sys.stderr.write("Done.\n")
	
	"""
fname_with_ID = '/tmp/seed_batch3_out.txt'
fname_with_phenotype = '/tmp/ft_LN_other phenotypes_10_batch_3_complete_7_29th.txt'
output_fname = '/tmp/batch_3_name2id.tsv'
AnalyzeSNPData.linkEcotypeIDFromSuziPhenotype(fname_with_ID, fname_with_phenotype, output_fname)
	"""
	
	
#2007-03-05 common codes to initiate database connection
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

import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment

"""
from codense.common import db_connect, form_schema_tables
hostname='dl324b-1'
dbname='yhdb'
schema = 'dbsnp'
hostname='zhoudb'
dbname='graphdb'
schema = 'dbsnp'
#conn, curs = db_connect(hostname, dbname, schema)

hostname='localhost'
dbname='stock20070829'
import MySQLdb
conn0 = MySQLdb.connect(db=dbname,host=hostname)
curs0 = conn0.cursor()

hostname='papaya.usc.edu'
dbname='stock_250k'
db_user='yh'
db_passwd = ''
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname, user=db_user, passwd=db_passwd)
curs = conn.cursor()

drivername='mysql'
schema = None
import Stock_250kDB
db_250k = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
db_250k.setup(create_tables=False)

drivername='mysql'
hostname='papaya.usc.edu'
dbname='stock'
db_user='yh'
db_passwd = ''
schema = None
from StockDB import StockDB
db_149 = StockDB(drivername=drivername, username=db_user,
				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
db_149.setup(create_tables=False)

"""

if __name__ == '__main__':
	ecotype_table = 'ecotype'
	calls_table = 'calls'
	#strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
	#strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
	#strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)
	
	#### 2009-5-2 code to check enrichment data by 3D plot. run by mayavi2 -x misc.py
	hostname='papaya.usc.edu'
	dbname='stock_250k'
	db_user='yh'
	db_passwd = 'yh324'

	drivername='mysql'
	schema = None
	import Stock_250kDB
	db_250k = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
					password=db_passwd, hostname=hostname, database=dbname, schema=schema)
	db_250k.setup(create_tables=False)
	
	#import pdb
	#pdb.set_trace()
	
	ler_blast_result_fname = '/Network/Data/250k/tmp-dazhe/ler_raw_CNV_QC.csv'
	max_delta_ratio = 0.4
	max_length_delta = 10000
	for max_length_delta in range(1,7):		# range(1,11)+[12,15,20]:
		max_length_delta = max_length_delta*10000
		output_fname = '/tmp/Ler-span-over-Col-mdr%s-mld%s.tsv'%(max_delta_ratio, max_length_delta)
		CNV.discoverLerContigSpanOverCol(ler_blast_result_fname, output_fname, min_no_of_matches=25,
										max_delta_ratio=max_delta_ratio, max_length_delta=max_length_delta)
	
	"""
	max_copy_number = 16
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_intensity_QNorm_sub_ref_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number)
	
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Col_intensity_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6909, data_source_id=9, cnv_type_id=None, max_copy_number=max_copy_number)
	
	input_fname_ls = [os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity.tsv')]
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number)
	input_fname_ls = []
	for i in range(1,6):
		input_fname_ls.append(os.path.expanduser('~/mnt2/panfs/250k/CNV/call_method_48_CNV_intensity_QNorm_sub_ref_chr%s.tsv'%i))
	output_fname_prefix = os.path.expanduser('~/tmp/call_48_Ler_intensity_QNorm_sub_ref_vs_true_copy_number_m%s'%max_copy_number)
	CNV.drawIntensityVsProbeTrueCopyNumber(db_250k, input_fname_ls, output_fname_prefix, \
										ecotype_id=6932, data_source_id=8, cnv_type_id=None, max_copy_number=max_copy_number)
	"""
	