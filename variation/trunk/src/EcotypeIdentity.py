#!/usr/bin/env python
"""
Usage: EcotypeIdentity.py [OPTIONS] -i

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-i ...,	input_fname
	-t ...,	identity table, 'ecotype_identity' (default)
	-p ...,	component2clique table, 'component2clique'(default)
	-q ...,	clique2ecotype table, 'clique2ecotype' (default)
	-c	commit the database submission
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	EcotypeIdentity.py -d stock20071008 -i stock20071008/data_d110_c0_5.tsv -c
	
Description:
	Partition strains into identity cliques. And store info into database.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/pymodule')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/pymodule')))
import psycopg, sys, getopt
from codense.common import db_connect, dict_map
from sets import Set
import networkx as nx

def get_continent_combo2clique_id_ls(curs, country2continent_table='at.country2continent', clique2ecotype_table='stock20071008.clique2ecotype', stock_db='stock20071008'):
	"""
	2007-10-16
	"""
	sys.stderr.write("Getting continent_combo2clique_id_ls ...")
	continent_combo2clique_id_ls = {}
	clique_id2continent_ls = {}
	curs.execute("select distinct ce.clique_id, cc.continent from %s cc, %s ce, %s.ecotype e, %s.country c, %s.site s, %s.address a where ce.ecotypeid=e.id and e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and c.id=cc.id"%(country2continent_table, clique2ecotype_table, stock_db, stock_db, stock_db, stock_db))
	rows = curs.fetchall()
	for row in rows:
		clique_id, continent = row
		if clique_id not in clique_id2continent_ls:
			clique_id2continent_ls[clique_id] = []
		clique_id2continent_ls[clique_id].append(continent)
	for clique_id, continent_ls in clique_id2continent_ls.iteritems():
		continent_ls.sort()
		continent_combo = tuple(continent_ls)
		if continent_combo not in continent_combo2clique_id_ls:
			continent_combo2clique_id_ls[continent_combo] = []
		continent_combo2clique_id_ls[continent_combo].append(clique_id)
	sys.stderr.write("Done.\n")
	return continent_combo2clique_id_ls


def get_clique_id2size(curs, clique2ecotype_table='stock20071008.clique2ecotype'):
	"""
	2007-10-16
	"""
	clique_id2size = {}
	curs.execute("select clique_id, count(ecotypeid) from %s group by clique_id"%(clique2ecotype_table))
	rows = curs.fetchall()
	for row in rows:
		clique_id, cnt = row
		clique_id2size[clique_id] = cnt
	return clique_id2size

def construct_site_graph_out_of_ecotype_id_ls(ecotype_id_ls, ecotypeid2pos):
	"""
	2007-10-18
		modified from construct_site_graph_out_of_strain_graph()
		node weighted by the number of ecotype_id's
	"""
	sys.stderr.write("Constructing site graph out of ecotype_id_ls...")
	site2pos = {}
	site2weight = {}
	import networkx as nx
	site_g = nx.XGraph()
	for e in ecotype_id_ls:
		pos = (round(ecotypeid2pos[e][0],2), round(ecotypeid2pos[e][1],2))
		site_g.add_node(pos)
		if pos not in site2weight:
			site2weight[pos] = 0
		site2weight[pos] += 1
		site2pos[pos] = pos
	sys.stderr.write("Done.\n")
	return site_g, site2weight, site2pos


def outputCliqueInLatexTable(curs, region, clique_id, table_label, fig_label, stock_db, clique2ecotype_table, outf):
	"""
	2007-10-18
	2007-10-21
		shorten some label names
			ecotypeid -> id
			stockparent -> stkparent
		add one more column, collector
	"""
	from pymodule.latex import outputMatrixInLatexTable
	data_matrix = []
	curs.execute("select e.id, e.name, e.nativename, e.stockparent, e.latitude, e.longitude, c.abbr, s.name, p.firstname, p.surname from %s.ecotype e, %s.address a, %s.site s, %s.country c, %s.person p, %s ce where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id and e.id=ce.ecotypeid and p.id=e.collectorid and ce.clique_id=%s order by nativename, stockparent"%(stock_db, stock_db, stock_db, stock_db, stock_db, clique2ecotype_table, clique_id))
	rows = curs.fetchall()
	ecotype_id_ls = []
	for row in rows:
		ecotype_id = row[0]
		ecotype_id_ls.append(ecotype_id)
		new_row = list(row[:-4])
		new_row.append('%s %s'%(row[-2], row[-1]))	#firstname surname
		new_row.append(row[-4])	#country
		new_row.append(row[-3])	#site
		data_matrix.append(new_row)
	caption = 'haplotype clique %s has %s ecotypes from %s. check Figure~\\ref{%s}.'%(clique_id, len(ecotype_id_ls), region, fig_label)
	header_ls = ['id', 'name', 'nativename', 'stkparent', 'lat', 'lon', 'collector',  'country', 'site']
	outf.write(outputMatrixInLatexTable(data_matrix, caption, table_label, header_ls))
	return ecotype_id_ls


def OutputCliquesGroupedByRegion(curs, region2clique_id_ls, clique_id2size, clique2ecotype_table, component2clique_table, output_fname, fig_output_dir, stock_db='stock20071008', min_clique_size=3, need_to_draw_figures=0):
	"""
	2007-10-18
	"""
	from variation.src.common import draw_graph_on_map, get_ecotypeid2pos
	from variation.src.common import get_pic_area
	from pymodule.latex import outputFigureInLatex
	import sys,os
	ecotypeid2pos = get_ecotypeid2pos(curs, '%s.ecotype'%stock_db)
	
	outf = open(output_fname, 'w')
	no_of_regions = 0
	clique_id2region = {}	#for the cc2clique latex table
	for region, clique_id_ls in region2clique_id_ls.iteritems():
		region = ' and '.join(region)	#region is tuple, convert it to string
		sys.stderr.write("%s:\n"%region)
		outf.write("\\subsection{%s}\n"%region)
		size_clique_id_ls = []
		for clique_id in clique_id_ls:
			if clique_id2size[clique_id]>=min_clique_size:
				size_clique_id_ls.append([clique_id2size[clique_id], clique_id])
		outf.write('region %s has %s cliques in total but %s cliques above min size %s.\n'%(region, len(clique_id_ls), len(size_clique_id_ls), min_clique_size))
		if len(size_clique_id_ls)==0:
			continue
		size_clique_id_ls.sort()
		caption = 'region %s has %s cliques in total but %s cliques above min size %s.'%(region, len(clique_id_ls), len(size_clique_id_ls), min_clique_size)
		no_of_regions += 1
		table_label = 't_region_%s'%no_of_regions
		clique_id_size_ls = [[row[1],row[0]] for row in size_clique_id_ls]
		header_ls = ['clique_id', 'size']
		outf.write(outputMatrixInLatexTable(clique_id_size_ls, caption, table_label, header_ls))
		
		for size, clique_id in size_clique_id_ls:
			sys.stderr.write("clique %s, size %s\n"%(clique_id, size))
			#output the table and get ecotype_id_ls
			fig_label = 'fclique%s'%clique_id
			table_label = 'tclique%s'%clique_id
			ecotype_id_ls = outputCliqueInLatexTable(curs, region, clique_id, table_label, fig_label, stock_db, clique2ecotype_table, outf)
			
			#draw the clique
			site_g, site2weight, site2pos = construct_site_graph_out_of_ecotype_id_ls(ecotype_id_ls, ecotypeid2pos)
			fig_fname = 'haplotype_%s_size%s_%ssites_map'%(clique_id,len(ecotype_id_ls), len(site2pos))	#just the name without filetype extension, like png, eps etc.
			if need_to_draw_figures:
				draw_graph_on_map(site_g, site2weight, site2pos, 'clique %s, size %s, %s sites of %s'%(clique_id, len(ecotype_id_ls), len(site2pos), region), pic_area=get_pic_area(site2pos.values(), 30), output_fname_prefix=os.path.join(fig_output_dir, fig_fname))
				
				#output the figure in latex
				fig_fname = os.path.join('figures/', '%s.png'%fig_fname)	#relative path of the figure to the latex file
				caption = 'map of haplotype clique %s has %s ecotypes from %s sites of %s. check Table~\\ref{%s}.'%(clique_id, len(ecotype_id_ls), len(site2pos), region, table_label)
				outf.write(outputFigureInLatex(fig_fname, caption, fig_label))
			clique_id2region[clique_id] = region
			del site_g, site2pos, site2weight, ecotype_id_ls
	outf.write("\\subsection{which component the clique is from}\n")
	cc_id_clique_id_region_ls = []
	cc_id2clique_id_ls = {}
	for clique_id, region in clique_id2region.iteritems():
		curs.execute("select cc_id from %s where clique_id=%s"%(component2clique_table, clique_id))
		rows = curs.fetchall()
		cc_id = rows[0][0]
		cc_id_clique_id_region_ls.append([cc_id, clique_id, region])
		if cc_id not in cc_id2clique_id_ls:
			cc_id2clique_id_ls[cc_id] = []
		cc_id2clique_id_ls[cc_id].append(clique_id)
	cc_id_clique_id_region_ls.sort()
	
	caption = '%s components versus %s cliques'%(len(cc_id2clique_id_ls), len(cc_id_clique_id_region_ls))
	table_label = 't_cc2clique_1'
	header_ls = ['component id', 'clique id', 'region']
	outf.write(outputMatrixInLatexTable(cc_id_clique_id_region_ls, caption, table_label, header_ls))
	del outf

"""
#2007-10-16
import os,sys
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule.latex import outputMatrixInLatexTable

hostname='localhost'
dbname='stock20071008'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()

country2continent_table='at.country2continent'
clique2ecotype_table='stock20071008.clique2ecotype'
stock_db='stock20071008'
continent_combo2clique_id_ls = get_continent_combo2clique_id_ls(curs, country2continent_table, clique2ecotype_table, stock_db)

clique2ecotype_table='stock20071008.clique2ecotype'
clique_id2size = get_clique_id2size(curs, clique2ecotype_table)
component2clique_table = 'stock20071008.component2clique'
stock_db = 'stock20071008'
output_fname = 'script/variation/doc/StockIdentityReport/tables_figures.tex'
fig_output_dir = 'script/variation/doc/StockIdentityReport/figures'
OutputCliquesGroupedByRegion(curs, continent_combo2clique_id_ls, clique_id2size, clique2ecotype_table, component2clique_table, output_fname, fig_output_dir, stock_db, min_clique_size=3, need_to_draw_figures=1)
"""
class EcotypeIdentity:
	"""
	2007-10-11
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', \
		input_fname=None, identity_table='ecotype_identity', component2clique_table='component2clique', clique2ecotype_table='clique2ecotype',\
		commit=0, debug=0, report=0):
		"""
		2007-10-11
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.identity_table = identity_table
		self.component2clique_table = component2clique_table
		self.clique2ecotype_table = clique2ecotype_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def construct_identity_pair_ls(self, strain_acc_list, header, data_matrix):
		"""
		2007-09-13
			input_fname is a data frame file with snps coded in integers. 1st row is the ecotype id in table ecotype
		2007-10-11 copied from misc.py
		"""
		sys.stderr.write("Constructing identity_pair_ls ...")
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
		sys.stderr.write("Done.\n")
		return identity_pair_ls
	
	def construct_graph_out_of_identity_pair(self, identity_pair_ls):
		"""
		2007-09-18
			construct a graph based on identity pairs in identity_pair_ls
			identity_pair_ls's entries will be converted to integers
		2007-10-11
			copied from misc.py
		"""
		sys.stderr.write("Constructing graph out of identity_pair_ls ...")
		g = nx.Graph()
		for identity_pair in identity_pair_ls:
			identity_pair = map(int, identity_pair)
			g.add_edge(identity_pair[0], identity_pair[1])
		sys.stderr.write("Done.\n")
		return g
	
	def expand_g_with_singleton_strain_id_ls(self, g, strain_acc_list):
		"""
		2007-10-11
			node in g is in ecotype id format
			so strain_acc (1st column) in input_fname should be integer (ecotype id)
		"""
		sys.stderr.write("Expand g with singleton_strain_id_ls ...")
		graph_node_set = Set(g.nodes())
		singleton_strain_id_ls = []
		for strain_acc in strain_acc_list:
			strain_acc = int(strain_acc)
			if strain_acc not in graph_node_set:
				g.add_node(strain_acc)
		sys.stderr.write("Done.\n")
		return g
	
	def compute_components_and_cliques(self, g):
		"""
		2007-10-11 each component is partitioned into cliques
		"""
		sys.stderr.write("Computing components and cliques from graph ...\n")
		cc_id2clique_id_ls = {}
		clique_id2ecotype_id_ls = {}
		from PartitionGraphIntoCliques import PartitionGraphIntoCliques
		PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(algorithm_type=2)
		g_cc_ls = nx.connected_components(g)
		clique_id = 0
		for i in range(len(g_cc_ls)):
			sys.stderr.write("%s\t%s"%('\x08'*20, i))
			cc_g = g.subgraph(g_cc_ls[i])
			PartitionGraphIntoCliques_ins.partition(cc_g.copy())
			clique_ls = PartitionGraphIntoCliques_ins.clique_ls
			cc_id2clique_id_ls[i+1] = []
			for j in range(len(clique_ls)):
				clique_id += 1
				cc_id2clique_id_ls[i+1].append(clique_id)
				clique_id2ecotype_id_ls[clique_id] = clique_ls[j]
		sys.stderr.write("Done.\n")
		return cc_id2clique_id_ls, clique_id2ecotype_id_ls
	
	def create_identity_table(self, curs, identity_table):
		"""
		2007-10-11
		"""
		sys.stderr.write("Creating identity_table %s ..."%identity_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			ecotypeid1	integer not null,\
			ecotypeid2	integer not null)"%identity_table)
		sys.stderr.write("Done.\n")
	
	def create_component2clique_table(self, curs, component2clique_table):
		"""
		2007-10-11
		"""
		sys.stderr.write("Creating component2clique_table %s ..."%component2clique_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			cc_id	integer	not null,\
			clique_id	integer not null)"%component2clique_table)
		sys.stderr.write("Done.\n")
		
	def create_clique2ecotype_table(self, curs, clique2ecotype_table):
		"""
		2007-10-11
		"""
		sys.stderr.write("Creating clique2ecotype_table %s ..."%clique2ecotype_table)
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			clique_id	integer	not null,\
			ecotypeid	integer not null)"%clique2ecotype_table)
		sys.stderr.write("Done.\n")
	
	def submit_identity_pairs(self, curs, g, identity_table):
		sys.stderr.write("Submitting identity pairs ...")
		for e in g.edges():
			ecotypeid1 = min(e)
			ecotypeid2 = max(e)
			curs.execute("insert into %s(ecotypeid1, ecotypeid2) values(%s, %s)"%\
			(identity_table, ecotypeid1, ecotypeid2))
		sys.stderr.write("Done.\n")
	
	def submit_cc_id2clique_id_ls(self, curs, cc_id2clique_id_ls, component2clique_table):
		sys.stderr.write("Submitting cc_id2clique_id_ls ...")
		for cc_id, clique_id_ls in cc_id2clique_id_ls.iteritems():
			for clique_id in clique_id_ls:
				curs.execute("insert into %s(cc_id, clique_id) values(%s, %s)"%\
				(component2clique_table, cc_id, clique_id))
		sys.stderr.write("Done.\n")
	
	def submit_clique_id2ecotype_id_ls(self, curs, clique_id2ecotype_id_ls, clique2ecotype_table):
		sys.stderr.write("Submitting clique_id2ecotype_id_ls ...")
		for clique_id, ecotype_id_ls in clique_id2ecotype_id_ls.iteritems():
			for ecotypeid in ecotype_id_ls:
				curs.execute("insert into %s(clique_id, ecotypeid) values(%s, %s)"%\
				(clique2ecotype_table, clique_id, ecotypeid))
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2007-10-11
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(self.input_fname)
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		identity_pair_ls = self.construct_identity_pair_ls(strain_acc_list, header, data_matrix)
		g = self.construct_graph_out_of_identity_pair(identity_pair_ls)
		g = self.expand_g_with_singleton_strain_id_ls(g, strain_acc_list)
		cc_id2clique_id_ls, clique_id2ecotype_id_ls = self.compute_components_and_cliques(g)
		
		if self.commit:
			self.create_identity_table(curs, self.identity_table)
			self.create_component2clique_table(curs, self.component2clique_table)
			self.create_clique2ecotype_table(curs, self.clique2ecotype_table)
			self.submit_identity_pairs(curs, g, self.identity_table)
			self.submit_cc_id2clique_id_ls(curs, cc_id2clique_id_ls, self.component2clique_table)
			self.submit_clique_id2ecotype_id_ls(curs, clique_id2ecotype_id_ls, self.clique2ecotype_table)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:t:p:q:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock'
	schema = 'dbsnp'
	input_fname = None
	identity_table = 'ecotype_identity'
	component2clique_table='component2clique'
	clique2ecotype_table='clique2ecotype'
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
		elif opt in ("-t",):
			identity_table = arg
		elif opt in ("-p",):
			component2clique_table = arg
		elif opt in ("-q",):
			clique2ecotype_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if input_fname and hostname and dbname and schema:
		instance = EcotypeIdentity(hostname, dbname, schema, input_fname, identity_table,\
			component2clique_table, clique2ecotype_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)