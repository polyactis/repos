#!/usr/bin/env python
"""
Usage: CreatePopulation.py [OPTIONS] -o 

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-o ...,	output table
	-s ...,	strain_info table, 'ecotype'(default)
	-g ...,	max_dist to connect 2 sites, '100'(km) (default)
	-c,	commit
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	CreatePopulation.py -o popid2ecotypeid_50 -g 50
	
Description:
	Divide strains into populations based on geographical distance. Two sites are connected
	if the great circle distance between them is less than max_dist. The resulting connected
	components from this network are populations.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import psycopg, sys, getopt, csv, re
from codense.common import db_connect
from common import nt2number, number2nt
from sets import Set
import networkx as nx

class CreatePopulation:
	"""
	2007-07-11
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', input_fname=None, \
		output_table=None, strain_info_table='ecotype', max_dist=100, commit=0, debug=0, report=0):
		"""
		2007-07-11
		2007-07-13
			input_fname is useless
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_table = output_table
		self.strain_info_table = strain_info_table
		self.max_dist = int(max_dist)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def cal_great_circle_distance(self, lat1, lon1, lat2, lon2, earth_radius=6372.795):
		"""
		2007-06-17 copied from 2007-07-11
		http://en.wikipedia.org/wiki/Great-circle_distance
		"""
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
	
	def get_pos2ecotypeid_ls(self, curs, strain_info_table):
		"""
		2007-07-11
		2007-07-13
			not restricted by strain_acc_list
			select all strains with latitude and longitude
		"""
		sys.stderr.write("Fetching ecotype, latitude, longitude ...")
		lat_lon_ls = []
		pos2ecotypeid_ls = {}
		curs.execute("select id, latitude, longitude from %s where latitude is not null and longitude is not null"%(strain_info_table))
		rows = curs.fetchall()
		for row in rows:
			ecotypeid, latitude, longitude = row
			pos = (latitude, longitude)
			if pos not in pos2ecotypeid_ls:
				pos2ecotypeid_ls[pos] = []
			pos2ecotypeid_ls[pos].append(ecotypeid)
			lat_lon_ls.append(pos)
		sys.stderr.write("Done.\n")
		return lat_lon_ls, pos2ecotypeid_ls
	
	def divide_data_by_geography(self, lat_lon_ls, max_dist=100):
		"""
		2007-07-11 copied from misc.py
		"""
		sys.stderr.write("Constructing data structure on lat_lon_ls ...")
		#each distinctive position (lat,lon) gets a unique label (an integer starting from 0)
		pos2node_label = {}
		node_label2pos_counts = {}	#records how many strains are from the same position
		for lat, lon in lat_lon_ls:
			pos = (lat, lon)
			if pos not in pos2node_label:
				pos2node_label[pos] = len(pos2node_label)
			node_label = pos2node_label[pos]
			if node_label not in node_label2pos_counts:
				node_label2pos_counts[node_label] = [pos, 0]
			node_label2pos_counts[node_label][1] += 1
		sys.stderr.write("Done.\n")
		sys.stderr.write("Constructing graph ...")
		g = nx.Graph()
		pos_ls = pos2node_label.keys()
		no_of_sites = len(pos_ls)
		for i in range(no_of_sites):
			g.add_node(pos2node_label[pos_ls[i]])
		
		for i in range(no_of_sites):
			for j in range(i+1, no_of_sites):
				dist = self.cal_great_circle_distance(pos_ls[i][0], pos_ls[i][1], pos_ls[j][0], pos_ls[j][1])
				if dist<=max_dist:
					g.add_edge(pos2node_label[pos_ls[i]], pos2node_label[pos_ls[j]])
		sys.stderr.write("Done.\n")
		return g, node_label2pos_counts
	
	def get_pop_center_pos2ecotypeid_ls(self, g, node_label2pos_counts, pos2ecotypeid_ls):
		"""
		2007-07-11
		"""
		sys.stderr.write("Computing weighted centers...")
		pop_id2center_pos_ecotypeid_ls = {}
		weighted_pos_ls = []
		count_sum_ls = []
		popid_ls = []
		no_of_populations = 0
		c_components = nx.connected_components(g)
		for component in c_components:
			#weighted average
			count_sum = 0
			lat_sum = 0.0
			lon_sum = 0.0
			ecotypeid_ls = []
			for node_label in component:
				count_sum += node_label2pos_counts[node_label][1]
				node_pos = node_label2pos_counts[node_label][0]
				lat_sum += node_pos[0]*node_label2pos_counts[node_label][1]
				lon_sum += node_pos[1]*node_label2pos_counts[node_label][1]
				ecotypeid_ls += pos2ecotypeid_ls[node_pos]
			pop_center_pos = (lat_sum/count_sum, lon_sum/count_sum)
			no_of_populations += 1
			pop_id2center_pos_ecotypeid_ls[no_of_populations] = [pop_center_pos, ecotypeid_ls]
			popid_ls.append(no_of_populations)
			weighted_pos_ls.append(pop_center_pos)
			count_sum_ls.append(count_sum)
		sys.stderr.write("%s populations. Done.\n"%(len(weighted_pos_ls)))
		return pop_id2center_pos_ecotypeid_ls, weighted_pos_ls, count_sum_ls, popid_ls
		

	def draw_clustered_strain_location(self, popid_ls, weighted_pos_ls, count_sum_ls, pic_area=[-180,-90,180,90], output_fname_prefix=None):
		"""
		2007-07-11
		draw populations derived from connected_components of the strain network
		#each pie denotes a population, with diameter proportional to the size of the population
		#each pie labeled with the number of strains in that population
		2007-07-13
			use popid as label
		"""
		sys.stderr.write("Drawing population map...")
		import pylab
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.1,0.1,0.9,0.9])
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,1,1])
		m.drawmeridians(pylab.arange(-180,180,60), labels=[1,1,1,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		euc_coord1_ls = []
		euc_coord2_ls = []
		ax=pylab.gca()
		for i in range(len(weighted_pos_ls)):
			lat, lon = weighted_pos_ls[i]
			euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
			ax.text(euc_coord1, euc_coord2, str(popid_ls[i]), size=3, alpha=0.5, horizontalalignment='center', verticalalignment='center', zorder=12)
		m.scatter(euc_coord1_ls, euc_coord2_ls, 5*count_sum_ls, marker='o', color='r', alpha=0.3, zorder=10, faceted=False)
		pylab.title("worldwide distribution of %s populations, labeled by popid"%(len(weighted_pos_ls)))
		if output_fname_prefix:
			pylab.savefig('%s_pop_map.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_pop_map.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_pop_map.png'%output_fname_prefix, dpi=300)
		del m, pylab, Basemap
		sys.stderr.write("Done.\n")
	
	#2007-07-09
	def DrawStrainNetwork(self, g, node_label2pos_counts,pic_area=[-180,-90,180,90], output_fname_prefix=None):
		sys.stderr.write("Drawing Strain Network...")
		import pylab
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.1,0.1,0.9,0.9])
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,1,1])
		m.drawmeridians(pylab.arange(-180,180,60), labels=[1,1,1,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		ax=pylab.gca()
		for e in g.edges():
			lat1, lon1 = node_label2pos_counts[e[0]][0]
			lat2, lon2 = node_label2pos_counts[e[1]][0]
			x1, y1 = m(lon1, lat1)
			x2, y2 = m(lon2, lat2)
			ax.plot([x1,x2],[y1,y2], alpha=0.5, zorder=12)
		pylab.title("Network of strains")
		if output_fname_prefix:
			pylab.savefig('%s_strain_network.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_strain_network.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_strain_network.png'%output_fname_prefix, dpi=300)
		del m, pylab, Basemap
		sys.stderr.write("Done.\n")
	
	def create_population_table(self, curs, output_table):
		"""
		2007-07-13
			add selected to the output_table
		"""
		sys.stderr.write("Creating table %s ..."%(output_table))
		curs.execute("create table %s(\
			id	integer primary key auto_increment,\
			popid	integer not null,\
			latitude	float,\
			longitude	float,\
			ecotypeid	integer not null,\
			selected	integer default 0)"%output_table)
		sys.stderr.write("Done.\n")
	
	def submit_pop_id2center_pos_ecotypeid_ls(self, curs, output_table, pop_id2center_pos_ecotypeid_ls):
		"""
		2007-07-11
		"""
		sys.stderr.write("Submitting to table %s ..."%(output_table))
		for pop_id, center_pos_ecotypeid_ls in pop_id2center_pos_ecotypeid_ls.iteritems():
			center_pos, ecotypeid_ls = center_pos_ecotypeid_ls
			for ecotypeid in ecotypeid_ls:
				curs.execute("insert into %s(popid, latitude, longitude, ecotypeid) values(%s, %s, %s, %s)"%\
				(output_table, pop_id, center_pos[0], center_pos[1], ecotypeid))
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2007-07-11
		2007-07-13
			input_fname is useless
		"""
		import MySQLdb
		conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		
		lat_lon_ls, pos2ecotypeid_ls = self.get_pos2ecotypeid_ls(curs, self.strain_info_table)
		g, node_label2pos_counts = self.divide_data_by_geography(lat_lon_ls, self.max_dist)
		pop_id2center_pos_ecotypeid_ls, weighted_pos_ls, count_sum_ls, popid_ls = self.get_pop_center_pos2ecotypeid_ls(g, node_label2pos_counts, pos2ecotypeid_ls)
		
		if self.commit:
			self.create_population_table(curs, self.output_table)
			self.submit_pop_id2center_pos_ecotypeid_ls(curs, self.output_table, pop_id2center_pos_ecotypeid_ls)
			#conn.commit()
		pic_area=[-120,10,100,80]
		self.draw_clustered_strain_location(popid_ls, weighted_pos_ls, count_sum_ls, pic_area, self.output_table)
		self.DrawStrainNetwork(g, node_label2pos_counts, pic_area, self.output_table)
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:o:s:g:cbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock'
	schema = 'dbsnp'
	input_fname = None
	output_table = None
	strain_info_table = 'ecotype'
	max_dist = 100
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
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-g",):
			max_dist = int(arg)
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if output_table and hostname and dbname and schema:
		instance = CreatePopulation(hostname, dbname, schema, input_fname, output_table, \
			strain_info_table, max_dist, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)