#!/usr/bin/env python
"""
Usage: DrawPopulation.py [OPTIONS] -o -p

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-p ...,	popid2ecotypeid_table
	-s ...,	strain_info table, 'ecotype'(default)
	-o ...,	output_fname_prefix
	-a ...,	picture area, lower left corner longitude, latitude,
		upper right corner longitude, latitude default: "-130,10,140,70"
	-l ...,	type of label, 1(popid), 2(pop size), 3(selfing rate)
	-g ...,	max_dist to connect 2 sites, for -w, '50'(km) (default)
	-y ...,	which method estimates selfing rate, 1(default)
	-f ...,	selfing_rate_table (needed for label type 3)
	-w	draw site network showing how nearby sites give rise to population
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	DrawPopulation.py -d stock20070829 -o s0829popid2ecotypeid_25 -p popid2ecotypeid_25
	
	DrawPopulation.py -d stock20070829 -o s0829popid2ecotypeid_25_Eur -p popid2ecotypeid_25 -a "-10,30,30,70" -w -g 25
	
	DrawPopulation.py -d stock20070829 -o s0829popid2ecotypeid_10_Eng__7_48_3_57_l3y1 -p popid2ecotypeid_10 -a "-7,48,3,57" -l3 -fpopid2s_10
	
	DrawPopulation.py -d stock20070829 -o s0829popid2ecotypeid_10_NorAm__95_35__65_52_l3y1 -p popid2ecotypeid_10 -a "-95,35,-65,52"  -l3 -fpopid2s_10
	
	DrawPopulation.py -d stock20070829 -o s0829popid2ecotypeid_10_Eur__10_30_30_70_l3y1 -p popid2ecotypeid_10 -a "-10,30,30,70" -l3 -fpopid2s_10
	
Description:
	Divide strains into populations based on geographical distance. Two sites are connected
	if the great circle distance between them is less than max_dist. The resulting connected
	components from this network are populations.
	Draw Population on the map.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import psycopg, sys, getopt
from codense.common import db_connect, dict_map
from sets import Set
import networkx as nx

class DrawPopulation:
	"""
	2007-08-29
	"""
	def __init__(self, hostname='localhost', dbname='stock', schema='dbsnp', \
		popid2ecotypeid_table=None, strain_info_table='ecotype', output_fname_prefix=None,\
		pic_area=[-130,10,140,70],\
		label_type=1, max_dist=100, which_method=1, selfing_rate_table=None,\
		draw_site_network=0, debug=0, report=0):
		"""
		2007-08-29
		2007-08-30
			add which_method and selfing_rate_table
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.popid2ecotypeid_table = popid2ecotypeid_table
		self.strain_info_table = strain_info_table
		self.output_fname_prefix = output_fname_prefix
		self.pic_area = pic_area
		self.label_type = int(label_type)
		self.max_dist = int(max_dist)
		self.which_method = int(which_method)
		self.selfing_rate_table = selfing_rate_table
		self.draw_site_network = int(draw_site_network)
		self.debug = int(debug)
		self.report = int(report)
		
		self.label_type2label_name = {1:'popid',
			2: 'pop size',
			3: 'selfing rate'}
	
	def get_popid2pos_size(self, curs, popid2ecotypeid_table):
		sys.stderr.write("Getting popid2pos_size ...")
		popid2pos_size = {}
		curs.execute("select popid, latitude, longitude, ecotypeid from %s where selected=1"%(popid2ecotypeid_table))
		rows = curs.fetchall()
		for row in rows:
			popid, latitude, longitude, ecotypeid = row
			if popid not in popid2pos_size:
				popid2pos_size[popid] = [None, 0]
			if popid2pos_size[popid][0]==None:
				popid2pos_size[popid][0] = (latitude, longitude)
			popid2pos_size[popid][1] += 1
		sys.stderr.write("Done.\n")
		return popid2pos_size
	
	def get_popid2selfing_rate(self, curs, selfing_rate_table, which_method):
		"""
		2007-08-30
		"""
		sys.stderr.write("Getting popid2selfing_rate ...")
		popid2selfing_rate = {}
		from EstimateSelfingRate import EstimateSelfingRate
		EstimateSelfingRate_i = EstimateSelfingRate()
		curs.execute("select popid, %s from %s"%(EstimateSelfingRate_i.method2table_entry[which_method][0], selfing_rate_table))
		rows = curs.fetchall()
		for row in rows:
			popid, avg_s = row
			popid2selfing_rate[popid] = avg_s
		sys.stderr.write("Done.\n")
		return popid2selfing_rate
	
	def draw_clustered_strain_location(self, label_ls, weighted_pos_ls, diameter_ls, label_type, label_type2label_name, pic_area=[-180,-90,180,90], output_fname_prefix=None, label_name=None):
		"""
		2007-07-11
		draw populations derived from connected_components of the strain network
		#each pie denotes a population, with diameter proportional to the size of the population
		#each pie labeled with the number of strains in that population
		2007-07-13
			use popid as label
		2007-07-17
			no parallels, no meridians
		2007-08-29
			copied from CreatePopulation.py
		2007-09-11
			add label name
		2007-10-14
			correct a bug in 5*diameter_ls. diameter_ls has to an array

		"""
		sys.stderr.write("Drawing population map...")
		import pylab
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		
		euc_coord1_ls = []
		euc_coord2_ls = []
		ax=pylab.gca()
		for i in range(len(weighted_pos_ls)):
			lat, lon = weighted_pos_ls[i]
			euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
			ax.text(euc_coord1, euc_coord2, str(label_ls[i]), size=5, alpha=0.5, horizontalalignment='center', verticalalignment='center', zorder=12)
		import numpy
		diameter_ls = numpy.array(diameter_ls)
		m.scatter(euc_coord1_ls, euc_coord2_ls, 5*diameter_ls, marker='o', color='r', alpha=0.3, zorder=10, faceted=False)
		
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])	#labels intersect the left, right, top bottom of the plot
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()
		pylab.title("worldwide distribution of %s populations, labeled by %s"%(len(weighted_pos_ls), label_type2label_name[label_type]))
		if output_fname_prefix:
			pylab.savefig('%s_pop_map.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_pop_map.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_pop_map.png'%output_fname_prefix, dpi=300)
		del m, pylab, Basemap
		sys.stderr.write("Done.\n")
	
	#2007-07-09
	def DrawSiteNetwork(self, g, node_label2pos_counts,pic_area=[-180,-90,180,90], output_fname_prefix=None):
		"""
		2007-07-17
			put ax.plot() right after Basemap() but after m.xxx() so that it'll zoom in
			use 'g' in ax.plot(), otherwise, ax.plot() alternates all colors.
			no parallels, no meridians
		2007-08-29 copied from CreatePopulation.py and renamed from DrawStrainNetwork
		"""
		sys.stderr.write("Drawing Site Network...")
		import pylab
		from matplotlib.toolkits.basemap import Basemap
		pylab.clf()
		fig = pylab.figure()
		fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill')
		
		ax=pylab.gca()
		for e in g.edges():
			lat1, lon1 = node_label2pos_counts[e[0]][0]
			lat2, lon2 = node_label2pos_counts[e[1]][0]
			x1, y1 = m(lon1, lat1)
			x2, y2 = m(lon2, lat2)
			ax.plot([x1,x2],[y1,y2], 'g', alpha=0.5, zorder=12)
		
		#m.drawcoastlines()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
		m.fillcontinents()
		m.drawcountries()
		m.drawstates()

		pylab.title("Network of strains")
		if output_fname_prefix:
			pylab.savefig('%s_site_network.eps'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_site_network.svg'%output_fname_prefix, dpi=300)
			pylab.savefig('%s_site_network.png'%output_fname_prefix, dpi=300)
		del m, pylab, Basemap
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		2007-08-30
			add label_type 3
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname)
		curs = conn.cursor()
		popid2pos_size = self.get_popid2pos_size(curs, self.popid2ecotypeid_table)
		popid_ls = popid2pos_size.keys()
		pos_size_ls = dict_map(popid2pos_size, popid_ls)
		weighted_pos_ls = [row[0] for row in pos_size_ls]
		diameter_ls = [row[1] for row in pos_size_ls]
		if self.label_type == 1:
			label_ls = popid_ls
		elif self.label_type == 2:
			label_ls = diameter_ls
		elif self.label_type == 3:
			if self.selfing_rate_table is None:
				sys.stderr.write("Label type is 3(selfing rate), but no selfing_rate_table specified\n")
				sys.exit(3)
			popid2selfing_rate = self.get_popid2selfing_rate(curs, self.selfing_rate_table, self.which_method)
			label_ls = []
			for popid in popid_ls:
				avg_s = '0'
				if popid in popid2selfing_rate:
					if popid2selfing_rate[popid]:	#not NULL
						avg_s = int(round(popid2selfing_rate[popid]*1000))
				label_ls.append(avg_s)
		self.draw_clustered_strain_location(label_ls, weighted_pos_ls, diameter_ls, self.label_type, self.label_type2label_name, pic_area=self.pic_area, output_fname_prefix=self.output_fname_prefix)
		if self.draw_site_network:
			from CreatePopulation import CreatePopulation
			CreatePopulation_instance = CreatePopulation()
			lat_lon_ls, pos2ecotypeid_ls = CreatePopulation_instance.get_pos2ecotypeid_ls(curs, self.strain_info_table)
			g, node_label2pos_counts = CreatePopulation_instance.divide_data_by_geography(lat_lon_ls, self.max_dist)
			self.DrawSiteNetwork(g, node_label2pos_counts, self.pic_area, self.output_fname_prefix)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "debug", "report", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:p:s:o:a:l:g:y:f:wbrh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'localhost'
	dbname = 'stock'
	schema = 'dbsnp'
	popid2ecotypeid_table = None
	strain_info_table = 'ecotype'
	output_fname_prefix = None
	pic_area = [-130,10,140,70]
	label_type = 1
	max_dist = 100
	which_method = 1
	selfing_rate_table = None
	draw_site_network = 0
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
		elif opt in ("-p",):
			popid2ecotypeid_table = arg
		elif opt in ("-s",):
			strain_info_table = arg
		elif opt in ("-o",):
			output_fname_prefix = arg
		elif opt in ("-a",):
			pic_area = arg.split(',')
			pic_area = map(float, pic_area)
		elif opt in ("-l",):
			label_type = int(arg)
		elif opt in ("-g",):
			max_dist = int(arg)
		elif opt in ("-y",):
			which_method = int(arg)
		elif opt in ("-f",):
			selfing_rate_table = arg
		elif opt in ("-w",):
			draw_site_network = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	
	if popid2ecotypeid_table and output_fname_prefix and hostname and dbname and schema:
		instance = DrawPopulation(hostname, dbname, schema, popid2ecotypeid_table,\
			strain_info_table, output_fname_prefix, pic_area, label_type, max_dist, \
			which_method, selfing_rate_table, draw_site_network, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
