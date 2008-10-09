"""
2007-02-19
	common stuff for variation/src
"""


"""
2007-03-29
	add mappings for '-', 'N', and ambiguous letters ('M', 'R', 'W', 'S', 'Y', 'K')
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule.SNP import ab2number, number2ab, NA_set, nt2number, \
	number2nt, number2color, nt_number_matching_matrix, get_nt_number2diff_matrix_index,\
	RawSnpsData_ls2SNPData, transposeSNPData, SNPData2RawSnpsData_ls

def get_chr_id2size(curs, chromosome_table='at.chromosome'):
	"""
	2008-10-07 curs could be elixirdb.metadata.bind
	2007-10-12
	"""
	sys.stderr.write("Getting chr_id2size ...")
	rows = curs.execute("select id, size from %s"%chromosome_table)
	chr_id2size = {}
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#2008-10-07 curs could be elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	for row in rows:
		if is_elixirdb:
			chr_id = row.id
			size = row.size
		else:
			chr_id, size = row
		chr_id2size[chr_id] = size
	sys.stderr.write("Done.\n")
	return chr_id2size

def get_chr_id2cumu_size(chr_id2size, chr_gap=None):
	"""
	2008-02-04
		add chr_id_ls
		turn chr_id all into 'str' form
	2008-02-01
		add chr_gap, copied from variation.src.misc
	2007-10-16
	"""
	sys.stderr.write("Getting chr_id2cumu_size ...")
	#if chr_gap not specified, take the 1/5th of the average chromosome size as its value
	if chr_gap==None:
		chr_size_ls = chr_id2size.values()
		chr_gap = int(sum(chr_size_ls)/(5.0*len(chr_size_ls)))
	
	chr_id_ls = chr_id2size.keys()
	chr_id_ls.sort()
	chr_id_ls = ['0'] + chr_id_ls	#chromosome 0 is a place holder. 
	chr_id2cumu_size = {'0':0}	#chr_id_ls might not be continuous integers. so dictionary is better
	for i in range(1,len(chr_id_ls)):
		chr_id = chr_id_ls[i]
		prev_chr_id = chr_id_ls[i-1]
		chr_id2cumu_size[chr_id] = chr_id2cumu_size[prev_chr_id] + chr_id2size[chr_id] + chr_gap
	sys.stderr.write("Done.\n")
	return chr_id2cumu_size, chr_gap, chr_id_ls

def get_pic_area(pos_ls, min_span=30):
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
	if min_span:
		latitude_margin = (min_span-(max_lat-min_lat))/2.0
		longitude_margin = (min_span-(max_lon-min_lon))/2.0
		if latitude_margin<0:
			latitude_margin=5
		if longitude_margin<0:
			longitude_margin=5
	else:
		latitude_margin = 5
		longitude_margin = 5
	if min_lon-longitude_margin>pic_area[0]:
		pic_area[0] = min_lon-longitude_margin
	if min_lat-latitude_margin>pic_area[1]:
		pic_area[1] = min_lat-latitude_margin
	if max_lon+longitude_margin<pic_area[2]:
		pic_area[2] = max_lon+longitude_margin
	if max_lat+latitude_margin<pic_area[3]:
		pic_area[3] = max_lat+latitude_margin
	return pic_area

def get_ecotypeid2pos(curs, ecotype_table, with_data_affiliated=0):
	"""
	2008-10-07 curs could be elixirdb.metadata.bind
	2007-09-18
	2007-10-09
		add with_data_affiliated
	"""
	import os, sys
	sys.stderr.write("Getting ecotypeid2pos from %s ..."%ecotype_table)
	ecotypeid2pos = {}
	if with_data_affiliated:
		rows = curs.execute("select distinct e.id, e.latitude, e.longitude from %s e, calls c where e.latitude is not null and e.longitude is not null and e.id=c.ecotypeid"%ecotype_table)
	else:
		rows = curs.execute("select id, latitude, longitude from %s where latitude is not null and longitude is not null"%ecotype_table)
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#2008-10-07 curs could be elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	for row in rows:
		if is_elixirdb:
			ecotypeid = row.id
			latitude = row.latitude
			longitude = row.longitude
		else:
			ecotypeid, latitude, longitude = row
		ecotypeid2pos[ecotypeid] = (latitude, longitude)
	sys.stderr.write("Done.\n")
	return ecotypeid2pos

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
		pylab.savefig('%s.eps'%output_fname_prefix, dpi=600)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=600)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=600)
	del fig, m, pylab
	sys.stderr.write("Done.\n")

def get_chr_pos_from_x_axis_pos(x_axis_pos, chr_gap, chr_id2cumu_size, chr_id_ls):
	"""
	2008-02-04
		split out from variation.src.GenomeBrowser.py
	2008-02-03
		get chromosome, position from the x axis position
		chr_id_ls is the sorted version of the keys of chr_id2cumu_size
	"""
	chr_id_chosen = None
	position = -1
	for i in range(1, len(chr_id_ls)):	#the 1st in chr_id_ls is fake chromosome 0
		prev_chr_id = chr_id_ls[i-1]
		chr_id = chr_id_ls[i]
		if chr_id2cumu_size[chr_id]>=x_axis_pos:
			chr_id_chosen = chr_id
			position = x_axis_pos - chr_id2cumu_size[prev_chr_id]
			break
	return chr_id_chosen, position



def get_accession_id2name(curs, accession_table='at.accession'):
	"""
	2008-02-07
	"""
	sys.stderr.write("Getting accession_id2name ...")
	accession_id2name = {}
	curs.execute("select id, name from %s"%(accession_table))
	rows = curs.fetchall()
	for row in rows:
		accession_id, name = row
		accession_id2name[accession_id] = name
	sys.stderr.write("Done.\n")
	return accession_id2name

def map_perlegen_ecotype_name2accession_id(curs, accession_table='at.accession', perlegen_table='chip.snp_combined_may_9_06_no_van'):
	"""
	2008-02-07
		in perlegen table, ecotype is a name. all but one('sha') matches accession.name due to case-insensitive match in mysql.
		In fact, perlegen_table.ecotype is all lower-case. accession.name is capitalized on the first-letter.
	"""
	sys.stderr.write("Getting perlegen_ecotype_name2accession_id ...")
	perlegen_ecotype_name2accession_id = {}
	curs.execute("select distinct a.id, s.ecotype from (select distinct ecotype from %s) as s , %s a where s.ecotype=a.name"%(perlegen_table, accession_table))
	rows = curs.fetchall()
	for row in rows:
		accession_id, ecotype_name = row
		perlegen_ecotype_name2accession_id[ecotype_name] = accession_id
	perlegen_ecotype_name2accession_id['sha'] = 89	#this exception
	sys.stderr.write("Done.\n")
	return perlegen_ecotype_name2accession_id

def map_accession_id2ecotype_id(curs, accession2ecotype_table='at.ecotype2accession'):
	"""
	2008-02-07
	
	"""
	sys.stderr.write("Getting  accession_id2ecotype_id ...")
	curs.execute("select ecotype_id, accession_id from %s"%(accession2ecotype_table))
	rows = curs.fetchall()
	accession_id2ecotype_id = {}
	for row in rows:
		ecotype_id, accession_id = row
		accession_id2ecotype_id[accession_id] = ecotype_id
	sys.stderr.write("Done.\n")
	return accession_id2ecotype_id

def get_ecotypeid2nativename(curs, ecotype_table='ecotype'):
	"""
	2008-10-07 curs could be elixirdb.metadata.bind
	2008-04-04
	
	"""
	sys.stderr.write("Getting ecotypeid2nativename ...")
	ecotypeid2nativename = {}
	rows = curs.execute("select id, nativename from %s"%(ecotype_table))
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#2008-10-07 curs could be elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	for row in rows:
		if is_elixirdb:
			ecotypeid = row.id
			nativename = row.nativename
		else:
			ecotypeid, nativename = row
		ecotypeid2nativename[ecotypeid] = nativename
	sys.stderr.write("Done.\n")
	return ecotypeid2nativename
