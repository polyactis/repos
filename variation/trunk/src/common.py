"""
2007-02-19
	common stuff for variation/src
"""


"""
2007-03-29
	add mappings for '-', 'N', and ambiguous letters ('M', 'R', 'W', 'S', 'Y', 'K')
"""
nt2number = {'-': -1,	#deletion
	'N': 0,
	'NA': 0,
	'A': 1,
	'C': 2,
	'G': 3,
	'T': 4,
	'AC':5,
	'CA':5,
	'M':5,
	'AG':6,
	'GA':6,
	'R':6,
	'AT':7,
	'TA':7,
	'W':7,
	'CG':8,
	'GC':8,
	'S':8,
	'CT':9,
	'TC':9,
	'Y':9,
	'GT':10,
	'TG':10,
	'K':10
	}

number2nt = {-1: '-',
	0: 'NA',
	1:'A',
	2:'C',
	3:'G',
	4:'T',
	5:'AC',
	6:'AG',
	7:'AT',
	8:'CG',
	9:'CT',
	10:'GT'
	}

#2007-04-16 entry[i,j] means whether nucleotide i and j matches. 0(NA) matches everything. singleton(1-4) matches itself and the doublet containing it. doublet(5-10) matches only itself.
nt_number_matching_matrix = [[1, 1,1,1,1,1, 1,1,1,1,1],
	[1, 1,0,0,0,1, 1,1,0,0,0],
	[1, 0,1,0,0,1, 0,0,1,1,0],
	[1, 0,0,1,0,0, 1,0,1,0,1],
	[1, 0,0,0,1,0, 0,1,0,1,1],
	[1, 1,1,0,0,1, 0,0,0,0,0],
	[1, 1,0,1,0,0, 1,0,0,0,0],
	[1, 1,0,0,1,0, 0,1,0,0,0],
	[1, 0,1,1,0,0, 0,0,1,0,0],
	[1, 0,1,0,1,0, 0,0,0,1,0],
	[1, 0,0,1,1,0, 0,0,0,0,1]]


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
	2007-09-18
	2007-10-09
		add with_data_affiliated
	"""
	import os, sys
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