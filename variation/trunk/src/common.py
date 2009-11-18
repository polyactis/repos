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
	resolution='l',projection='mill', ax=pylab.gca())
	
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
	2009-4-10
		curs could be cursor of MySQLdb connection or ElixirDB.metadata.bind
	2008-02-07
		in perlegen table, ecotype is a name. all but one('sha') matches accession.name due to case-insensitive match in mysql.
		In fact, perlegen_table.ecotype is all lower-case. accession.name is capitalized on the first-letter.
	"""
	sys.stderr.write("Getting perlegen_ecotype_name2accession_id ...")
	perlegen_ecotype_name2accession_id = {}
	rows = curs.execute("select distinct a.id, s.ecotype from (select distinct ecotype from %s) as s , %s a where s.ecotype=a.name"%(perlegen_table, accession_table))
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#2008-10-07 curs could be elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	
	for row in rows:
		if is_elixirdb:
			accession_id = row.id
			ecotype_name = row.ecotype
		else:
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

def get_ecotypeid2nativename(curs, ecotype_table='stock.ecotype'):
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

def getNativename2EcotypeIDLs(curs, ecotype_table='stock.ecotype', turnUpperCase=False):
	"""
	2009-7-23
		add argument turnUpperCase: whether to turn nativename into uppercase
	2009-5-16
		reverse version of get_ecotypeid2nativename()
	"""
	sys.stderr.write("Getting nativename2ecotypeid_ls ...")
	dc = {}
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
		if turnUpperCase:
			nativename = nativename.upper()
		if nativename not in dc:
			dc[nativename] = []
		dc[nativename].append(ecotypeid)
	sys.stderr.write("Done.\n")
	return dc

def getNativename2TgEcotypeIDSet(curs, table='stock.ecotypeid2tg_ecotypeid', turnUpperCase=False):
	"""
	2009-7-23
		add argument turnUpperCase: whether to turn nativename into uppercase
	2009-5-28
		this function =  getNativename2EcotypeIDLs() + get_ecotypeid2tg_ecotypeid()
	"""
	sys.stderr.write("Getting nativename2tg_ecotypeid_ls ...")
	dc = {}
	rows = curs.execute("select tg_ecotypeid, nativename from %s"%(table))
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#2008-10-07 curs could be elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	for row in rows:
		if is_elixirdb:
			ecotypeid = row.tg_ecotypeid
			nativename = row.nativename
		else:
			ecotypeid, nativename = row
		if turnUpperCase:
			nativename = nativename.upper()
		if nativename not in dc:
			dc[nativename] = set()
		dc[nativename].add(ecotypeid)
	sys.stderr.write("Done.\n")
	return dc

def getEcotypeInfo(db, country_order_type=1):
	"""
	2009-09-2
		add region into ecotype_obj
	2008-10-08
		use ecotype_id2ecotype_obj to summarize
			ecotypeid2pos
			ecotypeid2nativename
			ecotypeid2country
	2008-10-08
		add option order_by_type
		get country2order
		moved from PlotGroupOfSNPs.py
		the db handle is not restricted to the stock database. could be any database on the same server.
		BUT StockDB has to be imported in the program where db connection is established just so that StockDB.Ecotype.table is setup while StockDB.Ecotype.table.metadata is None.
	2008-10-07
	"""
	sys.stderr.write("Getting  Ecotype info ... ")
	import StockDB
	from pymodule import PassingData
	ecotype_info = PassingData()
	if country_order_type==1:
		order_seq_sentence = 'c.latitude, c.longitude'
	else:
		order_seq_sentence = 'c.longitude, c.latitude'
	rows = db.metadata.bind.execute("select e.id as ecotype_id, e.nativename, e.latitude, e.longitude, a.region, \
		c.abbr as country, c.latitude as country_latitude, \
		c.longitude as country_longitude \
		from stock.%s e, stock.%s s, stock.%s a, stock.%s c where e.siteid=s.id and s.addressid=a.id and \
		a.countryid=c.id order by %s "%(getattr(StockDB.Ecotype.table, 'name', 'ecotype'), \
		getattr(StockDB.Site.table, 'name', 'site'), getattr(StockDB.Address.table, 'name', 'address'), \
		getattr(StockDB.Country.table, 'name', 'country'), order_seq_sentence))
	ecotype_id2ecotype_obj = {}
	country2order = {}
	for row in rows:
		ecotype_obj = PassingData()
		for key,value in row.items():	#not iteritems() for RowProxy object
			setattr(ecotype_obj, key, value)
		ecotype_id2ecotype_obj[row.ecotype_id] = ecotype_obj
		if row.country not in country2order:
			country2order[row.country] = len(country2order)
	ecotype_info.ecotype_id2ecotype_obj = ecotype_id2ecotype_obj
	ecotype_info.country2order = country2order
	sys.stderr.write("Done.\n")
	return ecotype_info

def get_total_gene_ls(curs, gene_table='genome.gene', tax_id=3702, debug=False):
	"""
	2008-10-21
		get a list of gene given tax_id
			no mitochondrial, chromosome has to be known,
	"""
	if debug:
		sys.stderr.write("Getting gene_id_ls ... ")
	rows = curs.execute("select distinct gene_id from %s where tax_id=%s and chromosome is not null and chromosome!='MT'"%\
								(gene_table, tax_id))
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#2008-10-07 this curs is not elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	gene_id_ls = []
	for row in rows:
		if is_elixirdb:
			gene_id = row.gene_id
		else:
			gene_id = row[0]
		gene_id_ls.append(row.gene_id)
	if debug:
		sys.stderr.write("Done.\n")
	return gene_id_ls

def get_phenotype_method_id_lsFromPhenData(phenData):
	"""
	2009-2-16
		phenData is a SNPData-class data structure read in from the file outputted by OutputPhenotype.py.
			each col_id is sth like '1_LD'.
		this function is to extract the integer phenotype_method_id out of col_id, in the same order as col_id_ls.
	"""
	phenotype_method_id_ls = []
	for col_id in phenData.col_id_ls:
		col_id_ls = col_id.split('_')
		phenotype_method_id=int(col_id_ls[0])
		phenotype_method_id_ls.append(phenotype_method_id)
	return phenotype_method_id_ls

def get_ecotypeid2tg_ecotypeid(curs, table='stock.ecotypeid2tg_ecotypeid', debug=False):
	"""
	2009-4-4
		get the mapping between ecotypeid and tg_ecotypeid
	"""
	if debug:
		sys.stderr.write("Getting ecotypeid2tg_ecotypeid ...")
	ecotypeid2tg_ecotypeid = {}
	rows = curs.execute("select * from %s"%(table))
	is_elixirdb = 1
	if hasattr(curs, 'fetchall'):	#this curs is not elixirdb.metadata.bind
		rows = curs.fetchall()
		is_elixirdb = 0
	
	#rows = StockDB.EcotypeIDStrainID2TGEcotypeID.query()
	for row in rows:
		if is_elixirdb:
			ecotypeid = row.ecotypeid
			tg_ecotypeid = row.tg_ecotypeid
		else:
			ecotypeid = row[0]
			tg_ecotypeid = row[6]
		ecotypeid2tg_ecotypeid[ecotypeid] = tg_ecotypeid
	if debug:
		sys.stderr.write("Done.\n")
	return ecotypeid2tg_ecotypeid

def get_ecotype_id_set_250k_in_pipeline(ArrayInfo):
	"""
	2009-7-22
		function copied from helloworld/controllers/Accession.py with ArrayInfo as argument.
	"""
	ecotype_id_set_250k_in_pipeline = set()
	for row in ArrayInfo.query():
		if row.maternal_ecotype_id==row.paternal_ecotype_id:	#no crosses.
			ecotype_id_set_250k_in_pipeline.add(row.maternal_ecotype_id)
	return ecotype_id_set_250k_in_pipeline

def fillInPhenotypeMethodID2ecotype_id_set(PhenotypeAvgTableCalss=None):
	"""
	2009-11-17
		return PhenotypeMethodID2ecotype_id_set
	"""
	sys.stderr.write("Filling up PhenotypeMethodID2ecotype_id_set ...")
	PhenotypeMethodID2ecotype_id_set = {}
	rows = PhenotypeAvgTableCalss.query()
	for row in rows:
		phenotype_method_id = row.method_id
		if phenotype_method_id not in PhenotypeMethodID2ecotype_id_set:
			PhenotypeMethodID2ecotype_id_set[phenotype_method_id] = set()
		PhenotypeMethodID2ecotype_id_set[phenotype_method_id].add(row.ecotype_id)
	sys.stderr.write("Done.\n")
	return PhenotypeMethodID2ecotype_id_set


def fillInCallMethodID2ecotype_id_set(CallInfoTableClass=None):
	"""
	2009-11-17
		return CallMethodID2ecotype_id_set
	"""
	sys.stderr.write("Filling up CallMethodID2ecotype_id_set...")
	CallMethodID2ecotype_id_set = {}
	rows = CallInfoTableClass.query()
	for row in rows:
		call_method_id = row.method_id
		if call_method_id not in CallMethodID2ecotype_id_set:
			CallMethodID2ecotype_id_set[call_method_id] = set()
		CallMethodID2ecotype_id_set[call_method_id].add(row.array.maternal_ecotype_id)
	sys.stderr.write("Done.\n")
	return CallMethodID2ecotype_id_set
	