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

#2008-05-06 ab2number and number2ab is for 384-illumina data
ab2number = {'N': 0,
	'NA': 0,
	'A': 1,
	'B': 2,
	'H': 3}

number2ab = {0: 'NA',
	1: 'A',
	2: 'B',
	3: 'H'}

#2008-05-12	a common NA set
from sets import Set
NA_set = Set([0, 'NA', 'N', -2, '|'])

nt2number = {'|': -2,	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	'-': -1,	#deletion
	'N': 0,
	'NA': 0,
	None: 0,
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

number2nt = {-2: '|',	#2008-01-07 not even tried. 'N'/'NA' is tried but produces inconclusive result.
	-1: '-',
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

number2color = {-1:(0,0,0), 0:(255,255,255), 1:(0,0,255), 2:(0,255,0), 3:(255,0,0), 4:(122,0,122), 5:(122,122,0), 6:(122,255,255), 7:(122,122,122), 8:(255,122,0), 9:(255,255,122), 10:(122,122,255) }

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

def get_nt_number2diff_matrix_index(number2nt):
	"""
	2008-01-01 copied from CmpAccession2Ecotype.py
	2007-10-31/2008-01-07
		nucleotide number ranges from -2 to 10.
		the diff_matrix_index ranges from 0 to 12.
	"""
	sys.stderr.write("Getting nt_number2diff_matrix_index from nt2number ...")
	nt_number2diff_matrix_index = {}
	number_nt_ls = []
	for number, nt in number2nt.iteritems():
		number_nt_ls.append([number,nt])
	number_nt_ls.sort()
	for i in range(len(number_nt_ls)):
		nt_number2diff_matrix_index[number_nt_ls[i][0]] = i
	sys.stderr.write("Done.\n")
	return nt_number2diff_matrix_index

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
	2008-04-04
	
	"""
	sys.stderr.write("Getting ecotypeid2nativename ...")
	ecotypeid2nativename = {}
	curs.execute("select id, nativename from %s"%(ecotype_table))
	rows = curs.fetchall()
	for row in rows:
		ecotypeid, nativename = row
		ecotypeid2nativename[ecotypeid] = nativename
	sys.stderr.write("Done.\n")
	return ecotypeid2nativename

def RawSnpsData_ls2SNPData(rawSnpsData_ls, report=0, use_nt2number=0):
	"""
	2008-05-11
		adapts RawSnpsData(bjarni's SNP data structure) to SNPData
		
		combine all chromsomes together
	"""
	import sys
	if report:
		sys.stderr.write("Converting RawSnpsData_ls to SNPData ...")
	from QC_250k import SNPData
	nt_dict_map = lambda x: nt2number[x]
	snpData = SNPData(row_id_ls = [], col_id_ls=[], data_matrix=[])
	for i in range(len(rawSnpsData_ls)):
		rawSnpsData = rawSnpsData_ls[i]
		chromosome = rawSnpsData.chromosome
		for j in range(len(rawSnpsData.positions)):
			if use_nt2number:
				data_row = map(nt_dict_map, rawSnpsData.snps[j])
			else:
				data_row = rawSnpsData.snps[j]
			snpData.data_matrix.append(data_row)
			snpData.row_id_ls.append((chromosome, rawSnpsData.positions[j]))
		if i==0:	#only need once
			for j in range(len(rawSnpsData.accessions)):
				if rawSnpsData.arrayIds:
					col_id = (rawSnpsData.arrayIds[j], rawSnpsData.accessions[j])
				else:
					col_id = rawSnpsData.accessions[j]
				snpData.col_id_ls.append(col_id)
	if report:
		sys.stderr.write("Done.\n")
	return snpData

def transposeSNPData(snpData, report=0):
	"""
	2008-05-18
		no more copy.deepcopy(snpData), data_matrix takes too long and too much memory
	05/12/08 fix a bug (return snpData)
	2008-05-11
	"""
	import sys
	if report:
		sys.stderr.write("Transposing SNPData ...")
	from pymodule import importNumericArray, SNPData
	num = importNumericArray()
	#copy except data_matrix
	import copy
	newSnpData = SNPData()
	"""
	for option_tuple in SNPData.option_default_dict:
		var_name = option_tuple[0]
		if var_name!='data_matrix':
			setattr(newSnpData, var_name, copy.deepcopy(getattr(snpData, var_name)))
	"""
	newSnpData.row_id_ls = copy.deepcopy(snpData.col_id_ls)
	newSnpData.col_id_ls = copy.deepcopy(snpData.row_id_ls)
	if isinstance(snpData.data_matrix, list):
		newSnpData.data_matrix = num.transpose(num.array(snpData.data_matrix))
	else:	#assume it's array type already. Numeric/numarray has ArrayType, but numpy doesn't
		newSnpData.data_matrix = num.transpose(snpData.data_matrix)
	if report:
		sys.stderr.write("Done.\n")
	return newSnpData

def SNPData2RawSnpsData_ls(snpData, use_number2nt=1, need_transposeSNPData=0, report=0, mask_untouched_deleltion_as_NA=1):
	"""
	2008-05-18
		add mask_untouched_deleltion_as_NA. turned on by default because bjarni's RawSnpsData structure only recognizes NA, A, T, C, G
		if col_id in newSnpData.col_id_ls is tuple of size >1, the 2nd one  in the tuple is taken as array id.
	2008-05-12
		reverse of RawSnpsData_ls2SNPData
		
		adapts SNPData (Yu's SNP data structure) to RawSnpsData(bjarni's SNP data structure)
		
		split into different chromsomes
	"""
	import sys
	if report:
		sys.stderr.write("Converting SNPData to RawSnpsData_ls ...")
	from snpsdata import RawSnpsData
	
	if need_transposeSNPData:
		newSnpData = transposeSNPData(snpData, report=report)
	else:
		newSnpData = snpData
	
	accessions = []
	arrayIds = []
	accession_id = None	#default
	array_id = None	#default
	for col_id in newSnpData.col_id_ls:
		if isinstance(col_id, tuple):
			if len(col_id)>0:
				accession_id = col_id[0]
			if len(col_id)>1:
				array_id = col_id[1]
		else:
			accession_id = col_id
		accessions.append(accession_id)
		if array_id is not None:
			arrayIds.append(array_id)
	
	if mask_untouched_deleltion_as_NA:
		number2nt[-2] = 'NA'	#mask -2 (untouched) as 'NA'
		number2nt[-1] = 'NA'	#mask -1 (deletion) as 'NA'
	number2nt_dict_map = lambda x: number2nt[x]
	rawSnpsData_ls = []
	
	rawSnpsData = RawSnpsData(accessions=accessions, arrayIds=arrayIds)
	rawSnpsData.snps = []
	rawSnpsData.positions = []
	row_id0 = newSnpData.row_id_ls[0]
	old_chromosome = row_id0.split('_')[:1]
	rawSnpsData.chromosome = old_chromosome
	rawSnpsData_ls.append(rawSnpsData)
	rawSnpsData_ls_index = 0
	for i in range(len(newSnpData.row_id_ls)):
		row_id = newSnpData.row_id_ls[i]
		new_chromosome, position = row_id.split('_')[:2]
		position = int(position)
		if new_chromosome!=old_chromosome:
			rawSnpsData = RawSnpsData(accessions=accessions, arrayIds=arrayIds)
			rawSnpsData.snps = []
			rawSnpsData.positions = []
			rawSnpsData_ls.append(rawSnpsData)
			rawSnpsData_ls_index += 1
			rawSnpsData.chromosome = new_chromosome
			old_chromosome = new_chromosome
		rawSnpsData_ls[rawSnpsData_ls_index].positions.append(position)
		if use_number2nt:
			data_row = map(number2nt_dict_map, newSnpData.data_matrix[i])
		else:
			data_row = newSnpData.data_matrix[i]
		rawSnpsData_ls[rawSnpsData_ls_index].snps.append(data_row)
	if report:
		sys.stderr.write("Done.\n")
	return rawSnpsData_ls