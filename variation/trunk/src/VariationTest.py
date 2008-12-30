#!/usr/bin/env python
""""
Usage: VariationTest.py -y TestCaseType [OPTIONS]

Option:
	-y ..., --type=...	which test case should be invoked.
	-h, --help              show this help

Examples:
	VariationTest.py -y 2

2007-03-08 1: TestTrioInference
2007-04-17 2: Test_find_smallest_vertex_set_to_remove_all_edges
2008-10-08 3: TestFetchSNPRegionPlot
2008-10-08 4: TestGetEcotypeInfo
2008-10-11 5: TestDrawMap
2008-10-17 6: TestScoreRankHistogram
2008-12-30 7: TestMAFVsScorePlot
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import unittest, os, sys, getopt, csv

class TestTrioInference(unittest.TestCase):
	"""
	2007-03-08
	"""
	def test_trio_ancestry_inference_by_DP(self):
		"""
		2007-04-16
			add a heterozygous call to test
		"""
		from MpiTrioAncestryInference import MpiTrioAncestryInference
		import Numeric
		data_matrix = Numeric.array([   [0,1,3,1,2,2,2,1,1,1],
						[1,2,4,3,2,2,1,4,2,4],
						[2,1,4,3,2,2,1,1,5,4]])
		chr_start_ls = [0,4,10]
		trio_arrangement_ls = [[0,1,2], [1,2,0], [2,0,1]]
		MpiTrioAncestryInference_instance = MpiTrioAncestryInference(debug=1)
		for trio_arrangement in trio_arrangement_ls:
			ancestry_ls, no_of_jumps = MpiTrioAncestryInference_instance.identify_ancestry_with_min_jumps(data_matrix[trio_arrangement[0]], data_matrix[trio_arrangement[1]], data_matrix[trio_arrangement[2]], chr_start_ls)
			print trio_arrangement
			print ancestry_ls
			print no_of_jumps


class Test_find_smallest_vertex_set_to_remove_all_edges(unittest.TestCase):
	"""
	2007-04-17
	"""
	def setUp(self):
		print
	
	def test_find_smallest_vertex_set_to_remove_all_edges(self):
		from dbSNP2data import dbSNP2data
		identity_pair_ls = [[1,2],[2,3],[2,4],[4,5]]
		import networkx as nx
		g = nx.Graph()
		g.add_edges_from(identity_pair_ls)
		from dbSNP2data import dbSNP2data
		dbSNP2data_instance = dbSNP2data()
		#import pdb
		#pdb.set_trace()
		vertex_list_to_be_deleted = dbSNP2data_instance.find_smallest_vertex_set_to_remove_all_edges(g)
		print 'graph'
		print identity_pair_ls
		print 'vertex_list_to_be_deleted'
		print vertex_list_to_be_deleted

class TestFetchSNPRegionPlot(unittest.TestCase):
	"""
	2008-10-08
		to test whether able to restore binary data encoded in text database type back into binary file.
	"""
	def setUp(self):
		print
	
	def test_fetchOneImageOut(self):
		import Stock_250kDB
		hostname='papaya.usc.edu'
		dbname='stock_250k'
		db_user='yh'
		db_passwd = ''
		drivername='mysql'
		schema = None
		db = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
						password=db_passwd, hostname=hostname, database=dbname, schema=schema)
		db.setup(create_tables=False)
		
		"""
		#2008-10-19	convert 64bit-encoded img data to direct binary data and store in a binary field(snp_region_plot.png_data
		#find out Binary db type can be handled correctly (give a length argument to Binary to avoid data truncation).
		import base64
		counter = 0
		print
		for snp_region_plot in Stock_250kDB.SNPRegionPlot.query.all():
			snp_region_plot.png_data = base64.b64decode(snp_region_plot.old_img_data)
			db.session.save_or_update(snp_region_plot)
			db.session.flush()
			counter += 1
			sys.stderr.write("%s\t%s"%('\x08'*80, counter))
		"""
		
		snp_region_plot = Stock_250kDB.SNPRegionPlot.get(1)
		outf = open('/tmp/snp_region_plot_1.png', 'wb')
		snp_region_plot
		outf.write(snp_region_plot.png_data)
		outf.close()

class TestGetEcotypeInfo(unittest.TestCase):
	"""
	2008-10-08
		to test what are the properties of (rows = db.metadata.bind.execute())
	"""
	def setUp(self):
		print
	
	def test_getEcotypeInfo(self):
		from common import getEcotypeInfo
		import StockDB, Stock_250kDB	#StockDB has to be setup otherwise, StockDB.Ecotype.table is None in getEcotypeInfo()
		hostname='papaya.usc.edu'
		dbname='stock_250k'
		db_user='yh'
		db_passwd = ''
		drivername='mysql'
		schema = None
		db = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
						password=db_passwd, hostname=hostname, database=dbname, schema=schema)
		#doesn't matter which database to connect as far as StockDB is imported
		#db = StockDB.StockDB(drivername=drivername, username=db_user,
		#				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
		db.setup(create_tables=False)
		import pdb
		pdb.set_trace()
		getEcotypeInfo(db)

class TestDrawMap(unittest.TestCase):
	"""
	2008-10-14
		test to connect the matplotlib's map with another plot
	"""
	def setUp(self):
		print
	
	def test_drawMap(self):
		sys.stderr.write("Drawing graph on a map ...\n")
		#import pdb
		#pdb.set_trace()
		import StockDB, Stock_250kDB	#StockDB has to be setup otherwise, StockDB.Ecotype.table is None in getEcotypeInfo()
		hostname='papaya.usc.edu'
		dbname='stock_250k'
		db_user='yh'
		db_passwd = ''
		drivername='mysql'
		schema = None
		db = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
						password=db_passwd, hostname=hostname, database=dbname, schema=schema)
		#doesn't matter which database to connect as far as StockDB is imported
		#db = StockDB.StockDB(drivername=drivername, username=db_user,
		#				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
		db.setup(create_tables=False)
		from common import getEcotypeInfo
		ecotype_info = getEcotypeInfo(db)
			
		from matplotlib.toolkits.basemap import Basemap
		#from mpl_toolkits.basemap import Basemap
		import pylab
		from matplotlib import rcParams
		rcParams['font.size'] = 6
		rcParams['legend.fontsize'] = 6
		#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
		rcParams['axes.labelsize'] = 4
		rcParams['axes.titlesize'] = 8
		rcParams['xtick.labelsize'] = 4
		rcParams['ytick.labelsize'] = 4
		
		pylab.clf()
		
		#fig = pylab.figure()
		#fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
		axe_map = pylab.axes([0.5, 0.02, 0.4, 0.8], frameon=False)
		axe_map.set_title("Global Arabidopsis Ecotype Distribution")
		axe_map.set_xlabel('ecotype is drawn as circles.')
		axe_map.set_ylabel('strains')
		pic_area=[-140,-40,140,70]
		m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
		resolution='l',projection='mill', ax=axe_map)
		"""
		llcrnrx = -self.rmajor
		llcrnry = -self.rmajor
		urcrnrx = -llcrnrx
		urcrnry = -llcrnry
		"""
		#m.drawcoastlines()
		#m.bluemarble()
		m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1], size=4, linewidth=0.1)
		m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1], size=4, linewidth=0.1)
		m.fillcontinents()
		m.drawcountries(linewidth=0.1)
		#m.drawstates()
		#m.drawlsmask((0,255,0,255), (0,0,255,255), lakes=True)
		#m.drawlsmask('coral','aqua',lakes=True)

		print "xlim:", axe_map.get_xlim()
		print "ylim:", axe_map.get_ylim()
		
		xlim = axe_map.get_xlim()
		ylim = axe_map.get_ylim()
		
		"""
		for strain_id in StrainID2PCAPosInfo.strain_id_ls:
			img_y_pos = StrainID2PCAPosInfo.strain_id2img_y_pos[strain_id]
			phenotype_row_index = phenData.row_id2row_index[strain_id]
			phenotype = phenData.data_matrix[phenotype_row_index][phenotype_col_index]
			ecotype_id = int(strain_id)
			ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
			if ecotype_obj:
				lat, lon = ecotype_obj.latitude, ecotype_obj.longitude
			else:
				sys.stderr.write("Warning: Ecotype %s not in ecotype_info (fetched from stock db).\n"%ecotype_id)
				continue
			if lat and lon:
				x, y = m(lon, lat)
				color = cmap(norm(phenotype))
				ax.plot([0, x], [img_y_pos, y], linestyle='--', alpha=0.2, linewidth=0.2)
				ax.scatter([x],[y], s=10, linewidth=0, facecolor=color)	#, zorder=10)
		"""
		
		#pylab.title("Global Arabidopsis Ecotype Distribution")
		output_fname_prefix = '/tmp/map'
		
		print "ylim:", axe_map.get_ylim()
		axe_map.set_xlim(xlim)
		axe_map.set_ylim(ylim)
		
		
		axe_chromosome = pylab.axes([0.05, 0.02, 0.8, 0.8], frameon=False)
		axe_chromosome.set_title("chromosome")
		
		#fix the two transformations before doing cross-axe drawings
		axe_map.transData.freeze()  # eval the lazy objects
		axe_map.transAxes.freeze()
		axe_chromosome.transData.freeze()  # eval the lazy objects
		axe_chromosome.transAxes.freeze()
		no_of_ecotypes = 200
		ecotype_id_ls = ecotype_info.ecotype_id2ecotype_obj.keys()[:no_of_ecotypes]
		no_of_ecotypes_drawn = 0
		for i in range(no_of_ecotypes):
			ecotype_id = ecotype_id_ls[i]
			y_pos = i/float(no_of_ecotypes)
			#y_pos = i/float(no_of_ecotypes)*ylim[1]
			ecotype_obj = ecotype_info.ecotype_id2ecotype_obj.get(ecotype_id)
			if ecotype_obj:
				lat, lon = ecotype_obj.latitude, ecotype_obj.longitude
			else:
				sys.stderr.write("Warning: Ecotype %s not in ecotype_info (fetched from stock db).\n"%ecotype_id)
				continue
			if lat and lon:
				x, y = m(lon, lat)
				#axe_map.plot([0, x], [y_pos, y], linestyle='--', alpha=0.2, linewidth=0.2)
				axe_map.set_xlim(xlim)
				axe_map.set_ylim(ylim)
				axe_map.scatter([x],[y], s=5, linewidth=0, facecolor='r', alpha=0.2, zorder=10)
				canvas_x, canvas_y = axe_map.transData.xy_tup((x,y))
				axe_chromosome_xy = axe_chromosome.transData.inverse_xy_tup((canvas_x,canvas_y))
				axe_chromosome.plot([0,axe_chromosome_xy[0]], [y_pos, axe_chromosome_xy[1]], linestyle='--', alpha=0.2, linewidth=0.2)
				no_of_ecotypes_drawn += 1
		#release two transformations
		axe_map.transData.thaw()  # eval the lazy objects
		axe_map.transAxes.thaw()
		axe_chromosome.transData.thaw()  # eval the lazy objects
		axe_chromosome.transAxes.thaw()
		#set to the same x/y_lim before cross-axe drawing
		axe_map.set_xlim(xlim)
		axe_map.set_ylim(ylim)
		axe_chromosome.set_xlim([0,1])
		axe_chromosome.set_ylim([0,1])
		if output_fname_prefix:
			pylab.savefig('%s.png'%output_fname_prefix, dpi=600)
			pylab.savefig('%s.svg'%output_fname_prefix)
		sys.stderr.write("%s ecotypes drawn. Done.\n"%(no_of_ecotypes_drawn))


class TestScoreRankHistogram(unittest.TestCase):
	"""
	2008-10-17
		test to check-in/restore binary data in ScoreRankHistogram
	"""
	def setUp(self):
		print
	
	def storeImage(self, db, input_fname_prefix, hist_type):
		import Stock_250kDB
		input_fname1 = '%s.png'%input_fname_prefix
		inf1 = open(input_fname1, 'rb')
		input_fname2 = '%s.svg'%input_fname_prefix
		inf2 = open(input_fname2, 'rb')
		
		score_rank_hist = Stock_250kDB.ScoreRankHistogram(phenotype_method_id=1, list_type_id=1)
		score_rank_hist.original_filename = input_fname_prefix
		score_rank_hist.hist_type = hist_type
		score_rank_hist.score_hist = inf1.read()
		score_rank_hist.score_hist_svg = inf2.read()
		del inf1, inf2
		db.session.save(score_rank_hist)
		db.session.flush()
	
	def test_SaveAndFetchOneImage(self):
		import Stock_250kDB
		hostname='papaya.usc.edu'
		dbname='stock_250k'
		db_user='yh'
		db_passwd = ''
		drivername='mysql'
		schema = None
		db = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
						password=db_passwd, hostname=hostname, database=dbname, schema=schema)
		db.setup(create_tables=False)
		call_method_id= 17
		min_distance = 20000
		get_closest = 0
		min_MAF = 0.1
		rows = Stock_250kDB.ScoreRankHistogramType.query.filter_by(call_method_id=call_method_id).filter_by(min_distance=min_distance).\
										filter_by(get_closest =get_closest).\
										filter(Stock_250kDB.ScoreRankHistogramType.min_MAF>=min_MAF-0.0001).filter(Stock_250kDB.ScoreRankHistogramType.min_MAF<=min_MAF+0.0001).\
										filter_by(allow_two_sample_overlapping = 0)
		if rows.count()>0:
			hist_type = rows.first()
		else:
			hist_type = Stock_250kDB.ScoreRankHistogramType(call_method_id=call_method_id, min_distance=min_distance,\
										get_closest =get_closest,
										min_MAF = min_MAF,
										allow_two_sample_overlapping = 0)
		input_fname_prefix = '/mnt/tmp/tmp/hist_of_results_by_gene_candidate_score_rank/Phenotype 1 LD vs 1 FT_short_score'
		self.storeImage(db, input_fname_prefix, hist_type)
		fetch_score_rank_hist = Stock_250kDB.ScoreRankHistogram.query.first()
		
		#import base64
		outf = open('/tmp/plot_1.png', 'wb')
		outf.write(fetch_score_rank_hist.score_hist)
		#outf.write(base64.b64decode(snp_region_plot.img_data))
		outf.close()
		
		outf = open('/tmp/plot_1.svg', 'wb')
		outf.write(fetch_score_rank_hist.score_hist_svg)
		#outf.write(base64.b64decode(snp_region_plot.img_data))
		outf.close()

class TestMAFVsScorePlot(unittest.TestCase):
	"""
	2008-12-30
		test to retrieve binary data from MAFVsScorePlot and store it in a file.
	"""
	def setUp(self):
		print
	
	def test_SaveAndFetchOneImage(self):
		import Stock_250kDB
		hostname='papaya.usc.edu'
		dbname='stock_250k'
		db_user='yh'
		db_passwd = '***'
		drivername='mysql'
		schema = None
		db = Stock_250kDB.Stock_250kDB(drivername=drivername, username=db_user,
						password=db_passwd, hostname=hostname, database=dbname, schema=schema)
		db.setup(create_tables=False)
		
		one_entry = Stock_250kDB.MAFVsScorePlot.query.first()
		
		outf = open('/tmp/plot_1.png', 'wb')
		outf.write(one_entry.png_data)
		outf.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "type="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hy:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	TestCaseDict = {1:TestTrioInference,
		2:Test_find_smallest_vertex_set_to_remove_all_edges,
		3:TestFetchSNPRegionPlot,
		4:TestGetEcotypeInfo,
		5:TestDrawMap,
		6:TestScoreRankHistogram,
		7:TestMAFVsScorePlot}
	type = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-y", "--type"):
			type = int(arg)
			
	if type:
		suite = unittest.TestSuite()
		suite.addTest(unittest.makeSuite(TestCaseDict[type]))
		unittest.TextTestRunner(verbosity=2).run(suite)
		
		"""
		#try to find a fancy to pass options to test class, not yet
		from pymodule import ProcessOptions
		main_class = TestCaseDict[type]
		po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
		instance = main_class(**po.long_option2value)
		"""
	else:
		print __doc__
		sys.exit(2)		