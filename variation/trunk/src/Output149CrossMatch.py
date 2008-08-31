#!/usr/bin/env python
"""

Examples:
	Output149CrossMatch.py -m 3 -x /tmp/149CrossMatch_m3.png -u yh -o /tmp/149CrossMatch_m3.tsv -r -s 5
	
	#QC_method_id=4, self-cross-match results are too huge and stored in a file. Image rendition is relegated to pymodule/DrawMatrix.py.
	Output149CrossMatch.py -m 4 -i ~/panfs/149CrossMatch/149_cross_match.tsv -o /tmp/149CrossMatch_m4_a0.2.tsv -s 5 -r -a 0.2
	
Description:
	Output 149 QCCrossMatch results (output of QC_149_cross_match.py/MpiQC149CrossMatch.py) into matrix (and draw the matrix).
	
	All strains are grouped by country. Countries are in the order of longitude, latitude.
	Within each country, strains are in the order of longitude, latitude.
"""
import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import getopt, csv, math
import Numeric, cPickle
from pymodule import PassingData, importNumericArray, write_data_matrix, SNPData
from variation.src import StockDB
from sets import Set
from pymodule.DrawMatrix import drawMatrix, drawLegend, drawContinousLegend, get_font, combineTwoImages, Value2Color

num = importNumericArray()

class Output149CrossMatch(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('QC_method_id', 1, int): [None, 'm', 1, 'id in table QC_method'],\
							("input_fname", 0, ): [None, 'i', 1, 'get the mismatch rate from this file, rather than the db'],\
							("max_mismatch_rate", 0, float): [1, 'a', 1, 'to filter out any data above the maximum mismatch rate'],\
							('min_no_of_non_NAs', 0, int):[20, 'n', 1, 'to filter out data below the minimum no_of_non_NA_pairs'],\
							('font_path', 1, ):['/usr/share/fonts/truetype/freefont/FreeSerif.ttf', 'e', 1, 'path of the font used to draw labels'],\
							('font_size', 1, int):[20, 's', 1, 'size of font, which determines the size of the whole figure.'],\
							("output_fname", 0, ): [None, 'o', 1, 'Filename to store data matrix'],\
							("fig_fname", 0, ): [None, 'x', 1, 'File name for the figure'],\
							("no_of_ticks", 1, int): [5, 't', 1, 'Number of ticks on the legend'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation. this commit happens after every db operation, not wait till the end.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-29
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getStrainidTargetidFromFile(self, db, QC_method_id, input_fname, max_mismatch_rate, min_no_of_non_NAs=20):
		"""
		2008-08-29
			to get strain id and target id set from the qc_cross_match result file.
		"""
		sys.stderr.write("Getting set of strain_id & target_id ... \n")
		i = 0
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()
		strain_id_set = Set()
		target_id_set = Set()
		for row in reader:
			id, strainid, target_id, qc_method_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, readme_id =row
			strainid = int(strainid)
			target_id = int(target_id)
			qc_method_id = int(qc_method_id)
			mismatch_rate = float(mismatch_rate)
			no_of_mismatches = int(no_of_mismatches)
			no_of_non_NA_pairs = int(no_of_non_NA_pairs)
			if qc_method_id == QC_method_id and no_of_non_NA_pairs>=min_no_of_non_NAs and mismatch_rate<=max_mismatch_rate:
				if QC_method_id==4:	#strain_id_set = target_id_set
					strain_id_set.add(strainid)
					strain_id_set.add(target_id)
				else:
					strain_id_set.add(strainid)
					target_id_set.add(target_id)
			i +=1
			if self.report and i%100000==0:
				sys.stderr.write("%s\t%s"%('\x08'*40, i))
			if self.debug and i>1000000:
				break
		if self.report:
			sys.stderr.write("%s\t%s\n"%('\x08'*40, i))
		return_data = PassingData()
		return_data.strain_id_set = strain_id_set
		return_data.target_id_set = target_id_set
		del reader
		sys.stderr.write("%s strainids and %s target_ids. Done.\n"%(len(strain_id_set), len(target_id_set)))
		return return_data
	
	def getStrainIDInfo(self, db, strain_id_info_query, strain_id_set=None):
		"""
		2008-08-29
		"""
		sys.stderr.write("Getting strain id info ...")
		rows = db.metadata.bind.execute(strain_id_info_query)
		strain_id_ls = []
		strain_id2index = {}
		strain_label_ls = []
		prev_country_abbr = None
		no_of_separators = 0
		for row in rows:
			if strain_id_set and row.strainid not in strain_id_set:	#skip
				continue
			if prev_country_abbr == None:
				prev_country_abbr = row.abbr
			elif row.abbr!=prev_country_abbr:
				prev_country_abbr = row.abbr
				no_of_separators += 1
				strain_id2index[-no_of_separators] = len(strain_id_ls)
				strain_id_ls.append(-no_of_separators)
				strain_label_ls.append('')
			strain_id2index[row.strainid] = len(strain_id_ls)
			strain_id_ls.append(row.strainid)
			if len(row.sitename)>10:
				sitename = row.sitename[:10]
			else:
				sitename = row.sitename
			strain_label_ls.append('%s_%s_%s_%s'%(row.abbr, sitename, row.nativename, row.strainid))
		strain_id_info = PassingData()
		strain_id_info.strain_id_ls = strain_id_ls
		strain_id_info.strain_id2index = strain_id2index
		strain_id_info.strain_label_ls = strain_label_ls
		sys.stderr.write("Done.\n")
		return strain_id_info
	
	def get_data_matrix(self, db, strain_id_info, target_id_info,  QC_method_id, max_mismatch_rate, min_no_of_non_NAs=20):
		"""
		2008-08-29
		"""
		sys.stderr.write("Getting data matrix ... \n")
		data_matrix = num.zeros([len(strain_id_info.strain_id_ls), len(target_id_info.strain_id_ls)], num.float)
		data_matrix[:] = -1
		i = 0
		block_size = 10000
		query = StockDB.QCCrossMatch.query.filter_by(qc_method_id=QC_method_id).filter(StockDB.QCCrossMatch.no_of_non_NA_pairs>min_no_of_non_NAs).filter(StockDB.QCCrossMatch.mismatch_rate<=max_mismatch_rate)
		rows = query.offset(i).limit(block_size)
		min_value = None
		max_value = None
		while rows.count()!=0:
			for row in rows:
				row_index = strain_id_info.strain_id2index[row.strainid]
				col_index = target_id_info.strain_id2index[row.target_id]
				data_value = row.mismatch_rate
				if data_value>=0:
					if min_value==None:
						min_value = data_value
					elif data_value<min_value:
						min_value = data_value
				
				if max_value==None:
					max_value=data_value
				elif data_value>max_value:
					max_value =data_value
				data_matrix[row_index, col_index] = data_value
				i += 1
			if self.report:
				sys.stderr.write("%s\t%s"%('\x08'*40, i))
			rows = query.offset(i).limit(block_size)
		sys.stderr.write("Done.\n")
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.min_value = min_value
		return_data.max_value = max_value
		return return_data
	
	def get_data_matrixFromFile(self, db, strain_id_info, target_id_info,  QC_method_id, input_fname, max_mismatch_rate, min_no_of_non_NAs=20):
		sys.stderr.write("Getting data matrix from  %s ... \n"%input_fname)
		data_matrix = num.zeros([len(strain_id_info.strain_id_ls), len(target_id_info.strain_id_ls)], num.float)
		data_matrix[:] = -1
		reader = csv.reader(open(input_fname), delimiter='\t')
		reader.next()
		min_value = None
		max_value = None
		i = 0
		for row in reader:
			id, strainid, target_id, qc_method_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, readme_id =row
			strainid = int(strainid)
			target_id = int(target_id)
			qc_method_id = int(qc_method_id)
			mismatch_rate = float(mismatch_rate)
			no_of_mismatches = int(no_of_mismatches)
			no_of_non_NA_pairs = int(no_of_non_NA_pairs)
			if qc_method_id == QC_method_id and no_of_non_NA_pairs>=min_no_of_non_NAs and mismatch_rate<=max_mismatch_rate:
				row_index = strain_id_info.strain_id2index.get(strainid)
				col_index = target_id_info.strain_id2index.get(target_id)
				if row_index is None or col_index is None:
					continue
				data_value = mismatch_rate
				if data_value>=0:
					if min_value==None:
						min_value = data_value
					elif data_value<min_value:
						min_value = data_value
				
				if max_value==None:
					max_value=data_value
				elif data_value>max_value:
					max_value =data_value
				data_matrix[row_index, col_index] = data_value
				if QC_method_id==4 and row_index!=col_index:	#149 self-cross-match
					data_matrix[col_index, row_index] = data_value
			i += 1
			if self.report and i%100000==0:
				sys.stderr.write("%s\t%s"%('\x08'*40, i))
			if self.debug and i>1000000:
				break
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.min_value = min_value
		return_data.max_value = max_value
		del reader
		sys.stderr.write("Done.\n")
		return return_data
	
	def markDataMatrixBoundary(self, data_matrix, strain_id_info, target_id_info):
		"""
		2008-08-29
			all those separators have id <0 and mark them with a value different from NA_value.
		"""
		sys.stderr.write("Marking data matrix boundaries ...")
		for strain_id in strain_id_info.strain_id2index:
			if strain_id<0:
				row_index = strain_id_info.strain_id2index[strain_id]
				data_matrix[row_index,:] = -3
		for target_id in target_id_info.strain_id2index:
			if target_id<0:
				col_index = target_id_info.strain_id2index[target_id]
				data_matrix[:, col_index] = -3
		sys.stderr.write("Done.\n")
		return data_matrix
		
	def run(self):	
		if self.debug:
			import pdb
			pdb.set_trace()
		db = StockDB.StockDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		order_by_sentence = " order by c.longitude, c.latitude, e.longitude, e.latitude, e.nativename "	#how to order strains.
		if self.input_fname and self.QC_method_id==4:
			id_set_data = self.getStrainidTargetidFromFile(db, self.QC_method_id, self.input_fname, self.max_mismatch_rate, self.min_no_of_non_NAs)
			sql_table_str = "from %s e, %s s, %s a, %s c"%(StockDB.Ecotype.table.name, StockDB.Site.table.name, StockDB.Address.table.name,\
								StockDB.Country.table.name)
			common_where_condition = "where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id %s " + order_by_sentence
			
			strain_where_condition = common_where_condition%(" and e.id=st.ecotypeid")
			strain_id_info_query = "select distinct st.id as strainid, e.id as ecotypeid, e.nativename, s.name as sitename, c.abbr %s, %s st %s"%(sql_table_str, StockDB.Strain.table.name, strain_where_condition)
		else:
			id_set_data = PassingData()
			id_set_data.strain_id_set = None
			id_set_data.target_id_set = None
			
			sql_table_str = "from %s q, %s e, %s s, %s a, %s c"%(StockDB.QCCrossMatch.table.name, StockDB.Ecotype.table.name, StockDB.Site.table.name, StockDB.Address.table.name,\
									StockDB.Country.table.name)
			common_where_condition = "where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id %s"+ " and q.qc_method_id=%s and q.no_of_non_NA_pairs>=%s and q.mismatch_rate<=%s "%\
				(self.QC_method_id, self.min_no_of_non_NAs, self.max_mismatch_rate) + order_by_sentence
			
			strain_where_condition = common_where_condition%(" and e.id=st.ecotypeid and st.id=q.strainid")
			strain_id_info_query = "select distinct q.strainid, e.id as ecotypeid, e.nativename, s.name as sitename, c.abbr %s, %s st %s"%(sql_table_str, StockDB.Strain.table.name, strain_where_condition)
		strain_id_info = self.getStrainIDInfo(db, strain_id_info_query, id_set_data.strain_id_set)
		if self.QC_method_id==4:
			target_id_info = strain_id_info
		else:
			target_where_condition = common_where_condition%(" and e.id=q.target_id")
			target_id_info_query = "select distinct e.id as strainid, e.id as ecotypeid, e.nativename, s.name as sitename, c.abbr %s %s"%(sql_table_str, target_where_condition)
			target_id_info = self.getStrainIDInfo(db, target_id_info_query)
		
		if self.input_fname:
			rdata = self.get_data_matrixFromFile(db, strain_id_info, target_id_info,  self.QC_method_id, self.input_fname, self.max_mismatch_rate, self.min_no_of_non_NAs)
		else:
			rdata = self.get_data_matrix(db, strain_id_info, target_id_info, self.QC_method_id, self.max_mismatch_rate, self.min_no_of_non_NAs)
		
		rdata.data_matrix = self.markDataMatrixBoundary(rdata.data_matrix, strain_id_info, target_id_info)
		
		header = ['strain info', ''] + target_id_info.strain_label_ls
		strain_acc_list = strain_id_info.strain_label_ls
		category_list = [1]*len(strain_acc_list)
		if SNPData.isDataMatrixEmpty(rdata.data_matrix):
			sys.stderr.write("Nothing fetched from database.\n")
			sys.exit(3)
		if self.output_fname:
			write_data_matrix(rdata.data_matrix, self.output_fname, header, strain_acc_list, category_list)
		
		if self.fig_fname:
			font = get_font(self.font_path, font_size=self.font_size)	#2008-08-01
			value2color_func = lambda x: Value2Color.value2HSLcolor(x, rdata.min_value, rdata.max_value)
			im_legend = drawContinousLegend(rdata.min_value, rdata.max_value, self.no_of_ticks, value2color_func, font)
			#im.save('%s_legend.png'%self.fig_fname_prefix)
			im = drawMatrix(rdata.data_matrix, value2color_func, strain_id_info.strain_label_ls,\
						target_id_info.strain_label_ls, with_grid=1, font=font)
			im = combineTwoImages(im, im_legend, font=font)
			im.save(self.fig_fname)
		
if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = Output149CrossMatch
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()