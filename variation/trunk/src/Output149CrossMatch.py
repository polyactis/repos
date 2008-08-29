#!/usr/bin/env python
"""

Examples:
	Output149CrossMatch.py -m 3 -x /tmp/149CrossMatch_m3.png -u yh -o /tmp/149CrossMatch_m3.tsv -r -s 5
	
	Output149CrossMatch.py -m 4 -o /tmp/149CrossMatch_m4.tsv -x /tmp/149CrossMatch_m4.png -s 5 -r
	
Description:
	Output 149 QCCrossMatch results (output of QC_149_cross_match.py/MpiQC149CrossMatch.py) into matrix and/or draw the matrix.
	
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
	
	def getStrainIDInfo(self, db, strain_id_info_query):
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
			strain_label_ls.append('%s_%s_%s_%s'%(row.sitename, row.abbr, row.nativename, row.strainid))
		strain_id_info = PassingData()
		strain_id_info.strain_id_ls = strain_id_ls
		strain_id_info.strain_id2index = strain_id2index
		strain_id_info.strain_label_ls = strain_label_ls
		sys.stderr.write("Done.\n")
		return strain_id_info
	
	def get_data_matrix(self, db, strain_id_info, target_id_info,  QC_method_id, min_no_of_non_NAs=20):
		"""
		2008-08-29
		"""
		sys.stderr.write("Getting data matrix ... \n")
		data_matrix = num.zeros([len(strain_id_info.strain_id_ls), len(target_id_info.strain_id_ls)], num.float)
		data_matrix[:] = -1
		i = 0
		block_size = 10000
		query = StockDB.QCCrossMatch.query.filter_by(qc_method_id=QC_method_id).filter(StockDB.QCCrossMatch.no_of_non_NA_pairs>min_no_of_non_NAs)
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
	
	def get_data_matrixFromFile(self, db, strain_id_info, target_id_info,  QC_method_id, input_fname, min_no_of_non_NAs=20):
		sys.stderr.write("Getting data matrix from  ... \n"%input_fname)
		data_matrix = num.zeros([len(strain_id_info.strain_id_ls), len(target_id_info.strain_id_ls)], num.float)
		data_matrix[:] = -1
		reader = csv.reader(open(input_fname), delimiter='\t')
		min_value = None
		max_value = None
		for row in reader:
			id, strainid, target_id, qc_method_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, readme_id =row
			strainid = int(strainid)
			target_id = int(target_id)
			qc_method_id = int(qc_method_id)
			mismatch_rate = float(mismatch_rate)
			no_of_mismatches = int(no_of_mismatches)
			no_of_non_NA_pairs = int(no_of_non_NA_pairs)
			if qc_method_id != QC_method_id:
				continue
			if no_of_non_NA_pairs<min_no_of_non_NAs:
				continue
			row_index = strain_id_info.strain_id2index[strainid]
			col_index = target_id_info.strain_id2index[target_id]
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
			i += 1
			if self.report and i%10000==0:
				sys.stderr.write("%s\t%s"%('\x08'*40, i))
		sys.stderr.write("Done.\n")
		return_data = PassingData()
		return_data.data_matrix = data_matrix
		return_data.min_value = min_value
		return_data.max_value = max_value
		return return_data
	
	def run(self):	
		if self.debug:
			import pdb
			pdb.set_trace()
		db = StockDB.StockDB(drivername=self.drivername, username=self.db_user,
						password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		
		sql_table_str = "from %s q, %s e, %s s, %s a, %s c"%(StockDB.QCCrossMatch.table.name, StockDB.Ecotype.table.name, StockDB.Site.table.name, StockDB.Address.table.name,\
								StockDB.Country.table.name)
		common_where_condition = "where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id %s"+ " and q.qc_method_id=%s order by c.abbr, s.name, e.nativename"%self.QC_method_id
		
		strain_where_condition = common_where_condition%(" and e.id=st.ecotypeid and st.id=q.strainid")
		strain_id_info_query = "select distinct q.strainid, e.id as ecotypeid, e.nativename, s.name as sitename, c.abbr %s, %s st %s"%(sql_table_str, StockDB.Strain.table.name, strain_where_condition)
		strain_id_info = self.getStrainIDInfo(db, strain_id_info_query)
		if self.QC_method_id==4:
			target_id_info = strain_id_info
		else:
			target_where_condition = common_where_condition%(" and e.id=q.target_id")
			target_id_info_query = "select distinct e.id as strainid, e.id as ecotypeid, e.nativename, s.name as sitename, c.abbr %s %s"%(sql_table_str, target_where_condition)
			target_id_info = self.getStrainIDInfo(db, target_id_info_query)
		
		if self.input_fname and os.path.isfile(self.input_fname):
			rdata = self.get_data_matrixFromFile(db, strain_id_info, target_id_info,  self.QC_method_id, self.input_fname)
		else:
			rdata = self.get_data_matrix(db, strain_id_info, target_id_info, self.QC_method_id)
		
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
			value2color_func = lambda x: Value2Color.value2HSLcolor(x, rdata.min_value, rdata.max_value, NA_value=-1, super_value=-2)
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