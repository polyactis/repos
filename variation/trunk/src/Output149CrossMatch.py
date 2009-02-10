#!/usr/bin/env python
"""

Examples:
	Output149CrossMatch.py -m 3 -x /tmp/149CrossMatch_m3.png -u yh -o /tmp/149CrossMatch_m3.tsv -r -s 5
	
	#QC_method_id=4, self-cross-match results are too huge and stored in a file. Image rendition is relegated to pymodule/DrawMatrix.py.
	Output149CrossMatch.py -m 4 -i ~/panfs/149CrossMatch/149_cross_match.tsv -o /tmp/149CrossMatch_m4_a0.2.tsv -s 5 -r -a 0.2
	
	#pretty similar to above, but no drawing. and align strains according to sequenom plates.
	Output149CrossMatch.py -m 4 -i ~/panfs/149CrossMatch/149SNPSequenomBlock_0_cross_match.tsv -o ~/panfs/149CrossMatch/149SNPSequenomBlock_0_cross_match_matrix_a0.3_g.tsv -a 0.3 -g 2
	
	#
	pretty similar to above, no drawing. and align row-strains according to sequenom plates and column-strains by (country, longitude).
	Output149CrossMatch.py -m 4 -i ~/panfs/149CrossMatch/149SNPSequenomBlock_0_cross_match.tsv -o ~/panfs/149CrossMatch/149SNPSequenomBlock_0_cross_match_matrix_a0.3_g.tsv -a 0.3 -g 3
	
Description:
	Output 149 QCCrossMatch results (output of QC_149_cross_match.py/MpiQC149CrossMatch.py) into matrix (and draw the matrix).
	
	By default, all strains are grouped by country. Countries are in the order of longitude, latitude.
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
import StockDB
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
							("how_to_group_strains", 0, ): [1, 'g', 1, '1: this program order/group strains by (country, longitude).\
								2: order/group strains by (sequenom plate, country, longitude).\
								3: order/group rows by (sequenom plate, country, longitude) and columns by (country, longitude), QC_method=4 only.'],\
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
		2008-09-10
			column in input_fname is determined on the fly
		2008-08-29
			to get strain id and target id set from the qc_cross_match result file.
		"""
		sys.stderr.write("Getting set of strain_id & target_id ... \n")
		reader = csv.reader(open(input_fname), delimiter='\t')
		#figure out which variable is in which column
		header = reader.next()
		col_name2index = {}
		for i in range(len(header)):
			column_name = header[i]
			col_name2index[column_name] = i
		
		strain_id_set = Set()
		target_id_set = Set()
		i = 0
		for row in reader:
			#id, strainid, target_id, qc_method_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, readme_id =row
			strainid = int(row[col_name2index['strainid']])	#2008-09-10
			target_id = int(row[col_name2index['target_id']])
			qc_method_id = int(row[col_name2index['qc_method_id']])
			mismatch_rate = float(row[col_name2index['mismatch_rate']])
			no_of_mismatches = int(row[col_name2index['no_of_mismatches']])
			no_of_non_NA_pairs = int(row[col_name2index['no_of_non_NA_pairs']])
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
	
	
	def testAllPlateIDinPlateSet(self, plate_id_ls, plate_id2plate_set):
		"""
		2008-09-12
		"""
		plate_set = None
		all_plate_id_in_plate_set = 1	#test whether all plate ids in previous plate sets or not
		for plate_id in plate_id_ls:
			if plate_id!=0:
				if plate_id not in plate_id2plate_set:
					all_plate_id_in_plate_set = 0
					break
				else:
					if plate_set==None:
						plate_set = plate_id2plate_set[plate_id]
					elif plate_id2plate_set[plate_id]!=plate_set:
						sys.stderr.write("This plate_id_ls, %s, has >1 plate_sets: %s, %s.\n"%(repr(plate_id_ls), plate_set, plate_id2plate_set[plate_id]))
						all_plate_id_in_plate_set = 0
						break
		return_data = PassingData()
		return_data.all_plate_id_in_plate_set = all_plate_id_in_plate_set
		return_data.plate_set = plate_set
		return return_data
	
	def alignStrainsAccordingToSeqPlate(self, db):
		"""
		2008-09-12
			group strains by sequenom plate set(usually 4 plates comprise 1 set)
		"""
		sys.stderr.write("Aligning strainid according to sequenom plate ...")
		#rows = db.metadata.bind.execute("select id, group_concat(seqinfoid order by seqinfoid) as plate_set from (select distinct s.id, c.seqinfoid from calls_byseq c, strain s where s.ecotypeid=c.ecotypeid and (s.plateid=c.plateid or s.plateid is NULL ) and (s.wellid=c.wellid  or s.wellid is null) order by s.id, c.seqinfoid, c.snpid ) as newt group by id order by plate_set")
		rows = StockDB.Strain.query.all()
		plate_set2strain_id_ls = {}
		plate_set2index = {}
		plate_id2plate_set = {}
		unprocessed_plate_set2strain_id_ls = {}
		seqinfoid_name_ls = ['seqinfoid1', 'seqinfoid2', 'seqinfoid3', 'seqinfoid4']
		for row in rows:
			plate_id_ls = []
			four_plate_complete = True
			for i in range(len(seqinfoid_name_ls)):
				seqinfoid = getattr(row, seqinfoid_name_ls[i], None)
				if seqinfoid is not None:
					plate_id_ls.append(seqinfoid)
				else:
					plate_id_ls.append(0)
					four_plate_complete = False
			strain_id = row.id
			plate_set = tuple(plate_id_ls)
			if four_plate_complete:
				if plate_set not in plate_set2index:
					plate_set2index[plate_set] = len(plate_set2index)
					for plate_id in plate_id_ls:
						plate_id2plate_set[plate_id] = plate_set
					plate_set2strain_id_ls[plate_set] = []
				plate_set2strain_id_ls[plate_set].append(strain_id)
			else:	#incomplete plate set, defer them to handle later
				pdata = self.testAllPlateIDinPlateSet(plate_id_ls, plate_id2plate_set)
				if pdata.all_plate_id_in_plate_set:
					plate_set2strain_id_ls[pdata.plate_set].append(strain_id)
				else:	#new plate_set, there might be complet 4-plate set to cover it later, so defer this.
					if plate_set not in unprocessed_plate_set2strain_id_ls:
						unprocessed_plate_set2strain_id_ls[plate_set] = []
					unprocessed_plate_set2strain_id_ls[plate_set].append(strain_id)
		
		for unprocessed_plate_set, strain_id_ls in unprocessed_plate_set2strain_id_ls.iteritems():
			plate_id_ls = unprocessed_plate_set
			pdata = self.testAllPlateIDinPlateSet(plate_id_ls, plate_id2plate_set)
			if pdata.all_plate_id_in_plate_set:
				for strain_id in strain_id_ls:
					plate_set2strain_id_ls[pdata.plate_set].append(strain_id)
			else:	#new plate set
				plate_set2index[unprocessed_plate_set] = len(plate_set2index)
				for plate_id in plate_id_ls:
					plate_id2plate_set[plate_id] = unprocessed_plate_set
				plate_set2strain_id_ls[unprocessed_plate_set] = []
				for strain_id in strain_id_ls:
					plate_set2strain_id_ls[unprocessed_plate_set].append(strain_id)
		
		strain_id2plate_set = {}
		for plate_set, strain_id_ls in plate_set2strain_id_ls.iteritems():
			for strain_id in strain_id_ls:
				strain_id2plate_set[strain_id] = plate_set
		
		plate_info = PassingData()
		plate_info.plate_set2strain_id_ls = plate_set2strain_id_ls
		plate_info.plate_set2index = plate_set2index
		plate_info.plate_id2plate_set = plate_id2plate_set
		plate_info.strain_id2plate_set = strain_id2plate_set
		sys.stderr.write("Done.\n")
		return plate_info
	
	def getStrainInfoGivenPlateInfo(self, db, plate_info, strain_id_info_query, strain_id_set=None):
		"""
		2008-09-13
			order/group the strains according to plate_set, country, strain longitude
			fetch appropriate labels for each strain
		"""
		sys.stderr.write("Getting strain_info given plate_info ...")

		
		#fetch appropriate label, and put strain id in country_longitude order within each plate
		plate_set2strain_id_ls_in_GPS_order = {}
		strain_id2label = {}
		rows = db.metadata.bind.execute(strain_id_info_query)
		for row in rows:
			if strain_id_set and row.strainid not in strain_id_set:	#skip
				continue
			plate_set = plate_info.strain_id2plate_set[row.strainid]
			if plate_set not in plate_set2strain_id_ls_in_GPS_order:
				plate_set2strain_id_ls_in_GPS_order[plate_set] = []
			plate_set2strain_id_ls_in_GPS_order[plate_set].append(row.strainid)
			
			if len(row.sitename)>10:	#cut short on the site name
				sitename = row.sitename[:10]
			else:
				sitename = row.sitename
			strain_label = '%s_%s_%s_%s_%s'%(row.abbr, sitename, row.nativename, row.strainid, repr(plate_set)[1:-1])
			strain_id2label[row.strainid] = strain_label
		
		#put in plate_set order, assign row index
		plate_set_ls = plate_set2strain_id_ls_in_GPS_order.keys()
		plate_set_ls.sort()
		no_of_plates = len(plate_set_ls)
		strain_id_ls = []
		strain_id2index = {}
		strain_label_ls = []
		for i in range(no_of_plates):
			plate_set = plate_set_ls[i]
			plate_strain_id_ls = plate_set2strain_id_ls_in_GPS_order[plate_set]
			if i!=0:	#insert separator, but not before the first plate_set
				strain_id2index[-i] = len(strain_id2index)
				strain_id_ls.append(-i)
				strain_label_ls.append('')
			for strain_id in plate_strain_id_ls:
				strain_id2index[strain_id] = len(strain_id2index)
				strain_id_ls.append(strain_id)
				strain_label_ls.append(strain_id2label[strain_id])
		
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
		"""
		2008-09-10
			column in input_fname is determined on the fly
		"""
		sys.stderr.write("Getting data matrix from  %s ... \n"%input_fname)
		data_matrix = num.zeros([len(strain_id_info.strain_id_ls), len(target_id_info.strain_id_ls)], num.float)
		data_matrix[:] = -1
		reader = csv.reader(open(input_fname), delimiter='\t')
		#figure out which variable is in which column
		header = reader.next()
		col_name2index = {}
		for i in range(len(header)):
			column_name = header[i]
			col_name2index[column_name] = i
		min_value = None
		max_value = None
		i = 0
		for row in reader:
			"""
			id, strainid, target_id, qc_method_id, mismatch_rate, no_of_mismatches, no_of_non_NA_pairs, readme_id =row
			strainid = int(strainid)
			target_id = int(target_id)
			qc_method_id = int(qc_method_id)
			mismatch_rate = float(mismatch_rate)
			no_of_mismatches = int(no_of_mismatches)
			no_of_non_NA_pairs = int(no_of_non_NA_pairs)
			"""
			strainid = int(row[col_name2index['strainid']])
			target_id = int(row[col_name2index['target_id']])
			qc_method_id = int(row[col_name2index['qc_method_id']])
			mismatch_rate = float(row[col_name2index['mismatch_rate']])
			no_of_mismatches = int(row[col_name2index['no_of_mismatches']])
			no_of_non_NA_pairs = int(row[col_name2index['no_of_non_NA_pairs']])
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
				if QC_method_id==4:	#149 self-cross-match
					row_index = strain_id_info.strain_id2index.get(target_id)
					col_index = target_id_info.strain_id2index.get(strainid)
					data_matrix[row_index, col_index] = data_value
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
		if self.QC_method_id ==4:
			sql_table_str = "from %s e, %s s, %s a, %s c"%(StockDB.Ecotype.table.name, StockDB.Site.table.name, StockDB.Address.table.name,\
								StockDB.Country.table.name)
			common_where_condition = "where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id %s " + order_by_sentence
			
			strain_where_condition = common_where_condition%(" and e.id=st.ecotypeid")
			strain_id_info_query = "select distinct st.id as strainid, e.id as ecotypeid, e.nativename, s.name as sitename, c.abbr %s, %s st %s"%(sql_table_str, StockDB.Strain.table.name, strain_where_condition)
		else:
			sql_table_str = "from %s q, %s e, %s s, %s a, %s c"%(StockDB.QCCrossMatch.table.name, StockDB.Ecotype.table.name, StockDB.Site.table.name, StockDB.Address.table.name,\
									StockDB.Country.table.name)
			common_where_condition = "where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id %s"+ " and q.qc_method_id=%s and q.no_of_non_NA_pairs>=%s and q.mismatch_rate<=%s "%\
				(self.QC_method_id, self.min_no_of_non_NAs, self.max_mismatch_rate) + order_by_sentence
			
			strain_where_condition = common_where_condition%(" and e.id=st.ecotypeid and st.id=q.strainid")
			strain_id_info_query = "select distinct q.strainid, e.id as ecotypeid, e.nativename, s.name as sitename, c.abbr %s, %s st %s"%(sql_table_str, StockDB.Strain.table.name, strain_where_condition)
		
		if self.how_to_group_strains==2 or self.how_to_group_strains==3:
			plate_info = self.alignStrainsAccordingToSeqPlate(db)
			id_set_data = PassingData()
			id_set_data.strain_id_set = None
			id_set_data.target_id_set = None
		elif self.input_fname:
			id_set_data = self.getStrainidTargetidFromFile(db, self.QC_method_id, self.input_fname, self.max_mismatch_rate, self.min_no_of_non_NAs)
		else:
			id_set_data = PassingData()
			id_set_data.strain_id_set = None
			id_set_data.target_id_set = None
		
		if self.how_to_group_strains==2 or self.how_to_group_strains==3:
			strain_id_info = self.getStrainInfoGivenPlateInfo(db, plate_info, strain_id_info_query, strain_id_set=None)
		else:
			strain_id_info = self.getStrainIDInfo(db, strain_id_info_query, id_set_data.strain_id_set)
		
		if self.QC_method_id==4:
			if self.how_to_group_strains==3:
				#2008-09-15 column strain id is in country, strain-longitude order
				target_id_info = self.getStrainIDInfo(db, strain_id_info_query, id_set_data.strain_id_set)
			else:
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