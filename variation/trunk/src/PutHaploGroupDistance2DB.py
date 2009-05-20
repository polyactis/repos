#!/usr/bin/env python
"""

Examples:
	PutHaploGroupDistance2DB.py -u yh -r -c
	
Description:
	2009-4-18
		1. calculate & save pairwise genetic distance between all haplotype groups
		2. calculate & save PCA results on SNP matrix from all haplotype groups
		3. calculate & save geographic distances between all ecotypes
		4. calculate & save GPS positions for the haplotype groups
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import csv
import traceback, gc, subprocess, numpy
import StockDB
from pymodule import SNPData, nt2number, calGreatCircleDistance
from common import get_ecotypeid2pos

class PutHaploGroupDistance2DB(object):
	__doc__ = __doc__	#use documentation in the beginning of the file as this class's doc
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, '', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('commit',0, int): [0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode, test-run 10 haplotype groups.'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2009-4-18
		"""
		
		from pymodule import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getHaploGroupSNPMatrix(self):
		"""
		2009-4-18
		"""
		sys.stderr.write("Getting HaploGroup SNP matrix ...")
		
		col_id_ls = []
		row_id_ls = []
		if self.debug:
			no_of_rows = 10
		else:
			no_of_rows = StockDB.HaploGroup.query.count()
		
		col_id2col_index = {}
		for row in StockDB.SNPs.query.order_by(StockDB.SNPs.chromosome).order_by(StockDB.SNPs.position):
			col_id_ls.append(row.id)
			col_id2col_index[row.id] = len(col_id2col_index)
		
		no_of_cols = len(col_id2col_index)
		
		data_matrix = numpy.zeros([no_of_rows, no_of_cols], numpy.int8)
		rows = StockDB.HaploGroup.query()
		row_index = 0
		for row in rows:
			data_rows = StockDB.FilteredCalls.query.filter_by(ecotypeid=row.ref_ecotypeid)
			row_index = len(row_id_ls)
			for one_call in data_rows:
				nt_number = nt2number[one_call.allele]
				col_index = col_id2col_index[one_call.snpid]
				data_matrix[row_index][col_index] = nt_number
			row_id_ls.append(row.id)
			if self.debug and row_index==no_of_rows-1:
				break
		snpData = SNPData(col_id_ls=col_id_ls, row_id_ls=row_id_ls, data_matrix=data_matrix)
		sys.stderr.write("Done.\n")
		return snpData
		
	def saveHaploGroupPairwiseDist(self, session, snpData):
		"""
		2009-4-18
		"""
		row_id2pairwise_dist = snpData.calRowPairwiseDist()
		sys.stderr.write("Saving pairwise distances between haplotype groups ... \n")
		counter = 0
		for row_id, pairwise_dist in row_id2pairwise_dist.iteritems():
			for row in pairwise_dist:
				mismatch_rate, row_id2, no_of_mismatches, no_of_non_NA_pairs = row[:4]
				entry = StockDB.HaploGroupPairwiseGeneticDist(haplo_group_id1=row_id, haplo_group_id2=row_id2, mismatch_rate=mismatch_rate, \
											no_of_mismatches=no_of_mismatches, no_of_non_NA_pairs=no_of_non_NA_pairs)
				session.save(entry)
				session.flush()
				counter += 1
				if counter%5000==0:
					sys.stderr.write("%s%s"%('\x08'*40, counter))
		sys.stderr.write("%s%s\t"%('\x08'*40, counter))
		sys.stderr.write("Done.\n")
	
	def saveHaploGroupPCA(self, session, snpData):
		"""
		2009-4-18
		"""
		sys.stderr.write("Saving haplotype group PCA information ...")
		import pca_module
		T, P, explained_var = pca_module.PCA_svd(snpData.data_matrix, standardize=False)
		for i in range(len(snpData.row_id_ls)):
			haplo_group_id = snpData.row_id_ls[i]
			for j in range(T.shape[1]):
				pc_value = T[i][j]
				pc_db_entry = StockDB.HaploGroupPCValues(haplo_group_id=haplo_group_id, which_pc=j+1, pc_value=pc_value)
				session.save(pc_db_entry)
				session.flush()
		
		for i in range(len(explained_var)):
			db_entry = StockDB.HaploGroupEigenValues(which_eigen=i+1, variance_perc=explained_var[i])
			session.save(db_entry)
			session.flush()
		sys.stderr.write("Done.\n")
	
	def calAndSaveGeographicDistanceBetweenEcotypes(self, session, ecotypeid2pos):
		"""
		2009-4-18
		"""
		sys.stderr.write("Calculating and saving geographic distance between ecotypes ... \n")
		ecotypeid_ls = ecotypeid2pos.keys()
		counter = 0
		for i in range(len(ecotypeid_ls)):
			ecotypeid1 = ecotypeid_ls[i]
			for j in range(i+1, len(ecotypeid_ls)):
				ecotypeid2 = ecotypeid_ls[j]
				pos1 = ecotypeid2pos[ecotypeid1]
				pos2 = ecotypeid2pos[ecotypeid2]
				geoDist = calGreatCircleDistance(pos1[0], pos1[1], pos2[0], pos2[1])
				db_entry = StockDB.EcotypePairwiseGeographicDist(ecotypeid1=ecotypeid1, ecotypeid2=ecotypeid2, distance=geoDist)
				session.save(db_entry)
				session.flush()
				counter += 1
				if self.debug and counter==100:
					break
				if counter%5000==0:
					sys.stderr.write("%s%s"%('\x08'*30, counter))
			if self.debug and counter>=100:
				break
		sys.stderr.write("%s%s\t"%('\x08'*30, counter))
		sys.stderr.write("Done.\n")
	
	def saveHaploGroupGPS(self, session, ecotypeid2pos, powerFactor=2):
		"""
		2009-4-21
			calculate the GPS position of a haplotype group:
				a normalized weighted sum of GPS positions of ecotypes within that group
				the weight is math.pow(frequency, powerFactor)
				
		"""
		sys.stderr.write("Calculating and saving GPS positions for the haplotype groups ... \n")
		rows = StockDB.HaploGroup.query()
		counter = 0
		for row in rows:
			latlon2count = {}
			for ecotype in row.ecotypes:
				if ecotype.id in ecotypeid2pos:
					pos = ecotypeid2pos[ecotype.id]
					if pos not in latlon2count:
						latlon2count[pos] = 0
					latlon2count[pos] += 1
			count_sum = 0
			lat_sum = 0.0
			lon_sum = 0.0
			for latlon, count in latlon2count.iteritems():
				weight = math.pow(count, powerFactor)
				count_sum += weight
				lat_sum += latlon[0]*weight
				lon_sum += latlon[1]*weight
			if count_sum!=0:
				row.latitude = lat_sum/count_sum
				row.longitude = lon_sum/count_sum
				session.save_or_update(row)
				session.flush()
			counter += 1
			if self.debug and counter==50:
				break
			if counter%200==0:
				sys.stderr.write("%s%s"%('\x08'*40, counter))
			
		sys.stderr.write("%s%s"%('\x08'*40, counter))
		sys.stderr.write("\tDone.\n")
			
	
	def run(self):
		"""
		"""
		if self.debug:
			import pdb
			pdb.set_trace()
			
		db = StockDB.StockDB(drivername=self.drivername, username=self.db_user,
				 		password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup()
		session = db.session
		#session.begin()
		
		snpData = self.getHaploGroupSNPMatrix()
		self.saveHaploGroupPairwiseDist(session, snpData)
		self.saveHaploGroupPCA(session, snpData)
		ecotypeid2pos = get_ecotypeid2pos(db.metadata.bind, StockDB.Ecotype.table.name)
		self.saveHaploGroupGPS(session, ecotypeid2pos)
		self.calAndSaveGeographicDistanceBetweenEcotypes(session, ecotypeid2pos)
				
		if self.commit:
			session.commit()
			session.clear()
		else:	#default is rollback(). to demonstrate good programming
			session.rollback()

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = PutHaploGroupDistance2DB
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()