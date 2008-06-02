#!/usr/bin/env python
"""

Examples:
	AddPopulationInfo.py -i /Network/Data/250k/db/reference_dataset/149SNP -l popid2ecotypeid_50 -o bin/structure_test/149_popid2ecotypeid_50.csv
	
Description:
	Replace 2nd column with population info. If no population info, the row is removed.
	
"""
from __init__ import *

class AddPopulationInfo(object):
	__doc__ = __doc__
	option_default_dict = {('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock', 'd', 1, 'database name', ],\
							('user', 1, ): [None, 'u', 1, 'database username', ],\
							('passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('population_table', 1, ): [None, 'l', 1, 'Table storing population id versus ecotypeid'],\
							('input_fname',1, ): [None, 'i', 1, ''],\
							('output_fname', 1, ): [None, 'o', 1, '', ],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self, **keywords):
		"""
		2008-06-02
		"""
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def run(self):
		"""
		"""
		import MySQLdb
		conn = MySQLdb.connect(db=self.dbname,host=self.hostname, user = self.user, passwd = self.passwd)
		curs = conn.cursor()
		
		snpData = SNPData(input_fname=self.input_fname, turn_into_array=1, ignore_2nd_column=1)
		from OutputPopulation import OutputPopulation
		
		popid2ecotypeid_ls = OutputPopulation.get_popid2ecotypeid_ls(curs, self.population_table)
		
		ecotypeid2popid = {}
		for popid, ecotypeid_ls in popid2ecotypeid_ls.iteritems():
			for ecotypeid in ecotypeid_ls:
				ecotypeid2popid[ecotypeid] = popid
		pop_id_ls = []
		rows_to_be_tossed_out = Set()
		for i in range(len(snpData.row_id_ls)):
			ecotype_id = int(snpData.row_id_ls[i])
			if ecotype_id not in ecotypeid2popid:
				rows_to_be_tossed_out.add(i)
				pop_id_ls.append(None)	#dont' know population, a placeholder
			else:
				pop_id_ls.append(ecotypeid2popid[ecotype_id])
		
		snpData.strain_acc_list = snpData.row_id_ls
		snpData.category_list = pop_id_ls
		
		snpData.tofile(self.output_fname, rows_to_be_tossed_out = rows_to_be_tossed_out)



if __name__ == '__main__':
	main_class = AddPopulationInfo
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()