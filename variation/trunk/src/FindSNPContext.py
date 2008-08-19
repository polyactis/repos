#!/usr/bin/env python
"""

Examples:
	FindSNPContext.py -u yh -c

Description:
	program to find the context (nearby genes) of a snp. It fills results into db table snps_context.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions, PassingData
from variation.src.Stock_250kDB import Stock_250kDB, QCMethod, CallQC, Snps, README, SnpsContext
from pymodule.db import formReadmeObj
from transfac.src.TFBindingSiteParse import TFBindingSiteParse
from annot.bin.codense.common import get_entrezgene_annotated_anchor

class FindSNPContext(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('output_fname', 0, ): [None, 'o', 1, 'if given, QC results will be outputed into it.'],\
							('max_upstream_distance', 0, int): [50000, 'a', 1, "maximum distance allowed between a SNP and a gene's start position given the SNP is upstream.\
								a negative value results the program to take the closest upstream gene."],\
							('max_downstream_distance', 0, int): [50000, 'm', 1, "maximum distance allowed between a SNP and a gene's stop position given the SNP is downstream.\
								a negative value results the program to take the closest downstream gene."],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-08-19
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def find_SNP_context(self, db, curs, snp_locus_table, snp_locus_context_table, entrezgene_mapping_table='genome.entrezgene_mapping',\
						annot_assembly_table='genome.annot_assembly', tax_id=3702, max_upstream_distance=50000, max_downstream_distance=50000, need_commit=0, debug=0):
		"""
		2008-08-12
			add max_upstream_distance and max_downstream_distance
			take target genes, which the SNP touches, is upstream or is downstream (within distance allowed)
			if none found, will try grab genes without strand info (no idea if it's upstream or downstream) if they exist
		2008-07-02
			modify to deal with stock_250k (mysql db on papaya)
		2007-04-29
			to determine which SNP is in coding or non-coding region
		"""
		sys.stderr.write("Finding SNP context ... \n")
		chromosome2anchor_gene_tuple_ls, gene_id2coord = get_entrezgene_annotated_anchor(curs, tax_id, entrezgene_mapping_table, annot_assembly_table)
		
		#from variation.src.Stock_250kDB import SnpsContext, Snps
		#session = db.session
		#if not debug:
		#	session.begin()
		"""
		for the 3 deletions
		4,268809,0
		4,269962,8
		4,270712,0
		curs.execute("select id, chromosome, position from %s where end_position is null and chromosome=4 and (position=268809 or position=269962 or position=270712)"%(snp_locus_table))	#only SNPs, not segments
		"""
		curs.execute("select id, chromosome, position from %s where end_position is null"%(snp_locus_table))	#only SNPs, not segments
		rows = curs.fetchall()
		offset_index = 0
		block_size = 5000
		#rows = Snps.query.filter_by(end_position!=None).offset(offset_index).limit(block_size)
		counter = 0
		#while rows.count()!=0:
		no_of_contexts = 0
		for row in rows:
			snp_locus_id, chromosome, position = row
			chromosome = str(chromosome)
			regulatory_coord = (chromosome, position, position)
			pdata = TFBindingSiteParse.return_target_gene_ls(regulatory_coord, chromosome2anchor_gene_tuple_ls, gene_id2coord, \
															max_upstream_distance, max_downstream_distance)
			target_gene_ls = []
			if pdata.regulatory_touch_target_gene_ls:
				target_gene_ls += pdata.regulatory_touch_target_gene_ls
			if pdata.regulatory_is_left_upstream_target_gene_ls or pdata.regulatory_is_right_upstream_target_gene_ls:
				target_gene_ls += pdata.regulatory_is_left_upstream_target_gene_ls + pdata.regulatory_is_right_upstream_target_gene_ls
			if pdata.regulatory_is_left_downstream_target_gene_ls or pdata.regulatory_is_right_downstream_target_gene_ls:
				target_gene_ls += pdata.regulatory_is_left_downstream_target_gene_ls + pdata.regulatory_is_right_downstream_target_gene_ls
			
			if not target_gene_ls:	#getting desperate, take the genes with no strand info as well
				if pdata.regulatory_is_left_target_gene_ls:
					target_gene_ls += pdata.regulatory_is_left_target_gene_ls
				if pdata.regulatory_is_right_target_gene_ls:
					target_gene_ls += pdata.regulatory_is_right_target_gene_ls
			
			for target_gene_tuple in target_gene_ls:
				gene_id = target_gene_tuple[1]
				left_or_right = target_gene_tuple[2]
				upstream_or_downstream = target_gene_tuple[3]
				disp_pos = target_gene_tuple[4]
				gene_start, gene_stop, gene_strand, gene_genomic_gi = gene_id2coord[gene_id]
				#snps_context = SnpsContext(snps_id=snp_locus_id, disp_pos=disp_pos, gene_id=gene_id, gene_strand=gene_strand,\
				#						left_or_right=left_or_right, disp_pos_comment=upstream_or_downstream)
				#session.save(snps_context)
				curs.execute("insert into %s(snps_id, disp_pos, gene_id, gene_strand, left_or_right, disp_pos_comment) values (%s, %s, %s, '%s', '%s', '%s')"%\
					(snp_locus_context_table, snp_locus_id, disp_pos, gene_id, gene_strand, left_or_right, upstream_or_downstream))
				no_of_contexts += 1
			offset_index += 1
			counter += 1
			#session.flush()
			if counter%5000==0:
				sys.stderr.write("%s%s"%('\x08'*10, counter))
			#rows = Snps.query.filter_by(end_position!=None).offset(offset_index).limit(block_size)
		if need_commit:
			curs.execute("commit")
			#session.commit()
			#session.clear()
		sys.stderr.write("%s SNPs with %s contexts.\n"%(counter, no_of_contexts))
		
	def run(self):
		import MySQLdb
		mysql_conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user = self.db_user, passwd = self.db_passwd)
		mysql_curs = mysql_conn.cursor()
		
		#2008-08-13 2 elixir dbs are bad. it causes tables to be cross-created in the two databases.
		#from transfac.src.GenomeDB import GenomeDatabase
		#genome_db = GenomeDatabase(drivername='mysql', username=user, password=passwd, hostname=hostname, database='genome')
		#from transfac.src.GenomeDB import getEntrezgeneAnnotatedAnchor
		#chromosome2anchor_gene_tuple_ls, gene_id2coord = getEntrezgeneAnnotatedAnchor(genome_db, tax_id=3702)
		#del genome_db
		
		db = Stock_250kDB(drivername=self.drivername, username=self.db_user, password=self.db_passwd, hostname=self.hostname, database=self.dbname)
		mysql_conn.autocommit(True)
		self.find_SNP_context(db, mysql_curs, Snps.table.name, SnpsContext.table.name, max_upstream_distance=self.max_upstream_distance,\
							max_downstream_distance=self.max_downstream_distance, need_commit=self.commit, debug=self.debug)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = FindSNPContext
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()