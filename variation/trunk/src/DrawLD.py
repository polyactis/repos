#!/usr/bin/env python
"""
Examples:
	#program to draw LD
	DrawLD.py -L /Network/Data/250k/tmp-yh/call_method_17_LD_m0.2.tsv -o ./call_method_17_LD -u yh -p secret -f 1000000 -t 1000000

Description:
	2008-10-01 program to draw LD
		
"""

import sys, os, math
#bit_number = math.log(sys.maxint)/math.log(2)
#if bit_number>40:       #64bit
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, cPickle
import warnings, traceback
from pymodule import PassingData, figureOutDelimiter, getColName2IndexFromHeader, getListOutOfStr
from sets import Set
import matplotlib; matplotlib.use("Agg")	#to avoid popup and collapse in X11-disabled environment
from matplotlib import rcParams
rcParams['font.size'] = 6
rcParams['legend.fontsize'] = 6
#rcParams['text.fontsize'] = 6	#deprecated. use font.size instead
rcParams['axes.labelsize'] = 6
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 6
rcParams['ytick.labelsize'] = 6

from matplotlib.patches import Polygon
from DrawSNPRegion import DrawSNPRegion, LD_statistic
import Stock_250kDB
import pylab

class DrawLD(DrawSNPRegion):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('LD_fname', 1, ): [None, 'L', 1, 'the file containing LD info, output of MpiLD.py', ],\
							('min_MAF', 1, float): [0.1, 'n', 1, 'minimum Minor Allele Frequency for both SNPs to be included'],\
							('chr1', 0, int): [1, 'c', 1, 'chromsome 1 (row)', ],\
							('start1', 0, int): [1, 'e', 1, 'start position of chromsome 1 (row)', ],\
							('stop1', 0, int): [100000, 'f', 1, 'stop position of chromsome 1 (row)', ],\
							('chr2', 0, int): [1, 'a', 1, 'chromsome 2 (column)', ],\
							('start2', 0, int): [1, 's', 1, 'start position of chromsome 2 (column)', ],\
							('stop2', 0, int): [100000, 't', 1, 'stop position of chromsome 2 (column)', ],\
							('min_gap', 0, int): [100000, '', 1, 'two loci have to be this apart to be drawn LD', ],\
							("output_fname_prefix", 1, ): [None, 'o', 1, 'image filename prefix'],\
							("which_LD_statistic", 1, int): [2, 'w', 1, 'which LD_statistic to plot, 1=r2, 2=|D_prime|, 3=|D|'],\
							('debug', 0, int):[0, 'b', 1, 'debug mode. 1=level 1 (pdb mode). 2=level 2 (same as 1 except no pdb mode)'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2008-10-10
			add option min_gap to draw LD between loci at least this apart
		2008-09-24
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def addOnePolygon(self, axe_LD, x_window, y_window, facecolor, xylim_data):
		x1, x2 = x_window
		y1, y2 = y_window
		if xylim_data.xlim[0]==-1:
			xylim_data.xlim[0] = x1
		elif x1<xylim_data.xlim[0]:
			xylim_data.xlim[0] = x1
		
		if xylim_data.xlim[1]==-1:
			xylim_data.xlim[1] = x2
		elif x2>xylim_data.xlim[1]:
			xylim_data.xlim[1] = x2
		
		if xylim_data.ylim[0]==-1:
			xylim_data.ylim[0] = y1
		elif y1<xylim_data.ylim[0]:
			xylim_data.ylim[0] = y1
		
		if xylim_data.ylim[1]==-1:
			xylim_data.ylim[1] = y2
		elif y2>xylim_data.ylim[1]:
			xylim_data.ylim[1] = y2
		xs = [x1, x1, x2, x2]
		ys = [y1, y2, y2, y1]
		poly = Polygon(zip(xs, ys), facecolor=facecolor, linewidth=0)
		axe_LD.add_patch(poly)
	
	def drawLD(self, axe_LD, LD_fname, row_snp_region, col_snp_region, min_MAF, which_LD_statistic, min_gap=100000):
		"""
		2008-10-03
		"""
		sys.stderr.write("Drawing LD from %s ...\n"%(LD_fname))
		reader = csv.reader(open(LD_fname), delimiter='\t')
		col_name2index = getColName2IndexFromHeader(reader.next())
		counter = 0
		real_counter= 0
		xylim_data = PassingData(xlim = [-1,-1], ylim=[-1,-1])
		for row in reader:
			snp1 = row[col_name2index['snp1']].split('_')
			snp1 = tuple(map(int, snp1))
			snp2 = row[col_name2index['snp2']].split('_')
			snp2 = tuple(map(int, snp2))
			allele1_freq = float(row[col_name2index['allele1_freq']])
			allele2_freq = float(row[col_name2index['allele2_freq']])
			if allele1_freq>=min_MAF and allele2_freq>=min_MAF:	#meet the minimum minor-allele-frequency
				LD_stat = float(row[col_name2index[LD_statistic.get_name(which_LD_statistic)]])
				LD_stat = abs(LD_stat)
				fc = self.r2ToRGBColor(LD_stat)
				if snp1 in row_snp_region.chr_pos2adjacent_window and snp2 in col_snp_region.chr_pos2adjacent_window:
					if snp1[0]==snp2[0] and abs(snp1[1]-snp2[1])<=min_gap:	#too close, ignore
						continue
					self.addOnePolygon(axe_LD, row_snp_region.chr_pos2adjacent_window[snp1], col_snp_region.chr_pos2adjacent_window[snp2], \
									fc, xylim_data)
					real_counter += 1
				if snp2 in row_snp_region.chr_pos2adjacent_window and snp1 in col_snp_region.chr_pos2adjacent_window:
					if snp1[0]==snp2[0] and abs(snp1[1]-snp2[1])<=min_gap:	#too close, ignore
						continue
					self.addOnePolygon(axe_LD, row_snp_region.chr_pos2adjacent_window[snp2], col_snp_region.chr_pos2adjacent_window[snp1], \
									fc, xylim_data)
					real_counter += 1
			counter += 1
			if counter%100000==0:
				sys.stderr.write('%s\t%s'%('\x08'*100, counter))
				if counter%500000==0 and self.debug>0:
					break
		del reader
		sys.stderr.write("%s LD drawn. Done.\n"%real_counter)
		return xylim_data
	
	def run(self):
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		snp_info = self.getSNPInfo(db)
		row_snp_region = self.findSNPsInRegion(snp_info, self.chr1, self.start1, self.stop1)
		col_snp_region = self.findSNPsInRegion(snp_info, self.chr2, self.start2, self.stop2)
		axe_LD = pylab.axes([0.1, 0.1, 0.8, 0.8], frameon=False)
		axe_LD.grid(True, alpha=0.3)
		axe_LD.title.set_text('LD chr%s(%s-%s) vs chr%s(%s-%s)'%(self.chr1, self.start1, self.stop1, self.chr2, self.start2, self.stop2))
		axe_LD.set_ylabel('chromosome %s'%self.chr1)
		axe_LD.set_xlabel('chromosome %s'%self.chr2)
		xylim_data = self.drawLD(axe_LD, self.LD_fname, row_snp_region, col_snp_region, self.min_MAF, self.which_LD_statistic, self.min_gap)
		axe_LD.set_xlim(xylim_data.xlim)
		axe_LD.set_ylim(xylim_data.ylim)
		
		axe_LD_legend = pylab.axes([0.92, 0.8, 0.08, 0.1], frameon=False)	#axes for the legend of LD
		axe_LD_legend.set_xticks([])
		axe_LD_legend.set_yticks([])
		self.drawLDLegend(axe_LD_legend, self.which_LD_statistic)
		
		png_output_fname = '%s.png'%self.output_fname_prefix
		pylab.savefig(png_output_fname, dpi=400)
		#pylab.savefig('%s.svg'%self.output_fname_prefix, dpi=300)

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = DrawLD
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()