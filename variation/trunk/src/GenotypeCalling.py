#!/usr/bin/env python
"""
Usage: GenotypeCalling.py [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, localhost(default)
	-d ..., --dbname=...	the database name, stock(default)
	-k ..., --schema=...	which schema in the database, dbsnp(default)IGNORE
	-p ...,	popid2ecotypeid_table
	-s ...,	strain_info table, 'ecotype'(default)
	-o ...,	output_fname_prefix
	-a ...,	picture area, lower left corner longitude, latitude,
		upper right corner longitude, latitude default: "-130,10,140,70"
	-l ...,	type of label, 1(popid), 2(pop size), 3(selfing rate)
	-g ...,	max_dist to connect 2 sites, for -w, '50'(km) (default)
	-y ...,	which method estimates selfing rate, 1(default)
	-f ...,	selfing_rate_table (needed for label type 3)
	-w	draw site network showing how nearby sites give rise to population
	-b, --debug	enable debug
	-r, --report	enable more progress-related output
	-h, --help	show this help

Examples:
	
	
Description:
	
2007-10-29 so far just functions to draw cluster plots investigating intensity/peak data of 149snp genotype calls. intended to do genotype calling for 149snp data
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


from pymodule.latex import outputFigureInLatex

"""
2007-10-26
	look at the intensity data of 149-snp data
"""
def MA_plot_for_calls(curs):
	"""
	2007-10-26
	"""
	curs.execute("select peak1, peak2 from calls where callhet is not null")
	rows = curs.fetchall()
	diff_ls = []
	avg_ls = []
	for row in rows:
		peak1, peak2 = row
		if peak1 and peak2:
			diff_ls.append(abs(peak1-peak2))
			avg_ls.append((peak1+peak2)/2)
	return diff_ls, avg_ls


def get_snpid2allele(curs, snps_sequenom_info_table='dbsnp.snps_sequenom_info', snps_table='stock20071008.snps'):
	import sys
	sys.stderr.write("Getting snpid2allele ...")
	snpid2allele = {}
	curs.execute("select ss.id, sl.ext1_call from %s sl, %s ss where sl.snpid=ss.snpid"%(snps_sequenom_info_table, snps_table))
	rows = curs.fetchall()
	for row in rows:
		snpid, allele = row
		snpid2allele[snpid] = allele
	sys.stderr.write("Done.\n")
	return snpid2allele


def get_peak1_peak2_color_ls(curs, snpid2allele, calls_table, extractionid=5):
	"""
	2007-11-08
		add extractionid
	2007-12-27
		calls_table is calls_byseq
	"""
	import sys
	sys.stderr.write("Getting peak1_peak2_color_ls ...")
	peak1_peak2_color_ls = []
	#curs.execute("select snpid, call1, call2, ext1height, ext2height from %s where extractionid=%s and ext1height is not null and ext2height is not null"%(calls_table, extractionid))
	curs.execute("select snpid, call1, call2, ext1height, ext2height from %s where ext1height is not null and ext2height is not null"%(calls_table))
	rows = curs.fetchmany(5000)
	while rows:
		for row in rows:
			snpid, call1, callhet, peak1, peak2 = row
			call1 = call1.upper()
			if callhet:
				color = 1	#green, heterozygous
			elif call1 ==snpid2allele[snpid]:
				color = 0	#red, allele1
			elif call1 == 'N' or call1 =='NA':
				color = 3	#not decided
			else:
				color = 2	#blue, allele2
			peak1_peak2_color_ls.append([peak1, peak2, color])
		rows = curs.fetchmany(5000)
	sys.stderr.write("Done.\n")
	return peak1_peak2_color_ls

def contrast_transform_peak1_peak2(peak1, peak2):
	import math
	S = peak1 + peak2
	if S!=0:
		y = 2*peak2/S-1
		if abs(y)<1:
			contrast = math.sinh(2*y)/math.sinh(2)
		else:
			contrast = y
	else:
		contrast = 0
	return contrast, S

def contrast_transform_peak1_peak2_color_ls(peak1_peak2_color_ls):
	"""
	2007-11-19
		transformation method from Moorhead2006 and Plagnol2007
		both paper shows some typos. vincent sent me his c source code
	"""
	ls = []
	for peak1, peak2, color in peak1_peak2_color_ls:
		contrast, S = contrast_transform_peak1_peak2(peak1, peak2)
		ls.append([contrast, S, color])
	return ls

def get_snpid2peak1_peak2_color_ls(curs, snpid2allele, calls_table, extractionid=5, data_type=1, ecotype2accession_table='at.ecotype2accession', data_2010=[]):
	"""
	2007-11-08
		add extractionid
	2007-12-28
		add data_type
			1: height
			2: area
			3: only the strains matched to 2010 data, use ecotype2accession_table
			4: 3 + call determined by 2010 data
	"""
	import sys
	sys.stderr.write("Getting snpid2peak1_peak2_color_ls ...")
	snpid2peak1_peak2_color_ls = {}
	if data_type==1:
		curs.execute("select snpid, call1, call2, ext1height, ext2height from %s where ext1height is not null and ext2height is not null"%(calls_table))
	elif data_type==2:
		curs.execute("select snpid, call1, call2, ext1area, ext2area from %s where ext1height is not null and ext2height is not null"%(calls_table))
	elif data_type==3:
		curs.execute("select c.snpid, c.call1, c.call2, c.ext1area, c.ext2area from %s c, %s ea where c.ext1height is not null and c.ext2height is not null and c.ecotypeid=ea.ecotype_id"%(calls_table, ecotype2accession_table))
	elif data_type==4:
		snp_acc2col_index, accession_id2row_index, data_matrix_2010 = data_2010
		curs.execute("select ea.accession_id, s.snpid, c.snpid, c.call1, c.call2, c.ext1area, c.ext2area from %s c, %s ea, snps s where c.ext1height is not null and c.ext2height is not null and c.ecotypeid=ea.ecotype_id and s.id=c.snpid"%(calls_table, ecotype2accession_table))
	rows = curs.fetchmany(5000)
	while rows:
		for row in rows:
			if data_type==4:
				accession_id, snp_acc, snpid, call1, callhet, peak1, peak2 = row
				#replace the call with the 2010 call
				if accession_id in accession_id2row_index and snp_acc in snp_acc2col_index:
					call_2010 = data_matrix_2010[accession_id2row_index[accession_id]][snp_acc2col_index[snp_acc]]
					from variation.src.common import number2nt
					call_2010_nt = number2nt[call_2010]
					call1 = call_2010_nt[0]
					if len(call_2010_nt)==2 and call_2010>0: #non NA
						callhet = call_2010_nt[1]
					else:
						callhet = None
				else:
					continue
			else:
				snpid, call1, callhet, peak1, peak2 = row
			if snpid not in snpid2peak1_peak2_color_ls:
				snpid2peak1_peak2_color_ls[snpid] = []
			call1 = call1.upper()
			if callhet:
				color = 1	#green, heterozygous
			elif call1 ==snpid2allele[snpid]:
				color = 0	#red, allele1
			elif call1 == 'N' or call1 =='NA':
				color = 3	#not decided
			else:
				color = 2	#blue, allele2
			snpid2peak1_peak2_color_ls[snpid].append([peak1, peak2, color])
		rows = curs.fetchmany(5000)
	sys.stderr.write("Done.\n")
	return snpid2peak1_peak2_color_ls

def contrast_transform_snpid2peak1_peak2_color_ls(snpid2peak1_peak2_color_ls):
	"""
	2007-11-19
		transformation method from Moorhead2006 and Plagnol2007
	"""
	new_snpid2peak1_peak2_color_ls = {}
	for snpid, peak1_peak2_color_ls in snpid2peak1_peak2_color_ls.iteritems():
		ls = []
		for peak1, peak2, color in peak1_peak2_color_ls:
			contrast, S = contrast_transform_peak1_peak2(peak1, peak2)
			ls.append([contrast, S, color])
		new_snpid2peak1_peak2_color_ls[snpid] = ls
	return new_snpid2peak1_peak2_color_ls

def get_ordered_snpid_ls_and_snpid2name(curs, snps_table='stock20071008.snps'):
	ordered_snpid_ls = []
	snpid2name = {}
	curs.execute("select id, snpid from %s order by chromosome, position"%(snps_table))
	rows = curs.fetchall()
	for row in rows:
		snpid, name = row
		ordered_snpid_ls.append(snpid)
		snpid2name[snpid] = name
	return ordered_snpid_ls, snpid2name

def drawClusterPlot(peak1_peak2_color_ls, color_picked=0, color_to_be_drawn='r', output_fname_prefix=None):
	import pylab
	xls = []
	yls = []
	for peak1, peak2, color in peak1_peak2_color_ls:
		if color==color_picked:
			xls.append(peak1)
			yls.append(peak2)
	pylab.title('cluster plots')
	pylab.xlabel('peak1')
	pylab.ylabel('peak2')
	ps=pylab.scatter(xls, yls, c=color_to_be_drawn, marker='o', alpha=0.2, faceted=False)
	if output_fname_prefix:
		#pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=150)
	pylab.show()


def drawClusterPlotForOneSNP(peak1_peak2_color_ls, figure_title='', output_fname_prefix=None):
	import pylab, os, sys
	sys.stderr.write("Drawing cluster plot for %s ..."%figure_title)
	pylab.clf()
	xls = []
	yls = []
	color_number2xy_ls = {}
	color_number2color_code = {0: 'r',\
		1: 'g',\
		2: 'b',\
		3: 'y'}
	color_number2genotype = {0: 'allele 1',\
		1: 'heterozygote', \
		2: 'allele 2',\
		3: 'NA'}
	for peak1, peak2, color_number in peak1_peak2_color_ls:
		if color_number not in color_number2xy_ls:
			color_number2xy_ls[color_number] = [[],[]]
		color_number2xy_ls[color_number][0].append(peak1)
		color_number2xy_ls[color_number][1].append(peak2)
	if figure_title:
		pylab.title(figure_title)
	else:
		pylab.title('cluster plots')
	pylab.xlabel('peak1')
	pylab.ylabel('peak2')
	pscatter_ls = []
	color_ls = []
	for color_number, xy_ls in color_number2xy_ls.iteritems():
		if len(xy_ls[0])>0:
			if color_number==0 or color_number==2:	#homozygotes use smaller alpha, more transparent cuz too much dat
				pscatter = pylab.scatter(xy_ls[0], xy_ls[1], c=color_number2color_code[color_number], marker='o', alpha=0.2, faceted=False)
			else:	
				pscatter = pylab.scatter(xy_ls[0], xy_ls[1], c=color_number2color_code[color_number], marker='o', alpha=0.6, faceted=False)
			pscatter_ls.append(pscatter)
			color_ls.append(color_number2genotype[color_number])
	pylab.legend(pscatter_ls, color_ls, shadow = True)
	if output_fname_prefix:
		#pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=150)
	sys.stderr.write("Done.\n")
	pylab.show()


def outputClusterPlotsForAllSNPs(ordered_snpid_ls, snpid2name, snpid2peak1_peak2_color_ls, fig_output_dir, output_fname, need_draw_figure=0):
	#from pymodule.latex import outputFigureInLatex
	import os,sys
	outf = open(output_fname, 'w')
	
	for snpid in ordered_snpid_ls:
		snp_name = snpid2name[snpid]
		fig_fname = 'cluster_plot_%s'%(snp_name.replace(' ', '_'))	#replace space with _
		if need_draw_figure:
			output_fname_prefix=os.path.join(fig_output_dir, fig_fname)
			drawClusterPlotForOneSNP(snpid2peak1_peak2_color_ls[snpid], snp_name, output_fname_prefix)
		#output the figure in latex
		fig_fname = os.path.join('figures/', '%s.png'%fig_fname)	#relative path of the figure to the latex file
		fig_label = 'fl%s'%(snp_name.replace(' ', ''))
		caption = 'cluster plot for %s.'%(snp_name)
		outf.write(outputFigureInLatex(fig_fname, caption, fig_label, need_floatpackage=1, fig_width=0.5))
	del outf

def read_2010_x_149SNP(input_fname):
	"""
	2007-12-30
	"""
	from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	snp_acc_ls = header[2:]
	snp_acc2col_index = dict(zip(snp_acc_ls, range(len(snp_acc_ls))))
	accession_id_ls = map(int, strain_acc_list)
	accession_id2row_index = dict(zip(accession_id_ls, range(len(accession_id_ls))))
	return snp_acc2col_index, accession_id2row_index, data_matrix

"""
#2007-10-29
#cluster plots investigating intensity/peak data of 149snp genotype calls
dbname = 'stock20071227'
snps_table='%s.snps'%dbname
snps_sequenom_info_table='dbsnp.snps_sequenom_info'
snpid2allele = get_snpid2allele(curs, snps_sequenom_info_table, snps_table)
calls_table = '%s.calls_byseq'%dbname
peak1_peak2_color_ls = get_peak1_peak2_color_ls(curs, snpid2allele, calls_table)

#new_peak1_peak2_color_ls = contrast_transform_peak1_peak2_color_ls(peak1_peak2_color_ls)
#peak1_peak2_color_ls = new_peak1_peak2_color_ls

output_dir = os.path.expanduser('~/script/variation/doc/PeakDataReport/figures')

import pylab
pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_allele1')
drawClusterPlot(peak1_peak2_color_ls, color_picked=0, color_to_be_drawn='r', output_fname_prefix=output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_het')
drawClusterPlot(peak1_peak2_color_ls, 1, 'g', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_allele2')
drawClusterPlot(peak1_peak2_color_ls, 2, 'b', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_NA')
drawClusterPlot(peak1_peak2_color_ls, 3, 'y', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_allele1_het')
drawClusterPlot(peak1_peak2_color_ls, 0, 'r', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 1, 'g', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_allele2_het')
drawClusterPlot(peak1_peak2_color_ls, 2, 'b', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 1, 'g', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_allele1_allele2')
drawClusterPlot(peak1_peak2_color_ls, 0, 'r', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 2, 'b', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_allele1_NA')
drawClusterPlot(peak1_peak2_color_ls, 0, 'r', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 3, 'y', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_homo_and_hetero')
drawClusterPlot(peak1_peak2_color_ls, 0, 'r', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 2, 'b', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 1, 'g', output_fname_prefix)

pylab.clf()
output_fname_prefix = os.path.join(output_dir, 'cluster_plots_all')
drawClusterPlot(peak1_peak2_color_ls, 0, 'r', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 1, 'g', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 2, 'b', output_fname_prefix)
drawClusterPlot(peak1_peak2_color_ls, 3, 'y', output_fname_prefix)

snpid2peak1_peak2_color_ls = get_snpid2peak1_peak2_color_ls(curs, snpid2allele, calls_table)
ordered_snpid_ls, snpid2name = get_ordered_snpid_ls_and_snpid2name(curs, snps_table)
#new_snpid2peak1_peak2_color_ls = contrast_transform_snpid2peak1_peak2_color_ls(snpid2peak1_peak2_color_ls)
#snpid2peak1_peak2_color_ls = new_snpid2peak1_peak2_color_ls

output_fname = 'script/variation/doc/PeakDataReport/tables_figures.tex'
fig_output_dir = 'script/variation/doc/PeakDataReport/figures'
outputClusterPlotsForAllSNPs(ordered_snpid_ls, snpid2name, snpid2peak1_peak2_color_ls, fig_output_dir, output_fname)

#check only data overlapping with 2010 data
snpid2peak1_peak2_color_ls = get_snpid2peak1_peak2_color_ls(curs, snpid2allele, calls_table, data_type=3)
drawClusterPlotForOneSNP(snpid2peak1_peak2_color_ls[2], 2, '/tmp/2')

#check only data overlapping with 2010 data but also use 2010 calls to color data points
input_fname_2010 = os.path.expanduser('~/script/variation/data/2010/data_2010_x_149SNP_y1.tsv')
snp_acc2col_index, accession_id2row_index, data_matrix_2010 = read_2010_x_149SNP(input_fname_2010)
data_2010=[snp_acc2col_index, accession_id2row_index, data_matrix_2010]
snpid2peak1_peak2_color_ls = get_snpid2peak1_peak2_color_ls(curs, snpid2allele, calls_table, data_type=4, data_2010=data_2010)
drawClusterPlotForOneSNP(snpid2peak1_peak2_color_ls[2], 2, '/tmp/2')



hostname='localhost'
dbname='stock20071227'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()
	
"""

import getopt, csv
import traceback, gc
from pymodule import process_function_arguments

class SimpleCall(object):
	"""

	Examples:
		GenotypeCalling.py -i /Network/Data/250k/db/intensity/ -o /tmp/simplecalls

	Description:
		apply simple calling algorithm on intensity matrices in input_dir outputted by DB_250k2Array.
		
		The output file is named like 'array_id'_array_call.tsv. The 'array_id' is extracted from input filename.
		If a file with same name in the output directory already exists, no call or output file would be generated.
	"""
	option_default_dict = {	('i', 'input_dir',1, 'directory containing intensity matricies outputted by DB_250k2Array.py', 1, ): None,\
							('o', 'output_dir',1, 'directory to hold output calls.', 1, ): None,\
							('b', 'debug',0, 'toggle debug mode', 0, int):0,\
							('r', 'report',0, 'toggle report, more verbose stdout/stderr.', 0, int):0}
	
	"""
	2008-04-20
		option_default_dict is a dictionary for option handling, including argument_default_dict info
		the key is a tuple, ('short_option', 'long_option', has_argument, description_for_option, is_option_required, argument_type)
		argument_type is optional
	"""
	def __init__(self, **keywords):
		"""
		2008-04-08
		"""
		from pymodule import process_function_arguments, turn_option_default_dict2argument_default_dict
		argument_default_dict = turn_option_default_dict2argument_default_dict(self.option_default_dict)
		
		#argument dictionary
		self.ad = process_function_arguments(keywords, argument_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def _simpleCall(self, input_fname, array_id, output_fname):
		"""
		2008-04-08
		"""
		reader = csv.reader(open(input_fname), delimiter='\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = reader.next()
		writer_header = header[:2]
		writer_header[1] = array_id
		writer.writerow(writer_header)
		for row in reader:
			snpid = row[0]
			int_ls = row[1:]
			int_ls = map(float, int_ls)
			sense1_int, sense2_int, antisense1_int, antisense2_int = int_ls
			snpid_split = snpid.split('_')
			allele1 = snpid_split[-2]
			allele2 = snpid_split[-1]
			
			if sense1_int>sense2_int and antisense1_int>antisense2_int:
				call = allele1
			elif sense1_int<sense2_int and antisense1_int<antisense2_int:
				call = allele2
			else:
				call = 'NA'
			writer.writerow([snpid, call])
		del reader, writer
	
	def call_one_dir(self, input_dir, output_dir):
		"""
		2008-04-11
			insert 'array' into output filename to make the 'array_id' in the beginning of the filename more obvious.
			check if output_fname exists or not. if yes, skip.
		2008-04-08
		"""
		if not os.path.isdir(output_dir):
			os.makedirs(output_dir)
		file_ls = os.listdir(input_dir)
		sys.stderr.write("\n\tTotally, %d files to be processed.\n"%len(file_ls))
		for i in range(len(file_ls)):
			filename = file_ls[i]
			array_id = filename.split('_')[0]
			sys.stderr.write("%d/%d:\t%s\n"%(i+1,len(file_ls),filename))
			input_fname = os.path.join(input_dir, filename)
			output_fname = os.path.join(output_dir, '%s_array_call.tsv'%array_id)
			if os.path.isfile(output_fname):
				sys.stderr.write("Output %s already exists. Ignore.\n"%output_fname)
				continue
			self._simpleCall(input_fname, array_id, output_fname)
	
	def run(self):
		"""
		"""
		self.call_one_dir(self.input_dir, self.output_dir)
		
		
if __name__ == '__main__':
	from pymodule import process_options, generate_program_doc
	main_class = SimpleCall
	opts_dict = process_options(sys.argv, main_class.option_default_dict, error_doc=generate_program_doc(sys.argv[0], main_class.option_default_dict)+main_class.__doc__)
	
	instance = main_class(**opts_dict)
	instance.run()

