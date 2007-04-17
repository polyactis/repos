"""
2007-03-05
	add the position info given by Yan Li from Borevitz Lab
"""
def fill_snp_locus_table_with_position_info(curs, input_fname, output_table, need_commit):
	import csv
	reader = csv.reader(open(input_fname))
	snp_locus_acc_list = reader.next()[1:]
	snp_locus_position_list = reader.next()[1:]
	for i in range(len(snp_locus_acc_list)):
		snp_locus_acc = snp_locus_acc_list[i]
		snp_locus_acc = snp_locus_acc.replace('.', ' ')
		snp_locus_position = snp_locus_position_list[i]
		chromosome, position = snp_locus_position.split('-')
		curs.execute("update %s set chromosome='%s', position=%s where acc='%s'"%(output_table, chromosome, position, snp_locus_acc))
	if need_commit:
		curs.execute("end")


fill_snp_locus_table_with_position_info(curs, './script/variation/data/snp_position.csv', 'dbsnp.snp_locus', need_commit=1)

"""
2007-03-05
	show an image visualizing SNP data
"""
def display_snp_matrix(input_fname, output_fname=None, need_sort=0, need_savefig=0):
	import csv, Numeric, pylab
	reader = csv.reader(open(input_fname), delimiter='\t')
	header = reader.next()
	data_matrix = []
	for row in reader:
		data_row = row[2:]
		data_row = map(int, data_row)
		data_matrix.append(data_row)
	del reader
	data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
	if need_sort:
		data_matrix.sort()
	data_matrix = Numeric.array(data_matrix)
	
	pylab.clf()
	pylab.imshow(data_matrix, interpolation='nearest')
	pylab.colorbar()
	if need_savefig:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()

"""
2007-03-15
"""
def reformat_data_for_chris(curs, input_fname, output_fname, snp_locus_table):
	from common import number2nt
	
	sys.stderr.write("Getting snp_acc2pos ...")
	snp_acc2pos = {}
	curs.execute("select acc, chromosome, position from %s"%snp_locus_table)
	rows = curs.fetchall()
	for row in rows:
		snp_acc, chromosome, position = row
		snp_acc2pos[snp_acc] = [chromosome, position]
	sys.stderr.write("Done.\n")
	
	import csv
	reader = csv.reader(open(input_fname), delimiter='\t')
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	header = reader.next()
	new_header = [header[0]] + header[2:]
	writer.writerow(new_header)
	
	chr_ls = ['']
	pos_ls = ['']
	for snp_acc in new_header[1:]:
		chr_ls.append(snp_acc2pos[snp_acc][0])
		pos_ls.append(snp_acc2pos[snp_acc][1])
	writer.writerow(chr_ls)
	writer.writerow(pos_ls)
	
	for row in reader:
		new_row = [row[0]]
		for call in row[2:]:
			nt = number2nt[int(call)]
			if nt=='NA':	#chris wants 'N' instead of 'NA'
				nt = 'N'
			new_row.append(nt)
		writer.writerow(new_row)
	del writer

hostname='zhoudb'
dbname='graphdb'
schema = 'dbsnp'
conn, curs = db_connect(hostname, dbname, schema)
reformat_data_for_chris(curs, '../data/justin_data_filtered.csv', '../data/justin_data_filtered_for_chris.csv', 'snp_locus')

"""
2007-03-18
	the input_fname is output of MpiTrioAncestryInference.py
	for each strain, it outputs no_of_pairs, no_of_distinct_parents, avg_no_of_jumps
"""
def TrioAncestrySummary(input_fname, output_fname):
	import sys, os, csv
	from sets import Set
	sys.stderr.write("Reading trio ancestry information ...")
	reader = csv.reader(open(input_fname), delimiter='\t')
	strain_index2data = {}	#value is list of [no_of_pairs, Set_of_parents, no_of_jump_list]
	for row in reader:
		i,j,k = row[:3]
		no_of_jumps = int(row[-1])
		i = int(i)
		j = int(j)
		k = int(k)
		if k not in strain_index2data:
			strain_index2data[k] = [0, Set(), []]
		strain_index2data[k][0] += 1
		strain_index2data[k][1].add(i)
		strain_index2data[k][1].add(j)
		strain_index2data[k][2].append(no_of_jumps)
	del reader
	sys.stderr.write('done\n')
	
	sys.stderr.write("Outputting summary data ...")
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	for strain_index, data in strain_index2data.iteritems():
		no_of_distinct_parents = len(data[1])
		avg_no_of_jumps = sum(data[2])/float(data[0])
		writer.writerow([strain_index, data[0], no_of_distinct_parents, avg_no_of_jumps])
	del writer
	sys.stderr.write('done\n')

TrioAncestrySummary('./script/variation/data/justin_data_filtered.trio_ances', './script/variation/data/justin_data_filtered.trio_ances.summary')

"""
2007-03-18
"""
def draw_histogram_of_data_from_specified_column(input_fname, column=-1, no_of_bins=20):
	import pylab, csv
	data_ls = []
	reader = csv.reader(open(input_fname), delimiter='\t')
	for row in reader:
		data_ls.append(float(row[column]))
	pylab.clf()
	pylab.hist(data_ls, no_of_bins)
	pylab.show()

draw_histogram_of_data_from_specified_column('./script/variation/data/justin_data_filtered.trio_ances.summary', 1)

"""
2007-03-19
"""
def get_one_ancestry_line(trio_ances_fname, triplet):
	import csv, sys
	sys.stderr.write("Getting one ancestry line...")
	reader = csv.reader(open(trio_ances_fname), delimiter='\t')
	ancestry_ls = None
	for row in reader:
		this_triplet = row[:3]
		this_triplet = map(int, this_triplet)
		if this_triplet[0] == triplet[0] and this_triplet[1]==triplet[1] and this_triplet[2] == triplet[2]:
			ancestry_ls = row[3:-2]	#-2 to discard the last chromosome separator and no_of_jumps
			break
	sys.stderr.write("Done.\n")
	del reader
	return ancestry_ls

def DrawTrioAncestry(data_matrix_fname, trio_ances_fname, triplet=[]):
	"""
	2007-04-17 double the layers to accomodate the heterozygotes
	"""
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
	if not triplet:
		print 'Error: triplet is not specified. '
		return
	ancestry_ls = get_one_ancestry_line(trio_ances_fname, triplet)
	import Numeric
	sub_data_matrix = Numeric.zeros([12, len(ancestry_ls)])
	from common import nt2number, number2nt
	no_of_chrs = 0
	for i in range(len(ancestry_ls)):
		if ancestry_ls[i]=='||':
			no_of_chrs += 1
		else:
			#first fill up the last 6 encoding lines
			sub_data_matrix[6,i] = 1
			sub_data_matrix[7,i] = 1
			sub_data_matrix[8,i] = 2
			sub_data_matrix[9,i] = 2
			ancestry_bit = int(ancestry_ls[i])
			if ancestry_bit == 2:
				sub_data_matrix[10,i] = 1
				sub_data_matrix[11,i] = 2
			else:
				sub_data_matrix[10,i] = sub_data_matrix[11,i] = ancestry_bit+1	#the original coding is 0 or 1, but 0 is reserved for NA
			for j in range(len(triplet)):	#second deal with the 1st 3 lines with real data
				SNP_call = data_matrix[triplet[j]][i-no_of_chrs]
				if SNP_call>4:
					SNP_call_1 = nt2number[number2nt[SNP_call][0]]
					sub_data_matrix[2*j,i] = SNP_call_1
					SNP_call_2 = nt2number[number2nt[SNP_call][1]]
					sub_data_matrix[2*j+1,i] = SNP_call_2
				else:
					sub_data_matrix[2*j,i] = sub_data_matrix[2*j+1,i] = SNP_call
	import pylab
	pylab.clf()
	pylab.imshow(sub_data_matrix, interpolation='nearest')
	pylab.colorbar()
	ytick_label_ls = []
	for i in range(2):	#repeat the labels
		for j in range(len(triplet)):
			ytick_label_ls.append(strain_acc_list[triplet[j]])	#2007-04-17 double it
			ytick_label_ls.append(strain_acc_list[triplet[j]])
	ytick_loc_ls = range(12)
	ytick_loc_ls.reverse()
	pylab.yticks(ytick_loc_ls, ytick_label_ls)
	pylab.show()

"""
2007-03-20
"""
def DrawStrainSNP_NA_PercHist(data_matrix_fname, output_fname, need_savefig=0):
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
	import Numeric
	data_matrix = Numeric.array(data_matrix)
	strain_NA_perc_ls = []
	for i in range(data_matrix.shape[0]):
		strain_NA_perc_ls.append(sum(data_matrix[i,:]==0)/float(data_matrix.shape[1]))	#0 is NA
	SNP_NA_perc_ls = []
	for i in range(data_matrix.shape[1]):
		SNP_NA_perc_ls.append(sum(data_matrix[:,i]==0)/float(data_matrix.shape[0]))	#0 is NA
	import pylab,os
	base_fname = os.path.basename(data_matrix_fname)
	pylab.clf()
	pylab.hist(strain_NA_perc_ls, 20)
	pylab.title("%s Strain NA perc histogram"%base_fname)
	if need_savefig:
		pylab.savefig('%s_strain_NA_perc.eps'%output_fname, dpi=300)
		pylab.savefig('%s_strain_NA_perc.svg'%output_fname, dpi=300)
		pylab.savefig('%s_strain_NA_perc.png'%output_fname, dpi=300)
	pylab.show()
	
	pylab.clf()
	pylab.hist(SNP_NA_perc_ls, 20)
	pylab.title("%s SNP NA perc histogram"%base_fname)
	if need_savefig:
		pylab.savefig('%s_SNP_NA_perc.eps'%output_fname, dpi=300)
		pylab.savefig('%s_SNP_NA_perc.svg'%output_fname, dpi=300)
		pylab.savefig('%s_SNP_NA_perc.png'%output_fname, dpi=300)
	pylab.show()

DrawStrainSNP_NA_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)

def DrawStrain_Heterozygotes_PercHist(data_matrix_fname, output_fname, need_savefig=0):
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
	import Numeric
	data_matrix = Numeric.array(data_matrix)
	strain_Hetero_perc_ls = []
	for i in range(data_matrix.shape[0]):
		strain_Hetero_perc_ls.append(sum(data_matrix[i,:]>4)/float(data_matrix.shape[1]))	#bigger than 4 is heterozygotes
	
	import pylab,os
	base_fname = os.path.basename(data_matrix_fname)
	pylab.clf()
	pylab.hist(strain_Hetero_perc_ls, 20)
	pylab.title("%s Strain Heterozygote perc histogram"%base_fname)
	if need_savefig:
		pylab.savefig('%s_strain_hz_perc.eps'%output_fname, dpi=300)
		pylab.savefig('%s_strain_hz_perc.svg'%output_fname, dpi=300)
		pylab.savefig('%s_strain_hz_perc.png'%output_fname, dpi=300)
	pylab.show()

DrawStrain_Heterozygotes_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)

"""
2007-03-21
"""
def DrawDistanceHistogram(data_matrix_fname, output_fname, need_savefig=0):
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
	import Numeric
	distance_ls = []
	no_of_NA_pairs = 0
	no_of_strains = len(data_matrix)
	no_of_snps = len(data_matrix[0])
	distance_matrix = Numeric.zeros([no_of_strains, no_of_strains], Numeric.Float)
	for i in range(no_of_strains):
		for j in range(i+1, no_of_strains):
			no_of_valid_pairs = 0.0
			no_of_matching_pairs = 0.0
			for k in range(no_of_snps):
				if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
					no_of_valid_pairs += 1
					if data_matrix[i][k] == data_matrix[j][k]:
						no_of_matching_pairs += 1
			if no_of_valid_pairs!=0:
				distance = 1 - no_of_matching_pairs/no_of_valid_pairs
				distance_matrix[i,j] = distance_matrix[j,i] = distance
				distance_ls.append(distance)
			else:
				no_of_NA_pairs += 1
	print "out of %s pairs, %s are NA"%((no_of_strains*(no_of_strains-1))/2, no_of_NA_pairs)
	import pylab
	pylab.clf()
	pylab.hist(distance_ls, 10)
	pylab.title("Histogram of non-NA distances")
	if need_savefig:
		pylab.savefig('%s_distance_hist.eps'%output_fname, dpi=300)
		pylab.savefig('%s_distance_hist.svg'%output_fname, dpi=300)
		pylab.savefig('%s_distance_hist.png'%output_fname, dpi=300)
	pylab.show()
	return distance_matrix

distance_matrix = DrawDistanceHistogram('./script/variation/data/justin_data_filtered.csv', './script/variation/data/justin_data_filtered' , 1)

"""
2007-03-21
	The stat is d(i,j)-|d(i,k)-d(j,k)| for k is a child of i and j
	the bigger this value is, the more divergent between two parents and child is in equal distance to two parents
"""
def cal_trio_stat(trio_ances_fname, distance_matrix, output_fname, need_savefig=0, report=0):
	import csv
	import Numeric
	no_of_strains = distance_matrix.shape[0]
	triplet2trio_stat = {}
	trio_stat_ls = []
	reader = csv.reader(open(trio_ances_fname), delimiter='\t')
	counter = 0
	for row in reader:
		triplet = row[:3]
		triplet = map(int, triplet)
		i,j,k = triplet
		trio_stat = distance_matrix[i,j] - abs(distance_matrix[i,k] - distance_matrix[j,k])
		triplet2trio_stat[tuple(triplet)] = trio_stat
		trio_stat_ls.append(trio_stat)
		counter += 1
		if report and counter%3000 == 0:
			sys.stderr.write("%s%s"%('\x08'*20, counter))
	if report:
		sys.stderr.write("%s%s\n"%('\x08'*20, counter))
	del reader
	import pylab
	pylab.clf()
	pylab.hist(trio_stat_ls, 20)
	pylab.title("hist of d(i,j)-|d(i,k)-d(j,k)|")
	if need_savefig:
		pylab.savefig('%s_trio_stat_his.eps'%output_fname, dpi=300)
		pylab.savefig('%s_trio_stat_hist.svg'%output_fname, dpi=300)
		pylab.savefig('%s_trio_stat_hist.png'%output_fname, dpi=300)
	pylab.show()
	return triplet2trio_stat

triplet2trio_stat = cal_trio_stat('./script/variation/data/justin_data_filtered.trio_ances', distance_matrix, './script/variation/data/justin_data_filtered.trio_ances', need_savefig=1, report=1)

"""
2007-03-28
"""
def strip_2010_strain_info(input_fname, output_fname):
	import csv
	reader = csv.reader(open(input_fname), delimiter='\t')
	#skip the 1st two sentences
	reader.next()
	reader.next()
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	for row in reader:
		for i in range(len(row)):
			row[i] = row[i].strip()
		writer.writerow(row)
	del reader, writer

strip_2010_strain_info('./script/variation/data/2010/2010_strain_info.csv', './script/variation/data/2010/2010_strain_info_stripped.csv')

"""
2007-04-09
"""
def draw_SNP_gap_histogram(curs, snp_locus_table, output_fname, need_savefig=0):
	SNP_gap_ls = []
	prev_chromosome = None
	prev_position = None
	curs.execute("select chromosome, position from %s order by chromosome, position"%snp_locus_table)
	rows = curs.fetchall()
	for row in rows:
		chromosome, position = row
		if prev_chromosome==None or prev_position==None:
			prev_chromosome = chromosome
			prev_position = position
		elif chromosome==prev_chromosome:
			SNP_gap_ls.append(position-prev_position)
		prev_chromosome = chromosome
		prev_position = position
	import pylab
	pylab.clf()
	pylab.hist(SNP_gap_ls, 20)
	pylab.title("hist of SNP gap")
	if need_savefig:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()
	return SNP_gap_ls

#2007-03-05 common codes to initiate database connection
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))

sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
from codense.common import db_connect, form_schema_tables
hostname='dl324b-1'
dbname='yhdb'
schema = 'dbsnp'
hostname='zhoudb'
dbname='graphdb'
schema = 'dbsnp'
conn, curs = db_connect(hostname, dbname, schema)
