#!/usr/bin/env python

import os,sys

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


"""
fill_snp_locus_table_with_position_info(curs, './script/variation/data/snp_position.csv', 'dbsnp.snp_locus', need_commit=1)
"""

"""
2007-03-05
	show an image visualizing SNP data
2007-06-05
	set aspect='auto' in imshow(), the default (pylab.image.rcParams['image.aspect'])='equal', which is bad
"""
def display_snp_matrix(input_fname, output_fname=None, need_sort=0, need_savefig=0, xlabel='', ylabel=''):
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
	pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
	pylab.colorbar()
	if xlabel:
		pylab.xticks([data_matrix.shape[1]/2], [xlabel])
	if ylabel:
		pylab.yticks([data_matrix.shape[0]/2], [ylabel])
	if need_savefig:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()

def make_snp_matrix_legend(value_ls, label_ls, output_fname=None):
	"""
	2007-10-25
		to pair with display_snp_matrix()
	"""
	import numpy, pylab
	label_ls_copy = label_ls[:]
	data_matrix = numpy.zeros([len(value_ls), 1], numpy.int)
	for i in range(len(value_ls)):
		data_matrix[i,0] = value_ls[i]
	label_ls_copy.reverse()	#pylab put the label on starting from the bottom of the picture
	pylab.clf()
	pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')
	pylab.yticks(range(len(value_ls)), label_ls_copy, fontsize=60, verticalalignment='bottom', horizontalalignment='right')
	pylab.xticks([],[])
	if output_fname:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()

"""
2007-10-25
value_ls = range(3)
label_ls = ['NA', 'A', 'C']
make_snp_matrix_legend(value_ls, label_ls)

"""


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

"""
hostname='zhoudb'
dbname='graphdb'
schema = 'dbsnp'
conn, curs = db_connect(hostname, dbname, schema)
reformat_data_for_chris(curs, '../data/justin_data_filtered.csv', '../data/justin_data_filtered_for_chris.csv', 'snp_locus')
"""

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

"""
TrioAncestrySummary('./script/variation/data/justin_data_filtered.trio_ances', './script/variation/data/justin_data_filtered.trio_ances.summary')
"""

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

"""
draw_histogram_of_data_from_specified_column('./script/variation/data/justin_data_filtered.trio_ances.summary', 1)
"""

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

"""
DrawStrainSNP_NA_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)
"""

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
	return strain_Hetero_perc_ls

"""
strain_Hetero_perc_ls = DrawStrain_Heterozygotes_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)
"""

"""
2007-03-21
2007-09-24 increase the #bins of histogram to 40
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
	pylab.hist(distance_ls, 40)
	pylab.title("Histogram of non-NA distances")
	if need_savefig:
		pylab.savefig('%s_distance_hist.eps'%output_fname, dpi=300)
		pylab.savefig('%s_distance_hist.svg'%output_fname, dpi=300)
		pylab.savefig('%s_distance_hist.png'%output_fname, dpi=300)
	pylab.show()
	return distance_matrix

"""
distance_matrix = DrawDistanceHistogram('./script/variation/data/justin_data_filtered.csv', './script/variation/data/justin_data_filtered' , 1)
"""

"""
2007-03-21
	The stat is d(i,j)-|d(i,k)-d(j,k)| for k is a child of i and j
	the bigger this value is, the more divergent between two parents and child is in equal distance to two parents
"""
def cal_trio_stat(trio_ances_fname, distance_matrix, output_fname, pic_fname_prefix, need_savefig=0, report=0):
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
	sys.stderr.write("outputting trio_stat...")
	writer = csv.writer(open(output_fname,'w'), delimiter='\t')
	for triplet, trio_stat in triplet2trio_stat.iteritems():
		writer.writerow(list(triplet)+[trio_stat])
	del writer
	sys.stderr.write("Done\n")
	import pylab
	pylab.clf()
	pylab.hist(trio_stat_ls, 20)
	pylab.title("hist of d(i,j)-|d(i,k)-d(j,k)|")
	if need_savefig:
		pylab.savefig('%s_trio_stat_his.eps'%pic_fname_prefix, dpi=300)
		pylab.savefig('%s_trio_stat_hist.svg'%pic_fname_prefix, dpi=300)
		pylab.savefig('%s_trio_stat_hist.png'%pic_fname_prefix, dpi=300)
	pylab.show()
	return triplet2trio_stat

"""
triplet2trio_stat = cal_trio_stat('./script/variation/data/justin_data_filtered.trio_ances', distance_matrix, './script/variation/data/justin_data_filtered.trio_ances.trio_stat', './script/variation/data/justin_data_filtered.trio_ances', need_savefig=1, report=1)
"""

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

"""
strip_2010_strain_info('./script/variation/data/2010/2010_strain_info.csv', './script/variation/data/2010/2010_strain_info_stripped.csv')
"""

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


"""
2007-04-30
"""
def blast_snp_segment(curs, snp_locus_table, output_fname, database_fname, flanking_seq_length=12, max_no_of_hits_to_be_outputted=3, blast_bin_path=os.path.expanduser('~/bin/blast/bin/blastall'), annot_assembly_table='sequence.annot_assembly', \
	raw_sequence_table='sequence.raw_sequence', tmp_blast_infname='/tmp/blast_input'):
	import sys
	from Bio.Blast import NCBIXML, NCBIStandalone
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	from codense.common import get_sequence_segment
	curs.execute("select s.id, s.acc, s.chromosome, s.position, a.gi from %s s, %s a\
		where a.tax_id=s.tax_id and a.chromosome=s.chromosome"%(snp_locus_table, annot_assembly_table))
	rows = curs.fetchall()
	counter = 0
	outf = open(output_fname, 'w')
	for row in rows:
		snp_locus_id, acc, chromosome, position, genomic_gi = row
		genomic_start = position-flanking_seq_length
		genomic_stop = position + flanking_seq_length
		seq = get_sequence_segment(curs, genomic_gi, genomic_start, genomic_stop, annot_assembly_table, raw_sequence_table)
		if seq:
			for candidate_snp in ['A', 'C', 'G', 'T']:
				seq = seq[:flanking_seq_length] + candidate_snp + seq[flanking_seq_length+1:]	#str doesn't support item assignment
				inf = open(tmp_blast_infname, 'w')
				inf.write('>%s(snp_locus_id=%s) candidate_snp=%s\n'%(acc, snp_locus_id, candidate_snp))
				inf.write('%s\n'%seq)
				del inf
				result_handle, error_info = NCBIStandalone.blastall(blast_bin_path, "blastn", database_fname, tmp_blast_infname, align_view=7)	#align_view=7 toggles the xml output
				blast_parser = NCBIXML.BlastParser()	#the parser has to be re-instanced every time otherwise the blast records parsed in the previous round stay in the parser
				blast_records = blast_parser.parse(result_handle)
				del blast_parser
				no_of_hits = min(max_no_of_hits_to_be_outputted, len(blast_records.alignments))
				outf.write("%s (id=%s) candidate_snp=%s\n"%(acc, snp_locus_id, candidate_snp))
				outf.write("query sequence = %s\n"%seq)
				for i in range(no_of_hits):
					outf.write("\t%s\n"%(blast_records.alignments[i].title))
					for hsp in blast_records.alignments[i].hsps:
						outf.write("\t\tquery_start=%s, sbjct_start=%s, frame=%s, identities=%s/%s E=%s\n"%(hsp.query_start, hsp.sbjct_start, hsp.frame, hsp.identities, 2*flanking_seq_length+1, hsp.expect))
						outf.write("\t\t%s\n"%hsp.query)
						outf.write("\t\t%s\n"%hsp.match)
						outf.write("\t\t%s\n"%hsp.sbjct)
						outf.write("\n")
					outf.write("\n")
				outf.write("\n")
		sys.stderr.write("%s%s"%('\x08'*10, counter))
		counter += 1
	del outf

"""
my_blast_db = os.path.expanduser("~/bin/blast/db/Arabidopsis_thaliana.main_genome.fasta")
blast_snp_segment(curs, 'snp_locus', './blast_149snps_vs_thaliana_len_25.txt', my_blast_db)


my_blast_db = os.path.expanduser("~/bin/blast/db/Arabidopsis_lyrata.main_genome.scaffolds.fasta")
blast_snp_segment(curs, 'snp_locus', './blast_149snps_vs_lyrata_len_51.txt', my_blast_db, flanking_seq_length=25)
"""

"""
2007-04-30
"""
def fill_snp_locus_table_with_25mer_thaliana_call(curs, snp_locus_table, annot_assembly_table='sequence.annot_assembly', \
	raw_sequence_table='sequence.raw_sequence', need_commit=0):
	import sys
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	from codense.common import get_sequence_segment
	curs.execute("select s.id, s.acc, s.chromosome, s.position, a.gi from %s s, %s a\
		where a.tax_id=s.tax_id and a.chromosome=s.chromosome"%(snp_locus_table, annot_assembly_table))
	rows = curs.fetchall()
	counter = 0
	for row in rows:
		snp_locus_id, acc, chromosome, position, genomic_gi = row
		genomic_start = position - 12
		genomic_stop = position + 12
		seq = get_sequence_segment(curs, genomic_gi, genomic_start, genomic_stop, annot_assembly_table, raw_sequence_table)
		if seq:
			thaliana_call = seq[12]
			curs.execute("update %s set thaliana_call='%s', flanking_25mer='%s' where id=%s"%(snp_locus_table, thaliana_call, seq, snp_locus_id))
		sys.stderr.write("%s%s"%('\x08'*10, counter))
		counter += 1
	if need_commit:
		curs.execute("end")

"""
fill_snp_locus_table_with_25mer_thaliana_call(curs, 'dbsnp.snp_locus', need_commit=1)
"""

#2007-06-17 function to calculate the great circle distance of a sphere given latitude and longitude
# http://en.wikipedia.org/wiki/Great-circle_distance
def cal_great_circle_distance(lat1, lon1, lat2, lon2, earth_radius=6372.795):
	import math
	lat1_rad = lat1*math.pi/180
	lon1_rad = lon1*math.pi/180
	lat2_rad = lat2*math.pi/180
	lon2_rad = lon2*math.pi/180
	long_diff = abs(lon1_rad-lon2_rad)
	sin_lat1 = math.sin(lat1_rad)
	cos_lat1 = math.cos(lat1_rad)
	sin_lat2 = math.sin(lat2_rad)
	cos_lat2 = math.cos(lat2_rad)
	spheric_angular_diff = math.atan2(math.sqrt(math.pow(cos_lat2*math.sin(long_diff),2) + math.pow(cos_lat1*sin_lat2-sin_lat1*cos_lat2*math.cos(long_diff),2)), sin_lat1*sin_lat2+cos_lat1*cos_lat2*math.cos(long_diff))
	return earth_radius*spheric_angular_diff

#2007-06-17 draw location of all strains onto a map
def test_draw_strain_location(pic_area=[-180,-90,180,90]):
	import pylab
	from matplotlib.toolkits.basemap import Basemap
	pylab.clf()
	fig = pylab.figure()
	m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
            resolution='l',projection='mill')
	#m.drawcoastlines()
	m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,1,1])
	m.drawmeridians(pylab.arange(-180,180,60), labels=[1,1,1,1])
	m.fillcontinents()
	m.drawcountries()
	m.drawstates()
	
	import MySQLdb
	db = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
	c = db.cursor()
	c.execute("""select latitude, longitude from ecotype where latitude is not null and longitude is not null""")
	lat_lon_ls = c.fetchall()
	euc_coord1_ls = []
	euc_coord2_ls = []
	for lat, lon in lat_lon_ls:
		euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
		euc_coord1_ls.append(euc_coord1)
		euc_coord2_ls.append(euc_coord2)
	m.scatter(euc_coord1_ls, euc_coord2_ls, 25, marker='o', zorder=10)
	pylab.title("worldwide distribution of %s strains"%(len(lat_lon_ls)))
	pylab.show()

#2007-06-17 calculate pairwise distance among strains by calling cal_great_circle_distance()
def cal_pairwise_distance_of_strains():
	import MySQLdb, pylab
	db = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
	c = db.cursor()
	c.execute("""select latitude, longitude from ecotype where latitude is not null and longitude is not null""")
	lat_lon_ls = c.fetchall()
	distance_ls = []
	no_of_sites = len(lat_lon_ls)
	for i in range(no_of_sites):
		for j in range(i+1, no_of_sites):
			distance_ls.append(cal_great_circle_distance(lat_lon_ls[i][0], lat_lon_ls[i][1], lat_lon_ls[j][0], lat_lon_ls[j][1]))
	pylab.clf()
	pylab.hist(distance_ls, 20)
	pylab.title("histogram of distances between %s strains(km)"%no_of_sites)
	pylab.show()

"""
#2007-06-18, calll cal_great_circle_distance()
import MySQLdb
conn = MySQLdb.connect(db="stock",host='natural.uchicago.edu', user='iamhere', passwd='iamhereatusc')
cursor = conn.cursor()
cursor.execute("select latitude, longitude from ecotype where latitude is not null and longitude is not null")
lat_lon_ls = cursor.fetchall()
"""


"""
2007-09-13
following functions to investigate inter/intra-population identity
"""
def construct_identity_pair_ls(input_fname):
	"""
	2007-09-13
		input_fname is a data frame file with snps coded in integers. 1st row is the ecotype id in table ecotype
	"""
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	no_of_strains = len(strain_acc_list)
	no_of_snps = len(header)-2
	identity_pair_ls = []
	for i in range(no_of_strains):
		for j in range(i+1, no_of_strains):
			no_of_same_cols = 0
			for k in range(no_of_snps):
				if data_matrix[i][k] == data_matrix[j][k] or data_matrix[i][k]==0 or data_matrix[j][k]==0:
					no_of_same_cols += 1
			if no_of_same_cols == no_of_snps:
				identity_pair_ls.append([strain_acc_list[i], strain_acc_list[j]])
	return identity_pair_ls


def construct_graph_out_of_identity_pair(identity_pair_ls):
	"""
	2007-09-18
		construct a graph based on identity pairs in identity_pair_ls
		identity_pair_ls's entries will be converted to integers
	"""
	sys.stderr.write("Constructing graph out of identity_pair_ls ...")
	import networkx as nx
	g = nx.XGraph()
	for identity_pair in identity_pair_ls:
		identity_pair = map(int, identity_pair)
		g.add_edge(identity_pair[0], identity_pair[1], 1)
	sys.stderr.write("Done.\n")
	return g

def get_singleton_strain_id_ls(g, input_fname):
	"""
	2007-10-01
		node in g is in ecotype id format
		so strain_acc (1st column) in input_fname should be integer (ecotype id)
	"""
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	sys.stderr.write("Getting singleton_strain_id_ls ...")
	from sets import Set
	graph_node_set = Set(g.nodes())
	singleton_strain_id_ls = []
	for strain_acc in strain_acc_list:
		strain_acc = int(strain_acc)
		if strain_acc not in graph_node_set:
			singleton_strain_id_ls.append(strain_acc)
	sys.stderr.write("Done.\n")
	return singleton_strain_id_ls

def construct_site_graph_out_of_strain_graph(g, ecotypeid2pos, node_weight_by_no_of_nodes=0):
	"""
	2007-09-20
		#reduce the drawing by omitting identical edges (linking same GPS positions), construct a graph whose node is (lat,lon) and then plot this graph.
	2007-10-01
		add node_weight_by_no_of_nodes
	"""
	sys.stderr.write("Constructing site graph out of strain graph...")
	site2pos = {}
	site2weight = {}
	import networkx as nx
	site_g = nx.XGraph()
	for e in g.edges():
		pos1 = (round(ecotypeid2pos[e[0]][0],2), round(ecotypeid2pos[e[0]][1],2))
		pos2 = (round(ecotypeid2pos[e[1]][0],2), round(ecotypeid2pos[e[1]][1],2))
		if not node_weight_by_no_of_nodes:	#2007-10-01
			if pos1 not in site2weight:
				site2weight[pos1] = 0
			site2weight[pos1] += 1
			if pos2 not in site2weight:
				site2weight[pos2] = 0
			site2weight[pos2] += 1
		if pos1==pos2:
			site_g.add_node(pos1)
		else:
			if not site_g.has_edge(pos1, pos2):
				site_g.add_edge(pos1, pos2, 0)
			site_g.adj[pos1][pos2] += 1
			site_g.adj[pos2][pos1] += 1
	if node_weight_by_no_of_nodes:	#2007-10-01
		for n in g.nodes():
			pos = (round(ecotypeid2pos[n][0],2), round(ecotypeid2pos[n][1],2))
			if pos not in site2weight:
				site2weight[pos] = 0
			site2weight[pos] += 1
	for n in site_g:
		site2pos[n] = n
	sys.stderr.write("Done.\n")
	return site_g, site2weight, site2pos

def collapseStrainsWithSamePos(ecotypeid_ls, ecotypeid2pos):
	"""
	2007-10-09
	"""
	sys.stderr.write("Collapsing Strains with same GPS ...")
	site2pos = {}
	site2weight = {}
	for ecotypeid in ecotypeid_ls:
		pos = (round(ecotypeid2pos[ecotypeid][0],2), round(ecotypeid2pos[ecotypeid][1],2))
		if pos not in site2weight:
			site2weight[pos] = 0
			site2pos[pos] = pos
		site2weight[pos] += 1
	sys.stderr.write("Done.\n")
	return site2weight, site2pos


def draw_strains_on_map(ecotypeid_ls, ecotypeid2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None, label_ls=None, need_collapse_strains_with_same_pos=0):
	"""
	2007-10-09
	2007-11-08
		add label_ls, need_collapse_strains_with_same_pos
	"""
	import os, sys
	sys.stderr.write("Drawing strains on a map ...\n")
	import pylab, math
	from matplotlib.toolkits.basemap import Basemap
	pylab.clf()
	fig = pylab.figure()
	fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
	m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
	resolution='l',projection='mill')
	
	
	sys.stderr.write("\tDrawing nodes ...")
	euc_coord1_ls = []
	euc_coord2_ls = []
	diameter_ls = []
	ax=pylab.gca()
	if need_collapse_strains_with_same_pos:
		site2weight, site2pos = collapseStrainsWithSamePos(ecotypeid_ls, ecotypeid2pos)
		for n in site2pos:
			lat, lon = site2pos[n]
			euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
			diameter_ls.append(8*math.sqrt(site2weight[n]))
	else:
		for i in range(len(ecotypeid_ls)):
			ecotypeid = ecotypeid_ls[i]
			lat, lon = ecotypeid2pos[ecotypeid]
			euc_coord1, euc_coord2 = m(lon, lat)
			euc_coord1_ls.append(euc_coord1)
			euc_coord2_ls.append(euc_coord2)
			if label_ls:
				ax.text(euc_coord1, euc_coord2, str(label_ls[i]), size=2, alpha=0.5, horizontalalignment='center', verticalalignment='center', zorder=12)
			diameter_ls.append(3)
	import numpy
	#diameter_ls = numpy.array(diameter_ls)
	m.scatter(euc_coord1_ls, euc_coord2_ls, diameter_ls, marker='o', color='r', alpha=0.4, zorder=10, faceted=False)
	sys.stderr.write("Done.\n")
	
	#m.drawcoastlines()
	m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])
	m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
	m.fillcontinents()
	m.drawcountries()
	m.drawstates()
	
	pylab.title(pic_title)
	if output_fname_prefix:
		pylab.savefig('%s.eps'%output_fname_prefix, dpi=600)
		#pylab.savefig('%s.svg'%output_fname_prefix, dpi=600)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=600)
	sys.stderr.write("Done.\n")

def display_matrix_of_component(input_fname, ecotypeid_ls, ecotypeid2pos, output_fname=None, need_sort=0, need_savefig=0):
	"""
	2007-09-20
		display the data from that component
		
	"""
	import csv, Numeric, pylab
	cc_ecotypeid_pos = []
	for ecotypeid in ecotypeid_ls:
		cc_ecotypeid_pos.append(ecotypeid2pos[ecotypeid])
	import numpy
	argsort_index = numpy.argsort(cc_ecotypeid_pos, 0)	#watch it's two dimensional
	ecotypeid2row_index = {}
	ytick_label_ls = []
	cc_size = len(ecotypeid_ls)
	for i in range(cc_size):
		ecotypeid_ls_index = argsort_index[i][1]	#sort based on longitude
		ecotypeid = ecotypeid_ls[ecotypeid_ls_index]
		ecotypeid2row_index[ecotypeid] = i
		ytick_label_ls.append('%s (%.2f, %.2f)'%(ecotypeid, ecotypeid2pos[ecotypeid][0], ecotypeid2pos[ecotypeid][1]))
	reader = csv.reader(open(input_fname), delimiter='\t')
	header = reader.next()
	data_matrix = [0]*cc_size
	for row in reader:
		ecotypeid = int(row[0])
		if ecotypeid in ecotypeid2row_index:
			data_row = row[2:]
			data_row = map(int, data_row)
			data_matrix[ecotypeid2row_index[ecotypeid]] = data_row
	del reader
	data_matrix.reverse()	#2007-03-06 reverse() due to the imshow()'s y axis starting from bottom
	if need_sort:
		data_matrix.sort()
	data_matrix = Numeric.array(data_matrix)
	
	pylab.clf()
	pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
	pylab.colorbar()
	pylab.yticks(range(cc_size), ytick_label_ls)
	if need_savefig:
		pylab.savefig('%s.eps'%output_fname, dpi=300)
		pylab.savefig('%s.svg'%output_fname, dpi=300)
		pylab.savefig('%s.png'%output_fname, dpi=300)
	pylab.show()

def get_ecotypeid2pos(curs, ecotype_table, with_data_affiliated=0):
	"""
	2007-09-18
	2007-10-09
		add with_data_affiliated
	"""
	sys.stderr.write("Getting ecotypeid2pos from %s ..."%ecotype_table)
	ecotypeid2pos = {}
	if with_data_affiliated:
		curs.execute("select distinct e.id, e.latitude, e.longitude from %s e, calls c where e.latitude is not null and e.longitude is not null and e.id=c.ecotypeid"%ecotype_table)
	else:
		curs.execute("select id, latitude, longitude from %s where latitude is not null and longitude is not null"%ecotype_table)
	rows = curs.fetchall()
	for row in rows:
		ecotypeid, latitude, longitude = row
		ecotypeid2pos[ecotypeid] = (latitude, longitude)
	sys.stderr.write("Done.\n")
	return ecotypeid2pos

def get_popid2pos_size_and_ecotypeid2pop_id(curs, popid2ecotypeid_table, selected=0):
	"""
	2007-09-13
	2007-09-24
		doesn't have to be selected
	"""
	sys.stderr.write("Getting popid2pos_size & ecotypeid2pop_id ...")
	popid2pos_size = {}
	ecotypeid2pop_id = {}
	curs.execute("select popid, latitude, longitude, ecotypeid from %s where selected=%s"%(popid2ecotypeid_table, selected))
	rows = curs.fetchall()
	for row in rows:
		popid, latitude, longitude, ecotypeid = row
		if popid not in popid2pos_size:
			popid2pos_size[popid] = [None, 0]
		if popid2pos_size[popid][0]==None:
			popid2pos_size[popid][0] = (latitude, longitude)
		popid2pos_size[popid][1] += 1
		ecotypeid2pop_id[ecotypeid] = popid
	sys.stderr.write("Done.\n")
	return popid2pos_size, ecotypeid2pop_id

def construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls):
	"""
	2007-09-18
	"""
	import os, sys
	popid2pos_size, ecotypeid2pop_id = get_popid2pos_size_and_ecotypeid2pop_id(curs, popid2ecotypeid_table)
	
	no_of_inter_pop_identities = 0
	no_of_valid_identities = 0	#valid means both strains belong to some population
	import networkx as nx
	g = nx.XGraph()
	popid2intra_pop_connections = {}
	for e in identity_pair_ls:
		e = map(int, e)
		popid1 = ecotypeid2pop_id.get(e[0])
		popid2 = ecotypeid2pop_id.get(e[1])
		if popid1 and popid2:	#both are classified into populations
			no_of_valid_identities += 1
			if popid1 not in popid2intra_pop_connections:
				popid2intra_pop_connections[popid1] = 0
				g.add_node(popid1)
			if popid2 not in popid2intra_pop_connections:
				popid2intra_pop_connections[popid2] = 0
				g.add_node(popid2)
			if popid1!=popid2:
				if not g.has_edge(popid1, popid2):
					g.add_edge(popid1, popid2, 0)
				g.adj[popid1][popid2] += 1	#increase the edge weight
				g.adj[popid2][popid1] += 1
				no_of_inter_pop_identities += 1
			else:
				popid2intra_pop_connections[popid1] += 1
	print '%s/%s=%s inter population identities'%(no_of_inter_pop_identities, no_of_valid_identities, float(no_of_inter_pop_identities)/no_of_valid_identities)
	return g, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id

def get_pic_area(pos_ls):
	"""
	2007-10-02
		get a good pic_area
	"""
	min_lat = pos_ls[0][0]
	max_lat = pos_ls[0][0]
	min_lon = pos_ls[0][1]
	max_lon = pos_ls[0][1]
	for i in range(1, len(pos_ls)):
		pos = pos_ls[i]
		if pos[0]<min_lat:
			min_lat = pos[0]
		if pos[0]>max_lat:
			max_lat = pos[0]
		if pos[1]<min_lon:
			min_lon = pos[1]
		if pos[1]>max_lon:
			max_lon = pos[1]
	pic_area = [-180, -90, 180, 90]	#this is the boundary
	if min_lon-5>pic_area[0]:
		pic_area[0] = min_lon-5
	if min_lat-5>pic_area[1]:
		pic_area[1] = min_lat-5
	if max_lon+5<pic_area[2]:
		pic_area[2] = max_lon+5
	if max_lat+5<pic_area[3]:
		pic_area[3] = max_lat+5
	return pic_area

def draw_graph_on_map(g, node2weight, node2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None, need_draw_edge=0):
	"""
	2007-09-13
		identity_pair_ls is a list of pairs of strains (ecotype id as in table ecotype)
	2007-10-08
		correct a bug in 4*diameter_ls, diameter_ls has to be converted to array first.
		sqrt the node weight, 8 times the original weight
	"""
	import os, sys
	sys.stderr.write("Drawing graph on a map ...\n")
	import pylab, math
	from matplotlib.toolkits.basemap import Basemap
	pylab.clf()
	fig = pylab.figure()
	fig.add_axes([0.05,0.05,0.9,0.9])	#[left, bottom, width, height]
	m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
	resolution='l',projection='mill')
	
	sys.stderr.write("\tDrawing nodes ...")
	euc_coord1_ls = []
	euc_coord2_ls = []
	diameter_ls = []
	for n in g.nodes():
		lat, lon = node2pos[n]
		euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
		euc_coord1_ls.append(euc_coord1)
		euc_coord2_ls.append(euc_coord2)
		diameter_ls.append(math.sqrt(node2weight[n]))
	import numpy
	diameter_ls = numpy.array(diameter_ls)
	m.scatter(euc_coord1_ls, euc_coord2_ls, 8*diameter_ls, marker='o', color='r', alpha=0.4, zorder=12, faceted=False)
	sys.stderr.write("Done.\n")
	
	if need_draw_edge:
		sys.stderr.write("\tDrawing edges ...")
		ax=pylab.gca()
		for popid1, popid2, no_of_connections in g.edges():
			lat1, lon1 = node2pos[popid1]
			lat2, lon2 = node2pos[popid2]
			x1, y1 = m(lon1, lat1)
			x2, y2 = m(lon2, lat2)
			ax.plot([x1,x2],[y1,y2], 'g', linewidth=math.log(no_of_connections+1)/2, alpha=0.2, zorder=10)
		sys.stderr.write("Done.\n")
	
	#m.drawcoastlines()
	m.drawparallels(pylab.arange(-90,90,30), labels=[1,1,0,1])
	m.drawmeridians(pylab.arange(-180,180,30), labels=[1,1,0,1])
	m.fillcontinents()
	m.drawcountries()
	m.drawstates()
	
	pylab.title(pic_title)
	if output_fname_prefix:
		pylab.savefig('%s_site_network.eps'%output_fname_prefix, dpi=600)
		pylab.savefig('%s_site_network.svg'%output_fname_prefix, dpi=600)
		pylab.savefig('%s_site_network.png'%output_fname_prefix, dpi=600)
	sys.stderr.write("Done.\n")

def check_cross_ocean_components(g, cc, ecotypeid2pos, longi_watershed=-30):
	from sets  import Set
	cross_ocean_cc2lon_pair_set = {}
	for i in range(len(cc)):
		gs = g.subgraph(cc[i])
		for e in gs.edges():
			l1 = ecotypeid2pos[e[0]][1]
			l2 = ecotypeid2pos[e[1]][1]
			if (l1<longi_watershed and l2>longi_watershed) or (l1>longi_watershed and l2<longi_watershed):
				if i not in cross_ocean_cc2lon_pair_set:
					cross_ocean_cc2lon_pair_set[i] = Set()
				cross_ocean_cc2lon_pair_set[i].add((l1, l2))
	return cross_ocean_cc2lon_pair_set

def get_clique_ls_indexed_by_cc_number(g):
	"""
	2007-10-01 each component is partitioned into cliques
		the resulting clique_ls is put in the order of connected components
		like clique_ls_indexed_by_cc_number[0] is clique_ls from connected component 0
	"""
	sys.stderr.write("Getting clique_ls_indexed_by_cc_number ...\n")
	import networkx as nx
	from pymodule.PartitionGraphIntoCliques import PartitionGraphIntoCliques
	PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(algorithm_type=2)
	g_cc_ls = nx.connected_components(g)
	clique_ls_indexed_by_cc_number = []
	for i in range(len(g_cc_ls)):
		sys.stderr.write("%s\t%s"%('\x08'*20, i))
		cc_g = g.subgraph(g_cc_ls[i])
		PartitionGraphIntoCliques_ins.partition(cc_g.copy())
		clique_ls_indexed_by_cc_number.append(PartitionGraphIntoCliques_ins.clique_ls)
	sys.stderr.write("Done.\n")
	return g_cc_ls, clique_ls_indexed_by_cc_number

def get_haplotype_size2number(clique_ls_indexed_by_cc_number):
	"""
	2007-10-02
	"""
	sys.stderr.write("Getting haplotype_size2number ...")
	haplotype_size2number = {}
	for i in range(len(clique_ls_indexed_by_cc_number)):
		for clique in clique_ls_indexed_by_cc_number[i]:
			haplotype_size = len(clique)
			if haplotype_size not in haplotype_size2number:
				haplotype_size2number[haplotype_size] = 0
			haplotype_size2number[haplotype_size] += 1
	sys.stderr.write("Done.\n")
	return haplotype_size2number

def get_ordered_clique_ls(clique_ls_indexed_by_cc_number):
	"""
	2007-10-02
		return cliques with size in descending order.
	"""
	sys.stderr.write("Getting ordered_clique_ls ...")
	clique_ls = []
	ordered_clique_ls = []
	import numpy
	for i in range(len(clique_ls_indexed_by_cc_number)):
		for clique in clique_ls_indexed_by_cc_number[i]:
			clique_ls.append(clique)
	clique_ls_argsort = list(numpy.argsort(map(len, clique_ls)))
	clique_ls_argsort.reverse()	#from big to small
	for index in clique_ls_argsort:
		ordered_clique_ls.append(clique_ls[index])
	sys.stderr.write("Done.\n")
	return ordered_clique_ls

def get_strains_within_longitude_span(ecotypeid_ls, ecotypeid2pos, longitude_span=[-12, 50]):
	"""
	2007-10-04
		check the map roughly
		[-12,50] is europe
		[-130,-60] is north america
		[60,90] is central asia
		[130, 150] is japan
		
	"""
	no_of_strains = 0
	strain_ls = []
	for ecotypeid in ecotypeid_ls:
		pos = ecotypeid2pos[ecotypeid]
		if pos[1]>=longitude_span[0] and pos[1]<=longitude_span[1]:
			no_of_strains += 1
			strain_ls.append(ecotypeid)
	print "#strains with %s is %s"%(longitude_span, no_of_strains)
	return strain_ls

def draw_histogram_of_haplotype_size(haplotype_size2number, output_fname_prefix=None):
	"""
	2007-10-02
	"""
	import pylab
	pylab.clf()
	haplotype_size_ls = []
	count_ls = []
	for haplotype_size, count in haplotype_size2number.iteritems():
		haplotype_size_ls.append(haplotype_size)
		count_ls.append(count)
	pylab.title("haplotype_size bar chart")
	pylab.xlabel("haplotype_size")
	pylab.ylabel("count")
	pylab.bar(haplotype_size_ls, count_ls)
	if output_fname_prefix:
		pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
	pylab.show()

def calGenoDistAndGeoDist(data_matrix_fname, ecotypeid2pos, longitude_span=[-180, 180], dist_range=[0.0,1.0]):
	"""
	2008-01-30
		add dist_range to restrict the output pairs within that genotype distance
	2008-01-29
		return strain_pair_ls as well
		add a general cutoff
	2007-10-03
		plot to see relationship between genotype distance and geographic distance
		-cal_great_circle_distance()
	2007-10-04
		add longitude_span
	"""
	sys.stderr.write("Calculating geno, geo distance ...\n")
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
	import Numeric
	strain_pair_ls = []
	geno_dist_ls = []
	geo_dist_ls = []
	no_of_NA_pairs = 0
	no_of_strains = len(data_matrix)
	no_of_snps = len(data_matrix[0])
	no_of_pairs = 0
	no_of_included_pairs = 0
	for i in range(no_of_strains):
		for j in range(i+1, no_of_strains):
			pos1 = ecotypeid2pos[int(strain_acc_list[i])]
			pos2 = ecotypeid2pos[int(strain_acc_list[j])]
			if pos1[1]>=longitude_span[0] and pos1[1]<=longitude_span[1] and pos2[1]>=longitude_span[0] and pos2[1]<=longitude_span[1]:
				no_of_pairs += 1
				no_of_valid_pairs = 0.0
				no_of_matching_pairs = 0.0
				for k in range(no_of_snps):
					if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
						no_of_valid_pairs += 1
						if data_matrix[i][k] == data_matrix[j][k]:
							no_of_matching_pairs += 1
				if no_of_valid_pairs!=0:
					distance = 1 - no_of_matching_pairs/no_of_valid_pairs
					if distance>=dist_range[0] and distance<=dist_range[1]:
						geno_dist_ls.append(distance)
						geo_dist = cal_great_circle_distance(pos1[0], pos1[1], pos2[0], pos2[1])
						geo_dist_ls.append(geo_dist)
						strain_pair_ls.append([int(strain_acc_list[i]), int(strain_acc_list[j])])
						no_of_included_pairs += 1
				else:
					no_of_NA_pairs += 1
	print "out of %s pairs, %s are NA, %s pairs' distance is within range(%s)."%(no_of_pairs, no_of_NA_pairs, no_of_included_pairs, dist_range)
	return strain_pair_ls, geno_dist_ls, geo_dist_ls

def calGenoDistAndGeoDistBetweenTwoAreas(data_matrix_fname, ecotypeid2pos, longitude_span1=[-180, 180], longitude_span2=[-180, 180]):
	"""
	2007-10-05
		modified from calGenoDistAndGeoDist with two longitude spans
	"""
	sys.stderr.write("Calculating geno, geo distance ...\n")
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(data_matrix_fname)
	import Numeric
	geno_dist_ls = []
	geo_dist_ls = []
	no_of_NA_pairs = 0
	no_of_strains = len(data_matrix)
	no_of_snps = len(data_matrix[0])
	no_of_pairs = 0
	for i in range(no_of_strains):
		for j in range(i+1, no_of_strains):
			pos1 = ecotypeid2pos[int(strain_acc_list[i])]
			pos2 = ecotypeid2pos[int(strain_acc_list[j])]
			if (pos1[1]>=longitude_span1[0] and pos1[1]<=longitude_span1[1] and pos2[1]>=longitude_span2[0] and pos2[1]<=longitude_span2[1]) or (pos2[1]>=longitude_span1[0] and pos2[1]<=longitude_span1[1] and pos1[1]>=longitude_span2[0] and pos1[1]<=longitude_span2[1]):
				no_of_pairs += 1
				no_of_valid_pairs = 0.0
				no_of_matching_pairs = 0.0
				for k in range(no_of_snps):
					if data_matrix[i][k]!=0 and data_matrix[j][k]!=0:
						no_of_valid_pairs += 1
						if data_matrix[i][k] == data_matrix[j][k]:
							no_of_matching_pairs += 1
				if no_of_valid_pairs!=0:
					distance = 1 - no_of_matching_pairs/no_of_valid_pairs
					geno_dist_ls.append(distance)
					geo_dist = cal_great_circle_distance(pos1[0], pos1[1], pos2[0], pos2[1])
					geo_dist_ls.append(geo_dist)
				else:
					no_of_NA_pairs += 1
	print "out of %s pairs, %s are NA"%(no_of_pairs, no_of_NA_pairs)
	return geno_dist_ls, geo_dist_ls

def plotHist(ls, title=None, output_fname_prefix=None, no_of_bins=40):
	import pylab
	pylab.clf()
	pylab.hist(ls, no_of_bins)
	if title:
		pylab.title(title)
	if output_fname_prefix:
		pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
	pylab.show()

def plotXY(xls, yls, title=None, xlabel=None, ylabel=None, output_fname_prefix=None):
	import pylab
	pylab.clf()
	
	pylab.title(title)
	pylab.xlabel(xlabel)
	pylab.ylabel(ylabel)
	pylab.plot(xls, yls, '.')
	if output_fname_prefix:
		pylab.savefig('%s.eps'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=300)
		pylab.savefig('%s.png'%output_fname_prefix, dpi=300)
	pylab.show()

def outputGraphInCliquerFormat(graph, output_fname):
	"""
	2007-10-05
		to test if cliquer suffers the same problem
	"""
	of = open(output_fname, 'w')
	of.write("p edge %s %s\n"%(graph.number_of_nodes(), graph.number_of_edges()))
	node_name2index = {}
	for n in graph.nodes():
		node_name2index[n] = len(node_name2index)+1
	for e in graph.edges():
		of.write("e %s %s\n"%(node_name2index[e[0]], node_name2index[e[1]]))
	del of

def get_ecotypeid_ls_given_popid(curs, popid2ecotypeid_table, popid):
	ecotypeid_ls = []
	curs.execute("select ecotypeid from %s where popid=%s and selected=1"%(popid2ecotypeid_table, popid))
	rows = curs.fetchall()
	for row in rows:
		ecotypeid_ls.append(row[0])
	return ecotypeid_ls

"""
2008-01-25
	following functions to investigate how identity spreads against distance in different regions.
"""
def find_span_index_out_of_span_ls(pos, longitude_span_ls):
	"""
	2008-01-25
		find the index of the longitude_span in which pos drops in
	"""
	span_index = -1
	for i in range(len(longitude_span_ls)):
		if pos[1]>=longitude_span_ls[i][0] and pos[1]<=longitude_span_ls[i][1]:
			span_index = i
	return span_index

def get_span_pair_count_ls(input_fname, ecotypeid2pos, span_ls):
	"""
	2008-01-30
		calculate the number of possible pairs for a specific span type (either single span or two crossing-spans)
		for normalization purpose in group_identity_pair_according_to_region_distance()
	"""
	sys.stderr.write("Getting span_pair_count_ls ...")
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	strain_acc_list = map(int, strain_acc_list)
	#1st check how many strains in each span
	span_count_ls = [0]*len(span_ls)
	for ecotypeid in strain_acc_list:
		pos = ecotypeid2pos[ecotypeid]
		for i in range(len(span_ls)):
			if pos[1]>=span_ls[i][0] and pos[1]<=span_ls[i][1]:
				span_count_ls[i] += 1
	#2nd calculate the number of pairs in each span
	span_pair_count_ls = []
	for i in range(len(span_count_ls)):
		span_pair_count_ls.append(span_count_ls[i]*span_count_ls[i]/2.0)
	#3rd calculate the number of pairs crossing two spans
	#check all combinations of span_ls
	for i in range(len(span_ls)):
		for j in range(i+1, len(span_ls)):
			span_pair_count_ls.append(span_count_ls[i]*span_count_ls[j])
	sys.stderr.write("Done.\n")
	return span_pair_count_ls

def group_identity_pair_according_to_region_distance(strain_pair_ls, geno_dist_ls, geo_dist_ls, ecotypeid2pos, span_ls, span_label_ls, geo_distance_window_size=100, max_identity_geno_dist=0.0, span_pair_count_ls=[]):
	"""
	2008-01-30
		add span_pair_count_ls to normalize the counts
	2008-01-29
		strain_pair_ls, geno_dist_ls, geo_dist_ls replace identity_pair_ls
		use max_identity_geno_dist to determine whether to call it identity or not
	2008-01-25
		to group the identity_pair_ls into different regions and distance windows.
		geo_distance_window_size is in unit km
	"""
	sys.stderr.write("Grouping identity pairs according to region and geographic distance ...\n")
	import math
	list_of_window_no2count = []
	#1. get the row labels (label for the span type of each identity_pair)
	#2. link the span_index_pair to the index of list_of_window_no2count
	#3. initiate list_of_window_no2count
	row_label_ls =[] #region designation
	span_index_pair2list_index = {}
	for i in range(len(span_ls)):
		row_label_ls.append(span_label_ls[i])
		span_index_pair = (i,i)
		span_index_pair2list_index[span_index_pair] = i
		list_of_window_no2count.append({})
	for i in range(len(span_ls)):
		for j in range(i+1, len(span_ls)):
			row_label_ls.append('%s - %s'%(span_label_ls[i], span_label_ls[j]))
			span_index_pair2list_index[(i,j)] = len(row_label_ls)-1 #both order have the same index
			span_index_pair2list_index[(j,i)] = len(row_label_ls)-1
			list_of_window_no2count.append({})

	#record all the counts
	max_window_no = -1 #used to determine the final maxtrix shape
	for i in range(len(strain_pair_ls)):
		strain_pair = strain_pair_ls[i]
		geno_dist = geno_dist_ls[i]
		geo_dist = geo_dist_ls[i]
		if geno_dist<=max_identity_geno_dist:
			pos1 = ecotypeid2pos[strain_pair[0]]
			pos2 = ecotypeid2pos[strain_pair[1]]
			span_index1 = find_span_index_out_of_span_ls(pos1, span_ls)
			span_index2 = find_span_index_out_of_span_ls(pos2, span_ls)
			if span_index1>=0 and span_index2>=0:	#both data are in span_ls
				window_no = math.floor(geo_dist/geo_distance_window_size) #take the floor as the bin number
				if window_no>max_window_no:
					max_window_no = window_no
				list_index = span_index_pair2list_index[(span_index1, span_index2)]
				window_no2count = list_of_window_no2count[list_index]
				if window_no not in window_no2count:
					window_no2count[window_no] = 0
				window_no2count[window_no] += 1
	#transform list_of_window_no2count to matrix
	import numpy
	if max_window_no>=0: #-1 means nothing recorded
		matrix_of_counts_by_region_x_geo_dist = numpy.zeros([len(list_of_window_no2count), max_window_no+1], numpy.float) #+1 because window_no starts from 0
		for i in range(len(list_of_window_no2count)):
			for j in range(max_window_no+1):
				if j in list_of_window_no2count[i]:
					if span_pair_count_ls:
						matrix_of_counts_by_region_x_geo_dist[i][j] = list_of_window_no2count[i][j]/float(span_pair_count_ls[i])	#2008-01-30 normalize
					else:
						matrix_of_counts_by_region_x_geo_dist[i][j] = list_of_window_no2count[i][j]
	else:
		matrix_of_counts_by_region_x_geo_dist = None
	#generate the labels for the columns of matrix_of_counts_by_region_x_geo_dist
	col_label_ls = [] #geo distance window designation
	col_index_ls = [] #the X-axis value
	for i in range(max_window_no+1):
		col_index = (i+0.5)*geo_distance_window_size #right in the middle of the window/bin
		col_index_ls.append(col_index)
		col_label_ls.append('%skm'%int(col_index)) #km is unit
	sys.stderr.write("Done.\n")
	return matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls

def draw_table_as_bar_chart(matrix, row_label_ls, col_label_ls, col_index_ls=None, output_fname_prefix=None):
	"""
	2008-01-30
		bars from different rows no longer sit on one another. get them on the same ground and one after another.
	2008-01-25
		originated from http://matplotlib.sourceforge.net/screenshots/table_demo.py
	"""
	sys.stderr.write("Drawing Table as bar char ...")
	import matplotlib
	import pylab
	from pymodule.colours import get_colours
	pylab.clf()
	pylab.axes([0.1, 0.3, 0.8, 0.6])   # leave room below the axes for the table
	# Get some pastel shades for the colours
	colours = get_colours(len(row_label_ls))
	colours = list(colours)
	#colours.reverse()
	rows = len(matrix) #matrix.shape[0] will also do, but in case matrix is a 2D list
	ind = pylab.arange(len(col_label_ls))  # the x locations for the groups
	cellText = []
	width = 0.2    # the width of the bars
	#width = (col_index_ls[1]-col_index_ls[0])/2     # the width of the bars is half the window size
	yoff = pylab.array([0.0] * len(col_label_ls)) # the bottom values for stacked bar chart
	bar_ins_ls = []
	step = 0.32
	for row in xrange(rows):
	    bar_ins = pylab.bar(ind+step*row, matrix[row], width, bottom=yoff, color=colours[row])
	    bar_ins_ls.append(bar_ins)
	    #yoff = yoff + matrix[row]
	    cellText.append(matrix[row])
	# Add a table at the bottom of the axes
	colours.reverse()
	cellText.reverse()
	row_label_ls_reverse = row_label_ls[:]
	row_label_ls_reverse.reverse()
	the_table = pylab.table(cellText=cellText, rowLabels=row_label_ls_reverse, rowColours=colours, colLabels=col_label_ls)
	pylab.ylabel("Counts")
	pylab.xlabel("Geographic Distance")
	#vals = arange(0, 2500, 500)
	#pylab.yticks(vals*1000, ['%d' % val for val in vals])
	#pylab.xticks([])
	pylab.title('Counts of identity pairs by region X geo distance')
	
	pylab.legend([bar_ins[0] for bar_ins in bar_ins_ls], row_label_ls, shadow=True)
	if output_fname_prefix:
		pylab.savefig('%s_identity_bar_char_by_region_x_geo_dist.svg'%output_fname_prefix, dpi=600)
		pylab.savefig('%s_identity_bar_char_by_region_x_geo_dist.png'%output_fname_prefix, dpi=600)
	pylab.show()
	sys.stderr.write("Done.\n")

"""
from misc import *

input_fname = 'script/variation/stock20071008/data_d110_c0_5.tsv'
identity_pair_ls = construct_identity_pair_ls(input_fname)
#store it in a file as it takes a long time to generate
import cPickle
of = open('%s_identity_pair_ls'%os.path.splitext(input_fname)[0], 'w')
cPickle.dump(identity_pair_ls, of)
del of

import cPickle
inf = open('%s_identity_pair_ls'%os.path.splitext(input_fname)[0], 'r')
identity_pair_ls = cPickle.load(inf)
del inf

#codes to group strains into populations
popid2ecotypeid_table = 'popid2ecotypeid_25'
pop_iden_g, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls)

popid2pos = {}
for popid, pos_size in popid2pos_size.iteritems():
	popid2pos[popid] = pos_size[0]

import math
popid2weight = {}
for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
	popid2weight[popid] = 8*math.log(intra_pop_connections+1)

draw_graph_on_map(pop_iden_g, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map1')

##graph of strains
strain_iden_g = construct_graph_out_of_identity_pair(identity_pair_ls)
ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
ecotypeid2weight = {}
for n in strain_iden_g:
	ecotypeid2weight[n] = 1

import networkx as nx
strain_iden_g_cc = nx.connected_components(strain_iden_g)

#2007-10-09 draw all strains on map
output_fname_prefix = '%s'%os.path.splitext(input_fname)[0]
strain_map_fname_prefix = '%s_strain_map'%output_fname_prefix
draw_strains_on_map(ecotypeid2pos.keys(), ecotypeid2pos, 'map of strains',  get_pic_area(ecotypeid2pos.values()), strain_map_fname_prefix)

ecotypeid2pos_with_data = get_ecotypeid2pos(curs, 'ecotype', 1)
draw_strains_on_map(ecotypeid2pos_with_data.keys(), ecotypeid2pos_with_data, 'map of strains with snp data',  get_pic_area(ecotypeid2pos_with_data.values()), '%s_with_data'%strain_map_fname_prefix)

#draw cc of strain_iden_g
for i in range(len(strain_iden_g_cc)):
	gs = strain_iden_g.subgraph(strain_iden_g_cc[i])
	site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
	for site, weight in site2weight.iteritems():
		site2weight[site] =  8*math.log(weight+1)
	
	display_matrix_of_component(input_fname, strain_iden_g_cc[i], ecotypeid2pos, output_fname='ecotype_identity_map_cc%s'%i, need_sort=0, need_savefig=1)
	draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, pic_area=[-130,34,35,65], output_fname_prefix='ecotype_identity_cc%s'%i)


#2007-10-02 use clique to represent unique haplotype
strain_iden_g_cc, clique_ls_indexed_by_cc_number = get_clique_ls_indexed_by_cc_number(strain_iden_g)
haplotype_size2number = get_haplotype_size2number(clique_ls_indexed_by_cc_number)
singleton_strain_id_ls = get_singleton_strain_id_ls(strain_iden_g, input_fname)
haplotype_size2number[1] += len(singleton_strain_id_ls)	#need to make up the singletons lost in graph
hist_output_fname_prefix = '%s_haplotype_size_bar'%os.path.splitext(input_fname)[0]
draw_histogram_of_haplotype_size(haplotype_size2number, hist_output_fname_prefix)
ordered_clique_ls = get_ordered_clique_ls(clique_ls_indexed_by_cc_number)
for i in range(len(ordered_clique_ls)):
	gs = strain_iden_g.subgraph(ordered_clique_ls[i])
	site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos, node_weight_by_no_of_nodes=1)
	display_matrix_of_component(input_fname, ordered_clique_ls[i], ecotypeid2pos, output_fname='haplotype_%s_size%s'%(i,len(ordered_clique_ls[i])), need_sort=0, need_savefig=1)
	draw_graph_on_map(site_g, site2weight, site2pos, 'haplotype %s size %s'%(i,len(ordered_clique_ls[i])), pic_area=get_pic_area(site2pos.values()), output_fname_prefix='haplotype_%s_size%s_map'%(i,len(ordered_clique_ls[i])))

europe_lon_span = [-12,50]
norame_lon_span = [-130,-60]
centralasia_lon_span = [60,90]
japan_lon_span = [130, 150]
norame_strain_ls = get_strains_within_longitude_span(ordered_clique_ls[0], ecotypeid2pos, norame_lon_span)
japan_strain_ls = get_strains_within_longitude_span(ordered_clique_ls[0], ecotypeid2pos, japan_lon_span)


cross_atlantic_cc2lon_pair_set = check_cross_ocean_components(strain_iden_g, ordered_clique_ls, ecotypeid2pos, longi_watershed=-30)
cross_eurasia_cc2lon_pair_set = check_cross_ocean_components(strain_iden_g, ordered_clique_ls, ecotypeid2pos, longi_watershed=70)


#2007-10-03 to see whether geographic distance correlates with genotype distance
import numpy
strain_pair_ls, geno_dist_ls, geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos)

#store them into cPickle
import cPickle
of = open('%s_strain_pair_geno_dist_geo_dist_ls'%os.path.splitext(input_fname)[0], 'w')
cPickle.dump([strain_pair_ls, geno_dist_ls, geo_dist_ls], of)
del of

geo_dist_ls = numpy.array(geo_dist_ls)
geno_dist_ls = numpy.array(geno_dist_ls)

output_fname_prefix = '%s'%os.path.splitext(input_fname)[0]
geo_output_fname_prefix = '%s_geo_distance_hist'%output_fname_prefix
plotHist(geo_dist_ls, title="Histogram of non-NA geographic distances", output_fname_prefix=geo_output_fname_prefix)

import random
index_ls = range(len(geo_dist_ls))
index_selected_ls = random.sample(index_ls, 10000)
geo_dist_selected_ls = geo_dist_ls[index_selected_ls]
geno_dist_selected_ls = geno_dist_ls[index_selected_ls]

geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist'%output_fname_prefix
plotXY(geo_dist_selected_ls, geno_dist_selected_ls, title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=geno_vs_geo_output_fname_prefix)
plot_LD(geo_dist_selected_ls, geno_dist_selected_ls,  title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=geno_vs_geo_output_fname_prefix)

#2007-10-04 divide data into continents to see whether geographic distance correlates with genotype distance
europe_lon_span = [-12,50]
norame_lon_span = [-130,-60]
centralasia_lon_span = [60,90]
japan_lon_span = [130, 150]
def sample_geno_geo_correlation(geno_dist_ls, geo_dist_ls, output_fname_prefix):
	import numpy
	geo_dist_ls = numpy.array(geo_dist_ls)
	geno_dist_ls = numpy.array(geno_dist_ls)
	import random
	index_ls = range(len(geo_dist_ls))
	index_selected_ls = random.sample(index_ls, 10000)
	geo_dist_selected_ls = geo_dist_ls[index_selected_ls]
	geno_dist_selected_ls = geno_dist_ls[index_selected_ls]
	plot_LD(geo_dist_selected_ls, geno_dist_selected_ls,  title='genotype vs geographic distance', xlabel="geographic dist", ylabel='genotype dist', output_fname_prefix=output_fname_prefix)

eur_strain_pair_ls, eur_geno_dist_ls, eur_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, europe_lon_span)
eur_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_eur'%output_fname_prefix
sample_geno_geo_correlation(eur_geno_dist_ls, eur_geo_dist_ls, eur_geno_vs_geo_output_fname_prefix)

noramer_strain_pair_ls, noramer_geno_dist_ls, noramer_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, norame_lon_span)
noramer_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_noramer'%output_fname_prefix
sample_geno_geo_correlation(noramer_geno_dist_ls, noramer_geo_dist_ls, noramer_geno_vs_geo_output_fname_prefix)

eur_noramer_strain_pair_ls, eur_noramer_geno_dist_ls, eur_noramer_geo_dist_ls =  calGenoDistAndGeoDistBetweenTwoAreas(input_fname, ecotypeid2pos, europe_lon_span, norame_lon_span)
eur_noramer_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_eur_noramer'%output_fname_prefix
sample_geno_geo_correlation(eur_noramer_geno_dist_ls, eur_noramer_geo_dist_ls, eur_noramer_geno_vs_geo_output_fname_prefix)

#2007-10-01 parition cc into cliques
from PartitionGraphIntoCliques import PartitionGraphIntoCliques
PartitionGraphIntoCliques_ins = PartitionGraphIntoCliques(0)
for for i in range(len(strain_iden_g_cc)):
	g_cc = strain_iden_g.subgraph(strain_iden_g_cc[i])
	PartitionGraphIntoCliques_ins.partition(g_cc.copy())
	map(len, PartitionGraphIntoCliques_ins.clique_ls)
	for j in range(len(PartitionGraphIntoCliques_ins.clique_ls)):
		gs = g_cc.subgraph(PartitionGraphIntoCliques_ins.clique_ls[j])
		site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
		for site, weight in site2weight.iteritems():
			site2weight[site] =  8*math.log(weight+1)
		
		display_matrix_of_component(input_fname, PartitionGraphIntoCliques_ins.clique_ls[j], ecotypeid2pos, output_fname='ecotype_identity_map_cc%sc%s'%(i,j), need_sort=0, need_savefig=1)
		draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s clique%s'%(i,j), pic_area=[-130,34,120,65], output_fname_prefix='ecotype_identity_cc%sc%s'%(i,j))

#2008-01-25 to investigate how identity spreads against distance in different regions.
europe_lon_span = [-12,50]
norame_lon_span = [-130,-60]
centralasia_lon_span = [60,90]
japan_lon_span = [130, 150]

span_ls = [europe_lon_span, norame_lon_span]
span_label_ls = ['europe', 'nor america']

dist_range = [0.0, 0.4]
strain_pair_ls, geno_dist_ls, geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, dist_range=dist_range)

def float2str(input):
	input = repr(input)
	input = input.replace('.', '_')
	return input

import cPickle
output_fname = '%s_strain_pair_geno_dist_%s_%s_geo_dist_ls'%(os.path.splitext(input_fname)[0], dist_range[0], dist_range[1])
output_fname = output_fname.replace('.', '_')
of = open(output_fname, 'w')
cPickle.dump([strain_pair_ls, geno_dist_ls, geo_dist_ls], of)
del of

matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls = group_identity_pair_according_to_region_distance(strain_pair_ls, geno_dist_ls, geo_dist_ls, ecotypeid2pos, span_ls, span_label_ls, geo_distance_window_size=100, max_identity_geno_dist=0.0)
draw_table_as_bar_chart(matrix_of_counts_by_region_x_geo_dist, row_label_ls, col_label_ls, col_index_ls)
"""


"""
2007-09-17
	check allele frequency
"""
def cal_maf_vector(input_fname):
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	import Numeric
	data_matrix = Numeric.array(data_matrix)
	from EstimateSelfingGeneration import EstimateSelfingGeneration
	EstimateSelfingGeneration_instance = EstimateSelfingGeneration()
	locus_allele_prob_vector = EstimateSelfingGeneration_instance.cal_locus_allele_prob_vector(data_matrix)
	maf_vector = Numeric.zeros(locus_allele_prob_vector.shape[0], Numeric.Float)	#the minor allele frequency vector
	for i in range(locus_allele_prob_vector.shape[0]):
		maf_vector[i] = min(locus_allele_prob_vector[i])
	import pylab
	pylab.hist(maf_vector, 10)
	return maf_vector
"""
input_fname = 'script/variation/stock20070829/data_row_na_col_na_bad_snps.tsv'
input_fname = '/mnt/hpc-cmb/KW/input/250K_method_5_after_imputation_noRedundant_051908.tsv'
maf_vector = cal_maf_vector(input_fname)
import pylab
pylab.clf()
pylab.plot(range(len(maf_vector)), maf_vector)
"""

"""
2008-06-27 read in pvalues from a file
"""
pvalue_fname = '/Network/Data/250k/db/results/type_1/394_results.tsv'
def plot_maf_vs_pvalue(maf_vector, input_fname, do_log10_transformation=True):
	from GenomeBrowser import GenomeBrowser
	genome_wide_result = GenomeBrowser.getGenomeWideResultFromFile(input_fname, do_log10_transformation=do_log10_transformation)
	pvalue_ls = [genome_wide_result.data_obj_ls[i].value for i in range(len(genome_wide_result.data_obj_ls))]
	import pylab
	pylab.clf()
	pylab.plot(maf_vector, pvalue_ls, '.')
	pylab.show()

"""
plot_maf_vs_pvalue(maf_vector, pvalue_fname)
for i in range(389, 758):
	pvalue_fname = '/Network/Data/250k/db/results/type_1/%s_results.tsv'%i
	plot_maf_vs_pvalue(maf_vector, pvalue_fname)
"""


"""
2007-09-21
	check duplicate calls in the database
"""
def check_inconsistent_duplicate_calls(curs, ecotype_table, calls_table, debug=0):
	"""
	2007-09-22
		buggy, don't use it. use check_inconsistent_duplicate_calls2()
	"""
	sys.stderr.write("Checking inconsistent duplicate calls ...\n")
	from common import nt2number
	no_of_strain_snp_pairs = 0
	no_of_inconsistent = 0
	inconsistency = 0
	old_call_number = 0
	old_strain_snp_pair = None
	curs.execute("select e.nativename, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
	rows = curs.fetchall()
	counter = 0
	from sets import Set
	strain_snp_pair_set = Set()
	inconsistent_dup_strain_snp_pair_set = Set()
	if debug:
		import pdb
		pdb.set_trace()
	for row in rows:
		nativename, snpid, call, callhet = row
		nativename = nativename.upper()	#bug here. same name appears in >1 forms differing in cases
		call = call.upper()
		if callhet:
			callhet.upper()
			call = call+callhet
		call_number = nt2number[call]
		strain_snp_pair = (nativename, snpid)
		if strain_snp_pair != old_strain_snp_pair:
			no_of_strain_snp_pairs += 1
			no_of_inconsistent += inconsistency
			if inconsistency == 1:
				inconsistent_dup_strain_snp_pair_set.add(old_strain_snp_pair)
			old_strain_snp_pair = strain_snp_pair
			inconsistency = 0
		else:
			if inconsistency == 0 and old_call_number!=0 and call_number!=0 and old_call_number!=call_number:	#if inconsistency==1, don't change it. comparison only between same pairs
				inconsistency = 1
		old_call_number = call_number
		counter += 1
		strain_snp_pair_set.add(strain_snp_pair)
	#don't miss out the last one
	no_of_inconsistent += inconsistency
	old_strain_snp_pair = strain_snp_pair
	sys.stderr.write('%s%s'%('\x08'*20, counter))
	sys.stderr.write("\nDone.\n")
	print "%s strain snp pairs"%(len(strain_snp_pair_set))
	print 'inconsistent ratio: %s/%s=%s'%(no_of_inconsistent, no_of_strain_snp_pairs, float(no_of_inconsistent)/no_of_strain_snp_pairs)
	return strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set

def check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=1, debug=0):
	"""
	use_ecotypeid_as_strainname controls whether it's nativename or ecotypeid in strain_snp_pair
	"""
	sys.stderr.write("Checking inconsistent duplicate calls ...\n")
	if strainname_type==1:
		print 'distinct on (nativename,snpid) pair:'
	elif strainname_type == 2:
		print 'distinct on (ecotypeid,snpid) pair:'
	elif strainname_type == 3:
		print 'distinct on (nativename, stockparent ,snpid) pair:'
	else:
		sys.stderr.write("unsupported strainname_type %s\n"%strainname_type)
		return None, None
	from common import nt2number
	no_of_strain_snp_pairs = 0
	no_of_inconsistent = 0
	old_strain_snp_pair = None
	duplicated_times = 0
	if strainname_type == 1:
		curs.execute("select e.nativename, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
	elif strainname_type == 2:
		curs.execute("select e.id, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, snpid"%(ecotype_table, calls_table))
	elif strainname_type == 3:
		curs.execute("select e.nativename, e.stockparent, c.snpid, c.call1, c.callhet from %s e, %s c where e.id=c.ecotypeid order by nativename, stockparent, snpid"%(ecotype_table, calls_table))
	rows = curs.fetchall()
	counter = 0
	from sets import Set
	strain_snp_pair_set = Set()
	inconsistent_dup_strain_snp_pair_set = Set()
	duplicate_call_type2counter = {}
	if debug:
		import pdb
		pdb.set_trace()
	for row in rows:
		if strainname_type==3:
			nativename, stockparent, snpid, call, callhet = row
		elif strainname_type==1:
			nativename, snpid, call, callhet = row
		elif strainname_type==2:
			ecotypeid, snpid, call, callhet = row
		if strainname_type!=2:
			nativename = nativename.upper()	#bug here. same name appears in >1 forms differing in cases
		call = call.upper()
		if callhet:
			callhet.upper()
			call = call+callhet
		call_number = nt2number[call]
		if strainname_type==1:
			strain_snp_pair = (nativename, snpid)
		elif strainname_type==2:
			strain_snp_pair = (ecotypeid, snpid)
		elif strainname_type==3:
			strain_snp_pair = (nativename, stockparent, snpid)
		if old_strain_snp_pair==None:	#first time
			old_strain_snp_pair = strain_snp_pair
			duplicated_times += 1
		elif strain_snp_pair != old_strain_snp_pair:
			if duplicated_times>1:
				no_of_strain_snp_pairs += 1
				if len(duplicate_call_type2counter)>1:
					no_of_inconsistent += 1
					inconsistent_dup_strain_snp_pair_set.add(old_strain_snp_pair)
			old_strain_snp_pair = strain_snp_pair
			duplicated_times = 1	#the current position is counted as 1
			duplicate_call_type2counter = {}
		else:
			duplicated_times += 1
		if call_number!=0:	#have to be put last
			if call_number not in duplicate_call_type2counter:
				duplicate_call_type2counter[call_number] = 0
			duplicate_call_type2counter[call_number] += 1
		counter += 1
		strain_snp_pair_set.add(strain_snp_pair)
	#don't miss out the last one
	if duplicated_times>1:
		no_of_strain_snp_pairs += 1
		if len(duplicate_call_type2counter)>1:
			no_of_inconsistent += 1
	old_strain_snp_pair = strain_snp_pair
	print "%s strain snp pairs"%counter
	print "%s distinct strain snp pairs"%(len(strain_snp_pair_set))
	print 'inconsistent ratio: %s/%s=%s'%(no_of_inconsistent, no_of_strain_snp_pairs, float(no_of_inconsistent)/no_of_strain_snp_pairs)
	sys.stderr.write("Done.\n")
	return strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set


def get_nativename_snp_distinct_set(curs, ecotype_table, calls_table):
	"""
	2007-09-22 find out the bug caused by the case-insensitive problem of sql.
	"""
	curs.execute("select distinct e.nativename, c.snpid from %s e, %s c where e.id=c.ecotypeid" %(ecotype_table, calls_table))
	from sets import Set
	strain_snp_pair_set = Set()
	rows = curs.fetchall()
	for row in rows:
		nativename, snpid = row
		strain_snp_pair = (nativename, snpid)
		strain_snp_pair_set.add(strain_snp_pair)
	return strain_snp_pair_set
"""
ecotype_table = 'ecotype'
calls_table = 'calls'
strain_snp_pair_set, inconsistent_dup_strain_snp_pair_set = check_inconsistent_duplicate_calls(curs, ecotype_table, calls_table)
strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)

nativename_snp_distinct_set = get_nativename_snp_distinct_set(curs, ecotype_table, calls_table)

strain_snp_pair_set - nativename_snp_distinct_set gives 1937 different pairs. turns out that sql's "distinct e.nativename, c.snpid" ignores the cases , like Kz-9 and KZ-9 are same. however python doesn't do that.
"""

"""
2007-09-23
	calculate the percentage of NAs in a data_matrix
"""
def calculate_NA_perc(input_fname):
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	import Numeric
	data_matrix = Numeric.array(data_matrix)
	no_of_NAs = 0
	for i in range(data_matrix.shape[0]):
		for j in range(data_matrix.shape[1]):
			if data_matrix[i][j] == 0:
				no_of_NAs += 1
	no_of_totals = data_matrix.shape[0]*data_matrix.shape[1]
	print '%s/%s=%s'%(no_of_NAs, no_of_totals, float(no_of_NAs)/no_of_totals)

"""
calculate_NA_perc('.//script/variation/stock20070919/data.tsv')
"""


"""
2007-09-24
	magnus's idea to reduce fake identities
"""
def get_data_matrix_sorted_by_NA_perc(input_fname):
	"""
	2007-09-24
		input_fname is a data frame file with snps coded in integers. 1st row is the ecotype id in table ecotype
	"""
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	no_of_strains = len(strain_acc_list)
	no_of_snps = len(header)-2
	data_matrix_sorted_by_NA_perc = []
	for i in range(no_of_strains):
		no_of_NAs = 0
		for j in range(no_of_snps):
			if data_matrix[i][j] == 0:
				no_of_NAs += 1
		data_matrix_sorted_by_NA_perc.append([float(no_of_NAs)/no_of_snps, data_matrix[i], strain_acc_list[i]])
	data_matrix_sorted_by_NA_perc.sort()
	return data_matrix_sorted_by_NA_perc

def get_strain_id2index_from_data_matrix_sorted_by_NA_perc(data_matrix_sorted_by_NA_perc):
	sys.stderr.write("Getting strain_id2index from data_matrix_sorted_by_NA_perc ...")
	strain_id2index = {}
	no_of_strains = len(data_matrix_sorted_by_NA_perc)
	no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
	for i in range(no_of_strains):
		ecotypeid = data_matrix_sorted_by_NA_perc[i][2]
		strain_id2index[ecotypeid] = i
	sys.stderr.write("Done.\n")
	return strain_id2index


def get_unique_haplotype_set(data_matrix_sorted_by_NA_perc):
	sys.stderr.write("Getting unique haplotype set ...")
	no_of_strains = len(data_matrix_sorted_by_NA_perc)
	no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
	from sets import Set
	unique_haplotype_set = Set()
	for i in range(no_of_strains):
		no_of_prev_identity_strains = 0
		for j in range(i):
			no_of_same_cols = 0
			for k in range(no_of_snps):
				call1 = data_matrix_sorted_by_NA_perc[i][1][k]
				call2 = data_matrix_sorted_by_NA_perc[j][1][k]
				if call1 == call2 or call1==0 or call2==0:
					no_of_same_cols += 1
			if no_of_same_cols == no_of_snps:
				no_of_prev_identity_strains += 1
		if no_of_prev_identity_strains==0:
			unique_haplotype_set.add(data_matrix_sorted_by_NA_perc[i][2])
	sys.stderr.write("Done.\n")
	return unique_haplotype_set


def assign_remainder_strains_to_unique_haplotype_set(unique_haplotype_set, data_matrix_sorted_by_NA_perc, strain_id2index):
	sys.stderr.write("Assigning remainder strains to unique_haplotype_set ...")
	identity_pair_ls = []
	no_of_strains = len(data_matrix_sorted_by_NA_perc)
	no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
	for i in range(no_of_strains):
		ecotypeid = data_matrix_sorted_by_NA_perc[i][2]
		if ecotypeid not in unique_haplotype_set:	#skip strains in unique_haplotype_set
			prev_identity_strains = []
			for seedid in unique_haplotype_set:
				seed_data_row = data_matrix_sorted_by_NA_perc[strain_id2index[seedid]][1]
				no_of_same_cols = 0
				for k in range(no_of_snps):
					call1 = data_matrix_sorted_by_NA_perc[i][1][k]
					call2 = seed_data_row[k]
					if call1 == call2 or call1==0 or call2==0:
						no_of_same_cols += 1
				if no_of_same_cols == no_of_snps:
					prev_identity_strains.append(seedid)
			if len(prev_identity_strains)==1:	#uniquely assigned
				identity_pair_ls.append([ecotypeid, prev_identity_strains[0]])
	sys.stderr.write("Done.\n")
	return identity_pair_ls

def construct_identity_pair_ls_given_strain_id_list(strain_id_list, data_matrix_sorted_by_NA_perc, strain_id2index):
	sys.stderr.write("Constructing identity_pair_ls given strain_id_list ...")
	identity_pair_ls = []
	no_of_strains = len(strain_id_list)
	no_of_snps = len(data_matrix_sorted_by_NA_perc[0][1])
	for i in range(no_of_strains):
		ecotypeid1 = repr(strain_id_list[i])	#convert int to char
		index_of_ecotypeid1 = strain_id2index[ecotypeid1]
		for j in range(i+1, no_of_strains):
			ecotypeid2 = repr(strain_id_list[j])
			index_of_ecotypeid2 = strain_id2index[ecotypeid2]
			no_of_same_cols = 0
			for k in range(no_of_snps):
				call1 = data_matrix_sorted_by_NA_perc[index_of_ecotypeid1][1][k]
				call2 = data_matrix_sorted_by_NA_perc[index_of_ecotypeid2][1][k]
				if call1 == call2 or call1==0 or call2==0:
					no_of_same_cols += 1
			if no_of_same_cols == no_of_snps:
				identity_pair_ls.append([ecotypeid1, ecotypeid2])
	sys.stderr.write("Done.\n")
	return identity_pair_ls


"""
input_fname = 'script/variation/stock20070919/data_d110_c0_5.tsv'
data_matrix_sorted_by_NA_perc = get_data_matrix_sorted_by_NA_perc(input_fname)

unique_haplotype_set = get_unique_haplotype_set(data_matrix_sorted_by_NA_perc)

strain_id2index = get_strain_id2index_from_data_matrix_sorted_by_NA_perc(data_matrix_sorted_by_NA_perc)

identity_pair_ls2 = assign_remainder_strains_to_unique_haplotype_set(unique_haplotype_set, data_matrix_sorted_by_NA_perc, strain_id2index)

import cPickle
of = open('%s_identity_pair_ls2'%os.path.splitext(input_fname)[0], 'w')
cPickle.dump(identity_pair_ls2, of)
del of
popid2ecotypeid_table = 'popid2ecotypeid_25'
pop_iden_g2, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls2)

popid2pos = {}
for popid, pos_size in popid2pos_size.iteritems():
	popid2pos[popid] = pos_size[0]

import math
popid2weight = {}
for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
	popid2weight[popid] = 8*math.log(intra_pop_connections+1)

draw_graph_on_map(pop_iden_g2, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map2')


strain_iden_g2 = construct_graph_out_of_identity_pair(identity_pair_ls2)

ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
ecotypeid2weight = {}
for n in strain_iden_g:
	ecotypeid2weight[n] = 1

import networkx as nx
strain_iden_g_cc2 = nx.connected_components(strain_iden_g2)

i=11


for i in range(len(strain_iden_g_cc2)):
	gs = strain_iden_g2.subgraph(strain_iden_g_cc2[i])
	site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
	for site, weight in site2weight.iteritems():
		site2weight[site] =  8*math.log(weight+1)
	
	display_matrix_of_component(input_fname, strain_iden_g_cc2[i], ecotypeid2pos, output_fname='ecotype_identity_tcc%s'%i, need_sort=0, need_savefig=1)
	draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, output_fname_prefix='ecotype_identity_tcc%s'%i)


#2007-09-25 complete the identity graph
identity_pair_ls2_complete = construct_identity_pair_ls_given_strain_id_list(strain_iden_g2.nodes(), data_matrix_sorted_by_NA_perc, strain_id2index)

pop_iden_g2_complete, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls2_complete)

popid2pos = {}
for popid, pos_size in popid2pos_size.iteritems():
	popid2pos[popid] = pos_size[0]

import math
popid2weight = {}
for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
	popid2weight[popid] = 8*math.log(intra_pop_connections+1)

draw_graph_on_map(pop_iden_g2_complete, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, output_fname_prefix='identity_map2_complete')

strain_iden_g2_complete = construct_graph_out_of_identity_pair(identity_pair_ls2_complete)
strain_iden_g_cc2_complete = nx.connected_components(strain_iden_g2_complete)

for i in range(len(strain_iden_g_cc2_complete)):
	gs = strain_iden_g2_complete.subgraph(strain_iden_g_cc2_complete[i])
	site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
	for site, weight in site2weight.iteritems():
		site2weight[site] =  8*math.log(weight+1)
	
	display_matrix_of_component(input_fname, strain_iden_g_cc2_complete[i], ecotypeid2pos, output_fname='ecotype_identity_cc2_complete%s'%i, need_sort=0, need_savefig=1)
	draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, output_fname_prefix='ecotype_identity_cc2_complete%s'%i)


cross_ocean_cc_set2_atlantic = check_cross_ocean_components(strain_iden_g2_complete, strain_iden_g_cc2_complete, ecotypeid2pos)

cross_ocean_cc_set2_eurasia = check_cross_ocean_components(strain_iden_g2_complete, strain_iden_g_cc2_complete, ecotypeid2pos, 60)
"""

class LD:
	"""
	2007-09-26
		calculate LD (r^2)
	2008-01-15
		put everything into class
	"""
	def __init__(self):
		pass
	
	def group_ordered_snps_into_chr_snp_2layer_ls(self, curs, snp_acc_list, snp_locus_table='snps'):
		"""
		2007-09-26
			assume snps are already in chromosome, position order, just need to find out
			where to stop the chromosome
		"""
		sys.stderr.write("Grouping ordered snps into chr_snp_2layer_ls ...")
		old_chromosome = -1
		chr_snp_2layer_ls = []
		snp_position_ls = []	#no chromosome, just position.
		chr_ls = []
		for i in range(len(snp_acc_list)):
			curs.execute("select chromosome, position from %s where snpid='%s'"%(snp_locus_table, snp_acc_list[i]))
			rows = curs.fetchall()
			chromosome, position = rows[0]
			snp_position_ls.append(position)
			if old_chromosome == -1:	#1st encounter
				old_chromosome = chromosome
			elif chromosome!=old_chromosome:
				chr_snp_2layer_ls.append(chr_ls)
				chr_ls = []
				old_chromosome = chromosome
			chr_ls.append(i)
		chr_snp_2layer_ls.append(chr_ls)
		sys.stderr.write("Done.\n")
		return chr_snp_2layer_ls, snp_position_ls
	
	def fill_in_snp_allele2index(self, diploid_allele, allele2index):
		from common import number2nt, nt2number
		if diploid_allele>4:
			nt = number2nt[diploid_allele]
			allele1 = nt2number[nt[0]]
			allele2 = nt2number[nt[1]]
		else:
			allele1 = allele2 = diploid_allele
		if allele1 not in allele2index:
			allele2index[allele1] = len(allele2index)
		if allele2 not in allele2index:
			allele2index[allele2] = len(allele2index)
		return allele1, allele2
	
	def calculate_LD(self, input_fname, curs, snp_locus_table='snps', debug=0, check_bit_ls=[1,0,0,0,0,0,0,0,0,0], chr_ls=[]):
		"""
		exclude pairs with one or two NAs
		exclude pairs both of who are heterozygous calls (can't figure out the phase)
		(only one of pairs is heterozygous is all right)
		2008-01-23 add chr_ls
		"""
		from FilterStrainSNPMatrix import FilterStrainSNPMatrix
		FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
		header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
		if debug:
			import pdb
			pdb.set_trace()
		import random
		no_of_strains = len(strain_acc_list)
		snp_acc_list = header[2:]
		chr_snp_2layer_ls, snp_position_ls = self.group_ordered_snps_into_chr_snp_2layer_ls(curs, snp_acc_list, snp_locus_table)
		no_of_chrs = len(chr_snp_2layer_ls)
		if len(chr_ls)==0:	#2008-01-23 if there's nothing in chr_ls
			chr_ls = range(no_of_chrs)
		r_square_ls = []
		distance_ls = []
		allele_freq_ls = []
		snp_pair_ls = []
		D_ls = []
		D_prime_ls = []
		import Numeric
		counter = 0
		period = len(check_bit_ls)
		for c in chr_ls:
			no_of_snps = len(chr_snp_2layer_ls[c])
			for i in range(no_of_snps):
				if len(check_bit_ls)>1:	#2008-01-23 no shuffling if there's only 1 bit
					random.shuffle(check_bit_ls)
				for j in range(i+1, no_of_snps):
					snp1_index = chr_snp_2layer_ls[c][i]
					snp2_index = chr_snp_2layer_ls[c][j]
					distance = abs(snp_position_ls[snp1_index]-snp_position_ls[snp2_index])
					if distance>10000:	#2008-01-15 skip SNPs two far
						continue
					counter += 1
					if check_bit_ls[counter%period]==0:	#2008-01-15 skip some portion of them randomly
						continue
					counter_matrix = Numeric.zeros([2,2])
					snp1_allele2index = {}
					snp2_allele2index = {}
					for k in range(no_of_strains):
						snp1_allele = data_matrix[k][snp1_index]
						snp2_allele = data_matrix[k][snp2_index]
						if snp1_allele!=0 and snp2_allele!=0 and not (snp1_allele>4 and snp2_allele>4):
							snp1_allele1, snp1_allele2 = self.fill_in_snp_allele2index(snp1_allele, snp1_allele2index)
							snp2_allele1, snp2_allele2 = self.fill_in_snp_allele2index(snp2_allele, snp2_allele2index)
							counter_matrix[snp1_allele2index[snp1_allele1],snp2_allele2index[snp2_allele1]] += 1
							counter_matrix[snp1_allele2index[snp1_allele2],snp2_allele2index[snp2_allele2]] += 1
					PA = sum(counter_matrix[0,:])
					Pa = sum(counter_matrix[1,:])
					PB = sum(counter_matrix[:,0])
					Pb = sum(counter_matrix[:,1])
					total_num = float(PA+Pa)
					try:
						PA = PA/total_num
						Pa = Pa/total_num
						PB = PB/total_num
						Pb = Pb/total_num
						PAB = counter_matrix[0,0]/total_num
						D = PAB-PA*PB
						PAPB = PA*PB
						PAPb = PA*Pb
						PaPB = Pa*PB
						PaPb = Pa*Pb
						Dmin = max(-PAPB, -PaPb)
						Dmax = min(PAPb, PaPB)
						if D<0:
							D_prime = D/Dmin
						else:
							D_prime = D/Dmax
						r2 = D*D/(PA*Pa*PB*Pb)
					except:	#2008-01-23 exceptions.ZeroDivisionError, Dmin or Dmax could be 0 if one of(-PAPB, -PaPb)  is >0 or <0
						sys.stderr.write('Unknown except, ignore: %s\n'%repr(sys.exc_info()[0]))
						continue
					allele_freq = (min(PA, Pa),min(PB, Pb))
					D_ls.append(D)
					D_prime_ls.append(D_prime)
					r_square_ls.append(r2)
					distance_ls.append(distance)
					allele_freq_ls.append(allele_freq)
					snp_pair_ls.append((snp_acc_list[i], snp_acc_list[j]))
		return D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls
	
	
	def plot_LD(self, x_ls, y_ls, title, xlabel, ylabel, max_dist=0, output_fname_prefix=None):
		import pylab
		import numpy
		import scipy.interpolate
		#x_ls = numpy.array(x_ls)
		#y_ls = numpy.array(y_ls)
		x_argsort_ls = numpy.argsort(x_ls)
		new_x_ls = []
		new_y_ls = []
		for i in range(len(x_argsort_ls)):
			if max_dist>0:
				if x_ls[x_argsort_ls[i]]<=max_dist:
					new_x_ls.append(x_ls[x_argsort_ls[i]])
					new_y_ls.append(y_ls[x_argsort_ls[i]])
			else:
				new_x_ls.append(x_ls[x_argsort_ls[i]])
				new_y_ls.append(y_ls[x_argsort_ls[i]])
		
		sp = scipy.interpolate.UnivariateSpline(new_x_ls,new_y_ls)
		step = (new_x_ls[-1]-new_x_ls[0])/100
		n_x_ls = numpy.arange(new_x_ls[0], new_x_ls[-1], step)
		n_y_ls = map(sp, n_x_ls)
		pylab.clf()
		pylab.title(title)
		pylab.xlabel(xlabel)
		pylab.ylabel(ylabel)
		pylab.plot(new_x_ls, new_y_ls, '.')
		pylab.plot(n_x_ls, n_y_ls)
		if output_fname_prefix:
			if max_dist:
				output_fname_prefix = '%s_%s'%(output_fname_prefix, max_dist)
			pylab.savefig('%s.eps'%output_fname_prefix, dpi=450)
			pylab.savefig('%s.svg'%output_fname_prefix, dpi=450)
			pylab.savefig('%s.png'%output_fname_prefix, dpi=450)
	
	def get_snp_acc2LD_data(self, LD_ls, distance_ls, snp_pair_ls):
		"""
		2008-01-24
			the previus model doesn't give good results.
			curve and real data are quite far apart. it's easily disrupted by points with low LD. axis distance=0 (the separating line) sometimes sits in the middle, like 5kb.
			so try Alex Platt's suggestion:
				r^2 = ax^b => log(r^2) = log(a) + b*log(x)

		2008-01-23
			process data for linear model fitting
			r^2 = 1/(a+bx) (page 483 of Hartl2007) => 1/r^2 = a+bx
		"""
		sys.stderr.write("Getting snp_acc2LD_data ...")
		snp_acc2LD_data = {}
		import math
		for i in range(len(snp_pair_ls)):
			snp_pair = snp_pair_ls[i]
			for snp_acc in snp_pair:
				if snp_acc not in snp_acc2LD_data:
					snp_acc2LD_data[snp_acc] = [[], []]
				if LD_ls[i]>0 and distance_ls[i]>0:
					snp_acc2LD_data[snp_acc][0].append(math.log(LD_ls[i]))
					snp_acc2LD_data[snp_acc][1].append([1, math.log(distance_ls[i])])
		sys.stderr.write("Done.\n")
		return snp_acc2LD_data
	
	def LD_linear_fitting(self, snp_acc2LD_data):
		"""
		2008-01-23
			actual linear fitting
		"""
		sys.stderr.write("Linear fitting for LD ...")
		from annot.bin.module_cc.linear_model import linear_model
		linear_model_ins = linear_model()
		snp_acc2LD_linear_fitting = {}
		for snp_acc, LD_data in snp_acc2LD_data.iteritems():
			LD_ls, X_ls = LD_data
			if len(LD_ls)>5:
				linear_model_ins.prepare_data(LD_ls, X_ls)
				linear_model_ins.run()
				coeff_list = linear_model_ins.coefficients()
				chisq_tuple = linear_model_ins.chisq_return()
				snp_acc2LD_linear_fitting[snp_acc] = [coeff_list, chisq_tuple]
				linear_model_ins.cleanup()
		del linear_model_ins
		sys.stderr.write("Done.\n")
		return snp_acc2LD_linear_fitting
	
	def check_one_snp_acc_LD_decay(self, snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc):
		"""
		2008-01-24
			follow the model change in get_snp_acc2LD_data()
		2008-01-23
			unfold the data processed in snp_acc2LD_data
			to see whether the theoretical curve (blue line) fits real data (red dots).
		"""
		import math
		a = snp_acc2LD_linear_fitting[snp_acc][0][0]
		b = snp_acc2LD_linear_fitting[snp_acc][0][1]
		theoretical_curve_func = lambda x: math.exp(a)*math.pow(x,b)
		x_ls = range(500,10000,100)
		y_ls = map(theoretical_curve_func, x_ls)
		import pylab
		pylab.clf()
		pylab.plot(x_ls, y_ls)
		
		inverse_func = lambda x: math.exp(x)
		r2_ls = map(inverse_func, snp_acc2LD_data[snp_acc][0])
		distance_ls = [math.exp(row[1]) for row in snp_acc2LD_data[snp_acc][1]]
		pylab.plot(distance_ls, r2_ls, 'r.')
		pylab.show()

"""
from variation.src.misc import LD
input_fname = os.path.expanduser('~/script/variation/stock20070919/data_d110_c0_5.tsv')
snp_locus_table = 'snps'
LD_ins = LD()
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table)

title = 'LD decay'
xlabel = 'Distance'
ylabel = r'$r^2$'
output_fname_prefix = '%s_LD_r2'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = 'D'
output_fname_prefix = '%s_LD_D'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = "D'"
output_fname_prefix = '%s_LD_D_prime'%os.path.splitext(input_fname)[0]
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 0, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
LD_ins.plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)


#2008-01-15 250k data
input_fname = os.path.expanduser('~/script/variation/genotyping/250ksnp/data/data_250k.tsv')
snp_locus_table = 'snps_250k'
LD_ins = LD()
#2008-01-23 chromosome 1 and no random skipping
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table, check_bit_ls=[1], chr_ls=[1])
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = LD_ins.calculate_LD(input_fname, curs, snp_locus_table, check_bit_ls=[1]+range(100), chr_ls=[2])

pickle_fname = os.path.expanduser('~/.pickle/LD_ins')
of = open(pickle_fname, 'w')
cPickle.dump([D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls], of)
del of

#2008-01-23 load the computed data by cPickle
import cPickle
inf = open(pickle_fname, 'r')
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = cPickle.load(inf)
del inf

snp_acc2LD_data = LD_ins.get_snp_acc2LD_data(r_square_ls, distance_ls, snp_pair_ls)
snp_acc2LD_linear_fitting = LD_ins.LD_linear_fitting(snp_acc2LD_data)

def check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc):
	import math
	a = snp_acc2LD_linear_fitting[snp_acc][0][0]
	b = snp_acc2LD_linear_fitting[snp_acc][0][1]
	theoretical_curve_func = lambda x: math.exp(a)*math.pow(x,b)
	x_ls = range(500,10000,100)
	y_ls = map(theoretical_curve_func, x_ls)
	import pylab
	pylab.clf()
	pylab.plot(x_ls, y_ls)
	
	inverse_func = lambda x: math.exp(x)
	r2_ls = map(inverse_func, snp_acc2LD_data[snp_acc][0])
	distance_ls = [math.exp(row[1]) for row in snp_acc2LD_data[snp_acc][1]]
	pylab.plot(distance_ls, r2_ls, 'r.')
	pylab.show()


LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting,'1_10219_T_A')
LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting,'1_3102_A_G')

for snp_acc in snp_acc2LD_linear_fitting:
	LD_ins.check_one_snp_acc_LD_decay(snp_acc2LD_data, snp_acc2LD_linear_fitting, snp_acc)
	print snp_acc2LD_linear_fitting[snp_acc]
	raw_input(":")
"""

"""
2007-10-16
"""
def get_snp_index2pos(snp_acc_list, curs, snps_table='stock20071008.snps'):
	sys.stderr.write("Getting snp_index2pos ...")
	snp_index2pos = {}
	for i in range(len(snp_acc_list)):
		curs.execute("select chromosome, position from %s where snpid='%s'"%(snps_table, snp_acc_list[i]))
		rows = curs.fetchall()
		chromosome, position = rows[0]
		snp_index2pos[i] = (chromosome, position)
	sys.stderr.write("Done.\n")
	return snp_index2pos

def get_shared_block_ls(row1, row2, snp_index2pos):
	shared_block_ls = []
	shared_block = []
	for i in range(len(row1)):
		if row1[i] == row2[i]:	#same allele
			if i==0 or snp_index2pos[i-1][0]==snp_index2pos[i][0]:	#1st snp or previous snp and this one are on same chromosome
				shared_block.append(i)
			elif shared_block:	#same allele, but previous snp is on different chromosome
				shared_block_ls.append(shared_block)
				shared_block = []	#clear it up
				shared_block.append(i)
		elif shared_block:	#different allele but shared_block is not empty
			shared_block_ls.append(shared_block)	#store the shared_block
			shared_block = []	#clear it up
	if row1[-1]==row2[-1]:	#last snp is same. then the shared_block is not put in yet.
		shared_block_ls.append(shared_block)
	return shared_block_ls

def get_all_pairs_with_max_shared_block_length_ls(strain_acc_list, snp_index2pos, data_matrix):
	sys.stderr.write("Getting all_pairs_with_max_shared_block_length_ls ...")
	all_pairs_with_max_shared_block_length_ls = []
	no_of_strains = len(strain_acc_list)
	for i in range(no_of_strains):
		for j in range(i+1, no_of_strains):
			shared_block_ls = get_shared_block_ls(data_matrix[i], data_matrix[j], snp_index2pos)
			shared_block_length_ls = map(len, shared_block_ls)
			shared_block_length_ls.sort()
			all_pairs_with_max_shared_block_length_ls.append([shared_block_length_ls[-1], (i,j)])	#take the maximum shared block length
	all_pairs_with_max_shared_block_length_ls.sort()
	sys.stderr.write("Done.\n")
	return all_pairs_with_max_shared_block_length_ls

def DrawTwoRowAndSharedBlock(data_matrix, strain_index_pair, strain_acc_list, snp_index2pos):
	"""
	2007-10-16
		used to check get_shared_block_ls() is correct
	"""
	import Numeric
	from common import nt2number, number2nt
	data_row1 = []
	data_row2 = []
	shared_block_row = []
	no_of_chrs = 0
	row1 = data_matrix[strain_index_pair[0]]
	row2 = data_matrix[strain_index_pair[1]]
	for i in range(len(row1)):
		#for shared_block_row: -1 is chromosome separator, 1 is same, 0 is different.
		if i!=0 and snp_index2pos[i-1][0]!=snp_index2pos[i][0]:	#a different chromosome
			shared_block_row.append(-1)
			data_row1.append(-1)	#-1 is chromosome separator
			data_row2.append(-1)
		if row1[i]==row2[i]:
			shared_block_row.append(1)
		else:
			shared_block_row.append(0)
		data_row1.append(row1[i])
		data_row2.append(row2[i])
	sub_data_matrix = Numeric.array([data_row1, data_row2, shared_block_row])
	import pylab
	pylab.clf()
	pylab.imshow(sub_data_matrix, interpolation='nearest')
	pylab.colorbar()
	ytick_label_ls = []
	for strain_index in strain_index_pair:
		ytick_label_ls.append(strain_acc_list[strain_index])
	ytick_label_ls.append('block')
	pylab.yticks([2,1,0], ytick_label_ls)	#it's reversed.
	pylab.show()


def DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size):
	"""
	2008-02-01
		add comments
	2007-10-16
	"""
	import pylab
	pylab.clf()
	#draw the chromosome separator as green circles
	for chr_id, cumu_size in chr_id2cumu_size.iteritems():
		pylab.plot([cumu_size], [1], 'go')
	#draw the snp first as red sticks
	for snp_index,pos in snp_index2pos.iteritems():
		chr_id, chr_pos = pos
		cumu_chr_pos = chr_id2cumu_size[chr_id-1]+chr_pos
		pylab.plot([cumu_chr_pos], [1], 'r|')
	#draw the blocks as lines crossing those red sticks
	for shared_block in shared_block_ls:
		if len(shared_block)>1:
			starting_snp_index = shared_block[0]
			ending_snp_index = shared_block[-1]
			
			chr_id = snp_index2pos[starting_snp_index][0]
			starting_chr_pos = snp_index2pos[starting_snp_index][1]
			ending_chr_pos = snp_index2pos[ending_snp_index][1]
			#add the cumulative chromosome size
			cumu_starting_chr_pos = chr_id2cumu_size[chr_id-1]+starting_chr_pos
			cumu_ending_chr_pos = chr_id2cumu_size[chr_id-1]+ending_chr_pos
			#draw the block
			pylab.plot([cumu_starting_chr_pos, cumu_ending_chr_pos], [1,1], c='b')
	#chromosome separator on the x axis
	xtick_label_ls = []
	xtick_loc_ls = []
	chr_id_ls = chr_id2cumu_size.keys()
	chr_id_ls.sort()
	for i in range(1,len(chr_id_ls)):
		chr_id = chr_id_ls[i]
		xtick_loc_ls.append(chr_id2cumu_size[chr_id])
		xtick_label_ls.append('%s:%s'%(chr_id, chr_id2size[chr_id]))
	pylab.xticks(xtick_loc_ls, xtick_label_ls)	#it's reversed.
	pylab.show()


"""
#2007-10-16 check shared block between pairwise strains

input_fname = 'script/variation/stock20071008/data_d110_c0_5_d001.tsv'
from FilterStrainSNPMatrix import FilterStrainSNPMatrix
FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
snp_acc_list = header[2:]
snp_index2pos = get_snp_index2pos(snp_acc_list, curs, snps_table='stock20071008.snps')
all_pairs_with_max_shared_block_length_ls = get_all_pairs_with_max_shared_block_length_ls(strain_acc_list, snp_index2pos, data_matrix)

#2007-10-16 test to draw the diagrams
from common import get_chr_id2size, get_chr_id2cumu_size
chr_id2size = get_chr_id2size(curs)
data_ls = get_chr_id2cumu_size(chr_id2size)
chr_id2cumu_size, chr_gap = data_ls[:2]
shared_block_ls = get_shared_block_ls(data_matrix[1], data_matrix[2], snp_index2pos)
DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size)

"""

"""
2007-11-08
	get 192 strains for suzi and draw them on the map
"""
def get_ecotypeid_ls_nativename_ls_of_192_strains(curs):
	ecotypeid_ls = []
	nativename_ls = []
	curs.execute("select b.ecotypeid, e.nativename from ecotype e, batch_ecotype b where b.ecotypeid=e.id and b.batchid=2")
	rows = curs.fetchall()
	for row in rows:
		ecotypeid, nativename = row
		ecotypeid_ls.append(ecotypeid)
		nativename_ls.append(nativename)
	return ecotypeid_ls, nativename_ls


"""
ecotypeid_ls, nativename_ls = get_ecotypeid_ls_nativename_ls_of_192_strains(curs)
ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
draw_strains_on_map(ecotypeid_ls, ecotypeid2pos, pic_title='192 strains',  pic_area=[-130,10,140,70], output_fname_prefix='/tmp/suzi_192', label_ls=nativename_ls, need_collapse_strains_with_same_pos=0)

"""


def simulatePvalue(curs, snps_table, output_fname):
	"""
	2008-01-30
		simulate uniform p-values for snps selected from db
	"""
	import csv
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	curs.execute("select chromosome, position from %s order by chromosome, position"%snps_table)
	rows = curs.fetchall()
	import random
	for row in rows:
		chromosome, position = row
		writer.writerow([chromosome, position, random.expovariate(1)])
	del writer

"""
snps_table = 'snps_250k'
output_fname = '/tmp/simulate.pvalue'
simulatePvalue(curs, snps_table, output_fname)
"""

def minusLogPvalue(input_fname, output_fname):
	"""
	2008-02-14
		take log of pvalue in input_fname
	"""
	import csv
	reader = csv.reader(open(input_fname), delimiter='\t')
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	import math
	for row in reader:
		chromosome, position, pvalue = row
		pvalue = float(pvalue)
		if pvalue>0:
			pvalue = -math.log(float(pvalue))
			writer.writerow([chromosome, position, pvalue])
	del writer, reader

def exponentialMinusLogPvalue(input_fname, output_fname):
	"""
	2008-02-14
		take log of pvalue in input_fname
	"""
	import csv
	reader = csv.reader(open(input_fname), delimiter='\t')
	reader.next()	#skip the header
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	import math
	for row in reader:
		chromosome, position, pvalue = row
		pvalue = float(pvalue)
		if pvalue>0:
			pvalue = math.exp(-float(pvalue))
			writer.writerow([chromosome, position, pvalue])
	del writer, reader


"""
2008-03-18
	check chiamo output
"""
def convertChiamoOutput(chiamo_infname, chiamo_outfname, ref_250k_infname, output_fname, posterior_min=0.95):
	from variation.src.common import nt2number
	import csv
	import numpy
	
	reader = csv.reader(open(ref_250k_infname), delimiter='\t')
	SNP_acc_ls = reader.next()[2:]
	del reader
	
	SNP_id2allele_ls = {}
	SNP_id2info = {}
	for SNP_acc in SNP_acc_ls:
		chr, pos, alleleA, alleleB = SNP_acc.split('_')
		SNP_id = '%s_%s'%(chr, pos)
		SNP_id2allele_ls[SNP_id] = [alleleA, alleleB]
		SNP_id2info[SNP_id] = [SNP_acc, len(SNP_id2allele_ls)-1]	#[name, index]
	
	reader = csv.reader(open(chiamo_infname), delimiter='\t')
	strain_id_dup_ls = reader.next()[5:]
	strain_id_ls = []
	for i in range(0, len(strain_id_dup_ls), 2):
		strain_id_ls.append(strain_id_dup_ls[i][:-2])
	del reader
	
	sys.stderr.write("Reading chiamo output ...\n")
	reader = csv.reader(open(chiamo_outfname), delimiter=' ')
	snp_acc_ls = []
	data_matrix = numpy.zeros([len(strain_id_ls), len(SNP_id2allele_ls)], numpy.integer)
	counter = 0
	for row in reader:
		SNP_id = row[0]
		if SNP_id not in SNP_id2allele_ls:
			continue
		allele_ls = SNP_id2allele_ls[SNP_id]
		snp_acc, col_index = SNP_id2info[SNP_id]
		snp_acc_ls.append(snp_acc)
		for i in range(5, len(row)-1, 3):
			posterior_ls = [float(row[i]), float(row[i+1]), float(row[i+2])]	#posterior for 3 classes
			max_posterior_index = numpy.argmax(posterior_ls)
			if posterior_ls[max_posterior_index]>=posterior_min:
				if max_posterior_index==2:	#heterozygous
					call = '%s%s'%(allele_ls[0], allele_ls[1])
				else:
					call = allele_ls[max_posterior_index]
				row_index = (i-5)/3
				try:
					data_matrix[row_index][col_index] = nt2number[call]
				except:
					print SNP_id, snp_acc, col_index, allele_ls, row_index, strain_id_ls[row_index], call
					return
		counter += 1
		if counter%2000==0:
			sys.stderr.write("%s%s"%('\x08'*40, counter))
				
	del reader
	sys.stderr.write("Done.\n")
	
	from variation.src.FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header = ['ecotypeid', 'ecotypeid'] + snp_acc_ls
	FilterStrainSNPMatrix_instance.write_data_matrix(data_matrix, output_fname, header, strain_id_ls, [1]*len(strain_id_ls))
	sys.stderr.write("Finished.\n")

"""
chiamo_infname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo_14SNPs.in')
chiamo_outfname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo.out_0_mcmc')
ref_250k_infname = os.path.expanduser('~/script/variation/genotyping/250ksnp/data/data_250k.tsv')
output_fname = os.path.expanduser('~/script/affy/250k_test/yanli8-29-07_chiamo_out.tsv')

convertChiamoOutput(chiamo_infname, chiamo_outfname, ref_250k_infname, output_fname, posterior_min=0.95)
"""





"""
2008-04-11
	this is a one-time 250k pipeline fix to avoid memory-hefty intensity re-output.
	
	form 2 lists of intensity matrix filenames based on the new and old array_info_table
	output them as two-column (old_fname, new_fname) into a output_fname.
	shell/file_batch_move.py reads the output_fname and handle name changing.

2008-04-11
	turns out to be useless because the header in each intensity matrix file has array_id embedded.
	have to output each array into intensity matrix.
"""
def output_intensity_fname(curs, new_array_info_table, old_array_info_table, output_fname):
	"""
	"""
	import csv
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	curs.execute("select n.id, o.id from %s n, %s o where n.original_filename=o.filename"%(new_array_info_table, old_array_info_table))
	rows = curs.fetchall()
	for row in rows:
		new_array_id, old_array_id = row
		new_fname = '%s_array_intensity.tsv'%new_array_id
		old_fname = '%s_array_intensity.tsv'%old_array_id
		writer.writerow([old_fname, new_fname])
	del writer

"""
new_array_info_table='stock_250k.array_info'
old_array_info_table='stock_250k.array_info_2008_04_11'
output_fname = '/tmp/intensity_fname.rename.tsv'
output_intensity_fname(curs, new_array_info_table, old_array_info_table, output_fname)
"""


"""
2008-05-31
	output results from stock_250k.results, totally db based, not reading results.filename
"""
def outputResults(db, results_method_id, output_fname):
	import csv
	import sqlalchemy as sql
	
	conn = db.connection	#establish the connection before referring db.tables (in case it hasn't been setup)
	sys.stderr.write("Getting marker_id2position ... ")
	marker_id2position = {}
	snps_table = db.tables['snps'].alias()
	results = conn.execute(sql.select([snps_table.c.id, snps_table.c.chromosome, snps_table.c.position, snps_table.c.end_position]))
	for row in results:
		marker_id2position[row.id] = (row.chromosome, row.position, row.end_position)
	del results
	sys.stderr.write("Done.\n")
	
	sys.stderr.write("Outputting results_method_id=%s ... "%results_method_id)
	results_table = db.tables['results'].alias()
	results = conn.execute(sql.select([results_table.c.snps_id, results_table.c.score], results_table.c.results_method_id==results_method_id,\
									  order_by=[results_table.c.snps_id]))
	
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	writer.writerow(['chromosome', 'position', 'score'])
	for row in results:
		chromosome, position, end_position = marker_id2position[row.snps_id]
		writer.writerow([chromosome, position, row.score])
	del writer
	sys.stderr.write("Done.\n")

"""
from variation.src.db import Stock_250kDatabase
db = Stock_250kDatabase(username='nordborglab',
				   password='papaya', hostname='papaya.usc.edu', database='stock_250k')
conn = db.connection	#establish the connection before referring db.tables (it needs to be setup)
outputResults(db, 5, '/tmp/5_results.tsv')
outputResults(db, 6, '/tmp/6_results.tsv')
"""





"""
2008-07-21
"""
def drawScoreHistogram(curs, results_method_id, list_type_id=1, do_log10_transformation=True, min_or_max_func='min'):
	from Stock_250kDB import ResultsMethod, GeneList, ResultsByGene, CandidateGeneRankSumTestResult
	rm = ResultsMethod.get(results_method_id)
	db_rows = GeneList.query.filter_by(list_type_id=list_type_id)
	from sets import Set
	gene_set = Set()
	for gl in db_rows:
		gene_set.add(gl.gene_id)
	
	import math
	score_ls1 = []
	score_ls2 = []
	#db_rows = ResultsByGene.query.filter_by(results_method_id=results_method_id)
	curs.execute("select r.gene_id, %s(r.score) as score from results_by_gene r  where r.results_method_id=%s group by r.gene_id"%(min_or_max_func, results_method_id))
	db_rows = curs.fetchall()
	for row in db_rows:
		gene_id, score = row
		if do_log10_transformation:
			score = -math.log(score)
		if gene_id in gene_set:
			score_ls1.append(score)
		else:
			score_ls2.append(score)
	import pylab
	pylab.clf()
	pylab.title('results_method_id=%s, (%s on %s) by list type: %s.'%(rm.id, rm.analysis_method.short_name, rm.phenotype_method.short_name, gl.list_type.short_name))
	n1 = pylab.hist(score_ls1, 100, alpha=0.4, normed=1)
	n2 = pylab.hist(score_ls2, 100, alpha=0.4, normed=1, facecolor='r')
	pylab.legend(['candidate gene list', 'non-candidate gene list'])
	pylab.show()
	return score_ls1, score_ls2

"""
score_ls1, score_ls2 = drawScoreHistogram(curs, 23, 1)
"""



"""
2008-07-31
"""
def getFRIAlignment(output_fname, alignment_id=1843):
	from variation.src.AtDB import AtDB, Sequence, Alignment
	db = AtDB(hostname='localhost')
	rows = Sequence.query.filter_by(alignment=alignment_id).order_by(Sequence.accession).all()
	"""
	#2008-07-31 output in fasta format
	outf = open(output_fname, 'w')
	is_target_alignment_outputted = 0
	for row in rows:
		if not is_target_alignment_outputted:
			outf.write('>ref\n')
			outf.write('%s\n'%row.alignment_obj.target)
			is_target_alignment_outputted = 1
		outf.write('>%s %s\n'%(row.accession, row.accession_obj.name))
		outf.write('%s\n'%(row.bases))
	del outf
	"""
	#2008-08-01 output in a yh SNP matrix format for DrawSNPMatrix.py
	import csv
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	is_target_alignment_outputted = 0
	for row in rows:
		if not is_target_alignment_outputted:
			header_row = ['name', 1]
			one_row = ['ref', 1]
			for i in range(len(row.alignment_obj.target)):
				base = row.alignment_obj.target[i]
				one_row.append(base)
				header_row.append(i+1)
			writer.writerow(header_row)
			writer.writerow(one_row)
			is_target_alignment_outputted = 1
		one_row = ['%s %s'%(row.accession, row.accession_obj.name), 1]
		for base in row.bases:
			one_row.append(base)
		writer.writerow(one_row)
	del writer

"""		
getFRIAlignment('/tmp/alignment_1843.fasta')
"""

"""
2008-08-04 investigate whether throwing off some rows help to increase significance
"""
def removeRowsBasedOnSNPAllele(input_fname, output_fname, SNP_label, allele='-'):
	from pymodule import read_data, write_data_matrix, nt2number
	header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
	snp_acc_ls = header[2:]
	snp_index = -1
	for i in range(len(snp_acc_ls)):
		if snp_acc_ls[i]==SNP_label:
			snp_index = i
	
	allele = nt2number[allele]
	rows_to_be_tossed_out = set()
	for i in range(len(data_matrix)):
		if data_matrix[i][snp_index]==allele:
			rows_to_be_tossed_out.add(i)
	
	write_data_matrix(data_matrix, output_fname, header, strain_acc_list, category_list, rows_to_be_tossed_out=rows_to_be_tossed_out)

"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_2.tsv')
SNP_label = '4_268809_0'
removeRowsBasedOnSNPAllele(input_fname, output_fname, SNP_label, allele='-')
"""


"""
2008-08-05
	NPUTE can't work with SNPs with >2 alleles
"""
def removeSNPsWithMoreThan2Alleles(input_fname, output_fname):
	from pymodule import SNPData
	snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1)
	newSNPData = snpData.removeSNPsWithMoreThan2Alleles(snpData)
	newSNPData.tofile(output_fname)

"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_only_2_alleles.tsv')
removeSNPsWithMoreThan2Alleles(input_fname, output_fname)
"""


"""
2008-08-05
	NPUTE output format is SNPXStrain by and large.
		1st and 2nd column are same as input's 1st row. 1st row is input's 1st column. 2nd row is input's 2nd column.
	
"""
def turnNPUTEOutputIntoYuFormat(input_fname, output_fname):
	from pymodule import SNPData
	snpData = SNPData(input_fname=input_fname, turn_into_integer=1, turn_into_array=1, double_header=1, ignore_2nd_column=1)
	snpData.col_id_ls = snpData.header[0][2:]	#take the first header
	
	from pymodule.SNP import transposeSNPData
	newSNPData = transposeSNPData(snpData)
	newSNPData.strain_acc_list = newSNPData.row_id_ls
	newSNPData.category_list = snpData.header[1][2:]
	newSNPData.tofile(output_fname)
"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.NPUTE.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.NPUTE.yh.tsv')
turnNPUTEOutputIntoYuFormat(input_fname, output_fname)
"""

"""
2008-08-05
	test boolean relationship between SNPs
"""

def returnTop2Allele(snp_allele2count):
	"""
	2008-08-06 remove redundant argument snp_allele_ls
	2008-08-05 in the descending order of count for each allele, assign index ascending
	"""
	snp_allele_count_ls = []
	snp_allele_ls = snp_allele2count.keys()
	for snp_allele in snp_allele_ls:
		snp_allele_count_ls.append(snp_allele2count[snp_allele])
	import numpy
	argsort_ls = numpy.argsort(snp_allele_count_ls)
	new_snp_allele2index = {}
	for i in [-1, -2]:
		snp_index = argsort_ls[i]	#-1 is index for biggest, -2 is next biggest
		new_snp_allele2index[snp_allele_ls[snp_index]] = -i-1
	return new_snp_allele2index

def booleanMergeSNPs(input_fname, output_fname, SNP_label1, SNP_label2, operator_type=1):	#1 is and, 2 is or
	"""
	2008-08-05
		alleles not in the top 2 are taken as NA. major allele is coded as 0. minor allele is coded as 1.
	"""
	from pymodule import read_data, write_data_matrix, nt2number
	header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=1)
	
	snp_acc_ls = header[2:]
	snp_index1 = -1
	snp_index2 = -1
	for i in range(len(snp_acc_ls)):
		if snp_acc_ls[i]==SNP_label1:
			snp_index1 = i
		if snp_acc_ls[i]==SNP_label2:
			snp_index2 = i
	
	snp_allele2count1 = {}
	snp_allele2count2 = {}
	no_of_rows = len(data_matrix)
	
	for i in range(no_of_rows):
		snp1_allele = data_matrix[i][snp_index1]
		snp2_allele = data_matrix[i][snp_index2]
		if snp1_allele!=0:
			if snp1_allele not in snp_allele2count1:
				snp_allele2count1[snp1_allele] = 0
			snp_allele2count1[snp1_allele] += 1
		if snp2_allele!=0:
			if snp2_allele not in snp_allele2count2:
				snp_allele2count2[snp2_allele] = 0
			snp_allele2count2[snp2_allele] += 1
	print snp_allele2count1
	print snp_allele2count2
	snp_allele2index1 = returnTop2Allele(snp_allele2count1)
	
	snp_allele2index2 = returnTop2Allele(snp_allele2count2)
	print snp_allele2index1
	print snp_allele2index2
	
	no_of_cols = 1
	new_data_matrix = data_matrix	#replace the 1st SNP's data with the new boolean result
	for i in range(no_of_rows):
		snp1_allele = data_matrix[i][snp_index1]
		snp2_allele = data_matrix[i][snp_index2]
		if snp1_allele in snp_allele2index1:
			snp_code1 = snp_allele2index1[snp1_allele]
		else:
			snp_code1 = 0
		
		if snp2_allele in snp_allele2index2:
			snp_code2 = snp_allele2index2[snp2_allele]
		else:
			snp_code2 = 0
		if operator_type==1:
			if snp1_allele in snp_allele2index1 and snp2_allele in snp_allele2index2:
				new_data_matrix[i][snp_index1] = (snp_code1 and snp_code2) + 1
		elif operator_type==2:
			if snp1_allele in snp_allele2index1 or snp2_allele in snp_allele2index2:
				new_data_matrix[i][snp_index1] = (snp_code1 or snp_code2) + 1
	write_data_matrix(new_data_matrix, output_fname, header, strain_acc_list, category_list)

"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2.tsv')
SNP_label1 = '4_268809_0'
SNP_label2 = '4_269962_8'
booleanMergeSNPs(input_fname, output_fname, SNP_label1, SNP_label2, operator_type=2)
input_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2.tsv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/alignment_1843_matrix_1_or_2_or_3.tsv')
SNP_label3 = '4_270712_0'
booleanMergeSNPs(input_fname, output_fname, SNP_label1, SNP_label3, operator_type=2)
"""


"""
	in processing FRI deletion data from Shindo2005. check plone doc, /research/variation/log-2008-07.
2008-08-05
"""
def outputShindo2005(input_fname, output_fname, which_type_of_id_to_output=1):
	"""
	2008-08-06
		in output. if which_type_of_id_to_output==1, output in ecotypeid; else output in accession_id
	"""
	from variation.src.AtDB import AtDB, Sequence, Alignment
	db = AtDB(hostname='localhost')
	from pymodule import read_data, write_data_matrix
	header, strain_acc_list, category_list, data_matrix = read_data(input_fname, turn_into_integer=0)
	
	no_of_rows = len(strain_acc_list)
	import numpy
	new_data_matrix = numpy.ones([no_of_rows, 2], numpy.int8)
	new_strain_acc_list = []
	new_category_list = []
	for i in range(len(strain_acc_list)):
		strain_acc = strain_acc_list[i].upper()
		#2008-08-06 for magnus's data.csv	no swedish letters.
		rows = db.metadata.bind.execute("select a2e.ecotype_id, a2e.nativename, a2e.accession_id from accession2tg_ecotypeid a2e, magnus_192_vs_accession m2a where upper(m2a.linename)='%s' and m2a.accession_id=a2e.accession_id"%\
								(strain_acc))
		
		#2008-08-06 for stuff extracted from Supplemental Figure 4A, Figure 4B. this query doesn't help to recognize those swedish letters
		#rows = db.metadata.bind.execute("select a2e.ecotype_id, a2e.nativename, a2e.accession_id from accession2tg_ecotypeid a2e where upper(a2e.accession_name)='%s'"%\
		#						(strain_acc))
		try:
			row = rows.fetchone()
			if which_type_of_id_to_output==1:
				new_strain_acc_list.append(row.ecotype_id)
			else:
				new_strain_acc_list.append(row.accession_id)
			new_category_list.append(row.nativename)
			deletion_code = int(category_list[i])
			deletion_code_index = deletion_code-2
			if deletion_code_index>=0:
				new_data_matrix[i][deletion_code_index] = 2	#allele 2
		except:
			print i, strain_acc
			new_strain_acc_list.append(strain_acc)
			new_category_list.append(strain_acc)
		
	new_header = header[:2] + ['4_268809_0', '4_269962_8']
	write_data_matrix(new_data_matrix, output_fname, new_header, new_strain_acc_list, new_category_list)

"""
input_fname = os.path.expanduser('~/script/variation/doc/FRI/Shindo2005_data.csv')
output_fname = os.path.expanduser('~/script/variation/doc/FRI/Shindo2005_data_SNP.tsv')
outputShindo2005(input_fname, output_fname)
"""


"""
2008-08-16
	check if the array_ids and ecotype_ids in call files generated by bjarni match the ones in db
"""
def checkBjarniFile(input_fname, curs):
	import csv
	reader = csv.reader(open(input_fname))
	array_id_ls = reader.next()[2:]
	array_id_ls = map(int, array_id_ls)
	print "%s arrays"%(len(array_id_ls))
	ecotype_id_ls = reader.next()[2:]
	ecotype_id_ls = map(int, ecotype_id_ls)
	print "%s ecotypes"%(len(ecotype_id_ls))
	for i in range(len(array_id_ls)):
		array_id = array_id_ls[i]
		ecotype_id = ecotype_id_ls[i]
		curs.execute("select tg_ecotypeid from stock.ecotypeid2tg_ecotypeid where ecotypeid=%s"%ecotype_id)
		rows = curs.fetchall()
		if rows:
			bjarni_ecotype_id = rows[0][0]
		else:
			sys.stderr.write( "ecotype_id %s has no tg_ecotypeid in ecotypeid2tg_ecotypeid.\n"%(ecotype_id))
			bjarni_ecotype_id = ecotype_id
		curs.execute("select maternal_ecotype_id from stock_250k.array_info where id=%s"%array_id)
		rows = curs.fetchall()
		maternal_ecotype_id = rows[0][0]
		if bjarni_ecotype_id!=maternal_ecotype_id:
			print i, array_id, bjarni_ecotype_id, maternal_ecotype_id
	del reader

"""
input_fname = '/Network/Data/250k/dataFreeze_080608/250K_f8_080608.csv'
checkBjarniFile(input_fname, curs)
"""


#2007-03-05 common codes to initiate database connection
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/test/python')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/variation/src')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

"""
from codense.common import db_connect, form_schema_tables
hostname='dl324b-1'
dbname='yhdb'
schema = 'dbsnp'
hostname='zhoudb'
dbname='graphdb'
schema = 'dbsnp'
#conn, curs = db_connect(hostname, dbname, schema)

hostname='localhost'
dbname='stock20070829'
import MySQLdb
conn0 = MySQLdb.connect(db=dbname,host=hostname)
curs0 = conn0.cursor()

hostname='papaya.usc.edu'
dbname='stock_250k'
db_user='yh'
db_passwd = ''
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname, user=db_user, passwd=db_passwd)
curs = conn.cursor()

drivername='mysql'
schema = None
from Stock_250kDB import Stock_250kDB
db = Stock_250kDB(drivername=drivername, username=db_user,
				password=db_passwd, hostname=hostname, database=dbname, schema=schema)
"""

if __name__ == '__main__':
	ecotype_table = 'ecotype'
	calls_table = 'calls'
	strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
	strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
	strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)
