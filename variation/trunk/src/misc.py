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
2007-04-29
	to determine which SNP is in coding or non-coding region
"""
def find_SNP_context(curs, snp_locus_table, snp_locus_context_table, entrezgene_mapping_table='sequence.entrezgene_mapping',\
	annot_assembly_table='sequence.annot_assembly', tax_id=3702, need_commit=0):
	import os,sys
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/transfac/src')))
	from TFBindingSiteParse import TFBindingSiteParse
	TFBindingSiteParse_instance =TFBindingSiteParse()
	from codense.common import get_entrezgene_annotated_anchor
	chromosome2anchor_gene_tuple_ls, gene_id2coord = get_entrezgene_annotated_anchor(curs, tax_id, entrezgene_mapping_table,\
		annot_assembly_table)
	curs.execute("select id, chromosome, position from %s"%(snp_locus_table))
	rows = curs.fetchall()
	counter = 0
	for row in rows:
		snp_locus_id, chromosome, position = row
		regulatory_coord = (chromosome, position, position)
		target_gene_ls, target_gene_ls_type = TFBindingSiteParse_instance.return_target_gene_ls(regulatory_coord, \
			chromosome2anchor_gene_tuple_ls, gene_id2coord)
		for gene_id, gene_start, gene_stop, gene_strand, gene_genomic_gi in target_gene_ls:
			if gene_strand == '+':
				disp_pos = position - gene_start
			else:
				disp_pos = gene_stop - position
			curs.execute("insert into %s(snp_locus_id, disp_pos, gene_id, gene_strand, disp_pos_comment) values (%s, %s, %s, '%s', '%s')"%\
				(snp_locus_context_table, snp_locus_id, disp_pos, gene_id, gene_strand, target_gene_ls_type))
		sys.stderr.write("%s%s"%('\x08'*10, counter))
		counter += 1
	if need_commit:
		curs.execute("end")

"""
conn, curs = db_connect(hostname, dbname, schema)
find_SNP_context(curs, 'snp_locus', 'snp_locus_context', need_commit=1)
"""

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


def draw_strains_on_map(ecotypeid_ls, ecotypeid2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None):
	"""
	2007-10-09
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
	
	site2weight, site2pos = collapseStrainsWithSamePos(ecotypeid_ls, ecotypeid2pos)
	
	sys.stderr.write("\tDrawing nodes ...")
	euc_coord1_ls = []
	euc_coord2_ls = []
	diameter_ls = []
	for n in site2pos:
		lat, lon = site2pos[n]
		euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
		euc_coord1_ls.append(euc_coord1)
		euc_coord2_ls.append(euc_coord2)
		diameter_ls.append(math.sqrt(site2weight[n]))
	import numpy
	diameter_ls = numpy.array(diameter_ls)
	m.scatter(euc_coord1_ls, euc_coord2_ls, 8*diameter_ls, marker='o', color='r', alpha=0.4, zorder=12, faceted=False)
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
		pylab.savefig('%s.svg'%output_fname_prefix, dpi=600)
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
	from PartitionGraphIntoCliques import PartitionGraphIntoCliques
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

def calGenoDistAndGeoDist(data_matrix_fname, ecotypeid2pos, longitude_span=[-180, 180]):
	"""
	2007-10-03
		plot to see relationship between genotype distance and geographic distance
		-cal_great_circle_distance()
	2007-10-04
		add longitude_span
	"""
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
					geno_dist_ls.append(distance)
					geo_dist = cal_great_circle_distance(pos1[0], pos1[1], pos2[0], pos2[1])
					geo_dist_ls.append(geo_dist)
				else:
					no_of_NA_pairs += 1
	print "out of %s pairs, %s are NA"%(no_of_pairs, no_of_NA_pairs)
	return geno_dist_ls, geo_dist_ls

def calGenoDistAndGeoDistBetweenTwoAreas(data_matrix_fname, ecotypeid2pos, longitude_span1=[-180, 180], longitude_span2=[-180, 180]):
	"""
	2007-10-05
		modified from calGenoDistAndGeoDist with two longitude spans
	"""
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
geno_dist_ls, geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos)
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

eur_geno_dist_ls, eur_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, europe_lon_span)
eur_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_eur'%output_fname_prefix
sample_geno_geo_correlation(eur_geno_dist_ls, eur_geo_dist_ls, eur_geno_vs_geo_output_fname_prefix)

noramer_geno_dist_ls, noramer_geo_dist_ls =  calGenoDistAndGeoDist(input_fname, ecotypeid2pos, norame_lon_span)
noramer_geno_vs_geo_output_fname_prefix = '%s_geno_vs_geo_dist_noramer'%output_fname_prefix
sample_geno_geo_correlation(noramer_geno_dist_ls, noramer_geo_dist_ls, noramer_geno_vs_geo_output_fname_prefix)

eur_noramer_geno_dist_ls, eur_noramer_geo_dist_ls =  calGenoDistAndGeoDistBetweenTwoAreas(input_fname, ecotypeid2pos, europe_lon_span, norame_lon_span)
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
maf_vector = cal_maf_vector(input_fname)
import pylab
pylab.clf()
pylab.plot(range(len(maf_vector)), maf_vector)
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

"""
2007-09-26
	calculate LD (r^2)
"""
def group_ordered_snps_into_chr_snp_2layer_ls(curs, snp_acc_list, snp_locus_table='snps'):
	"""
	2007-09-26
		assume snps are already in chromosome, position order, just need to find out
		where to stop the chromosome
	"""
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
	return chr_snp_2layer_ls, snp_position_ls

def fill_in_snp_allele2index(diploid_allele, allele2index):
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

def calculate_LD(input_fname, curs, snp_locus_table='snps', debug=0):
	"""
	exclude pairs with one or two NAs
	exclude pairs both of who are heterozygous calls (can't figure out the phase)
	(only one of pairs is heterozygous is all right)
	"""
	from FilterStrainSNPMatrix import FilterStrainSNPMatrix
	FilterStrainSNPMatrix_instance = FilterStrainSNPMatrix()
	header, strain_acc_list, category_list, data_matrix = FilterStrainSNPMatrix_instance.read_data(input_fname)
	if debug:
		import pdb
		pdb.set_trace()
	no_of_strains = len(strain_acc_list)
	snp_acc_list = header[2:]
	chr_snp_2layer_ls, snp_position_ls = group_ordered_snps_into_chr_snp_2layer_ls(curs, snp_acc_list, snp_locus_table)
	no_of_chrs = len(chr_snp_2layer_ls)
	r_square_ls = []
	distance_ls = []
	allele_freq_ls = []
	snp_pair_ls = []
	D_ls = []
	D_prime_ls = []
	import Numeric
	for c in range(no_of_chrs):
		no_of_snps = len(chr_snp_2layer_ls[c])
		for i in range(no_of_snps):
			for j in range(i+1, no_of_snps):
				snp1_index = chr_snp_2layer_ls[c][i]
				snp2_index = chr_snp_2layer_ls[c][j]
				counter_matrix = Numeric.zeros([2,2])
				snp1_allele2index = {}
				snp2_allele2index = {}
				for k in range(no_of_strains):
					snp1_allele = data_matrix[k][snp1_index]
					snp2_allele = data_matrix[k][snp2_index]
					if snp1_allele!=0 and snp2_allele!=0 and not (snp1_allele>4 and snp2_allele>4):
						snp1_allele1, snp1_allele2 = fill_in_snp_allele2index(snp1_allele, snp1_allele2index)
						snp2_allele1, snp2_allele2 = fill_in_snp_allele2index(snp2_allele, snp2_allele2index)
						counter_matrix[snp1_allele2index[snp1_allele1],snp2_allele2index[snp2_allele1]] += 1
						counter_matrix[snp1_allele2index[snp1_allele2],snp2_allele2index[snp2_allele2]] += 1
				PA = sum(counter_matrix[0,:])
				Pa = sum(counter_matrix[1,:])
				PB = sum(counter_matrix[:,0])
				Pb = sum(counter_matrix[:,1])
				total_num = float(PA+Pa)
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
				distance = abs(snp_position_ls[snp1_index]-snp_position_ls[snp2_index])
				allele_freq = (min(PA, Pa),min(PB, Pb))
				D_ls.append(D)
				D_prime_ls.append(D_prime)
				r_square_ls.append(r2)
				distance_ls.append(distance)
				allele_freq_ls.append(allele_freq)
				snp_pair_ls.append((snp_acc_list[i], snp_acc_list[j]))
	return D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls


def plot_LD(x_ls, y_ls, title, xlabel, ylabel, max_dist=0, output_fname_prefix=None):
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


"""
input_fname = './script/variation/stock20070919/data_d110_c0_5.tsv'
D_ls, D_prime_ls, r_square_ls, distance_ls, allele_freq_ls, snp_pair_ls = calculate_LD(input_fname, curs, snp_locus_table='snps')

title = 'LD decay'
xlabel = 'Distance'
ylabel = r'$r^2$'
output_fname_prefix = '%s_LD_r2'%os.path.splitext(input_fname)[0]
plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 0, output_fname_prefix)
plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
plot_LD(distance_ls, r_square_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = 'D'
output_fname_prefix = '%s_LD_D'%os.path.splitext(input_fname)[0]
plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 0, output_fname_prefix)
plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
plot_LD(distance_ls, D_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

ylabel = "D'"
output_fname_prefix = '%s_LD_D_prime'%os.path.splitext(input_fname)[0]
plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 0, output_fname_prefix)
plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 500000, output_fname_prefix)
plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 1000000, output_fname_prefix)
plot_LD(distance_ls, D_prime_ls, title, xlabel, ylabel, 5000000, output_fname_prefix)

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

def get_chr_id2cumu_size(chr_id2size):
	"""
	2007-10-16
	"""
	chr_id_ls = chr_id2size.keys()
	chr_id_ls.sort()
	chr_id_ls = [0] + chr_id_ls	#chromosome 0 is a place holder. 
	chr_id2cumu_size = {0:0}	#chr_id_ls might not be continuous integers. so dictionary is better
	for i in range(1,len(chr_id_ls)):
		chr_id = chr_id_ls[i]
		prev_chr_id = chr_id_ls[i-1]
		chr_id2cumu_size[chr_id] = chr_id2cumu_size[prev_chr_id] + chr_id2size[chr_id]
	return chr_id2cumu_size

def DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size):
	"""
	2007-10-16
	"""
	import pylab
	pylab.clf()
	#draw the chromosome separator
	for chr_id, cumu_size in chr_id2cumu_size.iteritems():
		pylab.plot([cumu_size], [1], 'go')
	#draw the snp first as red circles
	for snp_index,pos in snp_index2pos.iteritems():
		chr_id, chr_pos = pos
		cumu_chr_pos = chr_id2cumu_size[chr_id-1]+chr_pos
		pylab.plot([cumu_chr_pos], [1], 'r|')
	#draw the blocks as lines crossing those red circles.
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
from common import get_chr_id2size
chr_id2size = get_chr_id2size(curs)
chr_id2cumu_size = get_chr_id2cumu_size(chr_id2size)
shared_block_ls = get_shared_block_ls(data_matrix[1], data_matrix[2], snp_index2pos)
DrawSharedBlock_ls(shared_block_ls, snp_index2pos, chr_id2cumu_size)

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

hostname='localhost'
dbname='stock20071008'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()

if __name__ == '__main__':
	ecotype_table = 'ecotype'
	calls_table = 'calls'
	strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
	strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
	strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)