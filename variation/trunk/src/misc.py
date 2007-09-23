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
	pylab.imshow(data_matrix, aspect='auto', interpolation='nearest')	#2007-06-05
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

"""
DrawStrain_Heterozygotes_PercHist('./script/variation/data/justin_data_y.csv', './script/variation/data/justin_data_y', 1)
"""

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

def divide_data_by_geography(cursor, lat_lon_ls, max_dist=100):
	import sys, os
	sys.stderr.write("Constructing data structure on lat_lon_ls ...")
	#each distinctive position (lat,lon) gets a unique label (an integer starting from 0)
	pos2node_label = {}
	node_label2pos_counts = {}	#records how many strains are from the same position
	for lat, lon in lat_lon_ls:
		pos = (lat, lon)
		if pos not in pos2node_label:
			pos2node_label[pos] = len(pos2node_label)
		node_label = pos2node_label[pos]
		if node_label not in node_label2pos_counts:
			node_label2pos_counts[node_label] = [pos, 0]
		node_label2pos_counts[node_label][1] += 1
	sys.stderr.write("Done.\n")
	sys.stderr.write("Constructing graph ...")
	import networkx as nx
	g = nx.Graph()
	pos_ls = pos2node_label.keys()
	no_of_sites = len(pos_ls)
	for i in range(no_of_sites):
		g.add_node(pos2node_label[pos_ls[i]])
	
	for i in range(no_of_sites):
		for j in range(i+1, no_of_sites):
			dist = cal_great_circle_distance(pos_ls[i][0], pos_ls[i][1], pos_ls[j][0], pos_ls[j][1])
			if dist<=max_dist:
				g.add_edge(pos2node_label[pos_ls[i]], pos2node_label[pos_ls[j]])
	c_components = nx.connected_components(g)
	sys.stderr.write("Done.\n")
	return g, c_components, node_label2pos_counts

#2007-06-18
#draw populations derived from connected_components of the strain network
#each pie denotes a population, with diameter proportional to the size of the population
#each pie labeled with the number of strains in that population
def draw_clustered_strain_location(c_components, node_label2pos_counts, pic_area=[-180,-90,180,90]):
	import sys, os
	sys.stderr.write("Computing weighted centers...")
	weighted_pos_ls = []
	count_sum_ls = []
	for component in c_components:
		#weighted average
		count_sum = 0
		lat_sum = 0.0
		lon_sum = 0.0
		for node_label in component:
			count_sum += node_label2pos_counts[node_label][1]
			lat_sum += node_label2pos_counts[node_label][0][0]*node_label2pos_counts[node_label][1]
			lon_sum += node_label2pos_counts[node_label][0][1]*node_label2pos_counts[node_label][1]
		weighted_pos_ls.append((lat_sum/count_sum, lon_sum/count_sum))
		count_sum_ls.append(count_sum)
	sys.stderr.write("Done.\n")
	sys.stderr.write("Drawing on map...")
	import pylab
	from matplotlib.toolkits.basemap import Basemap
	pylab.clf()
	fig = pylab.figure()
	fig.add_axes([0.1,0.1,0.9,0.9])
	m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
            resolution='l',projection='mill')
	#m.drawcoastlines()
	m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,1,1])
	m.drawmeridians(pylab.arange(-180,180,60), labels=[1,1,1,1])
	m.fillcontinents()
	m.drawcountries()
	m.drawstates()
	euc_coord1_ls = []
	euc_coord2_ls = []
	ax=pylab.gca()
	for i in range(len(weighted_pos_ls)):
		lat, lon = weighted_pos_ls[i]
		euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
		euc_coord1_ls.append(euc_coord1)
		euc_coord2_ls.append(euc_coord2)
		ax.text(euc_coord1, euc_coord2, str(count_sum_ls[i]), size=6, alpha=0.5, horizontalalignment='center', verticalalignment='center', zorder=12)
	m.scatter(euc_coord1_ls, euc_coord2_ls, 5*count_sum_ls, marker='o', color='r', alpha=0.3, zorder=10)
	pylab.title("worldwide distribution of %s locations"%(len(weighted_pos_ls)))
	pylab.show()
	sys.stderr.write("Done.\n")

#2007-07-09
def DrawStrainNetwork(g, node_label2pos_counts,pic_area=[-180,-90,180,90]):
	import sys, os
	sys.stderr.write("Drawing Strain Network...")
	import pylab
	from matplotlib.toolkits.basemap import Basemap
	pylab.clf()
	#fig = pylab.figure()
	#fig.add_axes([0.1,0.1,0.9,0.9])
	m = Basemap(llcrnrlon=pic_area[0],llcrnrlat=pic_area[1],urcrnrlon=pic_area[2],urcrnrlat=pic_area[3],\
            resolution='l',projection='mill')
	#m.drawcoastlines()
	m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,1,1])
	m.drawmeridians(pylab.arange(-180,180,60), labels=[1,1,1,1])
	m.fillcontinents()
	m.drawcountries()
	m.drawstates()
	ax=pylab.gca()
	for e in g.edges():
		lat1, lon1 = node_label2pos_counts[e[0]][0]
		lat2, lon2 = node_label2pos_counts[e[1]][0]
		x1, y1 = m(lon1, lat1)
		x2, y2 = m(lon2, lat2)
		ax.plot([x1,x2],[y1,y2], alpha=0.5, zorder=12)
	pylab.title("Network of strains")
	pylab.show()
	sys.stderr.write("Done.\n")


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

def construct_site_graph_out_of_strain_graph(g, ecotypeid2pos):
	"""
	2007-09-20
		#reduce the drawing by omitting identical edges (linking same GPS positions), construct a graph whose node is (lat,lon) and then plot this graph.
	"""
	sys.stderr.write("Constructing site graph out of strain graph...")
	site2pos = {}
	site2weight = {}
	import networkx as nx
	site_g = nx.XGraph()
	for e in g.edges():
		pos1 = (round(ecotypeid2pos[e[0]][0],2), round(ecotypeid2pos[e[0]][1],2))
		pos2 = (round(ecotypeid2pos[e[1]][0],2), round(ecotypeid2pos[e[1]][1],2))
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
	
	for n in site_g:
		site2pos[n] = n
	sys.stderr.write("Done.\n")
	return site_g, site2weight, site2pos

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

def get_ecotypeid2pos(curs, ecotype_table):
	"""
	2007-09-18
	"""
	sys.stderr.write("Getting ecotypeid2pos from %s ..."%ecotype_table)
	ecotypeid2pos = {}
	curs.execute("select id, latitude, longitude from %s where latitude is not null and longitude is not null"%ecotype_table)
	rows = curs.fetchall()
	for row in rows:
		ecotypeid, latitude, longitude = row
		ecotypeid2pos[ecotypeid] = (latitude, longitude)
	sys.stderr.write("Done.\n")
	return ecotypeid2pos

def get_popid2pos_size_and_ecotypeid2pop_id(curs, popid2ecotypeid_table):
	"""
	2007-09-13
	"""
	sys.stderr.write("Getting popid2pos_size & ecotypeid2pop_id ...")
	popid2pos_size = {}
	ecotypeid2pop_id = {}
	curs.execute("select popid, latitude, longitude, ecotypeid from %s where selected=1"%(popid2ecotypeid_table))
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

def draw_graph_on_map(g, node2weight, node2pos, pic_title,  pic_area=[-130,10,140,70], output_fname_prefix=None):
	"""
	2007-09-13
		identity_pair_ls is a list of pairs of strains (ecotype id as in table ecotype)
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
	
	ax=pylab.gca()
	
	sys.stderr.write("\tDrawing nodes ...")
	euc_coord1_ls = []
	euc_coord2_ls = []
	diameter_ls = []
	for n in g.nodes():
		lat, lon = node2pos[n]
		euc_coord1, euc_coord2 = m(lon, lat)	#longitude first, latitude 2nd
		euc_coord1_ls.append(euc_coord1)
		euc_coord2_ls.append(euc_coord2)
		diameter_ls.append(node2weight[n])
	m.scatter(euc_coord1_ls, euc_coord2_ls, diameter_ls, marker='o', color='r', alpha=0.4, zorder=12, faceted=False)
	sys.stderr.write("Done.\n")
	
	sys.stderr.write("\tDrawing edges ...")
	for popid1, popid2, no_of_connections in g.edges():
		lat1, lon1 = node2pos[popid1]
		lat2, lon2 = node2pos[popid2]
		x1, y1 = m(lon1, lat1)
		x2, y2 = m(lon2, lat2)
		ax.plot([x1,x2],[y1,y2], 'g', linewidth=math.log(no_of_connections+1)/2, alpha=0.2, zorder=10)
	sys.stderr.write("Done.\n")
	
	#m.drawcoastlines()
	m.drawparallels(pylab.arange(-9,90,30), labels=[1,1,0,1])
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
"""
input_fname = 'script/variation/stock20070829/data_row_na_col_na_bad_snps.tsv'
identity_pair_ls = construct_identity_pair_ls(input_fname)
popid2ecotypeid_table = 'popid2ecotypeid_25'
pop_iden_g, popid2intra_pop_connections, popid2pos_size, ecotypeid2pop_id = construct_pop_identity_graph_from_strain_identity_ls(curs, popid2ecotypeid_table, identity_pair_ls)

popid2pos = {}
for popid, pos_size in popid2pos_size.iteritems():
	popid2pos[popid] = pos_size[0]

import math
popid2weight = {}
for popid, intra_pop_connections in popid2intra_pop_connections.iteritems():
	popid2weight[popid] = 8*math.log(intra_pop_connections+1)

draw_graph_on_map(pop_iden_g, popid2weight, popid2pos, 'inter population identity map '+popid2ecotypeid_table, pic_area=[-130,34,35,65], output_fname_prefix='/tmp/identity.map')


strain_iden_g = construct_graph_out_of_identity_pair(identity_pair_ls)
ecotypeid2pos = get_ecotypeid2pos(curs, 'ecotype')
ecotypeid2weight = {}
for n in strain_iden_g:
	ecotypeid2weight[n] = 1

import networkx as nx
strain_iden_g_cc = nx.connected_components(strain_iden_g)

i=11


for i in range(len(strain_iden_g_cc)):
	gs = strain_iden_g.subgraph(strain_iden_g_cc[i])
	site_g, site2weight, site2pos = construct_site_graph_out_of_strain_graph(gs, ecotypeid2pos)
	for site, weight in site2weight.iteritems():
		site2weight[site] =  8*math.log(weight+1)
	display_matrix_of_component(input_fname, strain_iden_g_cc[i], ecotypeid2pos, output_fname='ecotype identity map cc%s'%i, need_sort=0, need_savefig=1)
	draw_graph_on_map(site_g, site2weight, site2pos, 'ecotype identity map cc%s'%i, pic_area=[-130,34,35,65], output_fname_prefix='ecotype_identity_cc%s'%i)
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
	else:
		sys.stderr.write("unsupported strainname_type %s\n"%strainname_type)
		return None, None
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
	sys.stderr.write('%s%s'%('\x08'*20, counter))
	sys.stderr.write("\nDone.\n")
	print "%s strain snp pairs"%(len(strain_snp_pair_set))
	print 'inconsistent ratio: %s/%s=%s'%(no_of_inconsistent, no_of_strain_snp_pairs, float(no_of_inconsistent)/no_of_strain_snp_pairs)
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
print '(nativename,snpid) pair:'
strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
print '(ecotypeid,snpid) pair:'
strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
print '(nativename, stockparent ,snpid) pair:'
strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)

nativename_snp_distinct_set = get_nativename_snp_distinct_set(curs, ecotype_table, calls_table)

strain_snp_pair_set - nativename_snp_distinct_set gives 1937 different pairs. turns out that sql's "distinct e.nativename, c.snpid" ignores the cases , like Kz-9 and KZ-9 are same. however python doesn't do that.
"""

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
#conn, curs = db_connect(hostname, dbname, schema)

hostname='localhost'
dbname='stock20070829'
import MySQLdb
conn0 = MySQLdb.connect(db=dbname,host=hostname)
curs0 = conn0.cursor()

hostname='localhost'
dbname='stock20070919'
import MySQLdb
conn = MySQLdb.connect(db=dbname,host=hostname)
curs = conn.cursor()

if __name__ == '__main__':
	ecotype_table = 'ecotype'
	calls_table = 'calls'
	print '(nativename,snpid) pair:'
	strain_snp_pair_set1, inconsistent_dup_strain_snp_pair_set1 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table)
	print '(ecotypeid,snpid) pair:'
	strain_snp_pair_set2, inconsistent_dup_strain_snp_pair_set2 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=2)
	print '(nativename, stockparent ,snpid) pair:'
	strain_snp_pair_set3, inconsistent_dup_strain_snp_pair_set3 = check_inconsistent_duplicate_calls2(curs, ecotype_table, calls_table, strainname_type=3)