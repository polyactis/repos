#!/usr/bin/python
"""
2008-05-15
	functions here largely superceded by more sophiscated MpiQCCall.py and PlotQCCall.py
2008-05-01
	impute-related functions
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))


def transform_oligo_prob_out2std_call_format(input_dir, output_dir, snpid_fname, snpallele_fname):
	"""
	2008-04-30
		transform oligo probability output files (from Rong's new_basecall.R) into standard call format:
			tab-delimited 3-column: SNP_ID, call, probability
	"""
	import os, sys
	inf = open(snpid_fname)
	snpid_ls = inf.readline().split()
	del inf
	from variation.src.common import number2nt
	inf = open(snpallele_fname)
	snpallele_ls = []
	for line in inf:
		allele1_number = int(line[0])
		line = inf.next()	#discard the 2nd line
		line = inf.next()
		allele2_number = int(line[0])
		snpallele_ls.append([number2nt[allele1_number], number2nt[allele2_number]])
		line = inf.next()	#discard the 4th line
	del inf
	sys.stderr.write("Get all alleles for %s SNPs\n"%(len(snpallele_ls)))
	
	import csv
	file_dir_ls = os.listdir(input_dir)
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	for  i in range(len(file_dir_ls)):
		file_dir = file_dir_ls[i]
		sys.stderr.write("%d/%d:\t%s"%(i+1,len(file_dir_ls),file_dir))
		if file_dir[-13:]!='genocalls.txt':
			sys.stderr.write("\tNot Genotype Call file. Ignored.\n")
			continue
		pathname = os.path.join(input_dir, file_dir)
		array_id = int(file_dir.split('_')[0])
		reader = csv.reader(open(pathname), delimiter=' ')
		output_fname = os.path.join(output_dir, '%s_call.tsv'%array_id)
		if os.path.isfile(output_fname):
			sys.stderr.write("\tFile already there. Ignored.\n")
			continue
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		writer.writerow(['SNP_ID', array_id, 'oligo_probability'])
		i = 0
		for row in reader:
			allele_index = int(row[0])
			snpallele = snpallele_ls[i][allele_index/2]
			writer.writerow([snpid_ls[i], snpallele, row[1]])
			i += 1
		del reader, writer
		sys.stderr.write("\n")

input_dir = os.path.expanduser('~/script/oligo/other_output/')
output_dir = os.path.expanduser('~/script/oligo/call_output/')
snpid_fname = os.path.expanduser('~/script/oligo/snp-IDs.txt')
snpallele_fname = os.path.expanduser('~/script/oligo/snp-alleles.txt')
transform_oligo_prob_out2std_call_format(input_dir, output_dir, snpid_fname, snpallele_fname)

sys.path.insert(0, os.path.join(os.path.expanduser('/home/crocea/script')))

def get_array_id2QC(call_method_id=1, **keywords):
	"""
	2008-04-30
		get the array_id 2 QC of calls from simple calling algorithm (call_method_id=1)
	"""
	import sys, os
	sys.stderr.write("Getting array_id2QC ... ")
	from variation.src.db import Stock_250kDatabase, Call_QC
	db = Stock_250kDatabase(username=keywords['user'],
				   password=keywords['passwd'], host=keywords['hostname'], database=keywords['dbname'])
	session = db.session
	array_id2QC = {}
	for row in session.query(Call_QC).list():
		if row.call_info_obj.method_id!=call_method_id:	#simple calling algorithm
			continue
		array_id = row.call_info_obj.array_id
		if array_id in array_id2QC:
			if row.no_of_non_NA_pairs>array_id2QC[array_id].no_of_non_NA_pairs:	#replace the QC if there are more non-NA pairs
				array_id2QC[array_id] = row
		else:
			array_id2QC[array_id] = row
	sys.stderr.write("Done.\n")
	del db
	return array_id2QC

keywords = {'user':'nordborglab',\
		'passwd':'papaya',\
		'hostname': 'papaya.usc.edu',\
		'dbname':'stock_250k'}
array_id2QC = get_array_id2QC(**keywords)

def selectCallFileBasedOnQC(input_dir, output_dir, array_id2QC, max_mismatch_rate=0.15):
	"""
	2008-04-30
	"""
	import sys, os, subprocess
	file_dir_ls = os.listdir(input_dir)
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	for  i in range(len(file_dir_ls)):
		file_dir = file_dir_ls[i]
		sys.stderr.write("%d/%d:\t%s"%(i+1,len(file_dir_ls), file_dir))
		array_id = int(file_dir.split('_')[0])
		if array_id not in array_id2QC:
			sys.stderr.write("\tNo QC available. Ignored.\n")
			continue
		mismatch_rate = array_id2QC[array_id].mismatch_rate
		if mismatch_rate<=max_mismatch_rate:
			pathname = os.path.join(input_dir, file_dir)
			output_fname = os.path.join(output_dir, file_dir)
			os.symlink(pathname, output_fname)
			"""
			cp_p = subprocess.Popen(['cp', pathname, output_fname], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
			cp_p_stdout_out = cp_p.stdout.read()
			cp_p_stderr_out = cp_p.stderr.read()
			if cp_p_stderr_out:
				sys.stderr.write(" Copy Error: %s.\n"%(cp_p_stderr_out))
			else:
				sys.stderr.write(" mismatch_rate=%s. Copied.\n"%(mismatch_rate))
			"""
			sys.stderr.write(" mismatch_rate=%s. Copied.\n"%(mismatch_rate))
		else:
			sys.stderr.write(" mismatch_rate=%s. No copy.\n"%(mismatch_rate))

input_dir = os.path.expanduser('~/script/oligo/call_output/')
output_dir = os.path.expanduser('~/script/oligo/call_output_max_mismatch_rate_less_0.15/')
selectCallFileBasedOnQC(input_dir, output_dir, array_id2QC, max_mismatch_rate=0.15)


def get_chr2no_of_snps(input_fname):
	"""
	2008-04-30
	"""
	import csv, os, sys
	sys.stderr.write("Getting chr2no_of_snps ... ")
	reader = csv.reader(open(input_fname), delimiter='\t')
	reader.next()	#skip the header
	chr2no_of_snps = {}
	for row in reader:
		SNP_id, call, probability = row[:3]
		tmp_ls = SNP_id.split('_')
		chr = int(tmp_ls[0])
		position = int(tmp_ls[1])
		if chr not in chr2no_of_snps:
			chr2no_of_snps[chr] = 0
		chr2no_of_snps[chr] += 1
	del reader
	sys.stderr.write("Done.\n")
	return chr2no_of_snps


def transformStdCalls2NPUTEInput(input_dir, output_dir, column_record_filename, min_probability=-1):
	"""
	2008-05-01
		prepare input files for NPUTE
		read call files from input_dir and split into different chromosomes for NPUTE input
		
		Output 1: SNP by Individual, genotype matrix
		Output 2: column_record_filename contains the input call filenames for each column.
	2008-04-30
	"""
	import sys, os
	import csv
	file_dir_ls = os.listdir(input_dir)
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	no_of_files = len(file_dir_ls)
	chr2no_of_snps = get_chr2no_of_snps(os.path.join(input_dir, file_dir_ls[0]))
	
	chr_ls = chr2no_of_snps.keys()
	chr_ls.sort()
	chr2index = {}
	i = 0
	for chr in chr_ls:
		chr2index[chr] = i
		i += 1
	#open an output file for each chromosome
	outf_ls = []
	for chr in chr_ls:	#store each chromosome into different files
		output_fname = os.path.join(output_dir, 'impute_input_chr%s'%chr)
		outf = open(output_fname, 'w')
		outf_ls.append(outf)
	
	#open all input files.
	sys.stderr.write("Open all input files ...\n")
	reader_ls = []
	column_record_f = open(column_record_filename, 'w')
	for  i in range(no_of_files):
		file_dir = file_dir_ls[i]
		sys.stderr.write("%d/%d:\t%s\n"%(i+1, no_of_files, file_dir))
		column_record_f.write('%s\n'%file_dir)
		pathname = os.path.join(input_dir, file_dir)
		reader = csv.reader(open(pathname), delimiter='\t')
		reader.next()	#skip the header (1st line) in each file
		reader_ls.append(reader)
	del column_record_f
	sys.stderr.write("Done.\n")
	
	
	#All input files are in same SNP order. For one SNP, read one line from each file and output them into corresponding chromosome outputfile
	sys.stderr.write("Outputting ...\n")
	no_of_snps = sum(chr2no_of_snps.values())
	for i in range(no_of_snps):
		index = None	#which file in outf_ls for this SNP
		if i%5000==0:
			sys.stderr.write("%s\tSNP %s/%s"%('\x08'*80, i+1, no_of_snps))
		for reader in reader_ls:
			row = reader.next()
			SNP_id, call, probability = row[:3]
			probability = float(probability)
			if probability<min_probability:
				call = '?'
			if call=='NA':
				call = '?'
			#figure out which output file based on chromosme
			if index==None:
				tmp_ls = SNP_id.split('_')
				chr = int(tmp_ls[0])
				#position = int(tmp_ls[1])
				index = chr2index[chr]
			outf_ls[index].write("%s"%call)
		#write the line separator
		outf_ls[index].write("\n")
	sys.stderr.write("%s\tSNP %s/%s"%('\x08'*80, i+1, no_of_snps))
	sys.stderr.write(".Done\n")
	del reader_ls
	del outf_ls


input_dir = os.path.expanduser('~/script/oligo/call_output_max_mismatch_rate_less_0.15/')
output_dir = os.path.expanduser('~/script/oligo/call_output_max_mismatch_rate_less_0.15_NPUTE_min_prob_0.85/')
column_record_filename = os.path.expanduser('~/script/oligo/call_output_max_mismatch_rate_less_0.15_NPUTE_min_prob_0.85.columns')
min_probability = 0.85
transformStdCalls2NPUTEInput(input_dir, output_dir, column_record_filename, min_probability)


def transformStdCalls2fastPHASEInput(input_dir, output_dir, min_probability=-1):
	"""
	2008-05-01
		for fastPHASE
	2008-04-30
	"""
	import sys, os
	import csv
	file_dir_ls = os.listdir(input_dir)
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	no_of_files = len(file_dir_ls)
	chr2no_of_snps = get_chr2no_of_snps(os.path.join(input_dir, file_dir_ls[0]))
	
	chr_ls = chr2no_of_snps.keys()
	chr_ls.sort()
	chr2index = {}
	i = 0
	for chr in chr_ls:
		chr2index[chr] = i
		i += 1
	outf_ls = []
	for chr in chr_ls:	#store each chromosome into different files
		output_fname = os.path.join(output_dir, 'impute_input_chr%s'%chr)
		outf = open(output_fname, 'w')
		outf.write('%s\n'%no_of_files)
		outf.write('%s\n'%chr2no_of_snps[chr])
		outf_ls.append(outf)
	
	for  i in range(len(file_dir_ls)):
		file_dir = file_dir_ls[i]
		sys.stderr.write("%d/%d:\t%s"%(i+1,len(file_dir_ls),file_dir))
		pathname = os.path.join(input_dir, file_dir)
		array_id = int(file_dir.split('_')[0])
		reader = csv.reader(open(pathname), delimiter='\t')
		
		SNP_per_chr_ls = []
		for chr in chr_ls:
			SNP_per_chr_ls.append([])
			outf_ls[chr2index[chr]].write('# id %s\n'%array_id)
		
		reader.next()	#skip the header
		for row in reader:
			SNP_id, call, probability = row[:3]
			probability = float(probability)
			if probability<min_probability:
				call = '?'
			if call=='NA':
				call = '?'
			tmp_ls = SNP_id.split('_')
			chr = int(tmp_ls[0])
			position = int(tmp_ls[1])
			index = chr2index[chr]
			SNP_per_chr_ls[index].append(call)
		
		#output into each file
		for chr in chr_ls:
			index = chr2index[chr]
			snp_str = ''.join(SNP_per_chr_ls[index])
			outf_ls[index].write('%s\n'%snp_str)
			#fastPHASE deals with diploid and solves the phasing issue as well.
			outf_ls[index].write('%s\n'%snp_str)
		del reader
		sys.stderr.write(".\n")
	del outf_ls
	sys.stderr.write("\n")

input_dir = os.path.expanduser('~/script/oligo/call_output_max_mismatch_rate_less_0.15/')
output_dir = os.path.expanduser('~/script/oligo/call_output_max_mismatch_rate_less_0.15_fastPHASE_min_prob_0.85/')
min_probability = 0.85
transformStdCalls2fastPHASEInput(input_dir, output_dir, min_probability)

def transformfastPhaseOutput2stdCallFormat(input_dir, output_dir, snpid_fname, snpallele_fname):
	"""
	2008-04-30
		just copied from transformStdCalls2fastPHASEInput, not working yet
	"""
	import sys, os
	import csv
	file_dir_ls = os.listdir(input_dir)
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	no_of_files = len(file_dir_ls)
	chr2no_of_snps = get_chr2no_of_snps(os.path.join(input_dir, file_dir_ls[0]))
	
	chr_ls = chr2no_of_snps.keys()
	chr_ls.sort()
	chr2index = {}
	i = 0
	for chr in chr_ls:
		chr2index[chr] = i
		i += 1
	outf_ls = []
	for chr in chr_ls:	#store each chromosome into different files
		output_fname = os.path.join(output_dir, 'fastPhase_input_chr%s'%chr)
		outf = open(output_fname, 'w')
		outf.write('%s\n'%no_of_files)
		outf.write('%s\n'%chr2no_of_snps[chr])
		outf_ls.append(outf)
	
	for  i in range(len(file_dir_ls)):
		file_dir = file_dir_ls[i]
		sys.stderr.write("%d/%d:\t%s"%(i+1,len(file_dir_ls),file_dir))
		pathname = os.path.join(input_dir, file_dir)
		array_id = int(file_dir.split('_')[0])
		reader = csv.reader(open(pathname), delimiter='\t')
		
		SNP_per_chr_ls = []
		for chr in chr_ls:
			SNP_per_chr_ls.append([])
			outf_ls[chr2index[chr]].write('# id %s\n'%array_id)
		
		reader.next()	#skip the header
		for row in reader:
			SNP_id, call, probability = row[:3]
			probability = float(probability)
			if probability<min_probability:
				call = '?'
			if call=='NA':
				call = '?'
			tmp_ls = SNP_id.split('_')
			chr = int(tmp_ls[0])
			position = int(tmp_ls[1])
			index = chr2index[chr]
			SNP_per_chr_ls[index].append(call)
		
		#output into each file
		for chr in chr_ls:
			index = chr2index[chr]
			snp_str = ''.join(SNP_per_chr_ls[index])
			outf_ls[index].write('%s\n'%snp_str)
			outf_ls[index].write('%s\n'%snp_str)
		del reader
		sys.stderr.write(".\n")
	del outf_ls
	sys.stderr.write("\n")


def get_min_prob2call_mismatch_rate_ls(input_dir, QC_method_id=1):
	"""
	2008-05-15 make up documentation for this function
		take the output of QC_250k.py (which takes the output of modified NPUTE.py), return stat for plot_min_prob2NA_mismatch_rate_ls()
		
	"""
	import os,sys, csv
	file_dir_ls = os.listdir(input_dir)
	min_prob2NA_mismatch_rate_ls = {}
	import re
	min_prob_p = re.compile(r'_y(0.\d+)_')
	for  i in range(len(file_dir_ls)):
		file_dir = file_dir_ls[i]
		sys.stderr.write("%d/%d:\t%s"%(i+1,len(file_dir_ls),file_dir))
		if file_dir.find('_m%s_QC.tsv'%QC_method_id)<0:
			sys.stderr.write("\tignore.\n")
			continue
		#tmp_ls = file_dir[:-4].split('_')
		#min_prob = tmp_ls[-1][1:]
		min_prob = min_prob_p.search(file_dir).groups()[0]
		pathname = os.path.join(input_dir, file_dir)
		reader = csv.reader(open(pathname), delimiter='\t')
		reader.next()
		_no_of_NAs = 0
		_no_of_totals = 0
		_no_of_mismatches = 0
		_no_of_non_NA_pairs = 0
		for row in reader:
			array_id, ecotypeid, NA_rate, mismatch_rate, no_of_NAs, no_of_totals, no_of_mismatches, no_of_non_NA_pairs =row
			_no_of_NAs += int(no_of_NAs)
			_no_of_totals += int(no_of_totals)
			_no_of_mismatches += int(no_of_mismatches)
			_no_of_non_NA_pairs += int(no_of_non_NA_pairs)
		NA_rate = float(_no_of_NAs)/_no_of_totals
		mismatch_rate = float(_no_of_mismatches)/_no_of_non_NA_pairs
		min_prob2NA_mismatch_rate_ls[min_prob] = [NA_rate, mismatch_rate]
		del reader
		sys.stderr.write(".\n")
	return min_prob2NA_mismatch_rate_ls

input_dir = os.path.expanduser('~/script/oligo/QC/')
input_dir = '/mnt/nfs/NPUTE_data/'
min_prob2NA_mismatch_rate_ls = get_min_prob2call_mismatch_rate_ls(input_dir)

def plot_min_prob2NA_mismatch_rate_ls(min_prob2NA_mismatch_rate_ls, output_fname, show_min_prob=1):
	import pylab
	pylab.clf()
	call_ls = []
	mismatch_ls = []
	min_prob_ls = min_prob2NA_mismatch_rate_ls.keys()
	min_prob_ls.sort()
	for min_prob, NA_mismatch in min_prob2NA_mismatch_rate_ls.iteritems():
		if show_min_prob:
			call_ls.append(min_prob)
		else:
			call_ls.append(1-NA_mismatch[0])
		mismatch_ls.append(1-NA_mismatch[1])
	pylab.plot(call_ls, mismatch_ls, '.')
	if show_min_prob:
		pylab.xlabel('min probability')
	else:
		pylab.xlabel('call rate')
	pylab.ylabel('accuracy')
	pylab.savefig('%s.png'%output_fname, dpi=80)
	pylab.savefig('%s.svg'%output_fname, dpi=80)


output_fname = os.path.expanduser('~/script/oligo/QC/call_rate_vs_mismatch_rate')
output_fname = os.path.expanduser('/tmp/call_rate_vs_mismatch_rate')
qc_method_id_
for i in [1,2,3,8]:
	output_fname = '/tmp/min_prob_vs_mismatch_rate_qc_%s'
	plot_min_prob2NA_mismatch_rate_ls(min_prob2NA_mismatch_rate_ls, output_fname)