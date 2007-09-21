"""
2007-09-11 compare calls between justin genotype array and perlegen data
"""
def get_SNP_probe_ls(snp_loci_fname):
	from common import nt2number, number2nt
	import Numeric
	import csv, sys, os
	sys.stderr.write("Reading snp loci data ...")
	SNP_probe_ls = []
	SNP_loci2index = {}
	reader = csv.reader(open(snp_loci_fname), delimiter='\t')
	reader.next()	#skip header
	for row in reader:
		chr, pos, allele = row[3:6]
		SNP_probe_ls.append([int(chr), int(pos), allele])
		SNP_loci = (int(chr), int(pos))
		if SNP_loci not in SNP_loci2index:	#every 4 probes, a new locus
			SNP_loci2index[SNP_loci] = len(SNP_loci2index)
	del reader
	sys.stderr.write("Done.\n")
	return SNP_probe_ls, SNP_loci2index

def get_strain2index(ordered_20_strain_fname):
	from common import nt2number, number2nt
	import Numeric
	import csv, sys, os
	sys.stderr.write("Getting strain2index ...")
	strain2index = {}
	reader = csv.reader(open(ordered_20_strain_fname), delimiter='\t')
	for row in reader:
		strain_name = row[0].strip()[:-1]
		strain2index[strain_name] = len(strain2index)
	del reader
	sys.stderr.write("Done.\n")
	return strain2index

def get_perlegen_call_matrix(SNP_loci2index, strain2index, perlegen_raw_fname):
	from common import nt2number, number2nt
	import Numeric
	import csv, sys, os
	sys.stderr.write("Reading perlegen data matrix ...")
	perlegen_call_matrix = Numeric.zeros([len(SNP_loci2index), len(strain2index)])
	reader = csv.reader(open(perlegen_raw_fname), delimiter=',')
	reader.next()	#skip the header
	for row in reader:
		chr, pos = row[1:3]
		calls = row[9]
		SNP_loci = (int(chr), int(pos))
		if SNP_loci in SNP_loci2index:
			index = SNP_loci2index[SNP_loci]
			for i in range(len(calls)):
				perlegen_call_matrix[index][i] = nt2number[calls[i]]
	del reader
	sys.stderr.write("Done.\n")
	return perlegen_call_matrix

def make_simple_justin_calls(SNP_loci2index, intensity_col_name_ls, SNP_intensity_matrix_fname, SNP_probe_ls):
	from common import nt2number, number2nt
	import Numeric
	import csv, sys, os
	sys.stderr.write("Making simple calls on justin data ...")
	justin_call_matrix = Numeric.zeros([len(SNP_loci2index), len(intensity_col_name_ls)])
	reader= csv.reader(open(SNP_intensity_matrix_fname), delimiter='\t')
	probe_count_looper = 0
	SNP_locus_counter = 0
	justin_intensity_4probe_matrix = Numeric.zeros([4, len(intensity_col_name_ls)])
	reader.next()	#skip header
	for row in reader:
		row = map(int, row[1:])
		justin_intensity_4probe_matrix[probe_count_looper] = row
		probe_count_looper += 1
		if probe_count_looper==4:	#one SNP locus every 4 probes
			probe_count_looper = 0
			for i in range(len(intensity_col_name_ls)):
				if justin_intensity_4probe_matrix[0][i]>justin_intensity_4probe_matrix[2][i] and justin_intensity_4probe_matrix[1][i]>justin_intensity_4probe_matrix[3][i]:
					justin_call_matrix[SNP_locus_counter][i] = nt2number[SNP_probe_ls[SNP_locus_counter*4][2]]
				elif justin_intensity_4probe_matrix[0][i]<justin_intensity_4probe_matrix[2][i] and justin_intensity_4probe_matrix[1][i]<justin_intensity_4probe_matrix[3][i]:
					justin_call_matrix[SNP_locus_counter][i] = nt2number[SNP_probe_ls[SNP_locus_counter*4+2][2]]
			SNP_locus_counter += 1
	del reader
	sys.stderr.write("Done.\n")
	return justin_call_matrix

snp_loci_fname = '/tmp/snp'
ordered_20_strain_fname = '/Network/Servers/oak.usc.edu/Volumes/RAID/Data/Perlegen_data/tmp'
perlegen_raw_fname = '/Network/Servers/oak.usc.edu/Volumes/RAID/Data/Perlegen_data/all.perlegen.snps'
intensity_col_name_ls = ['col', 'col', 'ler', 'ler', 'van', 'van']	#for yanli_8-8-07
SNP_intensity_matrix_fname = '/tmp/yanli_8-8-07.rawdata'	#this is the raw data read by readcel_raw() from readcel.R

SNP_probe_ls, SNP_loci2index = get_SNP_probe_ls(snp_loci_fname)
strain2index = get_strain2index(ordered_20_strain_fname)
perlegen_call_matrix = get_perlegen_call_matrix(SNP_loci2index, strain2index, perlegen_raw_fname)
justin_call_matrix = make_simple_justin_calls(SNP_loci2index, intensity_col_name_ls, SNP_intensity_matrix_fname, SNP_probe_ls)

SNP_intensity_matrix_fname = '/tmp/6atSNPtilxarrays8-29-07.rawdata'	#for 6 at SNPtilx arrays 8-29-07
intensity_col_name_ls = ['bay', 'bor', 'br', 'bur', 'c24', 'col']
justin_call_matrix2 = make_simple_justin_calls(SNP_loci2index, intensity_col_name_ls, SNP_intensity_matrix_fname, SNP_probe_ls)

def compare_SNPcalls(strain2index, perlegen_call_matrix, justin_call_matrix, intensity_col_name_ls):
	from common import nt2number, number2nt
	import Numeric
	import csv, sys, os
	
	for j in range(len(intensity_col_name_ls)):
		perlegen_col_index = strain2index[intensity_col_name_ls[j]]
		diff_tag2counter = {}
		diff_tag_dict = {
			'NA_vs_NA':2,
			'NA_vs_call':3,
			'call_vs_NA':4,
			'call_ineq_call':5,
			'call_eq_call':6}	#the value is not useful here
		for tag in diff_tag_dict:
			diff_tag2counter[tag] = 0
		for i in range(len(SNP_loci2index)):
			if perlegen_call_matrix[i,perlegen_col_index]==0 and justin_call_matrix[i,j]==0:
				tag = 'NA_vs_NA'
			elif perlegen_call_matrix[i,perlegen_col_index]==0 and justin_call_matrix[i,j]!=0:
				tag = 'NA_vs_call'
			elif perlegen_call_matrix[i,perlegen_col_index]!=0 and justin_call_matrix[i,j]==0:
				tag = 'call_vs_NA'
			elif perlegen_call_matrix[i,perlegen_col_index] != justin_call_matrix[i,j]:
				tag = 'call_ineq_call'
			elif perlegen_call_matrix[i,perlegen_col_index] == justin_call_matrix[i,j]:
				tag = 'call_eq_call'
			diff_tag2counter[tag] += 1
		print diff_tag2counter
		total_calls = len(SNP_loci2index)
		for tag in diff_tag2counter:
			print tag, ':', float(diff_tag2counter[tag])/total_calls

