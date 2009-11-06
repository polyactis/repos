"""
2009-10-31
	module for CNV (copy-number-variation)-related functions & classes
"""

import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
from utils import dict_map, importNumericArray, figureOutDelimiter, PassingData
from SNP import GenomeWideResult, DataObject
import fileinput
from utils import getColName2IndexFromHeader

def getCNVDataFromFileInGWA(input_fname_ls, array_id, max_amp=-0.33, min_amp=-0.33, min_size=50, min_no_of_probes=None, report=False):
	"""
	2009-10-31
		get deletion (below max_amp) or duplication (above min_amp) from files (output by RunGADA.py)
	"""
	sys.stderr.write("Getting CNV calls for array %s, min_size %s, min_no_of_probes %s from %s ..."%\
					(array_id, min_size, min_no_of_probes, repr(input_fname_ls)))
	
	gwr_name = "(a-id %s)"%(array_id)
	gwr = GenomeWideResult(name=gwr_name)
	gwr.data_obj_ls = []	#list and dictionary are crazy references.
	gwr.data_obj_id2index = {}
	genome_wide_result_id = id(gwr)
	
	amp_ls = []
	array_id2array = {}
	counter = 0
	real_counter = 0
	no_of_segments = 0
	input_handler = fileinput.input(input_fname_ls)
	header = input_handler.readline().strip().split('\t')
	col_name2index = getColName2IndexFromHeader(header)
	ecotype_id = None
	for line in input_handler:
		if line.find("array_id")!=-1:
			continue
		line = line.strip()
		row = line.split('\t')
		cnv_array_id = int(row[col_name2index['array_id']])
		cnv_ecotype_id = int(row[col_name2index.get('ecotype_id', col_name2index['array_id'])])
		counter += 1
		if cnv_array_id==array_id:
			no_of_segments += 1
			if ecotype_id is None:
				ecotype_id = cnv_ecotype_id
			start_probe = row[col_name2index['start_probe']].split('_')	# split chr_pos
			start_probe = map(int, start_probe)
			start_probe_id = row[col_name2index.get('start_probe_id', col_name2index['start_probe'])]
			
			stop_probe = row[col_name2index['end_probe']].split('_')
			stop_probe = map(int, stop_probe)
			end_probe_id = row[col_name2index.get('end_probe_id', col_name2index['end_probe'])]
			
			no_of_probes = int(row[col_name2index['length']])
			if min_no_of_probes is not None and no_of_probes<min_no_of_probes:
				continue
			amplitude = float(row[col_name2index['amplitude']])
			segment_chromosome = start_probe[0]
			segment_start_pos = start_probe[1]-12
			segment_stop_pos = stop_probe[1]+12
			segment_length = abs(segment_stop_pos-segment_start_pos)
			if min_size is not None and segment_length<min_size:
				continue
			if amplitude<=max_amp or amplitude>=min_amp:
				real_counter += 1
				data_obj = DataObject(chromosome=segment_chromosome, position=segment_start_pos, stop_position=segment_stop_pos, \
									value=amplitude)
				data_obj.comment = 'start probe-id %s, end probe-id %s, no of probes %s'%\
							(start_probe_id, end_probe_id, no_of_probes)
				data_obj.genome_wide_result_id = genome_wide_result_id
				gwr.add_one_data_obj(data_obj)
				
		if report and counter%10000==0:
			sys.stderr.write('%s%s\t%s\t%s'%('\x08'*80, counter, no_of_segments, real_counter))
	sys.stderr.write("\n")
	
	if gwr.max_value<3:	# insertion at y=3
		gwr.max_value=3
	if gwr.min_value>-1:	# deletion at y = -1
		gwr.min_value = -1
	gwr.name = '%s '%ecotype_id +  gwr.name
	setattr(gwr, 'ecotype_id', ecotype_id)
	sys.stderr.write(" %s segments. Done.\n"%(len(gwr.data_obj_ls)))
	return gwr