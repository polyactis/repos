#!/usr/bin/env python
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

def get_overlap_ratio(span1_ls, span2_ls):
	"""
	2009-12-13
		calculate the two overlap ratios for two segments
	"""
	segment_start_pos, segment_stop_pos = span1_ls
	qc_start, qc_stop = span2_ls
	overlap_length = max(0, segment_stop_pos - qc_start) - max(0, segment_stop_pos - qc_stop) - max(0, segment_start_pos - qc_start)	# accomodates 6 scenarios
	overlap_length = float(overlap_length)
	overlap1 = overlap_length/(qc_stop-qc_start)
	overlap2 = overlap_length/(segment_stop_pos-segment_start_pos)
	return overlap1, overlap2

def is_reciprocal_overlap(span1_ls, span2_ls, min_reciprocal_overlap=0.6):
	"""
	2009-12-12
		return True if both overlap ratios are above the min_reciprocal_overlap
	"""
	overlap1, overlap2 = get_overlap_ratio(span1_ls, span2_ls)
	if overlap1>=min_reciprocal_overlap and overlap2>=min_reciprocal_overlap:
		return True
	else:
		return False

class CNVSegmentBinarySearchTreeKey(object):
	"""
	2009-12-12
		a key designed to represent a CNV segment in a binary search tree (BinarySearchTree.py), which could be used to do == or >, or < operators
	"""
	def __init__(self, chromosome=None, span_ls=None, min_reciprocal_overlap=0.6):
		self.chromosome = chromosome
		self.span_ls = span_ls
		self.min_reciprocal_overlap = min_reciprocal_overlap
	
	def __lt__(self, other):
		"""
		2009-12-12
			whether self is less than other, in chromosomal order
		"""
		if self.chromosome==other.chromosome:
			if len(self.span_ls)==1:
				if len(other.span_ls)==1:
					return self.span_ls[0]<other.span_ls[0]
				elif len(other.span_ls)>1:
					return self.span_ls[0]<other.span_ls[0]
				else:
					return None
			elif len(self.span_ls)>1:
				if len(other.span_ls)==1:
					return self.span_ls[0]<other.span_ls[0]
				elif len(other.span_ls)>1:
					overlap1, overlap2 = get_overlap_ratio(self.span_ls, other.span_ls)
					if overlap1>0 and overlap1<=1 and  overlap2>0 and overlap2<=1:	# there's overlap between two segments 
						return False
					else:
						return self.span_ls[1]<other.span_ls[0]	# whether the stop of this segment is ahead of the start of other
				else:
					return None
		else:
			return self.chromosome<other.chromosome
	
	def __le__(self, other):
		"""
		2009-12-12
			whether self is less than or equal to other, in chromosomal order
		"""
		if self.chromosome==other.chromosome:
			if len(self.span_ls)==1:	# self is a point
				if len(other.span_ls)==1:
					return self.span_ls[0]<=other.span_ls[0]
				elif len(other.span_ls)>1:
					return self.span_ls[0]<=other.span_ls[0]
				else:
					return None
			elif len(self.span_ls)>1:
				if len(other.span_ls)==1:
					return self.span_ls[0]<=other.span_ls[0]
				elif len(other.span_ls)>1:
					is_overlap = is_reciprocal_overlap(self.span_ls, other.span_ls, min_reciprocal_overlap=self.min_reciprocal_overlap)
					if is_overlap:
						return True
					else:
						return self.span_ls[1]<=other.span_ls[0]	# whether the stop of this segment is ahead of the start of other
				else:
					return None
		else:
			return self.chromosome<other.chromosome
		
		
	def __eq__(self, other):
		"""
		2009-12-12
		"""
		if self.chromosome==other.chromosome:
			if len(self.span_ls)==1:
				if len(other.span_ls)==1:
					return self.span_ls[0]==other.span_ls[0]
				elif len(other.span_ls)>1:
					return self.span_ls[0]>=other.span_ls[0] and self.span_ls[0]<=other.span_ls[1]	# equal if self is within the "other" segment
				else:
					return None
			elif len(self.span_ls)>1:
				if len(other.span_ls)==1:	# self is a segment. other is a point position.
					return self.span_ls[0]>=other.span_ls[0] and self.span_ls[1]<=other.span_ls[0]	# if self includes the point position, yes it's equal
				elif len(other.span_ls)>1:
					# need to calculate min_reciprocal_overlap
					return is_reciprocal_overlap(self.span_ls, other.span_ls, min_reciprocal_overlap=self.min_reciprocal_overlap)
					#return self.span_ls[1]<other.span_ls[0]	# whether the stop of this segment is ahead of the start of other
				else:
					return None
		else:
			return False
				
	def __ne__(self, other):
		"""
		2009-12-12
			whether self is not equal to other
		"""
		return not self.__eq__(other)
		
	def __ge__(self, other):
		"""
		2009-12-12
		"""
		if self.chromosome==other.chromosome and len(self.span_ls)>1 and len(other.span_ls)>1:
			overlap1, overlap2 = get_overlap_ratio(self.span_ls, other.span_ls)
			if overlap1>=min_reciprocal_overlap and overlap2>=min_reciprocal_overlap:
				return True
			elif overlap1>0 and overlap1<=1 and  overlap2>0 and overlap2<=1:	# there's overlap between two segments 
				return False
			else:
				return not self.__lt__(other)
		else:
			return not self.__lt__(other)
		
	def __gt__(self, other):
		"""
		2009-12-12
		"""
		if self.chromosome==other.chromosome and len(self.span_ls)>1 and len(other.span_ls)>1:
			overlap1, overlap2 = get_overlap_ratio(self.span_ls, other.span_ls)
			if overlap1>0 and overlap1<=1 and  overlap2>0 and overlap2<=1:	# there's overlap between two segments 
				return False
			else:
				return not self.__le__(other)
		else:
			return not self.__le__(other)

	def __str__(self):
		"""
		2009-12-13
		"""
		return "chromosome: %s, span_ls: %s"%(self.chromosome, repr(self.span_ls))

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

# test program if this file is run
if __name__ == "__main__":
	import os, sys
	#import pdb
	#pdb.set_trace()
	cnv_ls = [[1, (2323,2600)], [2,(50000,)], [3,(43214,78788)], [5, (43242,)], [5,(144,566)], [5,(150,500)], [5,(500,950)], [5, (43241, 43242)]]
	no_of_cnvs = len(cnv_ls)
	min_reciprocal_overlap = 0.6
	for i in range(no_of_cnvs):
		for j in range(i+1, no_of_cnvs):
			cnv1 = cnv_ls[i]
			cnv2 = cnv_ls[j]
			segmentKey1 = CNVSegmentBinarySearchTreeKey(chromosome=cnv1[0], span_ls=cnv1[1], min_reciprocal_overlap=min_reciprocal_overlap)
			segmentKey2 = CNVSegmentBinarySearchTreeKey(chromosome=cnv2[0], span_ls=cnv2[1], min_reciprocal_overlap=min_reciprocal_overlap)
			print segmentKey1, "vs", segmentKey2
			print ">", segmentKey1>segmentKey2
			print ">=", segmentKey1>=segmentKey2
			print "<", segmentKey1<segmentKey2
			print "<=", segmentKey1<=segmentKey2
			print "==", segmentKey1==segmentKey2