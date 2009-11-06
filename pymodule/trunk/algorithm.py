#!/usr/bin/env python
"""
2009-11-2
	module for various algorithms
"""

import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from ProcessOptions import  ProcessOptions
import copy

def listSubsets(element_ls, subset_size=None):
	"""
	2009-11-2
		list all possible subsets of element_ls. If subset_size is given, only produce subsets of that size.
			Otherwise, all of them.
		It's backtracking algorithm at play.
	"""
	sys.stderr.write("Listing all subsets of a list of %s elements ... "%(len(element_ls)))
	no_of_elements = len(element_ls)
	candidate_ls_ls = [[0,1]]	# 0 or 1 in k-th candidate_ls_ls signifies whether element_ls[k] would be included or not.
	#for i in range(no_of_elements):
	#	candidate_ls_ls.append([0, 1])
	
	one_solution=[]
	solution_ls=[]
	k = 0
	while k>=0:
		while len(candidate_ls_ls[k])>0:
			next_element = candidate_ls_ls[k].pop()
			if len(one_solution)==no_of_elements:
				one_solution[k] = next_element
			else:
				one_solution.append(next_element)
			k += 1
			if k==no_of_elements:
				if subset_size is None or sum(one_solution)==subset_size:
					one_subset = []
					for i in range(no_of_elements):
						if one_solution[i]==1:
							one_subset.append(element_ls[i])
					solution_ls.append(one_subset)		# python data by reference. have to use copy to sever the tie.
				break
			if len(candidate_ls_ls)<=k:	# still growing
				candidate_ls_ls.append([0,1])
			else:	# fully grown, now just replace the candidate list
				candidate_ls_ls[k] = [0, 1]
		k -= 1
	sys.stderr.write("Done.\n")
	return solution_ls


if __name__ == '__main__':
	#import pdb
	#pdb.set_trace()
	
	element_ls = [1,6,7]
	print listSubsets(element_ls, subset_size=1)
	
	print listSubsets(element_ls, subset_size=2)
	
	print listSubsets(element_ls)