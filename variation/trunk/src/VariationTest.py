#!/usr/bin/env python
""""
Usage: VariationTest.py -y TestCaseType [OPTIONS]

Option:
	-y ..., --type=...	which test case should be invoked.
	-h, --help              show this help

Examples:
	VariationTest.py -y 2

2007-03-08 1: TestTrioInference
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import unittest, os, sys, getopt, csv

class TestTrioInference(unittest.TestCase):
	"""
	2007-03-08
	"""
	def test_trio_ancestry_inference_by_DP(self):
		"""
		2007-04-16
			add a heterozygous call to test
		"""
		from MpiTrioAncestryInference import MpiTrioAncestryInference
		import Numeric
		data_matrix = Numeric.array([   [0,1,3,1,2,2,2,1,1,1],
						[1,2,4,3,2,2,1,4,2,4],
						[2,1,4,3,2,2,1,1,5,4]])
		chr_start_ls = [0,4,10]
		trio_arrangement_ls = [[0,1,2], [1,2,0], [2,0,1]]
		MpiTrioAncestryInference_instance = MpiTrioAncestryInference(debug=1)
		for trio_arrangement in trio_arrangement_ls:
			ancestry_ls, no_of_jumps = MpiTrioAncestryInference_instance.identify_ancestry_with_min_jumps(data_matrix[trio_arrangement[0]], data_matrix[trio_arrangement[1]], data_matrix[trio_arrangement[2]], chr_start_ls)
			print trio_arrangement
			print ancestry_ls
			print no_of_jumps

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "type="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hy:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	TestCaseDict = {1:TestTrioInference}
	type = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-y", "--type"):
			type = int(arg)
			
	if type:
		suite = unittest.TestSuite()
		suite.addTest(unittest.makeSuite(TestCaseDict[type]))
		unittest.TextTestRunner(verbosity=2).run(suite)

	else:
		print __doc__
		sys.exit(2)		