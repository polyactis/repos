"""
2007-02-19
	common stuff for variation/src
"""


"""
2007-03-29
	add mappings for '-', 'N', and ambiguous letters ('M', 'R', 'W', 'S', 'Y', 'K')
"""
nt2number = {'-': -1,	#deletion
	'N': 0,
	'NA': 0,
	'A': 1,
	'C': 2,
	'G': 3,
	'T': 4,
	'AC':5,
	'CA':5,
	'M':5,
	'AG':6,
	'GA':6,
	'R':6,
	'AT':7,
	'TA':7,
	'W':7,
	'CG':8,
	'GC':8,
	'S':8,
	'CT':9,
	'TC':9,
	'Y':9,
	'GT':10,
	'TG':10,
	'K':10
	}

number2nt = {0: 'NA',
	1:'A',
	2:'C',
	3:'G',
	4:'T',
	5:'AC',
	6:'AG',
	7:'AT',
	8:'CG',
	9:'CT',
	10:'GT'
	}

#2007-04-16 entry[i,j] means whether nucleotide i and j matches. 0(NA) matches everything. singleton(1-4) matches itself and the doublet containing it. doublet(5-10) matches only itself.
nt_number_matching_matrix = [[1, 1,1,1,1,1, 1,1,1,1,1],
	[1, 1,0,0,0,1, 1,1,0,0,0],
	[1, 0,1,0,0,1, 0,0,1,1,0],
	[1, 0,0,1,0,0, 1,0,1,0,1],
	[1, 0,0,0,1,0, 0,1,0,1,1],
	[1, 1,1,0,0,1, 0,0,0,0,0],
	[1, 1,0,1,0,0, 1,0,0,0,0],
	[1, 1,0,0,1,0, 0,1,0,0,0],
	[1, 0,1,1,0,0, 0,0,1,0,0],
	[1, 0,1,0,1,0, 0,0,0,1,0],
	[1, 0,0,1,1,0, 0,0,0,0,1]]
