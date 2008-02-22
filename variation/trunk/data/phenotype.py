import Gnumeric
wb = Gnumeric.workbooks()[1]
s = wb.sheets()[0]
f = lambda x: s[x,75].get_value()
kas_2_data = [f(i) for i in range(15,19)]
f2 = lambda j : [s[i,j].get_value() for i in range(25,29)]
matrix = [f2(j) for j in range(1,187)]
import numpy
from rpy import r
dist_ls = [numpy.average((numpy.array(matrix[i])-kas_2_data)*(numpy.array(matrix[i])-kas_2_data)) for i in range(len(matrix))]
from sets import Set
good_index_set = Set(range(186))
bad_index_set = Set([1, 8, 9, 17, 102, 104, 108, 119, 139, 140, 154, 157, 177, 180])

good_index_set = good_index_set - bad_index_set
dist_ls = [numpy.average((numpy.array(matrix[i])-kas_2_data)*(numpy.array(matrix[i])-kas_2_data)) for i in list(good_index_set)]
dist_ls_argsort = numpy.argsort(dist_ls)
dist_ls_argsort[0]

numpy.average((numpy.array(matrix[51])-kas_2_data)*(numpy.array(matrix[41])-kas_2_data))
