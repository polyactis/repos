import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
    sys.path.insert(0, os.path.expanduser('~/lib64/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
    sys.path.insert(0, os.path.expanduser('~/lib/python'))
    sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import psycopg2 as psycopg
import sys, getopt, csv, re
from codense.common import db_connect, dict_map
from common import nt2number, number2nt
from sets import Set
import networkx as nx
from pymodule import importNumericArray, ProcessOptions, SNPData, write_data_matrix, read_data, TwoSNPData
num = importNumericArray()