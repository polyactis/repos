import os,sys
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
#import sqlalchemy
from variation.src import Stock_250kDB, AtDB, StockDB
from pymodule import GenomeDB
#from variation.src.Stock_250kDB import GeneListType
#from variation.src.Stock_250kDB import *


db = None
genome_db = None
at_db = None
gene_annotation = None
CandidateGeneTopSNPTestRMType_id_min_distance2ScoreRankHistogramType_id = {}
ecotype_name2accession_id = None

# 2009-11-17
CallMethodID2ecotype_id_set = {}
PhenotypeMethodID2ecotype_id_set = {}
PhenotypeAndCallMethodID2ecotype_id_set = {}