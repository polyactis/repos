"""
2007-10-15. __init__.py for pymodule. pymodule is a concatenation of all common functions/classes.
"""
import os, sys
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from ProcessOptions import ProcessOptions, generate_program_doc, process_options, process_function_arguments, turn_option_default_dict2argument_default_dict
from SNP import write_data_matrix, read_data, SNPData, GenomeWideResults, GenomeWideResult, DataObject, getGenomeWideResultFromFile,\
	nt2number, number2nt, SNPInfo, number2single_char_nt
from TwoSNPData import TwoSNPData, QualityControl
from utils import PassingData, dict_map, importNumericArray, figureOutDelimiter, get_gene_symbol2gene_id_set, \
	FigureOutTaxID, getColName2IndexFromHeader, getListOutOfStr, runLocalCommand, getGeneIDSetGivenAccVer, calGreatCircleDistance

from Genome import GeneModel

from yh_matplotlib import assignMatPlotlibHueColorToLs, drawName2FCLegend
from PCA import PCA