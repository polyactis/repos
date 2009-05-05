import logging

from helloworld.lib.base import *
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule import PassingData

log = logging.getLogger(__name__)
import helloworld.model as model
from HelpOtherControllers import HelpothercontrollersController as hc

class ScorerankhistogramController(BaseController):

	def index(self):
		"""
		2008-10-30
			fetch ScoreRankHistogramType in particular order
		"""
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		rows = ScoreRankHistogramType.query.order_by(ScoreRankHistogramType.results_type).\
			order_by(ScoreRankHistogramType.null_distribution_type_id).\
			order_by(ScoreRankHistogramType.min_distance).order_by(ScoreRankHistogramType.get_closest).\
			order_by(ScoreRankHistogramType.min_MAF).order_by(ScoreRankHistogramType.allow_two_sample_overlapping).all()
		c.score_rank_hist_type_ls = rows
		return render('/score_rank_histgram.html')
	
	def type(self, id=None):
		ScoreRankHistogram = model.Stock_250kDB.ScoreRankHistogram
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		query = ScoreRankHistogram.query
		if id is None:
			id = request.params.get('id', 1)
		hist_type_id = id
		
		extra_condition = 's.hist_type_id=%s'%hist_type_id
		
		c.phenotype_info  = hc.getPhenotypeInfo(ScoreRankHistogram.table.name, extra_condition)
		c.list_info = hc.getListTypeInfo(ScoreRankHistogram.table.name, extra_condition)
		
		c.hist_type = ScoreRankHistogramType.get(hist_type_id)
		rows = query.filter_by(hist_type_id=hist_type_id).all()
		data_matrix = numpy.zeros([len(c.list_info.list_type_id_ls), len(c.phenotype_info.phenotype_method_id_ls)], numpy.int)
		data_matrix[:] = -1
		counter = 0
		for row in rows:
			row_index = c.list_info.list_type_id2index[row.list_type_id]
			col_index = c.phenotype_info.phenotype_method_id2index[row.phenotype_method_id]
			data_matrix[row_index, col_index] = row.id
			counter +=1
		c.counter = counter
		c.data_matrix = data_matrix
		return render('/score_rank_histgram_type.html')
	
	def showHistogram(self, id=None):
		ScoreRankHistogram = model.Stock_250kDB.ScoreRankHistogram
		c.score_rank_histgram = ScoreRankHistogram.get(id)
		return render('/score_rank_histgram_single.html')
	
	def showHistogramByQuery(self, id=None, type_id=None, list_type_id=None):
		"""
		2008-11-05
			similar to showHistogram, but identifying ScoreRankHistogram by results_id, type_id, list_type_id
			
			to be called by DisplayTopSNPTestRM.py's template
		"""
		type_id = request.params.get('type_id', type_id)
		list_type_id = request.params.get('list_type_id', list_type_id)
		results_id = request.params.get('id', id)
		if results_id is None:
			results_id = request.params.get('results_id', None)
		if results_id is None:
			return 'Nothing'
		
		ScoreRankHistogram = model.Stock_250kDB.ScoreRankHistogram
		rm = model.Stock_250kDB.ResultsMethod.get(results_id)
		if rm:
			row = ScoreRankHistogram.query.filter_by(phenotype_method_id=rm.phenotype_method_id).\
					filter_by(list_type_id=list_type_id).filter_by(hist_type_id=type_id).first()
			c.score_rank_histgram = row
			return render('/score_rank_histgram_single.html')
		else:
			return 'Nothing'
		
	
	def getImage(self, id=None, img_type='score_hist'):
		img_type = request.params.get('img_type', img_type)
		if img_type[-3:]=='svg':
			response.headers['Content-type'] = 'image/svg'
		else:
			response.headers['Content-type'] = 'image/png'
		ScoreRankHistogram = model.Stock_250kDB.ScoreRankHistogram
		if id is None:
			id = request.params.get('id', None)
		get_img_data_success = 0
		if id:
			score_rank_histgram = ScoreRankHistogram.get(id)
			if score_rank_histgram:
				img_data = getattr(score_rank_histgram, img_type, None)
				if img_data:
					get_img_data_success = 1
					img_data = img_data.__str__()	#2008-12-26	img_data is a buffer. weird!
				else:
					get_img_data_success = 0
		if not get_img_data_success:
			response.headers['Content-type'] = 'text/plain'
			img_data = 'No Image'
		return img_data