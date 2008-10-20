import logging

from helloworld.lib.base import *
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule import PassingData

log = logging.getLogger(__name__)
import helloworld.model as model


class ScorerankhistogramController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		rows = ScoreRankHistogramType.query.order_by(ScoreRankHistogramType.id).all()
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
		
		c.phenotype_info  = h.getPhenotypeInfo(ScoreRankHistogram.table.name, extra_condition)
		c.list_info = h.getListTypeInfo(ScoreRankHistogram.table.name, extra_condition)
		
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
				else:
					get_img_data_success = 0
		if not get_img_data_success:
			response.headers['Content-type'] = 'text/plain'
			img_data = 'No Image'
		return img_data