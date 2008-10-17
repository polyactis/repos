import logging

from helloworld.lib.base import *
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
from pymodule import PassingData

log = logging.getLogger(__name__)
import helloworld.model as model
from sets import Set
import numpy

class ScorerankhistogramController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		rows = ScoreRankHistogramType.query.order_by(ScoreRankHistogramType.id).all()
		c.score_rank_hist_type_ls = rows
		return render('/score_rank_histgram.html')
	
	def getPhenotypeInfo(self, hist_type_id=1):
		"""
		"""
		ScoreRankHistogram = model.Stock_250kDB.ScoreRankHistogram
		PhenotypeMethod = model.Stock_250kDB.PhenotypeMethod
		rows = model.db.metadata.bind.execute("select distinct s.phenotype_method_id, p.biology_category_id, p.short_name from %s s, %s p where \
			p.id=s.phenotype_method_id and s.hist_type_id=%s order by p.biology_category_id, s.phenotype_method_id"\
			%(ScoreRankHistogram.table.name, PhenotypeMethod.table.name, hist_type_id))
		phenotype_method_id_ls = []
		phenotype_method_id2index = {}
		phenotype_method_label_ls = []
		prev_biology_category_id = None
		no_of_separators = 0
		for row in rows:
			if prev_biology_category_id == None:
				prev_biology_category_id = row.biology_category_id
			elif row.biology_category_id!=prev_biology_category_id:
				prev_biology_category_id = row.biology_category_id
				#add a blank phenotype id as separator
				no_of_separators += 1
				phenotype_method_id2index[-no_of_separators] = len(phenotype_method_id_ls)
				phenotype_method_id_ls.append(-no_of_separators)
				phenotype_method_label_ls.append('=====')
			phenotype_method_id2index[row.phenotype_method_id] = len(phenotype_method_id_ls)
			phenotype_method_id_ls.append(row.phenotype_method_id)
			phenotype_method_label_ls.append('%s %s'%(row.phenotype_method_id, row.short_name))
		phenotype_info = PassingData()
		phenotype_info.phenotype_method_id2index = phenotype_method_id2index
		phenotype_info.phenotype_method_id_ls = phenotype_method_id_ls
		phenotype_info.phenotype_method_label_ls = phenotype_method_label_ls
		return phenotype_info
	
	def getListTypeInfo(self, hist_type_id=1):
		"""
		2008-08-29
			add -1 as a separator into list_type_id_ls
		"""
		ScoreRankHistogram = model.Stock_250kDB.ScoreRankHistogram
		GeneListType = model.Stock_250kDB.GeneListType
		rows = model.db.metadata.bind.execute("select distinct s.list_type_id, p.biology_category_id, p.short_name from %s s, %s p where \
			p.id=s.list_type_id and s.hist_type_id=%s order by p.biology_category_id, s.list_type_id"\
			%(ScoreRankHistogram.table.name, GeneListType.table.name, hist_type_id))
		list_type_id_ls = []
		list_type_id2index = {}
		list_type_label_ls = []
		prev_biology_category_id = None
		no_of_separators = 0
		for row in rows:
			if prev_biology_category_id == None:
				prev_biology_category_id = row.biology_category_id
			elif row.biology_category_id!=prev_biology_category_id:
				prev_biology_category_id = row.biology_category_id
				no_of_separators += 1
				list_type_id2index[-no_of_separators] = len(list_type_id_ls)
				list_type_id_ls.append(-no_of_separators)
				list_type_label_ls.append('====\n====')
			list_type_id2index[row.list_type_id] = len(list_type_id_ls)
			list_type_id_ls.append(row.list_type_id)
			list_type_label_ls.append('%s %s'%(row.list_type_id, row.short_name))
		list_info = PassingData()
		list_info.list_type_id2index = list_type_id2index
		list_info.list_type_id_ls = list_type_id_ls
		list_info.list_type_label_ls = list_type_label_ls
		return list_info
	
	def type(self, id=None):
		ScoreRankHistogram = model.Stock_250kDB.ScoreRankHistogram
		ScoreRankHistogramType = model.Stock_250kDB.ScoreRankHistogramType
		query = ScoreRankHistogram.query
		if id is None:
			id = request.params.get('id', 1)
		hist_type_id = id
		c.phenotype_info  = self.getPhenotypeInfo(hist_type_id)
		c.list_info = self.getListTypeInfo(hist_type_id)
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