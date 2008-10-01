"""
2008-10-01
	module related to Genome
"""
class GeneModel(object):
	def __init__(self, **keywords):
		"""
		2008-10-01
			moved from transfac.src.GenomeDB in order for the return of GenomeDB.get_gene_id2model() to be pickled independent of transfac.src
		2008-10-01
			a class to hold all stuff related to a gene (Gene+EntrezgeneMapping)
			it's hierarchical. Its gene_commentaries contains also GeneModel.
		"""
		for argument_key, argument_value in keywords.iteritems():
			setattr(self, argument_key, argument_value)
		if not hasattr(self, 'gene_commentaries'):
			self.gene_commentaries = []
