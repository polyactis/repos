import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render
from helloworld import model
from pylons.decorators import jsonify
import simplejson
from helloworld.lib.base import h, config
import gviz_api, datetime, re
from pymodule import PassingData

log = logging.getLogger(__name__)

class AccessionController(BaseController):
	ecotype_central_view = 'stock.ecotype_info_with_haplogroup'	#could also be ecotype_info.
	
	def index(self):
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		# return 'Hello World'
		c.find250kAccessionsURL = h.url_for(controller="Accession", action="find250kAccessions", id=None)
		return render('/Accession.html')
	
	@jsonify
	def json(self):
		"""
		2009-4-3
			a toy function responding to the Pyjamas's jsonrpc example. useless now.
		"""
		#body = request.body[:int(request.environ['CONTENT_LENGTH'])]
		body = request.body
		
		#body = '{"params": ["anything interesting"], "id": 3, "method": "reverse"}'	#2009-3-31 synthesize one
		body = request.params.get('json')	#2009-3-31 request.body is nothing from asyncPost(). so now request comes from asyncGet() with the json structure in parameter 'json'.
		json = simplejson.loads(body)
		data = json['params'][0]	#2009-3-31 suddenly the json packed by pyjamas (0.5, previously it's newest) becomes a tuple (variable-name, variable-value), rather than just variable-value.
		#data = request.params.get('params', 'abc')
		method = json['method']
		#method = request.params.get('method', 'echo')
		if method == 'reverse':
			data = data[::-1]
		if method == 'uppercase':
			data = data.upper()
		if method == 'lowercase':
			data = data.lower()
		return dict(result=data)
	
	
	def processRegexpString(cls, p_str):
		"""
		2009-4-3
			if p_str includes '(' or ')', mysql complains: ERROR 1139 (42000): Got error 'parentheses not balanced' from regexp
		"""
		p_str = p_str.replace('(', '\(')	#python re module only needs one '\'
		p_str = p_str.replace(')', '\)')
		p_str_sql = p_str.replace('(', '\\\(')	#mysql needs more
		p_str_sql = p_str_sql.replace(')', '\\\)')
		return PassingData(p_str=p_str, p_str_sql=p_str_sql)
	processRegexpString = classmethod(processRegexpString)
	
	def findAccessionsByName(self):
		"""
		2009-4-14
		"""
		name = request.params.get('name', 'Algutsrum')	#default is 'Algutsrum'
		name = self.processRegexpString(name).p_str_sql
		condition = "name rlike '%s' or stockparent rlike '%s' or nativename rlike '%s' or alias rlike '%s'"%\
						(name, name, name, name)
		return self.findAccessions(condition)
	
	def findAccessionsByID(self):
		"""
		2009-4-14
		"""
		ecotype_id = request.params.get('id', '100')	#default is 100
		condition = "ecotypeid=%s"%(ecotype_id)
		return self.findAccessions(condition)
	
	@classmethod
	def findAccessions(cls, condition=None, extra_tables=None):
		"""
		2009-4-19
			add argument 'extra_tables'
		2009-4-14
			become classmethod
		2009-4-3
			find accession with name exact matching name, nativename, stockparent or alias
		"""
		if condition:
			condition = 'where %s'%condition
		else:
			condition = ""
		table_str = '%s v'%cls.ecotype_central_view
		if extra_tables:
			table_str += ', %s'%extra_tables
		rows = model.db.metadata.bind.execute("select v.* from %s %s"%\
													(table_str, condition))
		# dictionary to test whether one ecotype is in one dataset
		if getattr(model, 'ecotype_id_set_2010', None)==None:	#this lead to faster access next time. mysteriously, every call to self.ecotype_id_set_2010 (property) would lead to its execution everytime. no reference is stored.
			model.ecotype_id_set_2010 = cls.ecotype_id_set_2010()
			model.ecotype_id_set_384 = cls.ecotype_id_set_384()
			model.ecotype_id_set_perlegen = cls.ecotype_id_set_perlegen()
			model.ecotype_id_set_250k =  cls.ecotype_id_set_250k()
			model.ecotype_id_set_250k_in_pipeline = cls.ecotype_id_set_250k_in_pipeline()
		
		dataset_name2ecotype_id_set = {'2010': model.ecotype_id_set_2010,
							'384': model.ecotype_id_set_384,
							'perlegen': model.ecotype_id_set_perlegen,
							'250k': model.ecotype_id_set_250k}
		
		#3rd finally construct the full data and turn it into json
		column_name_type_ls = [("ecotypeid", ("number", "Ecotype ID")), ("tg_ecotypeid", ("number", "Unique Ecotype ID")), \
							("name", ("string","Name")), \
							("nativename", ("string","Native Name")), ("alias", ("string","Alias")), \
							("stockparent",("string", "Stock Parent")), ("haplo_group_id",("string", "149SNP Haplo-Group")), \
							("2010",("string", "2010")), \
							("384",("string", "384")), \
							("perlegen",("string", "perlegen")), \
							("250k",("string", "250k")), \
							("latitude",("number", "Latitude")), ("longitude",("number", "Longitude")), \
							("geographic_integrity",("string", "Geographic Info Quality")),\
							("site_name", ("string", "Collection Site")), ("country", ("string", "Country")),\
							("collector", ("string", "Collector")), ("collectiondate", ("date", "Collection Date"))]
		
		description = dict(column_name_type_ls)
		return_ls = []
		for row in rows:
			entry = dict()
			for column_name_type in column_name_type_ls:
				column_name = column_name_type[0]
				column_type = column_name_type[1][0]
				
				if column_type=='string':
					default_value = ''
				elif column_type =='number':
					default_value = -1
				elif column_type=='date':
					default_value = datetime.date(2050, 1, 1)
				else:
					default_value = None
				
				if column_name == 'collector':
					column_value = '%s %s'%(row.firstname, row.surname)
				elif column_name == '2010' or column_name=='384' or column_name=='perlegen' or column_name=='250k':
					if row.tg_ecotypeid in dataset_name2ecotype_id_set[column_name]:
						column_value = 'yes'
					elif column_name=='250k' and row.tg_ecotypeid in model.ecotype_id_set_250k_in_pipeline:
						column_value = 'in pipeline'
					else:
						column_value = 'no'
				elif column_name=='haplo_group_id':
					column_value = getattr(row, column_name, default_value)
					#haplo_group_name = model.Stock.HaploGroup.get(column_value).short_name
					if column_value:
						column_value = "<a href=%s target='_blank'>%s</a>"%(h.url_for(controller="Accession", action='haploGroup', id=column_value), column_value)
				else:
					column_value = getattr(row, column_name, default_value)
				entry[column_name] = column_value
			return_ls.append(entry)
		
		data_table = gviz_api.DataTable(description)
		data_table.LoadData(return_ls)
		column_ls = [row[0] for row in column_name_type_ls]
		json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		response.headers['Content-Type'] = 'application/json'
		return json_result
	
	def findAccessionsNameLike(self, namelike):
		"""
		2009-4-3
			find all accessions with name beginned with 'namelike'
		"""
		name_processed_ob = self.processRegexpString(namelike)
		namelike = name_processed_ob.p_str_sql
		rows = model.db.metadata.bind.execute("select * from %s where name rlike '^%s' or stockparent rlike '^%s' or nativename rlike '^%s' or alias rlike '^%s'"%\
											(self.ecotype_central_view, namelike, namelike, namelike, namelike))
		name_set = set()
		name_p = re.compile(r'^%s'%name_processed_ob.p_str, re.IGNORECASE)
		columns_to_be_added = ['name', 'nativename', 'stockparent', 'alias']
		for row in rows:
			for col in columns_to_be_added:
				col_value = getattr(row, col)
				if col_value and name_p.search(col_value):
					name_set.add(col_value)
		return list(name_set)
	
	#@jsonify	#2009-4-3, comment it out and using simplejson.dumps() directly because the default encoding 'utf-8' doesn't work due to some swedish letters.
	def autoComplete(self):
		"""
		2009-4-3
			autoComplete server end for AccessionByName.java
		"""
		namelike = request.params.get('namelike')
		name_ls = self.findAccessionsNameLike(namelike)
		name_ls.sort()
		if len(name_ls)>100:
			name_ls = name_ls[:100]
		response.headers['Content-Type'] = 'application/json'
		return simplejson.dumps(dict(result=name_ls), encoding='latin1')
	
	def haploGroup(self, id=None):
		"""
		2009-4-5
			server end to return haplotype group json data structure
		"""
		if id is None:
			id = 1
		c.haplo_group_id = id
		c.getHaploGroupURL = h.url_for(controller="Accession", action="getHaploGroup", id=id)
		return render('/HaploGroup.html')
	
	def getHaploGroup(self, id=None):
		"""
		2009-4-19
		"""
		if id is None:
			id = 1
		c.haplo_group_id = id
		extra_tables = 'stock.haplo_group2ecotype h'
		#rows = model.db.metadata.bind.execute("select ecotype_id from stock.haplo_group2ecotype where haplo_group_id=%s"%c.haplo_group_id)
		#ecotype_id_ls = [str(row.ecotype_id) for row in rows]
		condition = "v.ecotypeid=h.ecotype_id and h.haplo_group_id=%s"%(c.haplo_group_id)
		return self.findAccessions(condition, extra_tables=extra_tables)
	
	def find250kAccessions(self):
		"""
		2009-4-16
		"""
		if getattr(model, 'ecotype_id_set_250k', None)==None:	#this lead to faster access next time. mysteriously, every call to self.ecotype_id_set_2010 (property) would lead to its execution everytime. no reference is stored.
			model.ecotype_id_set_250k =  self.ecotype_id_set_250k()
		ecotype_id_ls = [str(ecotype_id) for ecotype_id in model.ecotype_id_set_250k]
		condition = "ecotypeid in (%s)"%(','.join(ecotype_id_ls))
		response.headers['Content-Type'] = 'application/json'
		return self.findAccessions(condition)
	
	"""
	2009-4-10
		sets of ecotype ids for each individual dataset
	"""
	@classmethod
	def ecotype_id_set_perlegen(cls, accession2ecotype_table='accession2tg_ecotypeid'):
		ecotype_id_set_perlegen = set()
		for ecotype_name, accession_id in model.ecotype_name2accession_id.iteritems():
			row = model.at_db.metadata.bind.execute("select * from %s where accession_id=%s"%(accession2ecotype_table, accession_id)).fetchone()
			ecotype_id_set_perlegen.add(row.ecotype_id)
		return ecotype_id_set_perlegen
	
	@classmethod
	def ecotype_id_set_384(cls):
		ecotype_id_set_384 = set()
		#rows = model.dbsnp.Accession.query()
		rows = model.db.metadata.bind.execute("select * from dbsnp.accession")
		for row in rows:
			ecotype_id_set_384.add(row.ecotype_id)
		return ecotype_id_set_384
	
	@classmethod
	def ecotype_id_set_2010(cls, accession2ecotype_table='accession2tg_ecotypeid'):
		ecotype_id_set_2010 = set()
		for row in model.at_db.metadata.bind.execute("select * from %s"%accession2ecotype_table):
			ecotype_id_set_2010.add(row.ecotype_id)
		return ecotype_id_set_2010
	
	@classmethod
	def ecotype_id_set_250k(cls):
		ecotype_id_set_250k = set()		#from good_call_method_id
		for row in model.Stock_250kDB.CallInfo.query.filter_by(method_id=int(config['app_conf']['good_call_method_id'])):
			ecotype_id_set_250k.add(row.array.maternal_ecotype_id)
		return ecotype_id_set_250k
	
	@classmethod
	def ecotype_id_set_250k_in_pipeline(cls):
		ecotype_id_set_250k_in_pipeline = set()
		for row in model.Stock_250kDB.ArrayInfo.query():
			if row.maternal_ecotype_id==row.paternal_ecotype_id:	#no crosses.
				ecotype_id_set_250k_in_pipeline.add(row.maternal_ecotype_id)
		return ecotype_id_set_250k_in_pipeline