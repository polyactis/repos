import logging

from pylons import request, response, session, tmpl_context as c
from pylons.controllers.util import abort, redirect_to

from helloworld.lib.base import BaseController, render, config, h, model
#from helloworld import model
from pymodule import figureOutDelimiter
import datetime, StringIO, csv, gviz_api, traceback, sys

log = logging.getLogger(__name__)

class UtilsController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/template.mako')
		# or, Return a response
		c.motionChartFormActionURL = h.url_for(controller="utils", action="uploadForMotionChart");
		return render('/Utils.html')
	
	
	vizDataType = set(['string', 'number', 'date'])
	
	def uploadForMotionChart(self, id=None):
		"""
		2009-8-31
			server end to handle the file uploaded from MotionChartAppMCPanel.
			read the data and return it in the google visualization data table json structure if no error
		"""
		myfile = request.POST.values()[-1]
		#myfile = StringIO.StringIO(myfile)
		
		try:
			reader = csv.reader(myfile.file, delimiter='\t')
			#create column_name_type_ls based on the header
			column_header_ls = reader.next()
			column_name_type_ls = []
			for i in range(len(column_header_ls)):
				column_header = column_header_ls[i]
				tmp_ls = column_header.split('|')
				column_header = tmp_ls[0]
				if i==0:
					column_type = 'string'	#2009-8-31 has to be string. ID for each row.
				elif len(tmp_ls)>1:
					column_type = tmp_ls[1]
					if column_type not in self.vizDataType:
						column_type = 'number'
				else:
					column_type = 'number'
				column_name_type_ls.append((column_header, (column_type, column_header)))
			
			#read each row and put it into a dictionary
			return_ls = []
			for row in reader:
				if len(row)<2:
					continue
				entry = {}
				for i in range(len(row)):
					column_name, type_desc = column_name_type_ls[i]
					column_type, column_desc = type_desc
					if row[i] == 'NA':
						column_value = None
					elif row[i] == '':
						column_value = None
					elif column_type == 'number':
						column_value = float(row[i])
					elif column_type == 'date':
						column_value = datetime.datetime.strptime(row[i], "%y-%m-%d")
					else:
						column_value = row[i]
					entry[column_name] = column_value
				entry["date"] = datetime.date(2009,2,3)
				return_ls.append(entry)
			
			#insert the date column right after the 1st column
			column_name_type_ls.insert(1, ("date", ("date", "Date")))
			description = dict(column_name_type_ls)
			data_table = gviz_api.DataTable(description)
			data_table.LoadData(return_ls)
			column_ls = [row[0] for row in column_name_type_ls]
			json_result = data_table.ToJSon(columns_order=column_ls)	#ToJSonResponse only works with google.visualization.Query
		except:
			loginfo = "Error encountered while converting data into json structure.\n"
			loginfo += repr(sys.exc_info())
			loginfo += repr(traceback.print_exc())
			log.error(loginfo)
			json_result = loginfo
		#response.headers['Content-Type'] = 'application/json'	#non-text Content-type renders the browser to prompt the user regarding how to deal with data
		#"text/plain"
		return json_result
		"""
		permanent_file = open(os.path.join(permanent_store, myfile.filename.lstrip(os.sep)), 'w')		
		shutil.copyfileobj(myfile.file, permanent_file)
		myfile.file.close()
		permanent_file.close()
		
		return 'Successfully uploaded: %s, description: %s' % \
			(myfile.filename, request.POST['description'])
		"""