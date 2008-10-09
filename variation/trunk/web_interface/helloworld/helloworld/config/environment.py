"""Pylons environment configuration"""
import os

from pylons import config

import helloworld.lib.app_globals as app_globals
import helloworld.lib.helpers
from helloworld.config.routing import make_map
from mako.lookup import TemplateLookup
import helloworld.model as model

def load_environment(global_conf, app_conf):
	"""Configure the Pylons environment via the ``pylons.config``
	object
	"""
	# Pylons paths
	root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	paths = dict(root=root,
				 controllers=os.path.join(root, 'controllers'),
				 static_files=os.path.join(root, 'public'),
				 templates=[os.path.join(root, 'templates')])

	# Initialize config with the basic options
	config.init_app(global_conf, app_conf, package='helloworld',
					template_engine='mako', paths=paths)

	config['routes.map'] = make_map()
	config['pylons.g'] = app_globals.Globals()
	config['pylons.h'] = helloworld.lib.helpers
	
	# Customize templating options via this variable
	tmpl_options = config['buffet.template_options']

	# CONFIGURATION OPTIONS HERE (note: all config options will override
	# any Pylons config options)
	
	#2008-10-05 Use the strict behaviour of the template context object
	config['pylons.strict_c'] = True
	
	#2008-10-05 Create the Mako TemplateLookup, with the default auto-escaping. it doesn't work though.
	config['pylons.g'].mako_lookup = TemplateLookup(
	directories=paths['templates'],
	module_directory=os.path.join(app_conf['cache_dir'], 'templates'),
	input_encoding='utf-8', output_encoding='utf-8',
	imports=['from webhelpers.html import escape'],
	default_filters=['escape'])
	
	#2008-10-05 setup the database connection
	model.drivername = config['drivername']
	model.hostname = config['hostname']
	model.dbname = config['dbname']
	model.schema = config['schema']
	model.db_user = config['db_user']
	model.db_passwd = config['db_passwd']
	model.pool_recycle = int(config['pool_recycle'])
	"""
	model.db = model.Stock_250kDB(drivername=model.drivername, username=model.db_user, password=model.db_passwd, \
							hostname=model.hostname, database=model.dbname, schema=model.schema)
	model.db.setup(create_tables=False)
	"""