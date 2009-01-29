"""Routes configuration

The more specific and detailed routes should be defined first so they
may take precedent over the more generic routes. For more information
refer to the routes manual at http://routes.groovie.org/docs/
"""
from pylons import config
from routes import Mapper

def make_map():
	"""Create, configure and return the routes Mapper"""
	map = Mapper(directory=config['pylons.paths']['controllers'],
				always_scan=config['debug'])
	map.minimization = False	#2008-12-26 if with minimization enabled, Routes assumes that you want to use the optional default variables automatically if they aren't specified in the URL.
	#like http://localhost/abc = http://localhost/abc/ although the latter is not specified here.
	
	# The ErrorController route (handles 404/500 error pages); it should
	# likely stay at the top, ensuring it can always be resolved
	map.connect('/error/{action}', controller='error')
	map.connect('/error/{action}/{id}', controller='error')
	#map.connect('error/:action/:id', controller='error')

	# CUSTOM ROUTES HERE
	map.connect('/', controller='hello', action='index')
	#map.connect('/', controller='greeting', action='index')	#2008-12-26 the root of the site. '/' is now required after map.minimization=False
	#map.connect('/', '/templates/index.html')	#2008-12-27 doesn't work
	
	#map.connect(':controller/:action/:id')
	map.connect('/{controller}')
	map.connect('/{controller}/')
	map.connect('/{controller}/{action}')
	map.connect('/{controller}/{action}/')
	map.connect('/{controller}/{action}/{id}')
	#map.connect('*url', controller='template', action='view')

	return map
