"""Pylons middleware initialization"""
from paste.cascade import Cascade
from paste.registry import RegistryManager
from paste.urlparser import StaticURLParser
from paste.deploy.converters import asbool

from pylons import config
from pylons.error import error_template
from pylons.middleware import error_mapper, ErrorDocuments, ErrorHandler, \
	StaticJavascripts, StatusCodeRedirect
from pylons.wsgiapp import PylonsApp
#2008-12-25. 0.9.6 => 0.9.7 Add these imports to the top
from beaker.middleware import CacheMiddleware, SessionMiddleware
from routes.middleware import RoutesMiddleware

from helloworld.config.environment import load_environment

def make_app(global_conf, full_stack=True, **app_conf):
	"""Create a Pylons WSGI application and return it

	``global_conf``
		The inherited configuration for this application. Normally from
		the [DEFAULT] section of the Paste ini file.

	``full_stack``
		Whether or not this application provides a full WSGI stack (by
		default, meaning it handles its own exceptions and errors).
		Disable full_stack when this application is "managed" by
		another WSGI middleware.

	``app_conf``
		The application's local configuration. Normally specified in the
		[app:<name>] section of the Paste ini file (where <name>
		defaults to main).
	"""
	# Configure the Pylons environment
	load_environment(global_conf, app_conf)

	# The Pylons WSGI app
	app = PylonsApp()

	# CUSTOM MIDDLEWARE HERE (filtered by error handling middlewares)
	
	#2008-12-25 0.9.6 => 0.9.7
	# Add these below the 'CUSTOM MIDDLEWARE HERE' line, or if you removed
	# that, add them immediately after the PylonsApp initialization
	app = RoutesMiddleware(app, config['routes.map'])
	app = SessionMiddleware(app, config)
	app = CacheMiddleware(app, config)
	
	if asbool(full_stack):
		# Handle Python exceptions
		app = ErrorHandler(app, global_conf,
						   **config['pylons.errorware'])	#2008-12-25 remove ", error_template=error_template" from ErrorHandler()

		# Display error documents for 401, 403, 404 status codes (and
		# 500 when debug is disabled)
		#if asbool(config['debug']):
		#	app = StatusCodeRedirect(app)
		#else:
		#	app = StatusCodeRedirect(app, [400, 401, 403, 404, 500])
		#app = ErrorDocuments(app, global_conf, mapper=error_mapper, **app_conf)

	# Establish the Registry for this application
	app = RegistryManager(app)

	# Static files
	#javascripts_app = StaticJavascripts()
	static_app = StaticURLParser(config['pylons.paths']['static_files'])
	app = Cascade([static_app, app])	#"javascripts_app, " is removed out of the list
	return app
