import logging

from helloworld.lib.base import *
#from pylons import app_globals

log = logging.getLogger(__name__)


class HelloController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		#return 'Hello World'
		#response.headers['Content-Type'] = 'text/plain'
		#return 'Hello from the index() action!'
		return render('/index.html')


	def serverinfo(self):
		import cgi
		import pprint
		c.pretty_environ = cgi.escape(pprint.pformat(request.environ))
		c.name = 'The Black Knight'
		session['name'] = 'mighty'
		session.save()
		return render('/serverinfo.mako')

	def app_globals_test(self):
		if g.message == 'Hello':
			content = g.message
			g.message = 'Hello World!'
			#return content
		else:
			pass
			#return g.message
		g.visits += 1
		return "You are visitor number %s." % g.visits

	def environ(self):
		result = '<html><body><h1>Environ</h1>'
		for key, value in request.environ.items():
			result += '%s: %r <br />'%(key, value)
		result += '</body></html>'
		return result

	def navigation(self):
		return render('/navigation.html')
