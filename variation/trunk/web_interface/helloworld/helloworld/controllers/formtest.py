import logging

from helloworld.lib.base import *

log = logging.getLogger(__name__)

import formencode
from formencode import htmlfill

class EmailForm(formencode.Schema):
	allow_extra_fields = True
	filter_extra_fields = True
	email = formencode.validators.Email(not_empty=True)
	date = formencode.validators.DateConverter(not_empty=True)


class FormtestController(BaseController):

	def index(self):
		# Return a rendered template
		#   return render('/some/template.mako')
		# or, Return a response
		return render('/simpleform.html')
	
	def form(self):
		return render('/simpleform.html')
	
	@validate(schema=EmailForm(), form='form', post_only=False, on_get=True)
	def submit(self):
		"""
		schema = EmailForm()
		try:
			form_result = schema.to_python(dict(request.params))
		except formencode.Invalid, error:
			c.form_result = error.value
			c.form_errors = error.error_dict or {}
			html = render('/simpleform.html')
			return htmlfill.render(
					html,
					defaults=c.form_result,
					errors=c.form_errors
					)
			
			#response.content_type = 'text/plain'
			#return 'Invalid: '+str(error)
		else:
		"""
		html = render('/simpleform.html')
		return htmlfill.render(
					html,
					defaults=self.form_result,
					#errors=self.form_errors
					)
		#return 'Your email is: %s, date selected was %r.'%(self.form_result.get('email'), self.form_result.get('date'))
		
		"""
		c.email_msg = ''
		email = request.params.get('email')
		if not email:
			c.email_msg = "Please enter a value"
		elif '@' not in email:
			c.email_msg = "An email address must contain at least one '@' character."
		else:
			domain = email.split('@')[1]
			if '.' not in domain:
				c.email_msg = "An email address domain must contain at least one '.' character."
			if not domain.split('.')[-1]:
				c.email_msg = "Please specify a domain type after the '.' character"
		if c.email_msg:
			c.email_value = email
			return render('/simpleform.html')
		return 'Your email is: %s' % request.params['email']
		"""
		#h.redirect_to(action='result')

	def result(self):
		#return 'Your email is: %s' % request.params['email']
		return 'Your data was successfully submitted.'