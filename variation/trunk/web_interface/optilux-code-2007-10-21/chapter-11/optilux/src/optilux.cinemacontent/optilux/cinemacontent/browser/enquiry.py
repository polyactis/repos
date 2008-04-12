"""An enquiry form using zope.formlib
"""

import re

from zope.interface import Interface

from zope import schema

from zope.formlib import form
from Products.Five.formlib import formbase
from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile

from Products.statusmessages.interfaces import IStatusMessage

from Acquisition import aq_inner
from Products.CMFCore.utils import getToolByName

from optilux.cinemacontent import CinemaMessageFactory as _

# Define a valiation method for email addresses
class NotAnEmailAddress(schema.ValidationError):
    __doc__ = _(u"Invalid email address")

check_email = re.compile(r"[a-zA-Z0-9._%-]+@([a-zA-Z0-9-]+\.)*[a-zA-Z]{2,4}").match
def validate_email(value):
    if not check_email(value):
        raise NotAnEmailAddress(value)
    return True
    
MESSAGE_TEMPLATE = """\
Enquiry from: %(name)s <%(email_address)s>

%(message)s
"""

class IEnquiryForm(Interface):
    """Define the fields of our form
    """
    
    subject = schema.TextLine(title=_(u"Subject"),
                              required=True)
                              
    name = schema.TextLine(title=_(u"Your name"),
                              required=True)
    
    email_address = schema.ASCIILine(title=_(u"Your email address"),
                                    description=_(u"We will use this to contact you if you request it"),
                                    required=True,
                                    constraint=validate_email)
    
    message = schema.Text(title=_(u"Message"),
                          description=_(u"Please keep to 1,000 characters"),
                          required=True,
                          max_length=1000)

class EnquiryForm(formbase.PageForm):
    form_fields = form.FormFields(IEnquiryForm)
    label = _(u"Make an enquiry")
    description = _(u"Got a question or comment? Please submit it using the form below!")
    
    # This trick hides the editable border and tabs in Plone
    def __call__(self):
        self.request.set('disable_border', True)
        return super(EnquiryForm, self).__call__()
    
    @form.action(_(u"Send"))
    def action_send(self, action, data):
        """Send the email to the site administrator and redirect to the
        front page, showing a status message to say the message was received.
        """
        
        context = aq_inner(self.context)
        
        mailhost = getToolByName(context, 'MailHost')
        urltool = getToolByName(context, 'portal_url')
        
        portal = urltool.getPortalObject()
        email_charset = portal.getProperty('email_charset')

        # Construct and send a message
        to_address = portal.getProperty('email_from_address')
        source = "%s <%s>" % (data['name'], data['email_address'])
        subject = data['subject']
        message = MESSAGE_TEMPLATE % data

        mailhost.secureSend(message, to_address, str(source),
                            subject=subject, subtype='plain',
                            charset=email_charset, debug=False,
                            From=source)
        
        # Issue a status message
        confirm = _(u"Thank you! Your enquiry has been received and we will respond as soon as possible")
        IStatusMessage(self.request).addStatusMessage(confirm, type='info')
        
        # Redirect to the portal front page. Return an empty string as the
        # page body - we are redirecting anyway!
        self.request.response.redirect(portal.absolute_url())
        return ''