"""The form used to reserve tickets. This demonstrates using zope.formlib
with a named template adapter, to get a custom template.
"""

from zope.component import getUtility
from zope.formlib import form

from Acquisition import aq_inner

from Products.Five.browser.pagetemplatefile import ViewPageTemplateFile
from Products.Five.formlib import formbase

from Products.statusmessages.interfaces import IStatusMessage

from plone.app.form import named_template_adapter
from plone.app.form.validators import null_validator

from plone.memoize.instance import memoize

from optilux.cinemacontent.interfaces import IReservation
from optilux.cinemacontent.interfaces import ReservationError
from optilux.cinemacontent.interfaces import IScreeningLocator
from optilux.cinemacontent.interfaces import ITicketReservations

from optilux.cinemacontent import CinemaMessageFactory as _
from optilux.cinemacontent.reservation import Reservation

# This is registered as an adapter factory. When formlib renders the page,
# it will look up this adapter to work out which template to use.
        
reserve_screening_formview = named_template_adapter(ViewPageTemplateFile('reserve.pt'))
        
class ReserveScreeningView(formbase.PageForm):
    """Reserve tickets for a screening
    """
    
    label = _(u"Reserve tickets")
    form_fields = form.FormFields(IReservation).omit('screening')
    
    @form.action(_(u"Reserve"))
    def action_reserve(self, action, data):
        """Reserve tickets
        """
        context = aq_inner(self.context)
        screening = self.screening()
        reservations = getUtility(ITicketReservations)

        try:
            reservations(Reservation(data['customer_name'], 
                                     data['num_tickets'],
                                     screening))
        except ReservationError, e:
            IStatusMessage(self.request).addStatusMessage(e.error_message, type='error')
        else:
                
            confirm = _(u"Thank you! Your tickets will be ready for collection at the front desk.")
            IStatusMessage(self.request).addStatusMessage(confirm, type='info')
            self.request.response.redirect(context.absolute_url())
            return ''
    
    # The null_validator ensures that we can cancel the form even if some
    # required fields are not entered
    @form.action(_(u"Cancel"), validator=null_validator)
    def action_cancel(self, action, data):
        """Cancel the reservation operation
        """
        context = aq_inner(self.context)
        
        confirm = _(u"Reservation cancelled.")
        IStatusMessage(self.request).addStatusMessage(confirm, type='info')
        
        self.request.response.redirect(context.absolute_url())
        return ''
        
    @memoize
    def screening(self):
        screening_id = self.request['screening_id']
        locator = getUtility(IScreeningLocator)
        return locator.screening_by_id(screening_id)