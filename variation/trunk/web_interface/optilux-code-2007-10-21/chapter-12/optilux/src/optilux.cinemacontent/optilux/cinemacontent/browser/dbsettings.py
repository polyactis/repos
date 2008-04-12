from zope.component import getUtility
from zope.formlib import form

from plone.app.controlpanel.form import ControlPanelForm

from collective.lead.interfaces import IDatabase

from optilux.cinemacontent.interfaces import IDatabaseSettings
from optilux.cinemacontent import CinemaMessageFactory as _

def reservations_database_settings(context):
    """When the control panel reads and writes form values, it will adapt
    its context (the Plone site root) to the form field interface 
    (IReservationDatabaseSettings), and read/write properties on this adapter.
    
    Since we want to write to a local/persistent utility which we know
    conforms to this interface, we can just use this. Therefore, we use
    this custom adapter factory to fetch this object in response to the
    adaptation.
    """
    return getUtility(IDatabaseSettings)
    
class ReservationsDatabaseControlPanel(ControlPanelForm):
    """Control panel form view for the Reservations Database settings.
    
    This uses zope.formlib to present a form from the interface above.
    """

    form_fields = form.FormFields(IDatabaseSettings)

    form_name = _(u"Reservations Database settings")
    label = _(u"Reservations Database settings")
    description = _(u"Please enter the appropriate connection settings for the database")
    
    def _on_save(self, data):
        """This method is called when the form is successfully saved. We use
        it to inform the database engine that the settings have changed.
        """
        db = getUtility(IDatabase, name='optilux.reservations')
        db.invalidate()
