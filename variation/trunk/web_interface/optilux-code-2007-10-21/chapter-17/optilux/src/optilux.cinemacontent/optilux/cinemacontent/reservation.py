from zope.interface import implements
from zope.component import getUtility

from optilux.cinemacontent.interfaces import IReservation
from optilux.cinemacontent.interfaces import ITicketReservations
from optilux.cinemacontent.interfaces import ReservationError

from optilux.cinemacontent.screening import Screening
from optilux.cinemacontent import CinemaMessageFactory as _

import sqlalchemy as sql
from collective.lead.interfaces import IDatabase

class Reservation(object):
    """A ticket reservation for a particular screening
    """
    
    implements(IReservation)
    
    customer_name = u""
    num_tickets = 0
    screening = None
    
    def __init__(self, customer_name, num_tickets, screening):
        self.customer_name = customer_name
        self.num_tickets = num_tickets
        self.screening = screening

class TicketReservations(object):
    """Make reservations in the reservations database
    """
    
    implements(ITicketReservations)
    
    def __call__(self, reservation):
        """Make a reservation
        """
        
        db = getUtility(IDatabase, name='optilux.reservations')
        session = db.session
        
        # Make sure there are still seats available
        screening = reservation.screening
        session.refresh(screening)
        
        if screening.remaining_tickets <= 0:
            raise ReservationError(_(u"This screening is sold out!"))
        elif screening.remaining_tickets < reservation.num_tickets:
            raise ReservationError(_(u"Not enough tickets remaining!"))
        
        # Otherwise, we can save the reservation
        screening.remaining_tickets -= reservation.num_tickets
        session.update(screening)
        session.save(reservation)
        session.flush()