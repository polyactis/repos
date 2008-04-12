"""Interfaces for view components
"""

from zope.interface import Interface

class IMyCinemasBox(Interface):
    """A box displaying "my cinemas"
    """
    
    def cinemas():
        """Yield a list of dicts, containing keys 'title', 'address'
        and 'url'.
        """