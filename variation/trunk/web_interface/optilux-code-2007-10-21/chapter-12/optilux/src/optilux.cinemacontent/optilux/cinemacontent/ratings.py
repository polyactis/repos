"""Ratings for films
"""

from BTrees.OOBTree import OOSet

from zope.interface import implements
from zope.component import adapts

from zope.annotation.interfaces import IAnnotations

from optilux.cinemacontent.interfaces import IFilm
from optilux.cinemacontent.interfaces import IRatings

POSITIVE_KEY = 'optilux.cinemacontent.ratings.positive'
NEGATIVE_KEY = 'optilux.cinemacontent.ratings.negative'

class FilmRatings(object):
    """Rate a film
    
    The user_token is a username or IP address or other string identifying
    users. We use this to try to avoid people voting more than once.
    
    Here is how it works. First, we create a dummy film. We need to make sure 
    it's annotatable. The standard Film content type is attribute annotatable 
    because all content in Plone is.
        
        >>> from zope.interface import implements
        >>> from zope.annotation.interfaces import IAttributeAnnotatable
        
        >>> from optilux.cinemacontent.interfaces import IFilm
        
        >>> class DummyFilm(object):
        ...     implements(IFilm, IAttributeAnnotatable)
        ...     film_code = u""
        ...     title = u""
        ...     summary = u""
        ...     teaser = u""
        ...     shown_from = None
        ...     shown_until = None
        
    We need to make sure the rating adapter is configured. Normally, this 
    would happen during ZCML processing. We also need the standard annotation
    adapter.
    
        >>> from zope.component import provideAdapter
        >>> from zope.annotation.attribute import AttributeAnnotations
        >>> provideAdapter(AttributeAnnotations)
    
        >>> from optilux.cinemacontent.ratings import FilmRatings
        >>> provideAdapter(FilmRatings)
    
    Now we can adapt a film to a IRatings and rate it.
    
        >>> test_film = DummyFilm()
        >>> ratings = IRatings(test_film)
        
        >>> ratings.score is None
        True
        
    Let's rate as different users. The score is calculated as the percentage
    of positive votes, returned as an integer. Once a user has voted once, he 
    cannot vote again.
        
        >>> ratings.available('user1')
        True
        
        >>> ratings.rate('user1', True)
        >>> ratings.available('user1')
        False
        >>> ratings.score
        100
        
        >>> ratings.rate('user1', True) # doctest: +ELLIPSIS
        Traceback (most recent call last):
        ...
        KeyError: 'Ratings not available for user1'
        
        >>> ratings.rate('user2', False)
        >>> ratings.score
        50
        
        >>> ratings.rate('user3', True)
        >>> ratings.score
        66
    """
    implements(IRatings)
    adapts(IFilm)
    
    def __init__(self, context):
        self.context = context
        
        # We assume IFilm is annotatable, in which case we can adapt to
        # IAnnotations and get a mapping-like object back. We manage all
        # ratings under a particular key. We store all the positive votes
        # in one key and all the negative votes in another. We can then
        # just count the number of items in each set to work out the score.
        
        annotations = IAnnotations(context)
        self.positive = annotations.setdefault(POSITIVE_KEY, OOSet())
        self.negative = annotations.setdefault(NEGATIVE_KEY, OOSet())
    
    @property
    def score(self):
        positives = len(self.positive)
        negatives = len(self.negative)
        total = positives + negatives
        
        if total == 0:
            return None
        else:
            return int((float(positives) / total) * 100)
        
    def available(self, user_token):
        return not (self.positive.has_key(user_token) or 
                    self.negative.has_key(user_token))
                    
    def rate(self, user_token, positive):
        if not self.available(user_token):
            raise KeyError("Ratings not available for %s" % user_token)
        if positive:
            self.positive.insert(user_token)
        else:
            self.negative.insert(user_token)