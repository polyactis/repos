# Roman numeral conversion
# Copyright (c) 2001 Mark Pilgrim
# This file is free software; you can redistribute it and/or modify it
# under the terms of the GPL-compatible Python 1.6.1 license:
#	 http://hdl.handle.net/1895.22/1013

# We'll be using regular expressions to match partial Roman numerals
# in fromRoman.  Read all about regular expressions in Python here:
#	 http://py-howto.sourceforge.net/regex/regex.html
import re

# What we really want here is an ordered structure of related elements,
# so we're storing them in a tuple of tuples.  Tuples are faster than
# lists but can be iterated through in the same manner, and they're
# perfect for static data structures like this.  We do not want to
# use dictionaries for this, because dictionaries do not retain their
# defined order (mapping have no order).
# The re.compile stuff pre-compiles a regular expression that we'll
# use in fromRoman to detect any valid number of occurrences of each
# symbol as we scan through the Roman numeral.
rom = (('M', 1000, re.compile('^MM?M?')), 
	('CM', 900, re.compile('^CM')), 
	('D', 500, re.compile('^D')), 
	('CD', 400, re.compile('^CD')), 
	('C', 100, re.compile('^CC?C?')), 
	('XC', 90, re.compile('^XC')), 
	('L', 50, re.compile('^L')), 
	('XL', 40, re.compile('^XL')), 
	('X', 10, re.compile('^XX?X?')), 
	('IX', 9, re.compile('^IX')), 
	('V', 5, re.compile('^V')), 
	('IV', 4, re.compile('^IV')), 
	('I', 1, re.compile('^II?I?')))

def toRoman(n):
	result = ""
	
	# assert is a Python built-in that takes 2 arguments: an expression
	# to evaluate, and an error message.
	# If the expression evaluates to 0, '', [], (), {}, None, or any other
	# false value, assert will raise an AssertionError with the given
	# error message.
	assert 0 < n < 4000, "number out of range (must be 1..3999)"
	
	# We can quickly iterate through our data structure (a tuple of tuples)
	# all pull out each value from each tuple in turn.  Called
	# "tuple unpacking" or "multi-variable assignment".  See
	#	 http://diveintopython.org/odbchelper_multiassign.html
	for roman, arabic, pattern in rom:
		while n >= arabic:
			# Python 2.0 supports augmented assignment, and it's
			# faster than doing result = result + roman.
			result += roman
			n -= arabic
	return result
	
def fromRoman(s):
	result = 0
	# We'll be munging s as we scan through it, so save the original version
	# for error reporting later.
	original = s
	
	# Python 2.0 supports string methods.  In earlier versions, you would
	# do this instead:
	#	 import string
	#	 s = string.upper(s)
	s = s.upper()
	
	for roman, arabic, pattern in rom:
		# pattern is the pre-compiled regular expression pattern,
		# which has methods to search a given string for the
		# regular expression the pattern was compiled with.
		match = pattern.search(s)

		# If the regular expression engine found a match,
		# the search method returns a MatchObject, otherwise
		# it returns None.  None is a false value in Python,
		# but a MatchObject (and almost any other class instance)
		# is true, so we can just say "if match:" here instead of
		# "if match <> None:" or some other more verbose incarnation.
		if match:
			# The magic starts here.  If we found a match,
			# what we found was 1 or more occurrences of the
			# same Roman numeral digit ("M", or "CM", or "D", or
			# whatever) in our data structure.  That means that
			# the arabic value of this found fragment is the
			# arabic value of a single digit * the number of
			# digits found.  match.group() gives us the found
			# fragment, but the length of that string may not
			# be the number of digits found because individual
			# Roman digits can be more than 1 character (e.g. "CM").
			# So the total arabic value is the arabic value of the
			# digit, times the character length of the found fragment,
			# divided by the character length of the digit.
			lenmatch = len(match.group())
			result += arabic * lenmatch / len(roman)
			
			# Now hack the processed part off the beginning of
			# the string and continue.
			s = s[lenmatch:]
	
	# 1. If s has anything left, that means the Roman numeral was mal-formed, 
	# because otherwise our regular expressions would have caught it.  This
	# depends on the fact that Roman numerals have to use the highest digits
	# possible, and our data structure was ordered from highest to lowest.
	# 2. locals() is a function that returns a dictionary of local variables, 
	# with keys being the variable names as strings and values being the
	# actual variables' values.  See
	#	 http://diveintopython.org/dialect_locals.html
	# 3. The %(original)s is dictionary-based string formatting;
	# it gets replaced by the value of the 'original' key in the given
	# dictionary; since the given dictionary is a dictionary of local
	# variables, this effectively inserts the value of <original> into
	# the string.  See
	#	 http://diveintopython.org/dialect_dictsub.html
	assert s == '', "invalid Roman numeral: %(original)s" % locals()
	
	return result
	
def test():
	# If you're serious about unit testing, check out PyUnit, which will be
	# part of the standard library in Python 2.1:
	#	 http://pyunit.sourceforge.net/

	# Test all the good values.  range is a built-in Python function that
	# returns a list of integers starting with the first parameter and
	# up to but not including the second parameter.  This is akin to
	# for (n=1; n<4000; n++) { ... } in C.
	for n in range(1, 4000):
		s = toRoman(n)
		fullcircle = fromRoman(s)
		assert n == fullcircle, 'FAILED: %(n)d --> %(s)s --> %(fullcircle)s'%locals()
		
	# Test some bad values for the arabic --> Roman conversion.
	# We expect the function to raise an AssertionError, so we
	# explicitly catch that; if toRoman doesn't raise an exception,
	# we treat that as an error and raise our own expection in response.
	# For more on exceptions, see
	#	 http://diveintopython.org/fileinfo_exception.html
	# Note that we're calling a function which is returning a result,
	# but we're never assigning the result to anything.  That's OK;
	# Python will just swallow the result and continue without error.
	for n in (4000, 0, -1):
		try:
			toRoman(n)
		except AssertionError:
			pass
		else:
			raise AssertionError, "FAILED: toRoman(%(n)d) should have failed but didn't" % locals()
	
	# Test some bad values for the Roman --> arabic conversion.
	for s in ('MMMM', 'CMCM', 'DD', 'CDCD', 'CCCC', 'XCXC', 'LL', 'XLXL', 
'XXXX', 
			'IXIX', 'VV', 'IVIV', 'IIII', 'IIMXCC', 'CMM'):
		try:
			fromRoman(s)
		except AssertionError:
			pass
		else:
			raise AssertionError, "FAILED: fromRoman('%(s)s') should have failed but didn't" % locals()
	
	# If we got this far, every test passed (otherwise we would have
	# raised an exception by now, or one of the functions we called would
	# have raised one and we intentionally wouldn't have caught it,
	# so we'd never get here).
	print "Test suite completed successfully."
	
# This is a standard trick in Python.  If this script is imported by another
# module, this module's __main__ attribute will be some descriptive value
# (generally the filename minus the .py extension).  But if the script is
# run by itself from the command line (or through an IDE), the module's
# __name__ attribute will be a special value, "__main__".  We can use this
# discrepancy to our advantage and define a test suite that will only run
# when the script is run by itself; the test suite will be completely ignored
# when the module is imported by another module.  See
#	 http://diveintopython.org/odbchelper_testing.html
if __name__ == "__main__":
	test()
	roman_str = 'VI'
	print '%s is %s.'%(roman_str, fromRoman(roman_str))

#---end of roman.py
