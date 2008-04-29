#!/usr/bin/python
"""
2007-10-17
	output (to stdout) abstracts into two parts, presentation and poster. each part is ordered by submitter/1st author's last name.
"""

import cgi
import cgitb
cgitb.enable()

import sys
import psycopg2 as psycopg
 
class abstract_output:
	"""
	2007-10-17
	"""
	def __init__(self):
		self.conn = psycopg.connect("host=localhost  dbname=yhdb user=yh password=123456")
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to retreat")
	
	def get_abstract_id_ls(self, curs, pref='Presentation', abstract_table='abstract'):
		curs.execute("select id, name from %s where pref='%s'"%(abstract_table, pref))
		rows = curs.fetchall()
		abstract_last_name_id_ls = []
		for row in rows:
			abstract_id, name = row
			last_name = name.split()[-1]
			abstract_last_name_id_ls.append((last_name, abstract_id))
		abstract_last_name_id_ls.sort()
		abstract_id_ls = []
		for last_name, abstract_id in abstract_last_name_id_ls:
			abstract_id_ls.append(abstract_id)
		return abstract_id_ls
	
	def output_abstract_in_order(self, curs, abstract_id_ls, abstract_table='abstract'):
		"""
		2008-04-29
			add name, email, pi, abstract_filename
		"""
		no_of_abstracts = 0
		for abstract_id in abstract_id_ls:
			curs.execute("select name, email, pi, pref, title, author_list, abstract, abstract_filename from %s where id=%s"%(abstract_table, abstract_id))
			rows = curs.fetchall()
			name, email, pi, pref, title, author_list, abstract, abstract_filename = rows[0]
			no_of_abstracts += 1
			print '%s %s'%(pref, no_of_abstracts)
			print "Name:", name
			print "Email:", email
			print "PI:", pi
			print "Title:", title
			print "Author List:", author_list
			if abstract:
				print "Abstract:", abstract
			elif abstract_filename:
				print "Abstract Filename:", abstract_filename
			else:
				print "Abstract:", None
			print
	
	def split_name_into_first_last_name(self, curs, register_table='register', commit=0):
		"""
		circa 2007-10-20
			split name into lastname and firstname for sorting purpose
		"""
		curs.execute("select id, name from %s"%register_table)
		rows = curs.fetchall()
		for row in rows:
			row_id, name = row
			name_ls = name.split()
			lastname = name_ls[-1]
			lastname = lastname[0].upper()+lastname[1:]	#turn 1st character into capital form
			firstname = ' '.join(name_ls[:-1])
			curs.execute("update %s set lastname='%s', firstname='%s' where id=%s"%(register_table, lastname, firstname, row_id))
		if commit:
			self.conn.commit()
	
	def output(self):
		print '<pre>'
		presentation_abstract_id_ls = self.get_abstract_id_ls(self.curs, pref='Presentation')
		self.output_abstract_in_order(self.curs, presentation_abstract_id_ls)
		poster_abstract_id_ls = self.get_abstract_id_ls(self.curs, pref='Poster')
		self.output_abstract_in_order(self.curs, poster_abstract_id_ls)
		print '</pre>'
		self.split_name_into_first_last_name(self.curs, register_table='register', commit=1)

if __name__ == '__main__':
	print "Content-Type: text/html\n\n"
	abstract_output_ins = abstract_output()
	abstract_output_ins.output()