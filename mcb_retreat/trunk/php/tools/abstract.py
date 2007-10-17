#!/usr/bin/python
"""
2007-10-17
	output abstracts into two parts, presentation and poster. each part is ordered by submitter/1st author's last name.
"""

import cgi,psycopg
import cgitb
cgitb.enable()

import psycopg,sys
 
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
		for abstract_id in abstract_id_ls:
			curs.execute("select pref,title,author_list, abstract from %s where id=%s"%(abstract_table, abstract_id))
			rows = curs.fetchall()
			pref, title, author_list, abstract = rows[0]
			print title
			print author_list
			print abstract
			print
	
	def output(self):
		presentation_abstract_id_ls = self.get_abstract_id_ls(self.curs, pref='Presentation')
		print '<p align=center ><h>Presentation</h></p>'
		print '<pre>'
		self.output_abstract_in_order(self.curs, presentation_abstract_id_ls)
		print '</pre>'
		poster_abstract_id_ls = self.get_abstract_id_ls(self.curs, pref='Poster')
		
		print '<p align=center ><h>Poster</h></p>'
		print '<pre>'
		self.output_abstract_in_order(self.curs, poster_abstract_id_ls)
		print '</pre>'

if __name__ == '__main__':
	print "Content-Type: text/html\n\n"
	abstract_output_ins = abstract_output()
	abstract_output_ins.output()