#!/usr/bin/python

#sys.argv[1] is the infname.
#sys.argv[2] is the table of SETdb to be submitted.
import sys,psycopg,string


class db_submit:
	
	#submit data from 'infname' to SETdb
	
	def __init__(self,infname):
		self.inf=open(infname,'r')
		self.conn=psycopg.connect("dbname=SETdb")
		self.curs=self.conn.cursor()
		
	def read_data(self):
		self.lines=[]
		self.line=self.inf.readline()
		while self.line:
			tmp=self.line[0:len(self.line)-1]
			self.line=string.split(tmp,'\t')
			self.lines.append(self.line)
			self.line=self.inf.readline()
	
	def submit(self,table):
		print "Submitting data to "+table+"...."
		if table=='dom_sdg':
			for i in range(len(self.lines)):
				self.curs.execute('insert into dom_sdg values (%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()

		if table=='dom_prot':
			for i in range(len(self.lines)):
				self.curs.execute('insert into dom_prot values (%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
		

		if table=='protein_new':
			for i in range(len(self.lines)):
				self.curs.execute('insert into protein_new values (%s,%s,%s,%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
			
		if table=='protein':
			for i in range(len(self.lines)):
				self.curs.execute('insert into protein values (%s,%s,%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
		
			
		if table=='temp_pr':
			for i in range(len(self.lines)):
				self.curs.execute('insert into temp_pr values (%s,%s,%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
		
			
		if table=='protein_pname':
			for i in range(len(self.lines)):
				self.curs.execute('insert into protein_pname values (%s,%s,%s,%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
		
		if table=='nt_chromdb':
			for i in range(len(self.lines)):
				self.curs.execute('insert into nt_chromdb values (%s,%s,%s,%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
		
		if table=='pr_chromdb':
			for i in range(len(self.lines)):
				self.curs.execute('insert into pr_chromdb values (%s,%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
		
		if table=='aa_chromdb':
			for i in range(len(self.lines)):
				self.curs.execute('insert into aa_chromdb values (%s,%s,%s,%s,%s)',self.lines[i])
				
			self.conn.commit()
		
		if table=='gp':
			for i in range(len(self.lines)):
				self.curs.execute('insert into gp values (%s,%s)',self.lines[i])
				
			self.conn.commit()
		
if __name__=='__main__':
	db_sub=db_submit(sys.argv[1])
	db_sub.read_data()
	db_sub.submit(sys.argv[2])
