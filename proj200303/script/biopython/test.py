#!/usr/bin/python

import sys,re,string,psycopg
from parsing import my_gb

def proteindb_check(inf,emblaccs_f):


	gp2gb_dict = my_gb.make_gp2gb_dict(emblaccs_f)
	
	line = inf.readline()
	accs = []
	while line:
		accs.append(string.split(line,'\t')[0])
		line = inf.readline()

	for item in gp2gb_dict:
		if item not in accs:
			sys.stderr.write('%s not found\n' % item)
	
def f2gb_compare(inf):

	"""
	>(Accession No.1)
	ACCESSION   (Accession No.2)
	....
	
	check whether the Accession No.1 is the same as Accession No.2.
	
	"""
	 
	p_fasta=re.compile(r'^>(\b\w*\b$)')
	p_accession=re.compile(r'^ACCESSION\s*(\b\w*\b)')
	p_end = re.compile(r'^\/\/')
	
	
	tmp = inf.readline()
	i=1
	print "No %d Record" % i
	
	while tmp:
		if p_fasta.search(tmp):
			fasta=p_fasta.search(tmp).group(1)
		if p_accession.search(tmp):
			accession=p_accession.search(tmp).group(1)
			sys.stdout.write( '%s\t%s\n' % (fasta, accession))
			if accession != fasta:
				
				sys.stderr.write('%s\t%s\n' % (fasta,accession))
			
	
		if p_end.search(tmp):
			i=i+1
			print "No %d Record" % i
		tmp=inf.readline()

def SETdb2sp_check():	
	conn = psycopg.connect("dbname=SETdb")
	curs = conn.cursor()

	curs.execute("select pacc from pname")

	rows = curs.fetchall()
	rows.sort()
	i=0
	for row in rows:
		for r in row:
			curs.execute("select pacc,sequence from protein where pacc=%s",(r,))
			row = curs.fetchall()
					
			pacc = row[0][0]
			sequence = row[0][1]
			
			curs.execute("select acc,seq from db_sp where fasta=%s",(r,))
			rows=curs.fetchall()
			
			if rows:
				for row in rows:
					sp_acc = row[0]
					sp_seq = row[1]
					if sp_seq != sequence:
						sys.stdout.write('>%dSETdb_%s\n%s\n>%ddb_sp_%s\n%s\n' % (2*i+1,pacc,sequence,2*i+2,sp_acc,sp_seq))
						i=i+1
						curs.execute('select nom from nomenclature where pacc=%s',(pacc,))
						row =curs.fetchall()
						if row:
							sys.stderr.write('%s:%s\n'%(pacc,row))

if __name__ == '__main__':
	
	#inf = open(sys.argv[1],'r')
	
	SETdb2sp_check()
	#emblaccs_f = open(sys.argv[2],'r')
	#proteindb_check(inf, emblaccs_f)
	
	#f2gb_compare(inf)
