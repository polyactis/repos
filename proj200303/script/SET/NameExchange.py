#!/usr/bin/python

import sys,re,psycopg,string
from biopython.parsing import my_gb

def value_pacc2embl_pr(pacc,embl2sp_dict,gp2gb_dict):

	
	dict_iterator = embl2sp_dict.iteritems()
	returnvalue = []
	while 1:
		try:
			item = dict_iterator.next()
			if item[1] == pacc and gp2gb_dict.has_key(item[0]):
				returnvalue.append( item[0])
		except:
			if returnvalue == []:
				return 'None'
			else:
				return '|'.join(returnvalue)

def value_pacc2embl_nt(pacc,embl2sp_dict,gp2gb_dict):

	
	dict_iterator = embl2sp_dict.iteritems()
	returnvalue = []
	while 1:
		try:
			item = dict_iterator.next()
			if item[1] == pacc and gp2gb_dict.has_key(item[0]):
				returnvalue.append( gp2gb_dict[item[0]])
		except:
			if returnvalue == []:
				return 'None'
			else:
				return '|'.join(returnvalue)

def value_pname2pacc(pname):
	conn = psycopg.connect("dbname=SETdb")
	curs = conn.cursor()

	curs.execute("select pacc from pname where pname= %s",(pname,))
	row = curs.fetchall()
	if row:
		return row[0][0]
	else:
		return "None"

def run_pacc2embl_pr(inf,outf,emblaccs_f):
	
	emblaccs_f.seek(0)
	embl2sp_dict = my_gb.make_embl2sp_dict(emblaccs_f)
	emblaccs_f.seek(0)
	gp2gb_dict = my_gb.make_gp2gb_dict(emblaccs_f)
	
	line = inf.readline()

	while line:
		tmp = string.split(line,'\t')[2]
		key = tmp[0:len(tmp)-1]
		sys.stderr.write('.')
		value = value_pacc2embl_pr(key,embl2sp_dict,gp2gb_dict)
		if value != "None":
			outf.write("%s\t%s\n" % (line[0:len(line)-1],value))
		else:
			sys.stderr.write("%s: Accession No. not found\n" % key)
			outf.write(line)

		line = inf.readline()
	
def run_pacc2embl_nt(inf,outf,emblaccs_f):
	
	emblaccs_f.seek(0)
	embl2sp_dict = my_gb.make_embl2sp_dict(emblaccs_f)
	emblaccs_f.seek(0)
	gp2gb_dict = my_gb.make_gp2gb_dict(emblaccs_f)
	
	line = inf.readline()

	while line:
		tmp = string.split(line,'\t')[2]
		key = tmp[0:len(tmp)-1]
		sys.stderr.write('.')
		value = value_pacc2embl_nt(key,embl2sp_dict,gp2gb_dict)
		if value != "None":
			outf.write("%s\t%s\n" % (line[0:len(line)-1],value))
		else:
			sys.stderr.write("%s: Accession No. not found\n" % key)
			outf.write(line)

		line = inf.readline()
	
def run_fasta_title_update(inf,outf,sdgnames_f):
	
	from Bio import Fasta
	parser = Fasta.RecordParser()
	fasta_iterator = Fasta.Iterator(inf,parser)

	
	sdgnames_f.seek(0)
	sdgnames_dict = my_gb.make_update_dict(sdgnames_f)
	
	while 1:
		record = fasta_iterator.next()

		if record is None:
			break
		
		if sdgnames_dict.has_key(record.title):
		
			value = sdgnames_dict[record.title]
		else:
			value = record.title

		outf.write('>%s\n%s\n' % (value,record.sequence))

def run_embl_pr2est(inf,outf,embl_pr2est_f):
	
	embl_pr2est_f.seek(0)
	embl_pr2est_dict = my_gb.make_update_dict(embl_pr2est_f)
	
	line = inf.readline()
	while line:
		tmp = string.split(line)
		key = tmp[len(tmp)-1]
		
		if embl_pr2est_dict.has_key(key):
			value = embl_pr2est_dict[key]
			line = line[0:len(line)-1] + '\t' + value + '\n'
		else:
			sys.stderr.write("%s: EST not found\n" % key)
			
		outf.write(line)

		line = inf.readline()
		
def run_sp_acc_update(inf,outf,updateaccs_f):
	
	updateaccs_f.seek(0)
	update_dict = my_gb.make_update_dict(updateaccs_f)
	
	line = inf.readline()
	while line:
		tmp = string.split(line)
		key = tmp[len(tmp)-1]
		
		if update_dict.has_key(key):
			value = update_dict[key]
			line = tmp[0] + '\t' + tmp[1] + '\t' + value + '\n'
		else:
			sys.stderr.write("%s: no update_acc found\n" % key)
			
		outf.write(line)

		line = inf.readline()
		
		

if __name__ == '__main__':
	inf = open(sys.argv[1],'r')
	outf = open(sys.argv[2],'w')
	
	#emblaccs_f = open('/home/awdong/script/biopython/emblaccs','r')
	#emblaccs_f = open('/home/awdong/script/biopython/tmp.emblaccs_mini','r')
	embl_pr2est_f = open('/home/awdong/script/embl_pr2est','r')
	#sdgnames_f = open('/home/awdong/tmp.csv','r')
	updateaccs_f = open('/home/awdong/script/biopython/updateaccs','r')
	
	#p_key = re.compile(r'\b\w*\b$')
	
	run_sp_acc_update(inf,outf,updateaccs_f)
	#run_embl_pr2est(inf,outf,embl_pr2est_f)
	#run_pacc2embl_nt(inf,outf,emblaccs_f)
	#run_fasta_title_update(inf,outf,sdgnames_f)
	
	inf.close()
	outf.close()
