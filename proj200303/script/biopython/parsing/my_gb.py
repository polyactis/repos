#!/usr/bin/python

import sys,re


def make_gp2gb_dict(inf):
	gp2gb_dict={}
	accs=[]
	tmp=inf.readline()
	
	while tmp:
		tmp=tmp.split()[1]
		tmp=tmp.split('|')
		accs=accs+tmp
		tmp = inf.readline()
	
	for i in range(0,len(accs)/2):
		gp2gb_dict[accs[2*i+1]]=accs[2*i]
	
	return gp2gb_dict

def make_embl2sp_dict(inf):
	
	"""
	file: emblaccs
	entry like this:
	Q9SUE7	AL035524|CAB36760.1|AL161572|CAB79593.1
	
	result:
	embl2sp_dict[CAB36760.1]=Q9SUE7
	embl2sp_dict[CAB79593.1]=Q9SUE7
	"""
	embl2sp_dict={}
	tmp=inf.readline()
	anti_content=re.compile(r'\s')

	while tmp:
		first_col=tmp.split()[0]
		second_col=tmp.split()[1]
		embls=second_col.split('|')
		
		for embl in embls:
			embl2sp_dict[embl]=first_col
			
		tmp=inf.readline()

	return embl2sp_dict

def make_update_dict(inf):

	"""
	file: updateaccs
	entry like this:
	AAN71912	Q84WW6
	
	result:
	update_dict[AAN71912]=Q84WW6
	"""
	
	update_dict={}
	tmp=inf.readline()

	while tmp:
		first_col=tmp.split()[0]
		second_col=tmp.split()[1]
		
		update_dict[first_col]=second_col

		tmp=inf.readline()

	return update_dict
	
if __name__ == '__main__':
	inf = open(sys.argv[1],'r')
	gp2gb_dict = make_gp2gb_dict(inf)

	print len(gp2gb_dict)
	inf.seek(0)
	embl2sp_dict = make_embl2sp_dict(inf)
	dict=embl2sp_dict.iteritems()
	while 1:
		try:
			item=dict.next()
			print item
		except:
			break

	print len(embl2sp_dict)
	inf.close()
