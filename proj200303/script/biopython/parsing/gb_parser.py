#!/usr/bin/python

from Bio import GenBank
from Bio import sequtils
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import sys

from my_gb import *


"""

This program is for extracting nucleotide sequences from GenBank files.
And also for extracting specific protein sequences corresponding with Swissprot's counterparts (stored in the emblaccs_handle), its CDS sequences and its locations.

input:
GenBank File

output:
GenBank.db
Protein.db
CDS.db
Location.db

"""

def my_gb(gb_handle, emblaccs_handle):

	gp2gb_dict = make_gp2gb_dict(emblaccs_handle)
	emblaccs_handle.seek(0)
	embl2sp_dict = make_embl2sp_dict(emblaccs_handle)

	GenBankdb_handle=open('GenBank.db','w')
	CDSdb_handle=open('CDS.db','w')
	Proteindb_handle=open('Protein.db','w')
	Locationdb_handle=open('Location.db','w')
	

	feature_parser = GenBank.FeatureParser()
	gb_iterator = GenBank.Iterator(gb_handle, feature_parser)
	gene_seq=Seq('',IUPAC.ambiguous_dna)
	
	while 1:
		cur_record = gb_iterator.next()
		if cur_record is None:
			break
		# now do something with the record
		#print cur_record.dbxrefs
		organism = cur_record.annotations['organism']
		gi = cur_record.annotations['gi']
		locus = cur_record.name
		gb_acc = cur_record.id
		sp_acc = embl2sp_dict[gb_acc]
	
		for feature in cur_record.features:
			if feature.type=='CDS' and feature.qualifiers.has_key('protein_id') and (feature.qualifiers['protein_id'][0] in gp2gb_dict):
				protein_id = feature.qualifiers['protein_id'][0]
				protein_seq =  feature.qualifiers['translation'][0]
				if feature.qualifiers.has_key('gene'):
					gene = feature.qualifiers['gene'][0]
				else:
					gene = 'gene_' + protein_id
				
				if feature.sub_features:
	
					for sub_feature in feature.sub_features:
						start=sub_feature.location.start.position
						end=sub_feature.location.end.position
						gene_seq=gene_seq+cur_record.seq[start:end]
						Locationdb_handle.write("%s\t%s\t%d\t%d\t%s\n" % (gene,sub_feature.strand,start,end,gb_acc))
				else:
					start=feature.location.start.position
					end=feature.location.end.position
					gene_seq=gene_seq+cur_record.seq[start:end]
					Locationdb_handle.write("%s\t%s\t%d\t%d\t%s\n" % (gene,feature.strand,start,end,gb_acc))
					
				if feature.strand==-1:
					gene_seq=Seq(sequtils.reverse(sequtils.complement(gene_seq)),IUPAC.ambiguous_dna)
			
				GenBankdb_handle.write("%s\t%s\t%s\t%s\t%s\n" % (locus,gb_acc,gi,organism,sp_acc))
				CDSdb_handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gene,organism,gene_seq.data,protein_id,gb_acc,sp_acc))
				Proteindb_handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (protein_id,organism,protein_seq,gene,gb_acc,sp_acc))
		
		gene_seq=Seq('',IUPAC.ambiguous_dna)
			
			

if __name__=='__main__':
	gb_handle = open(sys.argv[1], 'r')
	emblaccs_handle = open(sys.argv[2],'r')
	my_gb(gb_handle,emblaccs_handle)
