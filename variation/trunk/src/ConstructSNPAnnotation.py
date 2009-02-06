#!/usr/bin/env python
"""

Examples:
	ConstructSNPAnnotation.py -u yh -c
	
	ConstructSNPAnnotation.py -u yh -m 0 -g -s /mnt/panfs/250k/snps_context_g1_m0 -j /Network/Data/250k/tmp-yh/at_gene_model_pickelf_with_cds_mrna_seq

Description:
	2008-1-6 program to construct the annotatioin (intergenic, synonymous, non-synonymous) of a snp.
	It fills results into db table snp_annotation.
	
	Always set min_distance=0 and get_closest=1 and pick the corresponding snps_context_picklef. The algorithm depends on that.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

import time, csv, getopt
import warnings, traceback
from pymodule import ProcessOptions, PassingData
import Stock_250kDB

#from annot.bin.codense.common import get_entrezgene_annotated_anchor
from DrawSNPRegion import DrawSNPRegion	#to handle gene_annotation_picklef
from GeneListRankTest import GeneListRankTest, SnpsContextWrapper	#to handle snps_context_picklef
from Bio.Seq import Seq	#translate cds sequence, decide whether synonymous/non SNP
from Bio.Alphabet import IUPAC
from Bio import Translate
from pymodule.SNP import nt2complement	#to complement single nucleotide

class ConstructSNPAnnotation(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							('output_fname', 0, ): [None, 'o', 1, 'if given, QC results will be outputed into it.'],\
							("gene_annotation_picklef", 0, ): [None, 'j', 1, 'given the option, If the file does not exist yet, store a pickled gene_annotation into it. If the file exists, load gene_annotation out of it.'],\
							("snps_context_picklef", 0, ): [None, 's', 1, 'given the option, if the file does not exist yet, to store a pickled snps_context_wrapper into it, min_distance and flag get_closest will be attached to the filename. If the file exists, load snps_context_wrapper out of it.'],\
							("min_distance", 1, int): [5000, 'm', 1, 'minimum distance allowed from the SNP to gene'],\
							("get_closest", 0, int): [1, 'g', 0, 'only get genes closest to the SNP within that distance'],\
							('tax_id', 0, int): [3702, '', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							('commit', 0, int):[0, 'c', 0, 'commit db transaction'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	
	def __init__(self,  **keywords):
		"""
		2009-1-5
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
	
	def getSNPAnnotationShortName2id(self):
		"""
		2009-1-5
		"""
		sys.stderr.write("Getting snp_annotation_short_name2id ...")
		snp_annotation_short_name2id = {}
		rows = Stock_250kDB.SNPAnnotationType.query.all()
		for row in rows:
			snp_annotation_short_name2id[row.short_name] = row.id
		sys.stderr.write("Done.\n")
		return snp_annotation_short_name2id
	
	def _constructSNPAnnotation(self, session, snp_info, snps_context_wrapper, gene_annotation, snp_annotation_short_name2id):
		"""
		2009-2-5
			bug fixed. when adding a box as UTR, make sure after the forloop above, the current box is still the UTR.
		2009-1-5
		"""
		sys.stderr.write("Constructing SNPAnnotation ...\n")
		standard_translator = Translate.unambiguous_dna_by_id[1]
		counter = 0
		real_counter = 0
		for i in range(len(snp_info.chr_pos_ls)):
			counter += 1
			snps_id, chr, pos, allele1, allele2 = snp_info.data_ls[i]
			snps_context_matrix = snps_context_wrapper.returnGeneLs(chr, pos)
			snp_annotation_type_short_name_ls = []	#each element is (snp_annotation_type_short_name, gene_id, gene_commentary_id, which_exon_or_intron, pos_within_codon)
			if snps_context_matrix:
				for snps_context in snps_context_matrix:
					snps_id, disp_pos, gene_id = snps_context
					gene_model = gene_annotation.gene_id2model.get(gene_id)
					if gene_model is None:
						continue
					for gene_commentary in gene_model.gene_commentaries:
						if gene_commentary.protein_box_ls:
							which_intron = -1	#which intron the SNP resides, starting from 1
							which_coding_exon = -1	#which exon the SNP resides in terms of the CDS sequence, starting from 1
							accum_intron_len = 0	#the cumulative length of all introns from the one after the first coding exon up till the SNP's position (including the intron the SNP is in if it's under the rule).
							exon_5_end_pos = -1	#5' end position of the whole CDS sequence. bad name 'exon'
							UTR_2nd = 0	#whether this is the 2nd UTR
							box_type = None
							UTR_type = None
							is_SNP_within_box = 0	#protect against the possibility that SNP overshoots box_ls
							for i in range(len(gene_commentary.box_ls)):
								box = gene_commentary.box_ls[i]
								start, stop, box_type, is_translated, protein_box_index = box
								if box_type=='exon' and is_translated==0:
									if UTR_2nd==0:
										if gene_model.strand=='-1':
											UTR_type='3UTR'
										else:
											UTR_type='5UTR'
										UTR_2nd += 1
									else:
										if gene_model.strand=='-1':
											UTR_type='5UTR'
										else:
											UTR_type='3UTR'
								elif box_type=='exon' and is_translated==1:
									if exon_5_end_pos==-1:
										exon_5_end_pos=start
									which_coding_exon += 1
								elif box_type=='intron':
									if which_coding_exon!=-1:	#exclude introns that stand before the 1st coding exon
										accum_intron_len += abs(stop-start)+1	#+1 because stop-start is intron_length-1
									which_intron += 1
								
								if pos>=start and pos<=stop:	#with this box
									is_SNP_within_box = 1
									break
							
							if gene_model.strand=='-1':
								#continue to read the box_ls to count no_of_introns
								no_of_introns = which_intron+1
								for j in range(i+1, len(gene_commentary.box_ls)):
									box = gene_commentary.box_ls[i]
									if box[2]=='intron':
										no_of_introns += 1
							
							if box_type!=None and is_SNP_within_box==1:	#box_type is the type of the box where the SNP resides
								if UTR_type!=None and box_type=='exon' and is_translated==0:	#it's UTR. bug fixed. make sure after the forloop above, the current box is still the UTR.
									snp_annotation_type_short_name_ls.append((UTR_type, gene_id, gene_commentary.gene_commentary_id))
								else:
									if gene_model.strand=='-1':	#reverse the order of exon/intron
										no_of_coding_exons = len(gene_commentary.protein_box_ls)
										#no_of_introns = len(gene_commentary.mrna_box_ls)-no_of_coding_exons	#not right
										which_coding_exon = no_of_coding_exons-which_coding_exon
										which_intron = no_of_introns-which_intron
									else:
										which_coding_exon += 1
										which_intron += 1
									if box_type=='intron':
										snp_annotation_type_short_name_ls.append(('intron', gene_id, gene_commentary.gene_commentary_id, which_intron))
										if pos-start<=1:	#within the donor/acceptor two-nucleotide
											if gene_model.strand=='-1':
												snp_annotation_type_short_name = 'splice-acceptor'	#on the 3' of this intron
											else:
												snp_annotation_type_short_name = 'splice-donor'	#on the 5' of this intron
											snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary.gene_commentary_id, which_intron))
										elif stop-pos<=1:
											if gene_model.strand=='-1':
												snp_annotation_type_short_name = 'splice-donor'	#on the 5' of this intron
											else:
												snp_annotation_type_short_name = 'splice-acceptor'	#on the 3' of this intron
											snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary.gene_commentary_id, which_intron))
									elif box_type=='exon':	#must be translated
										try:
											SNP_index_in_CDS = pos - exon_5_end_pos - accum_intron_len
											if gene_model.strand=='-1':	#reverse
												SNP_index_in_CDS = len(gene_commentary.cds_sequence)-SNP_index_in_CDS-1	#-1 because SNP_index_in_CDS starts from 0
												gene_allele1 = nt2complement[allele1]
												gene_allele2 = nt2complement[allele2]
											else:
												gene_allele1 = allele1
												gene_allele2 = allele2
											SNP_index_in_CDS = int(SNP_index_in_CDS)	#SNP_index_in_CDS is type long. without int(), cds_seq[SNP_index_in_CDS] returns a Bio.Seq with one nucleotide, rather than a single-char string
											SNP_index_in_peptide = SNP_index_in_CDS/3
											
											SNP_index_in_peptide = int(SNP_index_in_peptide)	#ditto as int(SNP_index_in_CDS), not necessary
											
											pos_within_codon = SNP_index_in_CDS%3+1	#pos_within_codon starts from 1
											cds_seq = Seq(gene_commentary.cds_sequence, IUPAC.unambiguous_dna)
											if SNP_index_in_CDS>=len(cds_seq):
												sys.stderr.write("Warning: SNP (%s, %s), SNP_index_in_CDS=%s, is beyond any of the boxes from gene %s (chr=%s, %s-%s), gene_commentary_id %s (%s-%s), box_ls=%s, cds-length=%s. counted as intergenic.\n"%\
																(chr, pos, SNP_index_in_CDS, gene_id, gene_model.chromosome, gene_model.start, gene_model.stop, \
																gene_commentary.gene_commentary_id, gene_commentary.start, gene_commentary.stop, repr(gene_commentary.box_ls), len(cds_seq)))
												snp_annotation_type_short_name_ls.append(['intergenic'])
												continue
											if cds_seq[SNP_index_in_CDS]!=gene_allele1 and cds_seq[SNP_index_in_CDS]!=gene_allele2:
												sys.stderr.write("Error: Neither allele (%s, %s) from SNP (%s,%s) matches the nucleotide, %s, from the cds seq of gene %s (gene_commentary_id=%s).\n"%\
																	(gene_allele1, gene_allele2, chr, pos, cds_seq[SNP_index_in_CDS], gene_id, gene_commentary.gene_commentary_id))
												sys.exit(3)
											cds_mut_ar = cds_seq.tomutable()
											cds_mut_ar[SNP_index_in_CDS] = gene_allele1
											peptide = standard_translator.translate(cds_mut_ar.toseq())
											
											alt_cds_mut_ar = cds_seq.tomutable()
											alt_cds_mut_ar[SNP_index_in_CDS] = gene_allele2
											alt_peptide = standard_translator.translate(alt_cds_mut_ar.toseq())
											aa = peptide[SNP_index_in_peptide]
											alt_aa = alt_peptide[SNP_index_in_peptide]
											if aa != alt_aa:
												snp_annotation_type_short_name = 'non-synonymous'
												comment = '%s->%s'%(aa, alt_aa)
											else:
												snp_annotation_type_short_name = 'synonymous'
												comment = None
											snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary.gene_commentary_id, which_coding_exon, pos_within_codon, comment))
											
											if aa != alt_aa:
												if aa=='*' or alt_aa=='*':
													snp_annotation_type_short_name = 'premature-stop-codon'	#could also be the last stop codon changing to something else and thereby extending the cds
													snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary.gene_commentary_id, which_coding_exon, pos_within_codon, comment))
												if SNP_index_in_peptide==0:
													snp_annotation_type_short_name = 'init-Met'
													snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary.gene_commentary_id, which_coding_exon, pos_within_codon, comment))
										except:
											traceback.print_exc()
											sys.stderr.write('%s.\n'%repr(sys.exc_info()))
											sys.stderr.write("Except encountered for SNP (%s, %s), gene %s (chr=%s, %s-%s), gene_commentary_id %s (%s-%s), box_ls=%s.\n"%\
												(chr, pos, gene_id, gene_model.chromosome, gene_model.start, gene_model.stop, \
												gene_commentary.gene_commentary_id, gene_commentary.start, gene_commentary.stop, repr(gene_commentary.box_ls)))
							elif box_type!=None and is_SNP_within_box==0:	#SNP is over the range of the gene. it happens when one gene has multiple alternative splicing forms. snps_context_wrapper uses the largest span to represent the gene.
								snp_annotation_type_short_name_ls.append(['intergenic'])
							else:
								sys.stderr.write("Warning: SNP (%s, %s) not in any of the boxes from gene %s (chr=%s, %s-%s), gene_commentary_id %s (%s-%s), box_ls=%s.\n"%\
												(chr, pos, gene_id, gene_model.chromosome, gene_model.start, gene_model.stop, \
												gene_commentary.gene_commentary_id, gene_commentary.start, gene_commentary.stop, repr(gene_commentary.box_ls)))
						else:
							if gene_model.type_of_gene=='pseudo':
								snp_annotation_type_short_name = gene_model.type_of_gene
							else:
								snp_annotation_type_short_name = gene_commentary.gene_commentary_type
							snp_annotation_type_short_name_ls.append((snp_annotation_type_short_name, gene_id, gene_commentary.gene_commentary_id, None))
			else:	#integenic
				snp_annotation_type_short_name_ls.append(['intergenic'])
			
			#now save everything into db
			for snp_annotation_type_tup in snp_annotation_type_short_name_ls:
				snp_annotation_type_short_name = snp_annotation_type_tup[0]
				if snp_annotation_type_short_name not in snp_annotation_short_name2id:
					ty = Stock_250kDB.SNPAnnotationType(short_name=snp_annotation_type_short_name)
					session.save(ty)
					session.flush()
					snp_annotation_short_name2id[snp_annotation_type_short_name] = ty.id
				if len(snp_annotation_type_tup)>=3:
					gene_id = snp_annotation_type_tup[1]
					gene_commentary_id = snp_annotation_type_tup[2]
				else:
					gene_id = None
					gene_commentary_id = None
				if len(snp_annotation_type_tup)>=4:
					which_exon_or_intron = snp_annotation_type_tup[3]
				else:
					which_exon_or_intron = None
				if len(snp_annotation_type_tup)>=5:
					pos_within_codon = snp_annotation_type_tup[4]
				else:
					pos_within_codon = None
				if len(snp_annotation_type_tup)>=6:
					comment = snp_annotation_type_tup[5]
				else:
					comment = None
				entry = Stock_250kDB.SNPAnnotation(snps_id=snps_id, gene_id=gene_id, gene_commentary_id=gene_commentary_id, \
												which_exon_or_intron=which_exon_or_intron, pos_within_codon=pos_within_codon,\
												comment=comment)
				entry.snp_annotation_type_id = snp_annotation_short_name2id[snp_annotation_type_short_name]
				session.save(entry)
				session.flush()
				real_counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s%s\t%s"%('\x08'*40, counter, real_counter))
		if self.report:
			sys.stderr.write("%s%s\t%s\n"%('\x08'*40, counter, real_counter))
		sys.stderr.write("Done.\n")
	
	def run(self):
		if self.debug:
			import pdb
			pdb.set_trace()
		db = Stock_250kDB.Stock_250kDB(drivername=self.drivername, username=self.db_user,
				   password=self.db_passwd, hostname=self.hostname, database=self.dbname, schema=self.schema)
		db.setup(create_tables=False)
		session = db.session
		#session.begin()
		snps_context_wrapper = GeneListRankTest.dealWithSnpsContextWrapper(self.snps_context_picklef, self.min_distance, self.get_closest)
		gene_annotation = DrawSNPRegion.dealWithGeneAnnotation(self.gene_annotation_picklef)
		snp_info = DrawSNPRegion.getSNPInfo(db)
		
		snp_annotation_short_name2id = self.getSNPAnnotationShortName2id()
		self._constructSNPAnnotation(session, snp_info, snps_context_wrapper, gene_annotation, snp_annotation_short_name2id)
		if self.commit:
			session.flush()
			session.commit()
		

if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = ConstructSNPAnnotation
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()