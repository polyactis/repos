#!/usr/bin/env python
"""
Examples:
	#setup database in postgresql
	GenomeDB.py -u crocea -k genome
	
	#setup database in mysql
	GenomeDB.py -v mysql -u yh -z papaya -d genome -k ""
	
Description:
	2008-07-09
	This is a wrapper for the genome database, build on top of elixir. supercedes the table definitions in genomedb.sql.
"""
from sqlalchemy.engine.url import URL
from elixir import Unicode, DateTime, String, Integer, UnicodeText, Text
from elixir import Entity, Field, using_options, using_table_options
from elixir import OneToMany, ManyToOne, ManyToMany, OneToOne
from elixir import setup_all, session, metadata, entities
from elixir.options import using_table_options_handler	#using_table_options() can only work inside Entity-inherited class.
from datetime import datetime
from sqlalchemy.schema import ThreadLocalMetaData, MetaData
from sqlalchemy.orm import scoped_session, sessionmaker

from db import ElixirDB

__session__ = scoped_session(sessionmaker(autoflush=False, transactional=False))
#__metadata__ = ThreadLocalMetaData() #2008-11-04 not good for pylon

__metadata__ = MetaData()

class SequenceType(Entity):
	"""
	2008-07-27
		a table storing meta information to be referenced by other tables
	"""
	type = Field(String(256), unique=True)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='sequence_type')
	using_table_options(mysql_engine='InnoDB')

class RawSequence(Entity):
	"""
	2008-07-27
		to store chunks of sequences of entries from AnnotAssembly
	"""
	annot_assembly_gi = Field(Integer)
	start = Field(Integer)
	stop = Field(Integer)
	sequence = Field(String(10000))	#each fragment is 10kb
	using_options(tablename='raw_sequence')
	using_table_options(mysql_engine='InnoDB')

class AnnotAssembly(Entity):
	"""
	2008-07-27
		table to store meta info of chromosome sequences
	"""
	gi = Field(Integer, primary_key=True)
	acc_ver = Field(String(32), unique=True)
	accession = Field(String(32))
	version = Field(Integer)
	tax_id = Field(Integer)
	chromosome = Field(String(256))
	start = Field(Integer)
	stop = Field(Integer)
	orientation = Field(String(1))
	sequence = Field(String(10000))
	raw_sequence_start = ManyToOne('RawSequence', colname='raw_sequence_start_id', onupdate='CASCADE')
	sequence_type = ManyToOne('SequenceType', colname='sequence_type_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(Text)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='annot_assembly')
	using_table_options(mysql_engine='InnoDB')

class EntrezgeneType(Entity):
	"""
	2008-07-28
		store the entrez gene types
	"""
	type = Field(String(256), unique=True)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='entrezgene_type')
	using_table_options(mysql_engine='InnoDB')

class EntrezgeneMapping(Entity):
	"""
	2008-07-27
		table to store position info of genes
	"""
	gene = ManyToOne('Gene', colname='gene_id', primary_key=True, ondelete='CASCADE', onupdate='CASCADE')
	tax_id = Field(Integer)
	genomic_accession = Field(String(32))
	genomic_version = Field(Integer)
	genomic_annot_assembly = ManyToOne('AnnotAssembly', colname='genomic_gi', ondelete='CASCADE', onupdate='CASCADE')
	chromosome = Field(String(256))
	strand = Field(String(4))
	start = Field(Integer)
	stop = Field(Integer)
	entrezgene_type = ManyToOne('EntrezgeneType', colname='entrezgene_type_id', ondelete='CASCADE', onupdate='CASCADE')
	comment = Field(Text)
	gene_commentaries = OneToMany('GeneCommentary')
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='entrezgene_mapping')
	using_table_options(mysql_engine='InnoDB')

class GeneCommentaryType(Entity):
	"""
	2008-07-28
	"""
	type = Field(String(256), unique=True)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='gene_commentary_type')
	using_table_options(mysql_engine='InnoDB')

class GeneCommentary(Entity):
	"""
	2008-07-28
		store different mRNAs/Peptides from the same gene or mRNA
	"""
	accession = Field(String(32))
	version = Field(Integer)
	gi = Field(Integer)
	gene = ManyToOne('EntrezgeneMapping', colname='gene_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_commentary = ManyToOne('GeneCommentary', colname='gene_commentary_id', ondelete='CASCADE', onupdate='CASCADE')
	gene_commentaries = OneToMany('GeneCommentary')
	start = Field(Integer)
	stop = Field(Integer)
	gene_commentary_type = ManyToOne('GeneCommentaryType', colname='gene_commentary_type_id', ondelete='CASCADE', onupdate='CASCADE')
	label = Field(Text)
	text = Field(Text)
	comment = Field(Text)
	gene_segments = OneToMany('GeneSegment')
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='gene_commentary')
	using_table_options(mysql_engine='InnoDB')
	
	def getSequence(self, box_ls):
		"""
		2009-01-03
		"""
		from db import get_sequence_segment
		seq = ''
		curs = __metadata__.bind	#or self.table.bind or self.table.metadata.bind
		for box in box_ls:
			genomic_start, genomic_stop = box[:2]
			seq += get_sequence_segment(curs, self.gene.genomic_gi, genomic_start, genomic_stop, AnnotAssembly.table.name,\
									RawSequence.table.name)
		if self.gene.strand=='-1' and seq:	#need to reverse complement
			from Bio.Seq import Seq
			from Bio.Alphabet import IUPAC
			seq1 = Seq(seq, IUPAC.unambiguous_dna)
			seq2 = seq1.reverse_complement()
			seq = seq2.data
		return seq
	
	def getCDSsequence(self):
		"""
		2009-01-03 implemented
		2008-09-22 implement later
		"""
		if not hasattr(self, 'protein_box_ls'):
			self.construct_protein_box()
		
		self.cds_sequence = self.getSequence(self.protein_box_ls)
		self.mrna_sequence = self.getSequence(self.mrna_box_ls)
	
	def construct_mrna_box(self):
		"""
		2008-09-22
		"""
		self.mrna_box_ls = []
		for gene_segment in self.gene_segments:
			self.mrna_box_ls.append([gene_segment.start, gene_segment.stop])
		self.mrna_box_ls.sort()
	
	def construct_protein_box(self):
		"""
		2008-09-22
			find the corresponding protein gene_commentary if it's available. and construct a protein_box_ls
		"""
		self.protein_box_ls = []
		if len(self.gene_commentaries)==1:
			gene_commentary = self.gene_commentaries[0]
		elif len(self.gene_commentaries)>1:
			print 'Warning: more than 1 gene_commentaries for this commentary id=%s, gene_id=%s.'%(self.id, self.gene_id)
			gene_commentary = self.gene_commentaries[0]
		else:
			gene_commentary = None
		if gene_commentary:
			for gene_segment in gene_commentary.gene_segments:
				self.protein_box_ls.append([gene_segment.start, gene_segment.stop])
			self.protein_label = gene_commentary.label
			self.protein_comment = gene_commentary.comment
			self.protein_text = gene_commentary.text
		else:
			self.protein_label = None
			self.protein_comment = None
			self.protein_text = None
		self.protein_box_ls.sort()
	
	def construct_annotated_box(self):
		"""
		2008-10-01
			fix a bug that coordinates of a whole untranslated mrna block are replaced by that of a protein block
		2008-10-01
			fix a bug that a whole untranslated mrna block got totally omitted.
		2008-09-22
			combine mrna_box_ls and protein_box_ls to partition the whole gene into finer segments.
			box_ls = []	#each entry is a tuple, (start, stop, box_type, is_translated, protein_box_index)
			box_type = 'intron' or 'exon'. is_translated = 0 or 1. if it's translated, protein_box_index is index in protein_box_ls.
		"""
		self.box_ls = []	#each entry is a tuple, (start, stop, box_type, is_translated, protein_box_index)
		if not hasattr(self, 'mrna_box_ls'):
			self.construct_mrna_box()
		
		if not hasattr(self, 'protein_box_ls'):
			self.construct_protein_box()
		
		no_of_mrna_boxes = len(self.mrna_box_ls)
		no_of_prot_boxes = len(self.protein_box_ls)
		j = 0	#index in protein_box_ls
		for i in range(no_of_mrna_boxes):
			mrna_start, mrna_stop = self.mrna_box_ls[i]
			if i>0:	#add intron if this is not the first exon
				intron_start = self.mrna_box_ls[i-1][1]+1	#stop of previous exon + 1
				intron_stop = mrna_start-1	#start of current exon + 1
				self.box_ls.append((intron_start, intron_stop, 'intron', 0, None))
			if j<no_of_prot_boxes:
				prot_start, prot_stop = self.protein_box_ls[j]
				if prot_start>=mrna_start and prot_stop<=mrna_stop:
					if prot_start>mrna_start:	#one untranslated exon
						self.box_ls.append((mrna_start, prot_start-1, 'exon', 0, None))
					self.box_ls.append((prot_start, prot_stop, 'exon', 1, j))	#this is the only translated box
					if prot_stop<mrna_stop:	#one more untranslated exon
						self.box_ls.append((prot_stop+1, mrna_stop, 'exon', 0, None))
					j += 1	#push protein box index up
				elif prot_stop<mrna_start:	#not supposed to happen
					sys.stderr.write("Error: protein box: [%s, %s] of gene id=%s is ahead of mrna box [%s, %s].\n"%\
									(prot_start, prot_stop, self.gene_id, mrna_start, mrna_stop))
				elif prot_start>mrna_stop:
					self.box_ls.append((mrna_start, mrna_stop, 'exon', 0, None))	#2008-10-1
				elif prot_start<=mrna_stop and prot_stop>mrna_stop:	#not supposed to happen
					sys.stderr.write("Error: protein box: [%s, %s] of gene id=%s is partial overlapping of mrna box [%s, %s].\n"%\
									(prot_start, prot_stop, self.gene_id, mrna_start, mrna_stop))
			else:
				self.box_ls.append((mrna_start, mrna_stop, 'exon', 0, None))
		
	
class GeneSegment(Entity):
	"""
	2008-07-28
		table to store position info of segments (exon, intron, CDS, UTR ...) of gene-products (mRNA, peptide) in table GeneProduct
	"""
	gene_commentary = ManyToOne('GeneCommentary', colname='gene_commentary_id', ondelete='CASCADE', onupdate='CASCADE')
	gi = Field(Integer)
	start = Field(Integer)
	stop = Field(Integer)
	gene_commentary_type = ManyToOne('GeneCommentaryType', colname='gene_commentary_type_id', ondelete='CASCADE', onupdate='CASCADE')
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='gene_segment')
	using_table_options(mysql_engine='InnoDB')

class Gene(Entity):
	"""
	2008-07-27
		table to store meta info of genes
	"""
	tax_id = Field(Integer)
	gene_id = Field(Integer, primary_key=True)
	gene_symbol = Field(String(128))
	locustag = Field(String(128))
	synonyms = Field(Text)
	dbxrefs = Field(Text)
	chromosome = Field(String(256))
	map_location = Field(String(256))
	description = Field(Text)
	type_of_gene = Field(String(128))
	symbol_from_nomenclature_authority = Field(String(256))
	full_name_from_nomenclature_authority = Field(Text)
	nomenclature_status = Field(String(64))
	other_designations = Field(Text)
	modification_date = Field(DateTime, default=datetime.now)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='gene')
	using_table_options(mysql_engine='InnoDB')

class Gene2go(Entity):
	"""
	2008-07-27
		table to store mapping between gene and GO
	"""
	tax_id = Field(Integer)
	gene = ManyToOne('Gene', colname='gene_id', ondelete='CASCADE', onupdate='CASCADE')
	go_id = Field(String(32))
	evidence = Field(String(3))
	go_qualifier = Field(String(64))
	go_description = Field(Text)
	pubmed_ids = Field(Text)
	category = Field(String(16))
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='gene2go')
	using_table_options(mysql_engine='InnoDB')

class README(Entity):
	title = Field(Text)
	description = Field(Text)
	created_by = Field(String(256))
	updated_by = Field(String(256))
	date_created = Field(DateTime, default=datetime.now)
	date_updated = Field(DateTime)
	using_options(tablename='readme')
	using_table_options(mysql_engine='InnoDB')

class Gene_symbol2id(Entity):
	"""
	2008-07-27
		a derived table from Gene in order to map all available gene names (symbol, locustag, synonym) to gene-id
	"""
	tax_id = Field(Integer)
	gene_symbol = Field(String(256))
	gene = ManyToOne('Gene', colname='gene_id', ondelete='CASCADE', onupdate='CASCADE')
	symbol_type = Field(String(64))
	using_options(tablename='gene_symbol2id')
	using_table_options(mysql_engine='InnoDB')

def getEntrezgeneAnnotatedAnchor(db, tax_id):
	"""
	2008-08-13
		similar to annot.bin.codense.common.get_entrezgene_annotated_anchor, but use elixir db interface
	"""
	sys.stderr.write("Getting entrezgene_annotated_anchor ...")
	chromosome2anchor_gene_tuple_ls = {}
	gene_id2coord = {}
	offset_index = 0
	block_size = 5000
	rows = EntrezgeneMapping.query.filter_by(tax_id=tax_id).offset(offset_index).limit(block_size)
	while rows.count()!=0:
		for row in rows:
			genomic_gi = row.genomic_gi
			chromosome = row.chromosome
			gene_id = row.gene_id
			strand = row.strand
			start = row.start
			stop = row.stop
			if strand=='1' or strand=='+1' or strand=='+':
				strand = '+'
			elif strand=='-1' or strand=='-':
				strand = '-'
			else:
				strand = strand
			
			if chromosome not in chromosome2anchor_gene_tuple_ls:
				chromosome2anchor_gene_tuple_ls[chromosome] = []
			
			chromosome2anchor_gene_tuple_ls[chromosome].append((start, gene_id))
			gene_id2coord[gene_id] = (start, stop, strand, genomic_gi)
		
		rows = EntrezgeneMapping.query.filter_by(tax_id=tax_id).offset(offset_index).limit(block_size)
	for chromosome in chromosome2anchor_gene_tuple_ls:	#sort the list
		chromosome2anchor_gene_tuple_ls[chromosome].sort()
	sys.stderr.write("Done.\n")
	return chromosome2anchor_gene_tuple_ls, gene_id2coord

class GenomeDatabase(ElixirDB):
	__doc__ = __doc__
	option_default_dict = ElixirDB.option_default_dict.copy()
	option_default_dict[('drivername', 1,)][0] = 'mysql'
	option_default_dict[('database', 1,)][0] = 'genome'
	def __init__(self, **keywords):
		"""
		2008-10-08
			simplified further by moving db-common lines to ElixirDB
		2008-07-09
		"""
		from ProcessOptions import ProcessOptions
		ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		self.setup_engine(metadata=__metadata__, session=__session__, entities=entities)
	
	def get_gene_id2model(self, tax_id=3702):
		"""
		2009-1-3
			add cds_sequence, mrna_sequence to new_gene_commentary
		2008-10-01
			sort EntrezgeneMapping by chromosome, start
		2008-10-1
			similar to variation.src.GenomeBrowser.get_gene_id2model() but adapts to the new genome db schema
			construct a data structure to hold whole-genome annotation of genes
		"""
		sys.stderr.write("Getting gene_id2model and chr_id2gene_id_ls ...\n")
		from Genome import GeneModel
		gene_id2model = {}
		chr_id2gene_id_ls = {}
		i = 0
		block_size = 5000
		query = EntrezgeneMapping.query.filter_by(tax_id=tax_id).order_by(EntrezgeneMapping.chromosome).order_by(EntrezgeneMapping.start)
		rows = query.offset(i).limit(block_size)
		while rows.count()!=0:
			for row in rows:
				chromosome = row.chromosome
				gene_id = row.gene_id
				if chromosome not in chr_id2gene_id_ls:
					chr_id2gene_id_ls[chromosome] = []
				chr_id2gene_id_ls[chromosome].append(gene_id)
				if gene_id not in gene_id2model:
					gene_id2model[gene_id] = GeneModel(gene_id=gene_id, chromosome=chromosome, gene_symbol=row.gene.gene_symbol,\
													locustag=row.gene.locustag, map_location=row.gene.map_location,\
													type_of_gene=row.entrezgene_type.type, type_id=row.entrezgene_type_id,\
													start=row.start, stop=row.stop, strand=row.strand, tax_id=row.tax_id)
					for gene_commentary in row.gene_commentaries:
						if not gene_commentary.gene_commentary_id:	#ignore gene_commentary that are derived from other gene_commentaries. they'll be handled within the parental gene_commentary.
							gene_commentary.construct_annotated_box()
							gene_commentary.getCDSsequence()
							#pass everyting to a database-independent object, easy for pickling
							new_gene_commentary = GeneModel(gene_id=gene_id, gene_commentary_id=gene_commentary.id, \
														start=gene_commentary.start, stop=gene_commentary.stop, \
														label=gene_commentary.label, comment=gene_commentary.comment, \
														text=gene_commentary.text,\
														gene_commentary_type=gene_commentary.gene_commentary_type.type,\
														protein_label=gene_commentary.protein_label,\
														protein_comment=gene_commentary.protein_comment,\
														protein_text=gene_commentary.protein_text,\
														mrna_box_ls=gene_commentary.mrna_box_ls,\
														protein_box_ls=gene_commentary.protein_box_ls,\
														box_ls=gene_commentary.box_ls,
														cds_sequence=gene_commentary.cds_sequence,
														mrna_sequence=gene_commentary.mrna_sequence)
							gene_id2model[gene_id].gene_commentaries.append(new_gene_commentary)
				else:
					sys.stderr.write("Error: gene %s already exists in gene_id2model.\n"%(gene_id))
				i += 1
			if self.report:
				sys.stderr.write("%s\t%s\t"%('\x08'*80, i))
			rows = query.offset(i).limit(block_size)
		
		sys.stderr.write("Done.\n")
		return gene_id2model, chr_id2gene_id_ls
	
if __name__ == '__main__':
	import sys, os, math
	bit_number = math.log(sys.maxint)/math.log(2)
	if bit_number>40:       #64bit
		sys.path.insert(0, os.path.expanduser('~/lib64/python'))
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
	else:   #32bit
		sys.path.insert(0, os.path.expanduser('~/lib/python'))
		sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

	from pymodule import ProcessOptions
	main_class = GenomeDatabase
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.setup()
	
	#2008-10-01	get gene model and pickle it into a file
	gene_id2model, chr_id2gene_id_ls = instance.get_gene_id2model()
	from pymodule import PassingData
	gene_annotation = PassingData()
	gene_annotation.gene_id2model = gene_id2model
	gene_annotation.chr_id2gene_id_ls = chr_id2gene_id_ls
	import cPickle
	picklef = open('/tmp/at_gene_model_pickelf', 'w')
	cPickle.dump(gene_annotation, picklef, -1)
	picklef.close()
	
	import sqlalchemy as sql
	#print dir(Gene)
	#print Gene.table.c.keys()
	#results = instance.session.query(Gene)
	#results = instance.session.execute(sql.select([Gene.table]))
	#print dir(results)
	#print dir(Gene.query)
	
	import pdb
	pdb.set_trace()
	i = 0
	block_size = 10
	rows = Gene.query.offset(i).limit(block_size)
	print dir(rows)
	while rows.count()!=0:
		print rows.count()
		for row in rows:
			i += 1
			print row.gene_id, row.gene_symbol
		if i>=5*block_size:
			break
		rows = Gene.query.offset(i).limit(block_size)
