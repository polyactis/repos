--2007-10-10 for database at
use at;
create table magnus_192(
	donorname	varchar(200),
	linename	varchar(200),
	barcode	varchar(200),
	stockparent	varchar(200),
	nativename	varchar(200),
	collectiondate	varchar(200),
	latitude	float,
	longitude	float,
	status	integer,
	sitename	varchar(200),
	dnastatus	varchar(200),
	bulkstatus	varchar(200),
	bulkdate	varchar(200),
	country	varchar(200),
	region	varchar(200),
	contact	varchar(200),
	stateprovince	varchar(200),
	city	varchar(200),
	description	varchar(200)
	);

create table magnus_192_status(
	id	integer,
	description	varchar(200)
	);


create table readme(
	id	integer auto_increment primary key,
	name	varchar(2000),
	description	varchar(60000),
	username	varchar(200),
	date_created	timestamp default 0,	--this column can't have "default CURRENT_TIMESTAMP" because there can be only one TIMESTAMP column with CURRENT_TIMESTAMP in DEFAULT or ON UPDATE clause.
	date_modified	TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
	);

create table magnus_192_vs_accession(
	linename	varchar(200),
	accession_id	integer
	);

create table ecotype_192_vs_accession_192(
	ecotype_id	integer,
	accession_id    integer
	);

create table ecotype2accession(
	ecotype_id	integer,
	accession_id    integer
	);

create table country2continent(
	id	integer primary key,
	name	varchar(100),
	abbr	varchar(10),
	continent	varchar(30)
	);

--2007-10-21 for database stock or its replicates
use stock;
create table person(
	id	integer primary key auto_increment,
	title	varchar(4),
	surname	varchar(40) not null,
	firstname	varchar(40) not null,
	email	varchar(100),
	donor	integer not null
	);

--2007-10-29 a copy of postgresql graphdb's dbsnp.snp_locus
use dbsnp;
create table snp_locus(
	id	integer,
	acc	varchar(20) primary key,
	tax_id	integer,
	chromosome	integer,
	strand	varchar(1),
	align_2010_start_pos	integer,
	freq	float,
	name	varchar(200),
	description	varchar(2000),
	date_created	timestamp default CURRENT_TIMESTAMP,
	position	integer,
	UEP_DIR	varchar(1),
	adjacent_sequence	varchar(200),
	lyrata_call	varchar(5),
	thaliana_call	varchar(5),
	flanking_25mer	varchar(200)
	);

create table readme(
	id	integer auto_increment primary key,
	name	varchar(2000),
	description	varchar(60000),
	username	varchar(200),
	date_created	timestamp default 0,
	date_modified	TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
	);

--2007-10-29 a table copying http://naturalsystems.uchicago.edu/naturalvariation/149SNPposII.csv
use dbsnp;
create table snps_sequenom_info(
	snpfragpos	varchar(200),
	chromosome	varchar(200),
	exact_SNP_position	integer,
	well	varchar(2),
	term	varchar(5),
	snpid	varchar(20) primary key,
	x2nd_pcrp	varchar(200),
	x1st_pcrp	varchar(200),
	amp_len	integer,
	up_conf	float,
	mp_conf	float,
	tm_nn	float,
	pcgc	float,
	pwarn	varchar(2),
	uep_dir	varchar(1),
	uep_mass	float,
	uep_seq	varchar(200),
	ext1_call	varchar(1),
	ext1_mass	float,
	ext1_seq	varchar(200),
	ext2_call	varchar(1),
	ext2_mass	float,
	ext2_seq	varchar(200)
	);

--2007-12-10 table to store 250k snps and calls
use stock;
create table probes_250k(
	id	integer auto_increment primary key,
	snpid	integer,
	seq	varchar(200),
	chromosome	integer,
	position	integer,
	allele	varchar(1),
	strand	varchar(10),
	xpos	integer,
	ypos	integer
	);

create table snps_250k(
	id	integer auto_increment primary key,
	snpid	varchar(200),
	chromosome	integer,
	position	integer,
	allele1	varchar(1),
	allele2	varchar(2)
	);

create table intensity_250k(
	id	integer auto_increment primary key,
	probeid	integer,
	intensity	integer,
	duplicate	integer,
	date_created	timestamp default 0,
	date_modified	TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
	);

create table calls_250k(
	id	integer auto_increment primary key,
	ecotypeid	integer,
	snpid	integer,
	snpcall	varchar(2),
	duplicate	integer,
	date_created	timestamp default 0,
	date_modified	TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
	);

create table calls_250k_duplicate_comment(
	id	integer auto_increment primary key,
	ecotypeid	integer,
	duplicate	integer,
	comment	varchar(2000),
	date_created	timestamp default 0,
	date_modified	TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
	);

