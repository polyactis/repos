create schema dbsnp;

set search_path to dbsnp;

create table strain_info(
	id	integer,
	acc	varchar primary key,
	tax_id	integer,
	name	varchar,
	category	varchar,
	description	varchar,
	date_created	timestamp default current_timestamp,
	pop_location	varchar,	--below until -- is from Diane.xls
	pop	varchar,
	line	integer,
	barcode	varchar,
	distance	float,
	latitude	float,
	longitude	float,
	comments	varchar,
	abrc_stock_acc	varchar,	--below is from 850 Natural Accessions.csv
	parental_stock_acc	varchar,
	abbr_name	varchar,
	full_name	varchar,
	days_to_flower	integer,	--negative means truncated >75 => -75
	country	varchar
	);

create table snp_locus(
	id	integer,
	acc	varchar primary key,
	tax_id	integer,
	chromosome	varchar,
	strand	varchar(1),
	align_2010_start_pos	integer,	--2007-04-01
	freq	float,
	name	varchar,
	description	varchar,
	date_created	timestamp default current_timestamp,
	--2007-04-01	add UEP_DIR, adjacent_sequence
	position	integer,
	UEP_DIR	varchar,
	adjacent_sequence	varchar,
	lyrata_call	varchar,	--2007-04-30
	thaliana_call	varchar,
	flanking_25mer	varchar
	);

create table justin_data(
	id	serial primary key,
	strain_id	integer,
	snp_id	integer,
	call	varchar,
	comment	varchar
	);

create or replace view view_justin_data as SELECT s1.acc as strain_acc, s1.category as strain_category,
	s2.acc as snp_acc, j.call from strain_info s1, snp_locus s2, justin_data j  
	where s1.id = j.strain_id and s2.id=j.snp_id order by strain_category, strain_acc, snp_acc;

create table readme(
	id	serial primary key,
	name	varchar,
	description	varchar,
	date_created	timestamp default current_timestamp
	);

create table strain_info_2010(
	acc	varchar primary key,
	region	varchar,
	population	varchar,
	plate_pos_col	varchar,
	plate_pos_row	varchar,
	stock_center_number	varchar
	);

--2007-04-30
create table snp_locus_context(
	id	serial primary key,
	snp_locus_id	integer,
	disp_pos	integer,
	gene_id	integer,
	gene_strand	varchar(1),
	disp_pos_comment	varchar
	);

--2007-05-15
create table longlat_192(
	name	varchar,
	accession	integer,
	original_name	varchar,
	longgrads	float,
	latgrads	float,
	rlonggrads	float,
	rlatgrads	float,
	longmins	integer,
	latmins	integer,
	colorcoller	integer,
	starplus	integer,
	shapecoller	integer,
	fricoller	integer,
	daysLDnoV	float,
	LeavesLDnoV	float,
	LeavesLDwithV	float,
	FLNnoV	float,
	FLNwithV	float,
	accessionmaria	varchar,
	longit	varchar,	
	lat	varchar,
	id	integer,
	accname	varchar,
	popname	varchar,
	regionDB	varchar,
	longitude	float,
	latitude	float,
	difflong	float,
	difflat	float,
	location	varchar,
	misc1	varchar,
	misc2	float,
	misc3	varchar,
	misc4	float
	);

--2007-06-09 table accession of database at from walnut
create table at_accession(
	id	integer,
	population	integer,
	region	integer,
	description	varchar,
	name	varchar,
	origin	varchar,
	date_created	timestamp,
	date_modified	timestamp,
	number	varchar
	);
