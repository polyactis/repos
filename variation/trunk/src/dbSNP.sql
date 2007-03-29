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
	position	integer,
	freq	float,
	name	varchar,
	description	varchar,
	date_created	timestamp default current_timestamp
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
