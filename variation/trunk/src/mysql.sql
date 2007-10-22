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
	date_created	timestamp default current_timestamp
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
