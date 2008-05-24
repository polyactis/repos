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

-- 2008-05-15 view to link entries in at.accession to target ecotypeid thru another view complete_2010_strains_in_stock and table stock.ecotypeid2tg_ecotypeid. target ecotypeid is the final ecotypeid for a group of duplicated ecotypeids.
create or replace view complete_2010_strains_in_stock2tg_ecotypeid as select distinct e.tg_ecotypeid, c.*, a.origin, a.number from complete_2010_strains_in_stock c, stock.ecotypeid2tg_ecotypeid e , accession a where c.ecotypeid=e.ecotypeid and c.accession_id=a.id order by accession_id;


--2008-05-18 offering the final linking between accession id and ecotype id. a view linking each accession.id to a stock.ecotypeid2tg_ecotypeid.tg_ecotypeid thru at.ecotype2accession.
create or replace view accession2tg_ecotypeid as select distinct e1.accession_id, a.name as accession_name, a.origin, a.number, e2.tg_ecotypeid, e.name, e.nativename, e.stockparent, e1.ecotype_id as intermediate_ecotype_id from ecotype2accession e1, stock.ecotype e, stock.ecotypeid2tg_ecotypeid e2, accession a where e1.accession_id=a.id and e1.ecotype_id=e2.ecotypeid and e2.tg_ecotypeid=e.id order by accession_id;

/*
create table readme(
	id	integer auto_increment primary key,
	name	varchar(2000),
	description	varchar(60000),
	created_by	varchar(200) default current_user,
	modified_by	varchar(200),
	date_created	timestamp default CURRENT_TIMESTAMP,	--this column can't have "default CURRENT_TIMESTAMP" because there can be only one TIMESTAMP column with CURRENT_TIMESTAMP in DEFAULT or ON UPDATE clause.
	date_modified	TIMESTAMP ON UPDATE CURRENT_TIMESTAMP
	);
*/

--2008-02-03 new standard readme table with two triggers
create table readme(
	id      integer auto_increment primary key,
	title    varchar(2000),
	description     varchar(60000),
	created_by      varchar(200),
	updated_by     varchar(200),
	date_created    timestamp default CURRENT_TIMESTAMP,
	date_updated   TIMESTAMP default 0
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_readme BEFORE INSERT ON readme
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_readme BEFORE UPDATE ON readme
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;


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

--2008-05-15 add foreign key constraints for tables created by GroupDuplicateEcotype.py but can't run, get "ERROR 1005 (HY000): Can't create table './stock/#sql-99_17.frm' (errno: 150)"
alter table ecotype_duplicate2tg_ecotypeid add foreign key (ecotypeid) references ecotype(id) on delete restrict on update cascade;
alter table ecotype_duplicate2tg_ecotypeid add foreign key (tg_ecotypeid) references ecotype(id) on delete restrict on update cascade;

alter table genotyping_all_na_ecotype add foreign key (ecotypeid) references ecotype(id) on delete restrict on update cascade;

alter table nativename_stkparent2tg_ecotypeid add foreign key (tg_ecotypeid) references ecotype(id) on delete restrict on update cascade;

-- 2008-05-18 view to check other info (where and country) about ecotype
create or replace view ecotype_info as select e.id as ecotypeid, e.name, e.stockparent, e.nativename, e.alias, e.siteid, s.name as site_name, c.abbr as country from ecotype e, site s, address a, country c where e.siteid=s.id and s.addressid=a.id and a.countryid=c.id;


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

--2008-02-12 table for processed/averaged phenotype
use at;
create table phenotype_p(
	id  integer auto_increment primary key,
	accession_id  integer  unique,
	region  varchar(200),
	LD  float,
	LD_sample_size  integer,
	LD_stdev  float,
	avrPph3  integer,
	avrRpm1  integer,
	avrRpt2  integer,
	avrB  integer,
	FRI_1Ler_2Col  varchar(1),
	Rps5  varchar(1),
	Rpm1  varchar(1),
	Rps2  varchar(1),
	LDV  float,
	SD  float,
	SDV  float,
	JIC0W  float,
	JIC2W  float,
	JIC4W  float,
	JIC8W  float,
	FLC  float,
	FRI  float);

--2008-04-14 add auto-updated created, modified, created_by, modified_by columns to at.phenotype
use at;

alter table phenotype add created_by varchar(200)  AFTER modified;

alter table phenotype add modified_by varchar(200) AFTER created_by;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_phenotype BEFORE INSERT ON phenotype
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
        if NEW.created is null then
               set NEW.created = CURRENT_TIMESTAMP();
        end if;
  END;
|

CREATE TRIGGER before_update_phentoype BEFORE UPDATE ON phenotype
  FOR EACH ROW BEGIN
        if NEW.modified_by is null then
                set NEW.modified_by = USER();
        end if;
        if NEW.modified=0 then
                set NEW.modified = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;


--2008-02-18 database for 250k
use stock_250k;
SET storage_engine=INNODB;

create table snps(
	id integer auto_increment primary key,
	name varchar(200) not null unique,
	chromosome integer,
	position integer,
	allele1 varchar(1),
	allele2 varchar(2),
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_snps BEFORE INSERT ON snps
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_snps BEFORE UPDATE ON snps
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|
DELIMITER ;

--2008-02-25 expand it to include tiling probe info
create table probes(
	id integer auto_increment primary key,
	snps_id integer,
	foreign key (snps_id) references snps(id) on delete cascade on update cascade,
	seq varchar(25),
	chromosome integer,
	position integer,
	allele varchar(1),
	strand varchar(20),
	xpos integer,
	ypos integer,
	direction varchar(20),
	gene varchar(50),
	RNA varchar(50),
	tu varchar(50),
	flank varchar(50),
	expressedClones float,
	totalClones integer,
	multiTranscript varchar(50),
	LerDel varchar(50),
	LerCopy integer,
	LerSNPdelL integer,
	LerSNPdelR integer,
	LerSNPpos integer,
	promoter BOOL,
	utr5 BOOL,
	utr3 BOOL,
	intron BOOL,
	intergenic BOOL,
	downstream BOOL,
	cda BOOL,
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_probes BEFORE INSERT ON probes
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_probes BEFORE UPDATE ON probes
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|
DELIMITER ;


create table strain_info(
	id integer auto_increment primary key,
	name varchar(40),
	description varchar(2000),
	ecotype_id integer,
	maternal_ecotype_id integer,
	paternal_ecotype_id integer,
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_strain_info BEFORE INSERT ON strain_info
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_strain_info BEFORE UPDATE ON strain_info
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

create table array_info(
	id integer auto_increment primary key,
	name varchar(40),
	filename varchar(1000),
	original_filename varchar(1000),
	description varchar(2000),
	ecotype_id integer,
	maternal_ecotype_id integer,
	paternal_ecotype_id integer,
	strain_id integer,
	md5sum varchar(100),
	experimenter varchar(200),
	samples varchar(20),
	dna_amount varchar(20),
	S260_280 float,
	total_vol varchar(20),
	hyb_vol varchar(20),
	seed_source varchar(100),
	method_name varchar(250),
	readme_id integer,
	foreign key (readme_id) references readme(id) on update cascade,
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP,
	CONSTRAINT array_info_strain_id_fk_constraint foreign key (strain_id) references strain_info(id) on delete cascade on update cascade
	)engine=INNODB;



DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_array_info BEFORE INSERT ON array_info
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_array_info BEFORE UPDATE ON array_info
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

create table array_data(
	id integer auto_increment primary key,
	array_id integer,
	-- foreign key (array_id) references array_info(id) on delete cascade on update cascade,
	probes_id integer,
	-- foreign key (probes_id) references probes(id) on delete cascade on update cascade,
	intensity float,
	--when probes_id is null, xpos and ypos will be used.
	xpos integer,
	ypos integer
	)engine=INNODB;

create table call_method(
	id integer auto_increment primary key,
	short_name varchar(20),
	method_description varchar(8000),
	data_description varchar(8000),
	comment varchar(8000),
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_c_method BEFORE INSERT ON call_method
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_c_method BEFORE UPDATE ON call_method
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

create table call_info(
	id integer auto_increment primary key,
	filename varchar(1000),
	description varchar(2000),
	array_id integer,
	NA_rate float,
	readme_id integer,
	foreign key (readme_id) references readme(id) on update cascade,
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP,
	foreign key (array_id) references array_info(id) on delete cascade on update cascade,
	method_id integer,
	foreign key (method_id) references call_method(id) on delete cascade on update cascade
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_call_info BEFORE INSERT ON call_info
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_call_info BEFORE UPDATE ON call_info
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

create table call_data(
	id integer auto_increment primary key,
	call_info_id integer,
	-- foreign key (call_info_id) references call_info(id) on delete cascade on update cascade,
	snps_id integer,
	-- foreign key (snps_id) references snps(id) on delete cascade on update cascade,
	snpcall varchar(2)
	)engine=INNODB;


create table if not exists qc_method(
	id integer auto_increment primary key,
	short_name varchar(30) unique,
	data1_type varchar(30) not NULL,
	data2_type varchar(30) not NULL,
	method_description varchar(8000),
	data_description varchar(8000),
	comment varchar(8000),
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_qc_method BEFORE INSERT ON qc_method
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_qc_method BEFORE UPDATE ON qc_method
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

create table if not exists call_qc(
	id integer auto_increment primary key,
	call_info_id integer,
	min_probability float,
	ecotype_id integer,
	duplicate integer,
	tg_ecotype_id integer,
	tg_duplicate integer,
	NA_rate float,
	no_of_NAs integer,
	no_of_totals integer,
	mismatch_rate float,
	no_of_mismatches integer,
	no_of_non_NA_pairs integer,
	created_by varchar(200),
	updated_by varchar(200),
	call_method_id integer,
	foreign key (call_method_id) references call_method(id) on delete cascade on update cascade,
	readme_id integer,
	foreign key (readme_id) references readme(id) on delete cascade on update cascade,
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP,
	foreign key (call_info_id) references call_info(id) on delete cascade on update cascade,
	qc_method_id integer,
	foreign key (qc_method_id) references qc_method(id) on delete cascade on update cascade
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_call_qc BEFORE INSERT ON call_qc
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_call_qc BEFORE UPDATE ON call_qc
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

create table if not exists snps_qc(
	id integer auto_increment primary key,
	snps_id integer,
	min_probability float,
	max_call_info_mismatch_rate float,
	snps_name varchar(200),
	tg_snps_name varchar(200),
	NA_rate float,
	no_of_NAs integer,
	no_of_totals integer,
	relative_NA_rate float,
	relative_no_of_NAs integer,
	relative_no_of_totals integer,
	mismatch_rate float,
	no_of_mismatches integer,
	no_of_non_NA_pairs integer,
	foreign key (snps_id) references snps(id) on delete cascade on update cascade,
	qc_method_id integer,
	foreign key (qc_method_id) references qc_method(id) on delete cascade on update cascade,
	call_method_id integer,
	foreign key (call_method_id) references call_method(id) on delete cascade on update cascade,
	readme_id integer,
	foreign key (readme_id) references readme(id) on delete cascade on update cascade,
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_snps_qc BEFORE INSERT ON snps_qc
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_snps_qc BEFORE UPDATE ON snps_qc
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;


--store the method
create table results_method(
	id integer auto_increment primary key,
	short_name varchar(30) unique,
	method_description varchar(8000),
	data_description varchar(8000),
	phenotype_method_id integer,
	call_method_id integer,
	comment varchar(8000),
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP,
	foreign key (phenotype_method_id) references phenotype_method(id) on delete RESTRICT on update cascade,
	foreign key (call_method_id) references call_method(id) on delete RESTRICT on update cascade
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_r_method BEFORE INSERT ON results_method
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_r_method BEFORE UPDATE ON results_method
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

--store the results
create table if not exists results(
	id integer auto_increment primary key,
	chr integer,
	start_pos integer,
	stop_pos integer,
	results_method_id integer,
	foreign key (results_method_id) references results_method(id) on delete cascade on update cascade,
	score float
	)engine=INNODB;

create table phenotype_method(
	id integer auto_increment primary key,
	short_name varchar(20),
	method_description varchar(8000),
	data_description varchar(8000),
	comment varchar(8000),
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_p_method BEFORE INSERT ON phenotype_method
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_p_method BEFORE UPDATE ON phenotype_method
  FOR EACH ROW BEGIN
        if NEW.updated_by is null then
                set NEW.updated_by = USER();
        end if;
        if NEW.date_updated=0 then
                set NEW.date_updated = CURRENT_TIMESTAMP();
        end if;
  END;
|

DELIMITER ;

--store the phenotype
create table phenotype(
	id integer auto_increment primary key,
	ecotype_id integer not null,
	value float,
	replicate integer,
	method_id integer not null,
	foreign key (method_id) references phenotype_method(id) on delete cascade on update cascade
	)engine=INNODB;

create table phenotype_avg(
	id integer auto_increment primary key,
	ecotype_id integer not null,
	value float,
	stdev float,
	sample_size integer,
	method_id integer not null,
	foreign key (method_id) references phenotype_method(id) on delete cascade on update cascade,
	readme_id integer,
	foreign key (readme_id) references readme(id) on delete cascade on update cascade
	)engine=INNODB;


--2008-04-30 create a view to view qc results for arrays
create or replace view view_qc as select e.id as ecotype_id, e.nativename, q.tg_ecotype_id, q.call_info_id, q.call_method_id, a.id as array_id, a.original_filename as array_original_filename, a.date_created as array_created, c.NA_rate, group_concat(q.mismatch_rate order by q.qc_method_id) as mismatch_rate, group_concat(q.qc_method_id order by q.qc_method_id) as qc_method, group_concat(q.no_of_non_NA_pairs order by q.qc_method_id) as no_of_non_NA_pairs from call_info c, array_info a, call_qc q , stock.ecotype e where e.id=q.ecotype_id and q.call_info_id = c.id and a.id=c.array_id group by q.ecotype_id, q.call_info_id, q.call_method_id order by nativename,  call_method_id, NA_rate,array_created, original_filename;

--2008-05-20 view the calls, arrays, ecotypes all together

create or replace view view_call as select c.id as call_info_id, c.filename, c.method_id as call_method_id, a.id as array_id, a.original_filename, a.maternal_ecotype_id as ecotype_id, e.nativename, e.stockparent from call_info c, array_info a, stock.ecotype e where a.id=c.array_id and e.id=a.maternal_ecotype_id order by nativename;