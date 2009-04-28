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
create or replace view complete_2010_strains_in_stock2tg_ecotypeid as 
	select distinct e.tg_ecotypeid, c.*, a.origin, a.number from complete_2010_strains_in_stock c,
	stock.ecotypeid2tg_ecotypeid e , accession a where c.ecotypeid=e.ecotypeid and c.accession_id=a.id order by accession_id;


--2008-05-18 offering the final linking between accession id and ecotype id. a view linking each accession.id to a stock.ecotypeid2tg_ecotypeid.tg_ecotypeid thru at.ecotype2accession.
create or replace view accession2tg_ecotypeid as select distinct e1.accession_id,
	a.name as accession_name, a.origin, a.number, e2.tg_ecotypeid as ecotype_id, e.name,
	e.nativename, e.stockparent, e1.ecotype_id as intermediate_ecotype_id from ecotype2accession e1,
	stock.ecotype e, stock.ecotypeid2tg_ecotypeid e2, accession a where e1.accession_id=a.id and
	e1.ecotype_id=e2.ecotypeid and e2.tg_ecotypeid=e.id order by accession_id;

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



--2007-10-21 for database stock or its replicates
use stock;
create table person(
	id integer primary key auto_increment,
	title varchar(4),
	surname varchar(40) not null,
	firstname varchar(40) not null,
	email varchar(100),
	username varchar(10),
	password varchar(15),
	usertype integer,
	donor integer not null
	)engine=INNODB;

-- 2009-3-30 add column tg_ecotypeid to view stock.ecotype_info
-- 2008-05-18 view to check other info (where and country) about ecotype
create or replace view ecotype_info as select e.id as ecotypeid, e.name, e.stockparent,
	e.nativename, e2t.tg_ecotypeid, e.latitude, e.longitude, e.alias, e.siteid, s.name as site_name, 
	c.abbr as country, p.firstname, p.surname, e.collectiondate, e.geographic_integrity_id from ecotype e,
	site s, address a, country c , person p, ecotypeid_strainid2tg_ecotypeid e2t 
	where p.id=e.collectorid and e.siteid=s.id and s.addressid=a.id and a.countryid=c.id
	and e.id=e2t.ecotypeid;

-- 2009-4-28 split table geographic_integrity out of ecotype_info and do left outer join to avoid no-geographic integrity rows being left out
create or replace view ecotype_info_with_gps_integrity as select e.*, g.short_name as geographic_integrity 
	from ecotype_info e left outer join
	geographic_integrity g on e.geographic_integrity_id=g.id;

-- 2009-3-30 table site_trip is not maintained. so useless. 
-- view doesn't allow selection clause in FROM, so construct more views based on ecotype_info to generate more info
create or replace view ecotype_info_with_trip as select e.*, s2t.tripid from
	ecotype_info e left outer join site_trip s2t on e.siteid=s2t.siteid;

-- 2009-3-30 table trip is not maintained. so this view is useless.
create or replace view ecotype_info_with_collectiondate as select e.*, t.collectiondate from
	ecotype_info_with_trip e left outer join trip t on e.tripid=t.id;

-- 2009-4-5 combine the haplogroup info here
create or replace view ecotype_info_with_haplogroup as select e.*, h.haplo_group_id from
	ecotype_info_with_gps_integrity e left outer join haplo_group2ecotype h on e.ecotypeid=h.ecotype_id;

-- 2008-08-08 manual tg_ecotypeid linking after GroupDuplicateEcotype.py http://papaya.usc.edu/2010/149-snps/149SNP-data-introduction
-- for Kas-1 & Kas-2
update ecotypeid_strainid2tg_ecotypeid set tg_ecotypeid=7183 where ecotypeid in (7185, 7183);
update ecotypeid_strainid2tg_ecotypeid set tg_ecotypeid=8424 where ecotypeid in (6925, 7184, 8315, 8424);

-- 2008-08-18 create a view to view qc results for db stock
create or replace view view_qc as select q.strainid, e.id as ecotype_id, e.nativename, q.target_id, 
	q.qc_method_id, qm.short_name as QC_method_name, q.NA_rate as QC_NA_rate, 
	q.mismatch_rate , q.no_of_mismatches, q.no_of_non_NA_pairs, q.created_by, q.updated_by, q.date_created, q.date_updated
	from call_qc q , ecotype e, qc_method qm where e.id=q.ecotypeid 
	and qm.id=q.qc_method_id order by nativename, strainid, qc_method_id;


-- 2008-08-11 don't need these anymore. GroupDuplicateEcotype.py can create these on the fly
-- 2008-08-11 turns out that the type of id in table ecotype is 'integer unsigned'. after modifying those columns to be unsigned, succeed.
-- 2008-05-15 attempt to add foreign key constraints for tables created by GroupDuplicateEcotype.py but all failed. get "ERROR 1005 (HY000): Can't create table './stock/#sql-99_17.frm' (errno: 150)"
alter table ecotype_duplicate2tg_ecotypeid add foreign key (ecotypeid) references ecotype(id) on delete restrict on update cascade;
alter table ecotype_duplicate2tg_ecotypeid add foreign key (tg_ecotypeid) references ecotype(id) on delete restrict on update cascade;
alter table genotyping_all_na_ecotype add foreign key (ecotypeid) references ecotype(id) on delete restrict on update cascade;
alter table nativename_stkparent2tg_ecotypeid add foreign key (tg_ecotypeid) references ecotype(id) on delete restrict on update cascade;


-- 2008-12-17 this part is outdated and defunct. check plone doc, /log/sql/, for updates.
-- 2008-08-11 two new tables to split calls_byseq, in StockDB.py, but due to this f*** unsigned, have to manually create them
CREATE TABLE strain (
	id INTEGER NOT NULL AUTO_INCREMENT primary key, 
	ecotypeid INTEGER unsigned, 
	extractionid tinyINT unsigned, 
	seqinfoid smallINT unsigned, 
	plateid VARCHAR(25), 
	wellid VARCHAR(3), 
	replicate BOOL,
	contaminant_type_id integer,
	created_by VARCHAR(128), 
	updated_by VARCHAR(128), 
	date_created DATETIME, 
	date_updated DATETIME, 
	UNIQUE (ecotypeid, plateid),
	 CONSTRAINT strain_ecotypeid_fk FOREIGN KEY(ecotypeid) REFERENCES ecotype (id) ON DELETE CASCADE ON UPDATE CASCADE, 
	 CONSTRAINT strain_extractionid_fk FOREIGN KEY(extractionid) REFERENCES extraction (id) ON DELETE CASCADE ON UPDATE CASCADE, 
	 CONSTRAINT strain_seqinfoid_fk FOREIGN KEY(seqinfoid) REFERENCES seqinfo (id) ON DELETE CASCADE ON UPDATE CASCADE,
	 FOREIGN KEY(contaminant_type_id) REFERENCES contaminant_type (id) ON DELETE restrict ON UPDATE CASCADE
)ENGINE=InnoDB;

CREATE TABLE calls (
	id INTEGER NOT NULL AUTO_INCREMENT primary key, 
	strainid INTEGER, 
	snpid INTEGER unsigned, 
	allele VARCHAR(5),	--#'call' is mysql reserved keyword
	 CONSTRAINT calls_strainid_fk FOREIGN KEY(strainid) REFERENCES strain (id) ON DELETE CASCADE ON UPDATE CASCADE, 
	 CONSTRAINT calls_snpid_fk FOREIGN KEY(snpid) REFERENCES snps (id) ON DELETE CASCADE ON UPDATE CASCADE
)ENGINE=InnoDB;



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


--2008-06-02 table to store selfing rates
create table popid2s_100(id integer primary key auto_increment,
	popid integer not null,
	avg_s_Jarne2006 float,
	std_s_Jarne2006 float,
	avg_s_Robertson1984 float,
	std_s_Robertson1984 float,
	avg_s_Weir1984 float,
	std_s_Weir1984 float,
	weir1984_multi_loci_s float,
	avg_s_Nordborg1997 float,
	std_s_Nordborg1997 float,
	s_g2_David2007 float)engine=INNODB;


--2008-02-18 database for 250k
use stock_250k;
SET storage_engine=INNODB;

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

create table call_data(
	id integer auto_increment primary key,
	call_info_id integer,
	-- foreign key (call_info_id) references call_info(id) on delete cascade on update cascade,
	snps_id integer,
	-- foreign key (snps_id) references snps(id) on delete cascade on update cascade,
	snpcall varchar(2)
	)engine=INNODB;

-- 2008-05-31 add two indices
create index results_snps_id_results_method_id_idx on results(snps_id, results_method_id);
create index results_results_method_id_idx on results(results_method_id);

--2008-04-30 create a view to view qc results for arrays
create or replace view view_qc as select e.id as ecotype_id, e.nativename, \
	a.id as array_id, q.tg_ecotype_id, q.call_info_id, 
	q.call_method_id, c.NA_rate as call_NA_rate,  q.qc_method_id, qm.short_name as QC_method_name, q.NA_rate as QC_NA_rate, 
	q.mismatch_rate , q.no_of_mismatches, q.no_of_non_NA_pairs, a.median_intensity, 
	a.original_filename as array_original_filename, a.date_created as array_created, q.date_created as qc_date
	from call_info c, array_info a, call_qc q , stock.ecotype e, qc_method qm where \
	e.id=q.ecotype_id 
	and q.call_info_id = c.id and a.id=c.array_id and qm.id=q.qc_method_id order by nativename, array_id, call_method_id,
	array_created, qc_method_id;

--2008-05-20 view the calls, arrays, ecotypes all together

create or replace view view_call as select c.id as call_info_id, c.filename, 
	c.method_id as call_method_id, a.id as array_id, a.original_filename,
	a.maternal_ecotype_id as ecotype_id, e.latitude, e.longitude, e.nativename, e.stockparent, s.name as site, ad.region, co.abbr as country \
	from call_info c, array_info a, stock.ecotype e, stock.site s, stock.address ad, stock.country co where \
	a.id=c.array_id and e.id=a.maternal_ecotype_id and e.siteid=s.id and s.addressid=ad.id and ad.countryid=co.id \
	order by nativename;

-- 2008-05-27 view the arrays
-- 2009-4-5 add a date array_created to view stock_250k.view_array
create or replace view view_array as select a.id as array_id, a.filename, 
	a.original_filename as array_filename,  a.maternal_ecotype_id,
	e1.nativename as maternal_nativename, e1.stockparent as maternal_stockparent,
	a.paternal_ecotype_id, e2.nativename as paternal_nativename,
	e2.stockparent as paternal_stockparent, a.date_created as array_created from array_info a, stock.ecotype e1,
	stock.ecotype e2 where e1.id=a.maternal_ecotype_id and e2.id=a.paternal_ecotype_id
	order by maternal_nativename, paternal_nativename;

--2008-07-18 view the rank sum test results
create or replace view view_rank_sum_test as select cgr.results_method_id,
	r.short_name as results_short_name, a.short_name as analysis_short_name, p.short_name as pheno_short_name,
	r.call_method_id, group_concat(cgr.pvalue order by cgr.list_type_id),
	group_concat(cgr.list_type_id order by cgr.list_type_id)
	from candidate_gene_rank_sum_test_result  cgr, results_method r,
	phenotype_method p, analysis_method a where cgr.results_method_id=r.id
	and r.phenotype_method_id=p.id and r.analysis_method_id=a.id group by cgr.results_method_id
	order by pheno_short_name, analysis_short_name;

-- 2008-08-21 view the different parameter setting of top SNP HG test
create or replace view view_top_snp_test_param as select distinct c.min_distance, c.get_closest, c.min_MAF, 
c.no_of_top_snps, r.call_method_id from results_method r, results_by_gene rbg, candidate_gene_top_snp_test c where r.id=rbg.results_method_id and rbg.id=c.results_by_gene_id;
 
-- 2008-08-21 view the different parameter setting of rank sum test
create or replace view view_rank_test_param as select distinct c.min_distance, c.get_closest, c.min_MAF, 
c.max_pvalue_per_gene, r.call_method_id from results_method r, results_by_gene rbg, candidate_gene_rank_sum_test_result c where r.id=rbg.results_method_id and rbg.id=c.results_by_gene_id;
