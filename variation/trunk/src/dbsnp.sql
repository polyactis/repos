-- 2008-05-02 schema for generic SNP database. first is to hold 384-loci illumina-multiplex data for 93 individuals

create database if not exists dbsnp;
use dbsnp;
SET storage_engine=INNODB;

create table snpset(
	id integer auto_increment primary key,
	acc varchar(200),
	description varchar(4000),
	created_by varchar(200),
	updated_by varchar(200),
	date_created timestamp default CURRENT_TIMESTAMP,
	date_updated TIMESTAMP
	)engine=INNODB;

DELIMITER |     -- change the delimiter ';' to '|' because ';' is used as part of one statement.

CREATE TRIGGER before_insert_snpset BEFORE INSERT ON snpset
  FOR EACH ROW BEGIN
        if NEW.created_by is null then
               set NEW.created_by = USER();
        end if;
  END;
|

CREATE TRIGGER before_update_snpset BEFORE UPDATE ON snpset
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

create table if not exists snps(
	id integer auto_increment primary key,
	acc varchar(200),
	chromosome integer,
	position integer,
	probe_sequence varchar(4000),
	allele1 varchar(2),
	allele2 varchar(2),
	allele3 varchar(2),
	allele4 varchar(2),
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

create table snps2snpset(
	id integer auto_increment primary key,
	snps_id integer,
	snpset_id integer,
	foreign key (snps_id) references snps(id) on delete cascade on update cascade,
	foreign key (snpset_id) references snpset(id) on delete cascade on update cascade
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

create table calls(
	id integer auto_increment primary key,
	ecotype_id integer,
	snps_id integer,
	genotype varchar(2),
	call_method_id integer,
	foreign key (snps_id) references snps(id) on delete cascade on update cascade,
	foreign key (call_method_id) references call_method(id) on delete cascade on update cascade
	)engine=INNODB;


create table readme(
	id      integer auto_increment primary key,
	title    varchar(2000),
	description     varchar(60000),
	created_by      varchar(200),
	updated_by     varchar(200),
	date_created    timestamp default CURRENT_TIMESTAMP,
	date_updated   TIMESTAMP default 0
	);

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


