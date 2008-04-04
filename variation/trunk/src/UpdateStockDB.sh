#!/bin/sh
##########
#2008-04-03 a script to restore uchicago's stock db on papaya
#
#####################
db_dump_fname=/tmp/uchicago.mysql.stock.db
echo mysqldump -h natural.uchicago.edu -u iamhere -p stock --ignore-table=stock.person --skip-lock-tables \> $db_dump_fname
mysqldump -h natural.uchicago.edu -u iamhere -p stock --ignore-table=stock.person --skip-lock-tables > $db_dump_fname

echo mysql -h papaya -u yh -p -e "create database stock"
mysql -h papaya -u yh -p -e "create database stock"

echo mysql -h papaya -u yh -p stock \<$db_dump_fname
#mysql -h papaya -u yh -p stock <$db_dump_fname

person_data_fname=$db_dump_fname.person.table

echo mysql -h natural.uchicago.edu -u iamhere -p stock -e "select id, title, SURNAME, FIRSTNAME, EMAIL, DONOR from person" \> $person_data_fname
mysql -h natural.uchicago.edu -u iamhere -p stock -e "select id, title, SURNAME, FIRSTNAME, EMAIL, DONOR from person" > $person_data_fname

no_of_lines=`wc $person_data_fname|awk '{print $1}'`
no_of_lines_need=`echo $no_of_lines-1|bc`

new_person_data_fname=$person_data_fname.new
tail -n $no_of_lines_need $person_data_fname > $new_person_data_fname

echo create table person on papaya, db stock.
mysql -h papaya -u yh -p stock -e "create table person(id integer primary key auto_increment,\
	title varchar(4),\
	surname varchar(40) not null,\
	firstname varchar(40) not null,\
	email varchar(100),\
	donor integer not null\
	)"

echo copy $new_person_data_fname to papaya to ensure \"load data infile\" could find it.
scp $new_person_data_fname yh@papaya:/tmp/

echo load person data into stock.person on papaya.
mysql -h papaya -u yh -p stock -e "load data infile '$new_person_data_fname' into table stock.person"

