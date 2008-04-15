create schema retreat;

set search_path to retreat;

create table register(
	id	serial primary key,
	name	varchar,
	firstname	varchar,
	lastname	varchar,
	email	varchar,
	pi	varchar,
	special_treat	varchar,
	nights	varchar,
	tshirt	varchar,
	roommate1	varchar,
	roommate2	varchar,
	roommate3	varchar,
	roommate4	varchar,
	date_created	timestamp default current_timestamp);

create table abstract(
	id	serial primary key,
	name	varchar,
	email	varchar,
	pi	varchar,
	pref	varchar,
	title	varchar,
	author_list	varchar,
	abstract	varchar,
	abstract_file_data text,
	abstract_file_content_type text,
	abstract_filename text,
	abstract_file_size integer,
	date_created	timestamp default current_timestamp);

create table carpool(
	id	serial primary key,
	name	varchar,
	email	varchar,
	phone	varchar,
	address	varchar,
	ride_type varchar,
	no_of_people integer,
	severity varchar,
	date_created timestamp default current_timestamp);