create schema retreat;

set search_path to retreat;

create table register(
	id	serial primary key,
	name	varchar,
	email	varchar,
	pi	varchar,
	tshirt	varchar,
	roommate1	varchar,
	roommate2	varchar,
	roommate3	varchar,
	roommate4	varchar);

create table abstract(
	id	serial primary key,
	name	varchar,
	email	varchar,
	pi	varchar,
	pref	varchar,
	title	varchar,
	author_list	varchar,
	abstract	varchar);
