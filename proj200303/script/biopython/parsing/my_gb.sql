create table db_genbank (
locus varchar,
gb_acc varchar,
gi varchar,
organism varchar,
sp_acc varchar
);

create table db_cds (
gene varchar,
organism varchar,
seq text,
protein_id varchar,
gb_acc varchar,
sp_acc varchar
);


create table db_protein (
protein_id varchar,
organism varchar,
seq varchar,
gene varchar,
gb_acc varchar,
sp_acc varchar
);

create table db_location (
gene varchar,
strand varchar,
start integer,
tail integer,
gb_acc varchar
);


create table db_sp (
fasta varchar,
acc varchar,
seq varchar
);
