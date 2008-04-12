/* This file contains table definitions for the optilux
  cinemas example.
 */

create database if not exists optilux;
use optilux;

-- Screenings
create table if not exists screening (
    screening_id integer unsigned not null auto_increment primary key,
    cinema_code char(4) not null,
    film_code char(4) not null,
    show_time datetime not null,
    remaining_tickets integer unsigned not null,
    index showing_cinema_code(cinema_code),
    index showing_film_code(film_code),
    index showing_show_time(show_time),
    index showing_remaining_tickets(remaining_tickets)
) engine=InnoDB;

-- Reservations
create table if not exists reservation (
    reservation_id integer unsigned not null auto_increment primary key,
    screening_id integer unsigned not null,
    num_tickets tinyint unsigned not null,
    customer_name varchar(64) not null,
    index reservation_num_tickets(num_tickets),
    foreign key(screening_id)
        references screening(screening_id)
            on update restrict
            on delete restrict
) engine=InnoDB;

-- Now set up an identical test database

create database if not exists optilux_test;
use optilux_test;

-- Screenings
create table if not exists screening (
    screening_id integer unsigned not null auto_increment primary key,
    cinema_code char(4) not null,
    film_code char(4) not null,
    show_time datetime not null,
    remaining_tickets integer unsigned not null,
    index showing_cinema_code(cinema_code),
    index showing_film_code(film_code),
    index showing_show_time(show_time),
    index showing_remaining_tickets(remaining_tickets)
) engine=InnoDB;

-- Reservations
create table if not exists reservation (
    reservation_id integer unsigned not null auto_increment primary key,
    screening_id integer unsigned not null,
    num_tickets tinyint unsigned not null,
    customer_name varchar(64) not null,
    index reservation_num_tickets(num_tickets),
    foreign key(screening_id)
        references screening(screening_id)
            on update restrict
            on delete restrict
) engine=InnoDB;