/////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation; either version 2.1 of the License,
// or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
/////////////////////////////////////////////////////////////////

#ifndef HASHL_HEADER
#define HASHL_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HASH_SIZE 270000

/*
** Structure for linked list of name + int value pairs.
*/
struct nlist {
  struct nlist *next;
  char *name;
  unsigned long value;
};

unsigned hash(char *s);
int hashtab_free(struct nlist *hashtab[]);
int hashtab_init(struct nlist *hashtab[], unsigned long n);
struct nlist *hashtab_install(char *name, int value, struct nlist *hashtab[]);
struct nlist *hashtab_lookup(char *s, struct nlist *hashtab[]);
int nlist_free(struct nlist *np);

#endif
