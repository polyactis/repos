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

#include "hash.h"

char *strdup(char *s);

/* form hash value from string */
unsigned hash(char *s) {
  unsigned hashval;

  for(hashval=0; *s != '\0'; s++)
    hashval = *s + 31 * hashval;

  return(hashval % HASH_SIZE);
}

/* hashtab_free: free entries in hashtab */
int hashtab_free(struct nlist *hashtab[]) {
  unsigned long i;

  for(i=0; i<HASH_SIZE; i++)
    nlist_free(hashtab[i]);

  return(EXIT_SUCCESS);
}

/* hashtab_init: initialize hashtab */
int hashtab_init(struct nlist *hashtab[], unsigned long n) {
  unsigned long i;

  for(i=0; i<n; i++)
    hashtab[i] = NULL;

  return(EXIT_SUCCESS);
}

/* hashtab_install: put(name,value) in hashtab */
struct nlist *hashtab_install(char *name, int value, struct nlist *hashtab[]) {
  unsigned hashval;
  struct nlist *np;

  if((np = hashtab_lookup(name,hashtab))==NULL) {
    /* name not found, add to head of list */
    np = (struct nlist *) malloc(sizeof(struct nlist));
    if(np == NULL || (np->name = strdup(name)) == NULL)
      return(NULL);
    hashval = hash(name);
    np->next = hashtab[hashval];
    hashtab[hashval] = np;
  }

  np->value = value;

  return(np);
}

/* hashtab_lookup: look for s in hashtab */
struct nlist *hashtab_lookup(char *s, struct nlist *hashtab[]) {
  struct nlist *np;

  for(np=hashtab[hash(s)]; np != NULL; np=np->next)
    if(!strcmp(s,np->name))
      return(np);

  return(NULL);
}

/* nlist_free: free entries in linked list */
int nlist_free(struct nlist *np) {

  if(np == NULL)
    return(EXIT_SUCCESS);

  if(np->next != NULL)
    nlist_free(np->next);

  if(np->name != NULL)
    free(np->name);

  return(EXIT_SUCCESS);
}

char *strdup(char *s) {
  char *c;

  if((c = (char *) malloc((1+strlen(s)) * sizeof(char)))==NULL)
    return(NULL);

  if(strcpy(c,s)==NULL)
    return(NULL);

  return(c);
}
