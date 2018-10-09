#ifndef __VECTOR_H
#define __VECTOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct vector_ {
	void **data;
	int size;
	int count;
	int valid;
} Vector;

#include "grouping/grouping.h"

int vector_get_valid(Vector *v);

void vector_set_valid(Vector *v, int valid);

void vector_init(Vector*);

int vector_count(Vector*);

void vector_add(Vector*, void*);

void vector_set(Vector*, int, void*);

void *vector_get(Vector*, int);

void vector_delete(Vector*, int);

void vector_free(Vector*);

void chain_free(Vector **chains, int size);

#endif /* __VECTOR_H */
