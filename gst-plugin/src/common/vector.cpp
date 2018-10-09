#include "vector.h"

/* */ 
void vector_init (Vector *v) {
   v->data = NULL;
   v->size = 0;
   v->count = 0;
   v->valid = 0;
}
 
/* */ 
int vector_count (Vector *v) {
   return v->count;
}

int vector_get_valid (Vector *v) {
   return v->valid;
}

void vector_set_valid (Vector *v, int valid) {
   v->valid = valid;
}

void vector_add(Vector *v, void *e)
{
	if (v->size == 0) {
		v->size = 10;
		v->data = (void **) malloc(sizeof(void *) * v->size);
		memset(v->data, '\0', sizeof(void*) * v->size);
	}
	/* Condition to increase v->data, last slot exhausted: */
	if (v->size == v->count) {
		v->size *= 2;
		v->data = (void **) realloc(v->data, sizeof(void *) * v->size);
	}
	v->data[v->count] = e;
	v->count++;
}
 
void vector_set (Vector *v, int index, void *e) {
   if (index >= v->count) {
      return;
   }
   v->data[index] = e;
}
 
void *vector_get(Vector *v, int index)
{
	if (index >= v->count) {
		/* FIXME: we must handle this. */
		//return;
		//return NULL;
	}

	return v->data[index];
}
 
void vector_delete (Vector *v, int index) {
   
   if (index >= v->count) {
      return;
   }
 
   v->data[index] = NULL;
 
   int i, j;
   void **newarr = (void**)malloc(sizeof(void*) * v->count * 2);
   for (i = 0, j = 0; i < v->count; i++) {
      if (v->data[i] != NULL) {
         newarr[j] = v->data[i];
         j++;
      }       
   }
   free(v->data);
   v->data = newarr;
   v->count--;
}
 
void vector_free(Vector *v)
{
	free(v->data);
}

void chain_free (Vector **chains, int size) {
   int i;
   for (i = 0; i < size; i++) {
      int j;
      for (j = 0; j < vector_count(chains[i]); j++) {
          region *r = (region *)vector_get(chains[i], j);
          vector_free(&(r->xp));
          vector_free(&(r->yp));
          free(r);
      }
      vector_free(chains[i]);
      free(chains[i]);
   }
   free(chains);
}

