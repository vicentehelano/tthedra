/*
 * Copyright (c) 2010 TTHedra project and any individual authors listed
 * elsewhere in this file.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License version 3 for more details.
 *
 * You should have received a copy of the GNU General Public License
 * version 3 along with this package (see COPYING file).
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * $URL$
 * $Id$
 *
 * Author(s): Vicente H. F. Batista
 */
#include <stdlib.h>
#include <assert.h>
#include <heap.h>

Heap *front_heap;
Heap *reject_heap;

void HEAP_Initialize_front_heap(Surface_mesh *front, Volume_mesh *domain)
{
  /* Definitions of local variables */
  double heapkey;
  unsigned long i;

  /* Allocate a block of memory to the main heap structure */
  front_heap = (Heap *) malloc(sizeof(Heap)); 
  assert(front_heap != NULL);     /*  Check if memory was allocated   */

  front_heap->nnds = 0;
  /* Calculate initial maximum number of nodes */
  i = front->nels;
  front_heap->maxnnds = 2;
  while ((i /= 2) != 0) {
    front_heap->maxnnds *= 2;
  }
  front_heap->maxnnds *= 2;
  front_heap->maxnnds--;

  /* Allocate the heap array (nodes + SENTINEL) */
  front_heap->node = (Heap_node *)
    malloc((front_heap->maxnnds + 1)*sizeof(Heap_node));
  assert(front_heap->node != NULL);   /* Check if memory was allocated    */

  /* Initialize the sentinel node */
  front_heap->node[0].id  = 0;
  front_heap->node[0].key = (double) SENTINEL;

  /* Insert the boundary faces into the heap  */
  for (i = 1; i <= front->nels; i++) {
    heapkey = HEAP_Key(&(front->element[i-1]), domain);
    HEAP_Insert(i, heapkey, front_heap);   /* Initial front of faces  */    
  }
}

double HEAP_Key(Triangle *t, Volume_mesh *domain)
{
  double sqrdperim;
  Point3D *v1, *v2, *v3;
  Point3D vec12, vec23, vec31;

  v1 = &(domain->node[t->vertex[0] - 1]);
  v2 = &(domain->node[t->vertex[1] - 1]);
  v3 = &(domain->node[t->vertex[2] - 1]);

  GEOM_SUB3D(vec12,(*v2),(*v1))
  GEOM_SUB3D(vec23,(*v3),(*v2))
  GEOM_SUB3D(vec31,(*v1),(*v3))

  sqrdperim = 0.0e+000;
  sqrdperim += GEOM_DOT3D(vec12,vec12);
  sqrdperim += GEOM_DOT3D(vec23,vec23);
  sqrdperim += GEOM_DOT3D(vec31,vec31);

  return (sqrdperim);
}

void HEAP_Insert(unsigned long  triangleid,  double trianglekey,
         Heap     *heap_ptr)
{
  /* Increase by one the heap number of nodes */
  (heap_ptr->nnds)++;

  /* Check if the available memory block if sufficient */
  HEAP_Check_memory(heap_ptr);

  /* Insert the new node at the last position of the heap array */
  heap_ptr->node[heap_ptr->nnds].id = triangleid;

/*  newtri = &(front->element[triangleid - 1]);*/

  heap_ptr->node[heap_ptr->nnds].key  = trianglekey;

  /* Traverse the heap and restore the heap condition */
  HEAP_Upwards(heap_ptr);
}

void HEAP_Check_memory(Heap *heap_ptr)
{
  if (heap_ptr->nnds > heap_ptr->maxnnds) {
    /* Update the maximum number of nodes */
    heap_ptr->maxnnds = 2*(heap_ptr->nnds) - 1;
    /* Reallocate the heap array (nodes + SENTINEL) */
    heap_ptr->node = (Heap_node *)
      realloc(heap_ptr->node, (heap_ptr->maxnnds + 1)*sizeof(Heap_node));
    assert(heap_ptr->node != NULL);   /* Check if memory was allocated    */
  }
}

void HEAP_Upwards(Heap *heap_ptr)
{
  /* Definitions of local variables */
  unsigned long father, son;
  Heap_node   newnode;

  /* Get the new node at the last position of the heap node */
  newnode = heap_ptr->node[heap_ptr->nnds];
  son   = heap_ptr->nnds;         /* 'newnode' position         */
  father  = son/2;              /* position of the 'newnode' father   */

  /* Restore the heap condition from the bottom upwards */
  while (heap_ptr->node[father].key >= newnode.key) {
      /* Move node down */
      heap_ptr->node[son] = heap_ptr->node[father];
      /* Move index up */
      son   =  father;
      father  /= 2;
  }
  /* Correct position of the new node */
  heap_ptr->node[son] = newnode;
}

void HEAP_Remove_root(Heap *heap_ptr)
{
  /* Put the last heap node at the root position and decrease 'nnds' */
  heap_ptr->node[1] = heap_ptr->node[heap_ptr->nnds];
  heap_ptr->node[heap_ptr->nnds].id  = 0;
  heap_ptr->node[heap_ptr->nnds].key = 0.0e+000;
  (heap_ptr->nnds)--;

  /* Traverse the heap and restore the heap condition */
  HEAP_Downwards(1, heap_ptr);
}

unsigned long HEAP_Get_root(Heap *heap_ptr)
{
  return (heap_ptr->node[1].id);
}

Boolean HEAP_Change_node_id(unsigned long oldid, unsigned long newid, Heap *heap_ptr)
{
  unsigned long nodepos;

  for (nodepos = 1; nodepos <= heap_ptr->nnds; nodepos++) {
    if (heap_ptr->node[nodepos].id == oldid) {
      heap_ptr->node[nodepos].id = newid;
      return (TRUE);
    }
  }
  return (FALSE);
}

void HEAP_Change_node_key(double newkey, unsigned long ndid, Heap *heap_ptr)
{
  unsigned long nodepos;

  for (nodepos = 1; heap_ptr->node[nodepos].id != ndid; nodepos++)
    continue;

  heap_ptr->node[nodepos].key = newkey;
  HEAP_Downwards(1, heap_ptr);
}

Boolean HEAP_Remove_node(unsigned long triangleid, Heap *heap_ptr)
{
  unsigned long nodepos;
  
  /* Search the position of the node with id number equal to 'triangleid' */
  for (nodepos = 1; nodepos <= heap_ptr->nnds; nodepos++) {
    if (heap_ptr->node[nodepos].id == triangleid) {
      /* Put the last heap node at the root position and decrease 'nnds' */
      heap_ptr->node[nodepos] = heap_ptr->node[(heap_ptr->nnds)--];
      /* Traverse the heap and restore the heap condition */
      HEAP_Downwards(nodepos, heap_ptr);
      return (TRUE);
    }
  }
  return (FALSE);
}

void HEAP_Downwards(unsigned long nodepos, Heap *heap_ptr)
{
  /* Definitions of local variables */
  unsigned long father, son;
  Heap_node   fthnode;

  /* Store the heap node at position 'nodepos' */
  fthnode = heap_ptr->node[nodepos];

  /* Restore the heap condition from the top downwards */
  father  = nodepos;
  while (father <= ((heap_ptr->nnds)/2)) {
    son = father + father;          /* Left son of 'father'         */
    /* Get the id number of the son with the smaller key */
    if (son < (heap_ptr->nnds) &&
        (heap_ptr->node[son].key > heap_ptr->node[son + 1].key))
      son++;                /* Right son of 'father'        */

    /* Compare the father and the son's keys */
    if (fthnode.key <= heap_ptr->node[son].key)
      break;
    /* Move node up */
    heap_ptr->node[father] = heap_ptr->node[son];
    /* Move index down */
    father = son;             
  }
  /* Correct position of the old top heap node */
  heap_ptr->node[father] = fthnode;   
}

void HEAP_Finalize(Heap *heap_ptr)
{
  /* Deallocate the heap array and set it to NULL */
  free(heap_ptr->node);
  heap_ptr->node = NULL;
  assert(heap_ptr->node == NULL);     /* Check if memory was deallocated    */

  /* Deallocate the main heap structure and set it to NULL */
  free(heap_ptr);
  heap_ptr     = NULL;
  assert(heap_ptr == NULL);         /* Check if memory was deallocated    */
}
