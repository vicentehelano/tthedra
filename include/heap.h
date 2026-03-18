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
#ifndef SAFT3D_HEAP_H
#define SAFT3D_HEAP_H

#include <geom.h>

enum Heap_const{
  SENTINEL = -1,
};

/* Each heap node contains the id number of a triangle and a key
   that is a double value variable. The key can be the perimeter,
   the area or other double value variable.
   The current key is the perimeter of a triangle.*/
typedef struct  heap_node {
  unsigned long id;
  double      key;
/*  Boolean operator<(const heap_node& node) const {
    if ((this->key < node.key)) { return true; }
    return false;
  }*/
} Heap_node;

/*
 * Main heap structure
 *
 * Comments:
 *    The main heap structures contains the current number of nodes
 *    ('nnds'), the maximum number of nodes ('maxnode') and a pointer
 *    to a heap node. The maximum number of nodes will be changed if
 *    more memory block is necessary.
 *
 * Last modified: 18/06/2004.
 */
typedef struct  heap {
  unsigned long nnds;
  unsigned long maxnnds;
  Heap_node   *node;
} Heap;

extern Heap *front_heap;      /* Pointer to a main heap structure   */
extern Heap *reject_heap;     /* Pointer to a main heap structure   */

/*
 * It allocates the main heap structure and respective sub-structures. The
 * memory allocated to the heap array is sufficient to store 14 levels of
 * the binary tree corresponding to the heap and the sentinel node.
 */
void HEAP_Initialize_front_heap(Surface_mesh *, Volume_mesh *);
void HEAP_Initialize(Heap *);

/*
 * It can return any double value variable. Actually, it returns
 * the perimeter of the triangle with id number 'triangleid' defined
 * in the boundary structure 'front'.
 */
double HEAP_Key(Triangle *, Volume_mesh *);

/*
 * It inserts a new node at the end of the heap array and
 * calls HEAP_Upwards() to restore the heap condition.
 */
void HEAP_Insert(unsigned long, double, Heap *);

/*
 * It checks if the currently available memory block is enough to
 * insert a new node in the heap array. If necessary, it allocates
 * a new memory block capable to store the old binary tree corresponding
 * to the heap plus one level.
 */
void HEAP_Check_memory(Heap *);

/*
 * Called by HEAP_Insert(). It traverses the heap from the bottom
 * upwards. If necessary the internal order is reestablished by
 * comparing father and son pairs.
 */
void HEAP_Upwards(Heap *);

/*
 * HEAP_Get_root()
 *
 * Prototype:
 *    unsigned long HEAP_Get_root(Heap *);
 *
 * Required Header:
 *    heap.h.
 *
 * Return Value:
 *    It returns the id number of the root node of the heap.
 *
 * Actual Parameters:
 *    It receives a pointer to a heap structure.
 *
 * Comments:
 *    It just returns the id number of the root node of the heap and
 *    do not remove the root from the heap.
 *
 * Last modified: 18/06/2004.
 */
unsigned long HEAP_Get_root(Heap *);

/*
 * It removes the root of the heap and calls HEAP_Downwards()
 * to restore the heap condition.
 */
void      HEAP_Remove_root(Heap *);

unsigned long HEAP_Get_root(Heap *);

Boolean     HEAP_Change_node_id(unsigned long, unsigned long, Heap *);

void      HEAP_Change_node_key(double, unsigned long, Heap *);

/*
 * It removes the node with id number equal to 'triangleid' and
 * calls HEAP_Downwards() to restore the heap condition.
 */
Boolean HEAP_Remove_node(unsigned long, Heap *);

/*
 * It traverses the heap from the top downwards and, if necessary, the
 * heap condition is restored by comparing father and sons pairs.
 */
void HEAP_Downwards(unsigned long, Heap *);

/*
 * It deallocates the main heap structure and respective sub-structures.
 */
void HEAP_Finalize(Heap *);

#endif /* SAFT3D_HEAP_H */
