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
#include <string.h>
#include <assert.h>

#include <geom.h>
#include <adjac.h>

/* List of faces adjacent to points */
Adjacency *fasupt;

void ADJAC_Initialize_fasupt(Surface_mesh *boundary)
{
  unsigned long i, elemid;
  unsigned long nd1, nd2, nd3;

  fasupt = (Adjacency *) malloc(sizeof(Adjacency));

  fasupt->nnds    = boundary->nnds;
  fasupt->maxnnds = boundary->maxnnds;

  fasupt->node = (Adjac_node *)
    malloc(fasupt->maxnnds*sizeof(Adjac_node));
  assert(fasupt->node != NULL);

  for (i = 0; i < fasupt->maxnnds; i++) {
    fasupt->node[i].nnbs   = 0;
    fasupt->node[i].maxnnbs  = INIT_NNEIGHBORS_FACE;
    fasupt->node[i].neighbor = (unsigned long *)
      malloc(INIT_NNEIGHBORS_FACE*sizeof(unsigned long));
    assert(fasupt->node[i].neighbor != NULL);

    memset(fasupt->node[i].neighbor, 0,
      INIT_NNEIGHBORS_FACE*sizeof(unsigned long));
  }

  for (i = 0; i < boundary->nels; i++) {
    nd1 = boundary->element[i].vertex[0];
    nd2 = boundary->element[i].vertex[1];
    nd3 = boundary->element[i].vertex[2];

    elemid = i + 1;
    ADJAC_Insert_neighbor(nd1, elemid, fasupt);
    ADJAC_Insert_neighbor(nd2, elemid, fasupt);
    ADJAC_Insert_neighbor(nd3, elemid, fasupt);
  }
}

void ADJAC_Insert_neighbor(unsigned long ptid, unsigned long ngbrid, Adjacency *fasupt)
{
  unsigned long ptpos, *luptr;
  int       new_nnbs, ngbrpos;

/*  if (ptpos > (fasupt->maxnnds - 1)) {
    fasupt->nnds++;
    fasupt->maxnnds *= 2;
    fasupt->node = (Adjac_node *)
      realloc(fasupt->node, fasupt->maxnnds*sizeof(Adjac_node));

    for (i = fasupt->nnds; i < fasupt->maxnnds; i++) {
      fasupt->node[i].nnbs   = 0;
      fasupt->node[i].neighbor = (unsigned long *)
        malloc(INIT_NNEIGHBORS_FACE*sizeof(unsigned long));
    }
  }
*/
  ptpos = ptid - 1;
  ngbrpos   = fasupt->node[ptpos].nnbs;
  new_nnbs = fasupt->node[ptpos].nnbs + 1;

  if (new_nnbs > fasupt->node[ptpos].maxnnbs) {
    fasupt->node[ptpos].maxnnbs *= 2;

    luptr = (unsigned long *) realloc(fasupt->node[ptpos].neighbor,
      fasupt->node[ptpos].maxnnbs*sizeof(unsigned long));
    assert(luptr != NULL);
    fasupt->node[ptpos].neighbor = luptr;

    luptr = (fasupt->node[ptpos].neighbor + ngbrpos);
    memset(luptr, 0, ngbrpos*sizeof(unsigned long));
  }
  
  fasupt->node[ptpos].neighbor[ngbrpos] = ngbrid;
  fasupt->node[ptpos].nnbs++;
}


void ADJAC_Remove_neighbor(unsigned long ngbrid, unsigned long ptid, Adjacency *adjc)
{
  int i, nnbs;
  unsigned long ptpos;

  ptpos = ptid - 1;

  nnbs = adjc->node[ptpos].nnbs;
  assert(adjc->node[ptpos].nnbs > 0);
  if (nnbs > 0) {
    for (i = 0; adjc->node[ptpos].neighbor[i] != ngbrid; i++)
      continue;

    if (i < (nnbs - 1)) {
      adjc->node[ptpos].neighbor[i] = adjc->node[ptpos].neighbor[nnbs - 1];
      adjc->node[ptpos].neighbor[nnbs - 1] = 0;
    }
    else {
      adjc->node[ptpos].neighbor[i] = 0;
    }

    adjc->node[ptpos].nnbs--;
  }
}

void ADJAC_Change_neighbor_id(unsigned long oldid, unsigned long newid,
                unsigned long ptid, Adjacency *adjc)
{
  int i;
  unsigned long ptpos;

  ptpos = ptid - 1;
  assert(adjc->node[ptpos].nnbs > 0);
  for (i = 0; adjc->node[ptpos].neighbor[i] != oldid; i++)
    continue;

  adjc->node[ptpos].neighbor[i] = newid;
}
