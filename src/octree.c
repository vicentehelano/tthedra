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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <octree.h>
#include <saft2d.h>

/*
 * Global variables
 *
 * Last modified: 19/07/2004.
 */
Octree *front_octree;     /* Pointer to a main octree structure */
/* Multiplicative factors:       OCT1  OCT2  OCT3  OCT4  OCT5  OCT6  OCT7    OCT8*/
static double xf[]  =   { -0.50,  +0.50,  -0.50,  +0.50,  -0.50,  +0.50,  -0.50,  +0.50};
static double yf[]  =   { -0.50,  -0.50,  +0.50,  +0.50,  -0.50,  -0.50,  +0.50,  +0.50};
static double zf[]  =   { -0.50,  -0.50,  -0.50,  -0.50,  +0.50,  +0.50,  +0.50,  +0.50};
/* Array with number of nodes per level (ten levels) */
static unsigned long level[] =  {1, 8, 64, 512, 4096, 32768, 262144, 2097152, 16777216, 134217728};

/*
 * OCTREE_Initialize()
 *
 * Prototype:
 *    void OCTREE_Initialize(Octree **, Boundary *);
 *
 * Required header:
 *    octree.h and stdlib.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives the address of a pointer to a octree structure yet not allocated and
 *    a pointer to a boundary structure.
 *
 * Comments:
 *    It allocates the main octree structure and respective sub-structures. The
 *    memory allocated to the heap array is sufficient to store 14 levels of
 *    the binary tree corresponding to the heap and the sentinel node.
 *
 * Last modified: 19/07/2004.
 */
void OCTREE_Initialize_front_octree(Volume_mesh *domain)
{
  /* Definitions of local variables */
  unsigned long j;
  int       i;

  /* Allocate a block of memory to the main octree structure */
  front_octree = (Octree *) malloc(sizeof(Octree));
  assert(front_octree != NULL);     /* Check if memory was allocated    */

  /* Initialize the maximum number of nodes */
  front_octree->maxnnds   = 0;
  /* Initial maximum number of levels (six levels) */
  front_octree->nlevels = INIT_NLEVELS;

  /* Calculate how many nodes are in an octree with six levels */
  for (i = 0; i < front_octree->nlevels; i++)
    front_octree->maxnnds += level[i];

  /* Allocate memory to the node array */
  front_octree->node    = (Octree_node *)
    malloc(front_octree->maxnnds*sizeof(Octree_node));
  assert(front_octree->node != NULL);   /* Check if memory was allocated    */

  /* Get the bounding volume of the domain */
  OCTREE_Get_bounds(front_octree, domain);

  /* Set initial status of node 1 to empty */
  front_octree->node[0].status  = 0;
  /* Set pointer 'point' of node 1 to NULL */
  front_octree->node[0].point   = NULL;

  /* Insert the domain points into the octree */
  for (j = 1; j <= domain->nnds; j++) {
    OCTREE_Insert(j, front_octree, domain);
  }
}

void OCTREE_Get_bounds(Octree *front_octree, Volume_mesh *domain)
{
  /* Definitions of local variables */
  unsigned long i;

  /* Initialize min and max front_octree variables */
  front_octree->bounds.min = domain->node[0];
  front_octree->bounds.max = domain->node[0];

  /* Get the min and max coordinates of the front_node array */
  for (i = 1; i < domain->nnds; i++) {
    if (domain->node[i].x < front_octree->bounds.min.x)
      front_octree->bounds.min.x = domain->node[i].x;

    if (domain->node[i].y < front_octree->bounds.min.y)
      front_octree->bounds.min.y = domain->node[i].y;

    if (domain->node[i].z < front_octree->bounds.min.z) 
      front_octree->bounds.min.z = domain->node[i].z;

    if (domain->node[i].x > front_octree->bounds.max.x)
      front_octree->bounds.max.x = domain->node[i].x;

    if (domain->node[i].y > front_octree->bounds.max.y)
      front_octree->bounds.max.y = domain->node[i].y;

    if (domain->node[i].z > front_octree->bounds.max.z)
      front_octree->bounds.max.z = domain->node[i].z;
  }

  /* Calculate the front_octree bounding box center */
  front_octree->center.x = front_octree->bounds.min.x +
    ((front_octree->bounds.max.x - front_octree->bounds.min.x)/2);

  front_octree->center.y = front_octree->bounds.min.y +
    ((front_octree->bounds.max.y - front_octree->bounds.min.y)/2);

  front_octree->center.z = front_octree->bounds.min.z +
    ((front_octree->bounds.max.z - front_octree->bounds.min.z)/2);
}



void OCTREE_Insert(unsigned long ptid, Octree *front_octree, Volume_mesh *domain)
{
  /* Definitions of local variables */
  int       i;
  unsigned long leafpos, ptpos;

  ptpos = ptid - 1;
  /* Get the leaf node that contains point 'ptid' */
  leafpos = OCTREE_Get_leaf(&(domain->node[ptpos]), front_octree);

  /* Check if the octree node 'leafid' is not full */
  if (front_octree->node[leafpos].status < BUCKETSIZE) {
    /* Check if the point id number array must be allocated */
    if (front_octree->node[leafpos].status == 0) {
      /* Check if node is really empty */
      assert(front_octree->node[leafpos].point == NULL);
      /* Allocate memory */
      front_octree->node[leafpos].point = (unsigned long *)
        malloc(BUCKETSIZE*sizeof(unsigned long));
      /* Check if memory was allocated */
      assert(front_octree->node[leafpos].point != NULL);
    }
    /* Increase node status by one */
    front_octree->node[leafpos].status++;
    /* Get the position of the point in the point id number array */
    ptpos = front_octree->node[leafpos].status - 1;
    /* Store the new point 'ptid' */
    front_octree->node[leafpos].point[ptpos] = ptid;
  }
  else {                    /* current node is full. Then...    */
    /* ... split node */
    OCTREE_Split_node(leafpos, front_octree);
    /* and redistribute the old points into its sons */
    for (i = 0; i < BUCKETSIZE; i++) {
      OCTREE_Insert(front_octree->node[leafpos].point[i], front_octree, domain);
    }
    /* Deallocate the point array of the new non-leaf node */
    free(front_octree->node[leafpos].point);
    /* Set the pointer 'point' of the new non-leaf node to NULL */
    front_octree->node[leafpos].point = NULL;
    /* Finally, insert the new point 'ptid' to the new octree arrangement */
    OCTREE_Insert(ptid, front_octree, domain);
  }
}



void OCTREE_Split_node(unsigned long leafpos, Octree *front_octree)
{
  /* Definitions of local variables */
  unsigned long i, firstson, lastson;

  /* Set status of node at 'leafpos' to non-leaf */
  front_octree->node[leafpos].status  = -1;

  /* Positions of the first and last sons of node at 'leafpos' */
  firstson = 8*leafpos + OCT1;
  lastson  = 8*leafpos + OCT8;

  /* Check if memory */
  OCTREE_Check_memory(firstson, front_octree);

  /* Set status of sons to empty and pointer 'point' to NULL */
  for (i = firstson; i <= lastson; i++) {
    front_octree->node[i].status  = 0;
    front_octree->node[i].point   = NULL;
  }
}


Octant OCTREE_Compare(Point3D *p, Point3D *center)
{
  if (p->x < center->x) {           /* 'p' can be in octants no. 1, 3, 5, 7 */
    if (p->y < center->y) {         /* 'p' can be in octants no. 1, 5   */
      if (p->z < center->z) {       /* Point 'p' is in octant no. 1     */
        return (OCT1);
      }
      else {                /* Point 'p' is in octant no. 5     */
        return (OCT5);
      }
    }
    else if (p->z < center->z) {      /* Point 'p' is in octant no. 3     */
      return (OCT3);
    }
    else {
      return (OCT7);            /* Point 'p' is in octant no. 7     */
    }
  }
  else if (p->y < center->y) {        /* 'p' can be in octants no. 2, 6   */
    if (p->z < center->z) {         /* Point 'p' is in octant no. 2     */
      return (OCT2);
    }
    else {
      return (OCT6);
    }
  }
  else if (p->z < center->z)  {       /* Point 'p' is in octant no. 4     */
    return (OCT4);
  }
  else {                    /* Point 'p' is in octant no. 8     */
    return (OCT8);
  }
}


void OCTREE_Check_memory(unsigned long firstson, Octree *front_octree)
{
  if (firstson > front_octree->maxnnds) {
    /* Increase one more level to the tree */
    front_octree->nlevels++;
    assert(front_octree->nlevels <= MAX_NLEVELS);
    /* Update the maximum number of nodes */
    front_octree->maxnnds += level[front_octree->nlevels - 1];
    /* Reallocate the octree node array */
    front_octree->node = (Octree_node *)
      realloc(front_octree->node, (front_octree->maxnnds)*sizeof(Octree_node));
    assert(front_octree->node != NULL);   /* Check if memory was allocated    */
  }
}



unsigned long OCTREE_Get_leaf(Point3D *p, Octree *front_octree)
{
  /* Definitions of local variables */
  unsigned int  octid, mfpos;
  unsigned long leafpos;
  double      dx, dy, dz;
  Point3D     center;
  
  /* Get the bounding box center */
  center  = front_octree->center;
  /* Get the bounding box edge lengths */
  dx    = front_octree->bounds.max.x - front_octree->bounds.min.x;
  dy    = front_octree->bounds.max.y - front_octree->bounds.min.y;
  dz    = front_octree->bounds.max.z - front_octree->bounds.min.z;

  /* Set 'leafid' to octree root id */
  leafpos = 0;
  /* Search for the leaf node where point 'p' is located */
  while (front_octree->node[leafpos].status < 0) {
    /* Get what octant 'p' is located */
    octid = OCTREE_Compare(p, &center);
    /* Get octant position for the multiplicative factors */
    mfpos = octid - 1;

    /* Calculate the edge lengths of this octant */
    dx /= 2;
    dy /= 2;
    dz /= 2;

    /* Calculate the center of this octant */
    center.x = center.x + xf[mfpos]*dx;
    center.y = center.y + yf[mfpos]*dy;
    center.z = center.z + zf[mfpos]*dz;

    /* Calculate the position in the octree array of this octant */
    leafpos = 8*leafpos + octid;
  }
  return (leafpos);
}




Boolean OCTREE_Remove(unsigned long ptid, Octree *front_octree, Volume_mesh *domain)
{
  /* Definitions of local variables */
  unsigned long i, leafpos;
  unsigned long ptpos, npts;
  unsigned long lastptpos, lastptid;

  ptpos = ptid - 1;
  /* Get the leaf node that contains point 'ptid' */
  leafpos = OCTREE_Get_leaf(&(domain->node[ptpos]), front_octree);

  /* Get the number of points inside node 'leafid' */
  npts = front_octree->node[leafpos].status;

  /* Get position and id number of the last point stored at 'leafid' */
  lastptpos = front_octree->node[leafpos].status - 1;
  lastptid  = front_octree->node[leafpos].point[lastptpos];

  /* Search 'ptid' inside node 'leaf' */
  for (i = 0; i < npts; i++) {
    if (ptid == front_octree->node[leafpos].point[i]) {
      if (i != lastptpos) {       /* If 'ptid' is not at last position */
        /* Put the last point at position of the removed point */
        front_octree->node[leafpos].point[i] = lastptid;
        /* And set last point to store the zero id */
        front_octree->node[leafpos].point[lastptpos] = 0;
      }
      else {
        /* Just set last point to store the zero id */
        front_octree->node[leafpos].point[lastptpos] = 0;
      }
      /* Then decrease by one the status of 'leafid' node */
      front_octree->node[leafpos].status--;

      /* If this node is going to be empty... */
      if (front_octree->node[leafpos].status == 0) {
        /* Deallocate the point array and... */
        free(front_octree->node[leafpos].point);
        /* Set pointer 'point' of node 'leafid' to NULL */
        front_octree->node[leafpos].point = NULL;
      }
      /* The remove operation was sucessful */
      return (TRUE);
    }
  }
  /* The remove operation was not sucessful */
  return (FALSE);
}

void OCTREE_Range_search(Sphere *s,       Candidate_points  *cdtpts,
             Octree *front_octree,  Volume_mesh     *domain)
{
  /* Definitions of local variables */
  Candidate_nodes *cdtnds;
  Bounding_box  *b;

  /* Allocate a Bounding_box structure pointed by 'b' */
  b = (Bounding_box *) malloc(sizeof(Bounding_box));
  assert(b != NULL);              /* Check if memory was allocated    */
  /* Get bounds of the front_octree */
  b->min = front_octree->bounds.min;
  b->max = front_octree->bounds.max;

  /* Allocate a Candidate_nodes structure pointed by 'cdtnds' */
  cdtnds = (Candidate_nodes *) malloc(sizeof(Candidate_nodes));
  assert(cdtnds != NULL);           /* Check if memory was allocated    */
  /* Initialize the number of candidate nodes */
  cdtnds->nnds = 0;
  /* Initialize the initial maximum number of candidate nodes */
  cdtnds->maxnnds = INIT_MAXCDTNDS;
  /* Allocate the unsigned long array pointed by 'cdtnds->node' */
  cdtnds->node = (unsigned long *) malloc(cdtnds->maxnnds*sizeof(unsigned long));
  assert(cdtnds->node != NULL);       /* Check if memory was allocated    */

  /* Get candidate nodes */
  OCTREE_Get_candidate_nodes(s, b, 0, cdtnds, front_octree);

  /* Initialize the number of candidate points */
  cdtpts->npts = 0;
  /* Initialize the maximum number of candidate points */
  cdtpts->maxnpts = BUCKETSIZE*cdtnds->nnds;  /* somei mais um, pois irei adicionar o ponto ideal no final */
  /* Allocate the unsigned long array pointed by 'cdtpts->point' */
  cdtpts->point = (unsigned long *)
    malloc((cdtpts->maxnpts)*sizeof(unsigned long));
  assert(cdtpts->point != NULL);        /* Check if memory was allocated    */

  /* Get candidate points */
  OCTREE_Get_candidate_points(s, cdtpts, cdtnds, front_octree, domain);
  
  /* Deallocate candidate nodes array */
  free(cdtnds->node);
  /* Deallocate the main Candidate_nodes structure */
  free(cdtnds);
  /* Deallocate the Bounding_box structure */
  free(b);
}


void OCTREE_Get_candidate_nodes(Sphere *s, Bounding_box *b,
                     unsigned long ndid,
                     Candidate_nodes *cdtnds,
                     Octree *front_octree)
{
  /* Definitions of local variables */
  int       i, mfpos, status;
  unsigned long *ulptr;
  double      dx, dy, dz;
  Point3D     fathctr, sonctr;
  Bounding_box  oct;

  /* Filtrating process */
  status = front_octree->node[ndid].status;
  if (OCTREE_Box_sphere_intersect(b, s) == TRUE) {
    /* node is a non-empty leaf node */
    if (status > 0 && status <= BUCKETSIZE) {
      (cdtnds->nnds)++;
      /* Check memory */
      if (cdtnds->nnds > cdtnds->maxnnds) {
        cdtnds->maxnnds *= 2;
        ulptr = (unsigned long *)
          realloc(cdtnds->node, cdtnds->maxnnds*sizeof(unsigned long));
        assert(ulptr != NULL);  /* Check if memory was allocated    */
        cdtnds->node = ulptr;
      }
      /* Store the node id inside the candidate nodes array */
      cdtnds->node[cdtnds->nnds - 1] = ndid;
    }/* node is a non-leaf node */
    else if (status < 0) {
      /* New edge length of sons */
      dx  = (b->max.x - b->min.x)/2.0;
      dy  = (b->max.y - b->min.y)/2.0;
      dz  = (b->max.z - b->min.z)/2.0;

      /* Center of father node */
      fathctr.x = b->min.x + dx;
      fathctr.y = b->min.y + dy;
      fathctr.z = b->min.z + dz;

      /* Traverse the sons and get candidate nodes */
      for (i = OCT1; i <= OCT8; i++) {
        /* Calculate center of the currently octree node */
        mfpos  = i - 1;
        sonctr.x = fathctr.x + xf[mfpos]*dx;
        sonctr.y = fathctr.y + yf[mfpos]*dy;
        sonctr.z = fathctr.z + zf[mfpos]*dz;

        /* Get the bounding box of the currently octree node */
        oct.min.x = sonctr.x + xf[0]*dx;
        oct.min.y = sonctr.y + yf[0]*dy;
        oct.min.z = sonctr.z + zf[0]*dz;

        oct.max.x = sonctr.x + xf[7]*dx;
        oct.max.y = sonctr.y + yf[7]*dy;
        oct.max.z = sonctr.z + zf[7]*dz;

        /* Get candidate nodes of sons */
        OCTREE_Get_candidate_nodes(s, &oct, 8*ndid + i, cdtnds, front_octree);
      }
    }
  }
}


void OCTREE_Get_candidate_points(Sphere      *s,    Candidate_points *cdtpts,
                 Candidate_nodes *cdtnds, Octree       *front_octree,
                 Volume_mesh   *domain)
{
  /* Definitions of local variables */
  int       j, npts;
  unsigned long i, ndid, pt2id;
  double      d2, r2;
  Point3D     pt1, pt2;

  /* Get point 1 coordinates and square of 's' radius */
  pt1 = s->center;
  r2 = (s->radius)*(s->radius);

  /* Check points inside each candidate node */
  for (i = 0; i < cdtnds->nnds; i++) {
    ndid = cdtnds->node[i];
    npts = front_octree->node[ndid].status;
    for (j = 0; j < npts; j++) {
      pt2id = front_octree->node[ndid].point[j];
      pt2   = domain->node[pt2id - 1];
      /* Calculate square of distance between 'pt1' and 'pt2' */
      d2    = (pt2.x - pt1.x)*(pt2.x - pt1.x) +
              (pt2.y - pt1.y)*(pt2.y - pt1.y) +
                (pt2.z - pt1.z)*(pt2.z - pt1.z);
      if (d2 <= r2) {           /* 'pt2' is inside Sphere 's'     */
        cdtpts->point[cdtpts->npts] = pt2id;
        (cdtpts->npts)++;
      }
    }
  }
}



Boolean OCTREE_Box_sphere_intersect(Bounding_box *b, Sphere *s)
{
  /* Definitions of local variables */
  double  dmin2, r2;
  
  /* Initialize local variables */
  dmin2 = 0.0;
  r2   = (s->radius)*(s->radius);

  /* Calculate the squared minimum distance (dmin2) to the sphere center */
  if (s->center.x < b->min.x)
    dmin2 += (s->center.x - b->min.x)*(s->center.x - b->min.x);
  else if(s->center.x > b->max.x)
    dmin2 += (s->center.x - b->max.x)*(s->center.x - b->max.x);

  if (s->center.y < b->min.y)
    dmin2 += (s->center.y - b->min.y)*(s->center.y - b->min.y);
  else if(s->center.y > b->max.y)
    dmin2 += (s->center.y - b->max.y)*(s->center.y - b->max.y);

  if (s->center.z < b->min.z)
    dmin2 += (s->center.z - b->min.z)*(s->center.z - b->min.z);
  else if(s->center.z > b->max.z)
    dmin2 += (s->center.z - b->max.z)*(s->center.z - b->max.z);

  if (dmin2 <= r2)
    return (TRUE);              /* Intersection             */

  return (FALSE);               /* No intersection            */
}

void OCTREE_Print(Octree *front_octree)
{
  /* Definitions of local variables */
  unsigned long i, id;
  unsigned long v1, v2, v3, v4, v5, v6, v7, v8;
  unsigned long ndid, npts, nels;
  Point3D     pt;
  Bounding_box  *b;
  FILE      *foct, *fpts, *fels;
  /* Output data file path          */
  #ifdef NDEBUG
    char fpath1[150] = "../data/output/";
    char fpath2[150] = "../data/output/";
    char fpath3[150] = "../data/output/";
  #else
    char fpath1[150] = "./../../data/output/";
    char fpath2[150] = "./../../data/output/";
    char fpath3[150] = "./../../data/output/";
  #endif

  /* Generate complete file path name */
  strcat(fpath1, projname);
  strcat(fpath1, "_oct.msh");
  strcat(fpath2, projname);
  strcat(fpath2, "_octpts.tmp");
  strcat(fpath3, projname);
  strcat(fpath3, "_octels.tmp");

  /* Opening files */
  foct = fopen(fpath1, "w");  /* Global octree file       */
  fpts = fopen(fpath2, "w");  /* Octree points file       */
  fels = fopen(fpath3, "w");  /* Octree elements file       */

  /* Initialize local variables */
  npts = 0;                 /* Number of octree points        */
  nels = 0;                 /* Number of octree elements      */
  ndid = 0;                 /* Octree node id number        */

  /* Allocate a Bounding_box structure pointed by 'b' */
  b = (Bounding_box *) malloc(sizeof(Bounding_box));
  assert(b != NULL);              /* Check if memory was allocated    */
  /* Get bounds of the front_octree */
  b->min = front_octree->bounds.min;
  b->max = front_octree->bounds.max;

  /* Print temporary files with points and elements separeted */
  OCTREE_Print_temps(b, ndid, &npts, &nels, fpts, fels, front_octree);
  
  /* Close files in mode 'w' */
  fclose(fpts);
  fclose(fels);

  /* Reopen closed files in new mode 'r' */
  fpts = freopen(fpath2, "r", fpts);
  fels = freopen(fpath3, "r", fels);
  
  
  fprintf(foct, "MESH    dimension 3 ElemType Hexahedra  Nnode 8\n\
    Coordinates");  
  for (i = 1; i <= npts; i++) {
    fscanf(fpts,  "\t%lu\t%lf\t%lf\t%lf\n", &id, &pt.x, &pt.y, &pt.z);
    fprintf(foct, "\t%lu\t%lf\t%lf\t%lf\n",  id,  pt.x,  pt.y,  pt.z);
  }
  fclose(fpts);
  remove(fpath2);
  fprintf(foct,"end coordinates\nElements");
  for (i = 1; i <= nels; i++) {
    fscanf(fels,  "\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",
      &id, &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8);
    fprintf(foct, "\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",
       id, v1,  v2,  v3,  v4,  v5,  v6,  v7,  v8);
  }
  fclose(fels);
  remove(fpath3);
  fprintf(foct,"end elements\n");
  fclose(foct);
}

void OCTREE_Print_temps(Bounding_box  *b,   unsigned long ndid,
            unsigned long *npts,  unsigned long *nels,
            FILE  *fpts, FILE *fels, Octree   *front_octree)
{
  /* Definitions of local variables */
  int       i, mfpos, status;
  double      dx, dy, dz;
  Point3D     sonctr, fathctr, corner;
  Bounding_box  oct;

  dx  = b->max.x - b->min.x;
  dy  = b->max.y - b->min.y;
  dz  = b->max.z - b->min.z;

  /* Center of father node */
  fathctr.x = b->min.x + dx/2.0;
  fathctr.y = b->min.y + dy/2.0;
  fathctr.z = b->min.z + dz/2.0;

  status = front_octree->node[ndid].status;
  if (status > 0 && status <= BUCKETSIZE) {
    (*nels)++;
    for (i = OCT1; i <= OCT8; i++) {
      mfpos = i - 1;
      corner.x = fathctr.x + xf[mfpos]*dx;
      corner.y = fathctr.y + yf[mfpos]*dy;
      corner.z = fathctr.z + zf[mfpos]*dz;
      (*npts)++;
      fprintf(fpts, "\t%lu\t%lf\t%lf\t%lf\n", (*npts), corner.x, corner.y, corner.z);
    }
    fprintf(fels, "\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\n",
      (*nels), (*npts - 1), (*npts - 3), (*npts - 2), (*npts - 0), (*npts - 5), (*npts - 7),
      (*npts - 6), (*npts - 4));

  }/* node is a non-leaf node */
  else if (front_octree->node[ndid].status == -1) {
    /* New edge length of sons */
    dx /= 2.0;
    dy /= 2.0;
    dz /= 2.0;

    /* Traverse the sons and get candidate nodes */
    for (i = OCT1; i <= OCT8; i++) {
      /* Calculate center of the currently octree node */
      mfpos  = i - 1;
      sonctr.x = fathctr.x + xf[mfpos]*dx;
      sonctr.y = fathctr.y + yf[mfpos]*dy;
      sonctr.z = fathctr.z + zf[mfpos]*dz;

      /* Get the bounding box of the currently octree node */
      oct.min.x = sonctr.x + xf[0]*dx;
      oct.min.y = sonctr.y + yf[0]*dy;
      oct.min.z = sonctr.z + zf[0]*dz;

      oct.max.x = sonctr.x + xf[7]*dx;
      oct.max.y = sonctr.y + yf[7]*dy;
      oct.max.z = sonctr.z + zf[7]*dz;

      /* Get candidate nodes of sons */
      OCTREE_Print_temps(&oct, 8*ndid + i, npts, nels, fpts, fels, front_octree);
    }
  }
}
