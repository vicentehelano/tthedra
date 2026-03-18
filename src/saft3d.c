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
#include <time.h>
#include <assert.h>
#include <saft3d.h>
#include <saft2d.h>
#include <octree.h>
#include <heap.h>
#include <adjac.h>
#include <quali.h>
#include <geom.h>

static Volume_mesh *mesh;
static Surface_mesh *front;
unsigned long     fase1_nels0, fase1_nels1;
unsigned long     fase2_nels0, fase2_nels1;
double cpu_time;
double cpu_time_fase1;
double cpu_time_fase2;
double cpu_time_intersec;
double clockbefore_intersec;
double clockafter_intersec;
char   status[8];

void SAFT3D_Initialize_mesh()
{
  unsigned long i;

  mesh = (Volume_mesh *) malloc(sizeof(Volume_mesh));
  assert(mesh != NULL);

  mesh->nnds = boundary->nnds;
  mesh->maxnnds = boundary->maxnnds;

  mesh->node = (Point3D *) malloc(mesh->maxnnds*sizeof(Point3D));
  assert(mesh->node != NULL);

  /* Get boundary nodes */
  for (i = 0; i < boundary->nnds; i++) {
    mesh->node[i] = boundary->node[i];
  }

  mesh->element = (Tetrahedron *)
    malloc(mesh->maxnnds*sizeof(Tetrahedron));
  assert(mesh->element != NULL);
}

void SAFT3D_Initialize_front()
{
  unsigned long i;

  front = (Surface_mesh *) malloc(sizeof(Surface_mesh));
  assert(front != NULL);

  front->nels   = boundary->nels;
  assert(front->nels >= 4);

  front->maxnels  = boundary->maxnels;

  front->element  = (Triangle *)
    malloc(front->maxnels*sizeof(Triangle));
  assert(front->element != NULL);

  /* Set the boundary as the initial front of faces */
  for (i = 0; i < boundary->nels; i++) {
    front->element[i] = boundary->element[i];
  }
}

/* based on:
 *   Rypl, D.,
 *   Sequential and Parallel Generation of Unstructured 3D Meshes
 *   PhD thesis
 *   Czech Technical University in Prague
 *   Prague, 1998
 */
void get_ideal_point_3(double elemsize, Point3D *nd1,
                 Point3D *nd2, Point3D *nd3,
                 Point3D *idealpt)
{
  double  height, nmagn, invnmagn;
  Point3D n, centroid;
  Point3D vec12, vec13;

  /* Calculate triangle centroid */
  GEOM_TRIANGLE_CENTROID3D(centroid,(*nd1),(*nd2),(*nd3));

  GEOM_SUB3D(vec12,(*nd2),(*nd1))
  GEOM_SUB3D(vec13,(*nd3),(*nd1))
  
  /* Calculate an unity vector normal to this triangle through the centroid */
  GEOM_CROSS3D(n,vec12,vec13)
  nmagn = sqrt(GEOM_DOT3D(n,n));
  invnmagn = (1.0e+000)/nmagn;
  n.x = n.x*invnmagn;
  n.y = n.y*invnmagn;
  n.z = n.z*invnmagn;

  /* Calculate ideal point */
  height     = (8.1649658092772e-001)*elemsize;
  idealpt->x = height*(n.x) + centroid.x;
  idealpt->y = height*(n.y) + centroid.y;
  idealpt->z = height*(n.z) + centroid.z;
}

void SAFT3D_Get_neighbor_faces(Candidate_points *cdtpts,
    Candidate_triangles *cdttri, Adjacency *fasupt)
{
  unsigned long i, k;
  int j;

  cdttri->maxntri = 0;
  for (i = 0; i < cdtpts->npts; i++) {
    cdttri->maxntri += fasupt->node[cdtpts->point[i] - 1].nnbs;
  }

  cdttri->triangle = (unsigned long *)
    malloc(cdttri->maxntri*sizeof(unsigned long));
  memset(cdttri->triangle, 0, cdttri->maxntri*sizeof(unsigned long));

  cdttri->ntri = 0;
  for (i = 0; i < cdtpts->npts; i++) {
    for (j = 0; j < (fasupt->node[cdtpts->point[i] - 1].nnbs); j++) {
      k = 0;
      while (fasupt->node[cdtpts->point[i] - 1].neighbor[j] != cdttri->triangle[k] && (k < cdttri->ntri)) {
        k++;
      }

      if (k == cdttri->ntri) {
        cdttri->triangle[k] = fasupt->node[cdtpts->point[i] - 1].neighbor[j];
        cdttri->ntri++;
      }
      else {
        continue;
      }
    }
  }
}

void SAFT3D_Get_range(Point3D *pt, Point3D *nd1, Point3D *nd2,
            Point3D *nd3, Sphere *range)
{
  double  sqrd1, sqrd2, sqrd3;
  Point3D v1;

  GEOM_SUB3D(v1,(*pt),(*nd1))
  sqrd1 = GEOM_DOT3D(v1,v1);

  GEOM_SUB3D(v1,(*pt),(*nd2))
  sqrd2 = GEOM_DOT3D(v1,v1);

  GEOM_SUB3D(v1,(*pt),(*nd3))
  sqrd3 = GEOM_DOT3D(v1,v1);

  range->center = *pt;
  range->radius = GEOM_MAX3(sqrd1,sqrd2,sqrd3);
  range->radius = (1.1e+000) * sqrt(range->radius);
}

void SAFT3D_Remove_selected_face(unsigned long triid, Candidate_triangles *cdttri)
{
  unsigned long i;
  
  assert(triid <= front->nels);
  for (i = 0; i < cdttri->ntri; i++) {
    if (cdttri->triangle[i] == triid)
      break;
  }
  assert(i <= cdttri->ntri);
  if (i < (cdttri->ntri - 1) ) {
    cdttri->triangle[i] = cdttri->triangle[cdttri->ntri - 1];
    cdttri->triangle[cdttri->ntri - 1] = 0;
    cdttri->ntri--;
  }
  else {
    cdttri->triangle[cdttri->ntri - 1] = 0;
    cdttri->ntri--;
  }
}


void SAFT3D_Remove_face_vertices(Triangle *t, Candidate_points *cdtpts)
{
  unsigned long i;

  for (i = 0; i < cdtpts->npts; i++) {
    if (cdtpts->point[i] == t->vertex[0])
      break;
  }
  assert(i <= cdtpts->npts);
  if (i < (cdtpts->npts - 1) ) {
    cdtpts->point[i] = cdtpts->point[cdtpts->npts - 1];
    cdtpts->point[cdtpts->npts - 1] = 0;
    cdtpts->npts--;
  }
  else {
    cdtpts->point[cdtpts->npts - 1] = 0;
    cdtpts->npts--;
  }

  for (i = 0; i < cdtpts->npts; i++) {
    if (cdtpts->point[i] == t->vertex[1])
      break;
  }
  assert(i <= cdtpts->npts);
  if (i < (cdtpts->npts - 1) ) {
    cdtpts->point[i] = cdtpts->point[cdtpts->npts - 1];
    cdtpts->point[cdtpts->npts - 1] = 0;
    cdtpts->npts--;
  }
  else {
    cdtpts->point[cdtpts->npts - 1] = 0;
    cdtpts->npts--;
  }

  for (i = 0; i < cdtpts->npts; i++) {
    if (cdtpts->point[i] == t->vertex[2])
      break;
  }
  assert(i <= cdtpts->npts);
  if (i < (cdtpts->npts - 1) ) {
    cdtpts->point[i] = cdtpts->point[cdtpts->npts - 1];
    cdtpts->point[cdtpts->npts - 1] = 0;
    cdtpts->npts--;
  }
  else {
    cdtpts->point[cdtpts->npts - 1] = 0;
    cdtpts->npts--;
  }
}

Boolean SAFT3D_Empty_tetrahedron_test(unsigned long nd4pos,
                    unsigned long idealptpos,
                    Point3D *nd1, Point3D *nd2,
                    Point3D *nd3, Point3D *nd4,
                    Candidate_points *cdtpts)
{
  unsigned long i, result;
  Point3D     *pt;

  result = TRUE;
  if (cdtpts->npts > 1) {
    if (nd4pos != idealptpos) {
      GEOM_SWAP(cdtpts->point[nd4pos],
        cdtpts->point[idealptpos - 1],i)
    }

    for (i = 0; i < idealptpos - 1; i++) {
      /* Obtem as coordenadas do ponto candidato a ser verificado */
      pt = &(mesh->node[cdtpts->point[i] - 1]);
      if (GEOM_In_tetrahedron(pt, nd1, nd2, nd3, nd4) == TRUE) {
        result = FALSE;
        break;
      }
    }

    if (nd4pos != idealptpos) {
      GEOM_SWAP(cdtpts->point[nd4pos],
        cdtpts->point[idealptpos - 1],i)
    }
  }
  return (result);
}


int SAFT3D_Triangle_nshared_nodes(Triangle *t1, Triangle *t2)
{
  int nshrdnds;

  nshrdnds = 0;
  if (t1->vertex[0] == t2->vertex[0])
    nshrdnds++;
  if (t1->vertex[0] == t2->vertex[1])
    nshrdnds++;
  if (t1->vertex[0] == t2->vertex[2])
    nshrdnds++;

  if (t1->vertex[1] == t2->vertex[0])
    nshrdnds++;
  if (t1->vertex[1] == t2->vertex[1])
    nshrdnds++;
  if (t1->vertex[1] == t2->vertex[2])
    nshrdnds++;

  if (t1->vertex[2] == t2->vertex[0])
    nshrdnds++;
  if (t1->vertex[2] == t2->vertex[1])
    nshrdnds++;
  if (t1->vertex[2] == t2->vertex[2])
    nshrdnds++;

  return (nshrdnds);
}


Boolean SAFT3D_Intersection_tests(Tetrahedron *newtet,
    Candidate_triangles *cdttri, unsigned long nshrdnds[3][8])
{
  unsigned long i, j, k;
  Point3D     *p1, *q1, *r1;
  Point3D     *p2, *q2, *r2;
  Triangle    newtri[3], oldtri;

  /* Face 01: */
  newtri[0].vertex[0] = newtet->vertex[1];
  newtri[0].vertex[1] = newtet->vertex[3];
  newtri[0].vertex[2] = newtet->vertex[2];

  /* Face 02: */
  newtri[1].vertex[0] = newtet->vertex[0];
  newtri[1].vertex[1] = newtet->vertex[2];
  newtri[1].vertex[2] = newtet->vertex[3];

  /* Face 03: */
  newtri[2].vertex[0] = newtet->vertex[0];
  newtri[2].vertex[1] = newtet->vertex[3];
  newtri[2].vertex[2] = newtet->vertex[1];

  /* Sorting: put smallest index at the first node to preserve orientation */
  for (i = 0; i < 3; i++) {
    while (newtri[i].vertex[0] > newtri[i].vertex[1] ||
      newtri[i].vertex[0] > newtri[i].vertex[2]) {
      k         = newtri[i].vertex[0];
      newtri[i].vertex[0] = newtri[i].vertex[1];
      newtri[i].vertex[1] = newtri[i].vertex[2];
      newtri[i].vertex[2] = k;
    }
  }

  for (i = 0; i < 3; i++) {
    p1 = &(mesh->node[newtri[i].vertex[0] - 1]);
    q1 = &(mesh->node[newtri[i].vertex[1] - 1]);
    r1 = &(mesh->node[newtri[i].vertex[2] - 1]);

    for (j = 0; j < cdttri->ntri; j++) {
      nshrdnds[i][0] = 0, nshrdnds[i][1] = 0, nshrdnds[i][2] = 0;
      nshrdnds[i][3] = 0, nshrdnds[i][4] = 0, nshrdnds[i][5] = 0;
      nshrdnds[i][6] = 0, nshrdnds[i][7] = 0;

      nshrdnds[i][7] = cdttri->triangle[j];
      oldtri   = front->element[cdttri->triangle[j] - 1];

      if (newtri[i].vertex[0] == oldtri.vertex[0]) {
        nshrdnds[i][0]++, nshrdnds[i][1]++, nshrdnds[i][4]++;
      }
      else if (newtri[i].vertex[0] == oldtri.vertex[1]) {
        nshrdnds[i][0]++, nshrdnds[i][1]++, nshrdnds[i][5]++;
      }
      else {
        if (newtri[i].vertex[0] == oldtri.vertex[2])
          nshrdnds[i][0]++, nshrdnds[i][1]++, nshrdnds[i][6]++;
      }

      if (newtri[i].vertex[1] == oldtri.vertex[0]) {
        nshrdnds[i][0]++, nshrdnds[i][2]++, nshrdnds[i][4]++;
      }
      else if (newtri[i].vertex[1] == oldtri.vertex[1]) {
        nshrdnds[i][0]++, nshrdnds[i][2]++, nshrdnds[i][5]++;
      }
      else {
        if (newtri[i].vertex[1] == oldtri.vertex[2])
          nshrdnds[i][0]++, nshrdnds[i][2]++, nshrdnds[i][6]++;
      }

      if (newtri[i].vertex[2] == oldtri.vertex[0]) {
        nshrdnds[i][0]++, nshrdnds[i][3]++, nshrdnds[i][4]++;
      }
      else if (newtri[i].vertex[2] == oldtri.vertex[1]) {
        nshrdnds[i][0]++, nshrdnds[i][3]++, nshrdnds[i][5]++;
      }
      else {
        if (newtri[i].vertex[2] == oldtri.vertex[2])
          nshrdnds[i][0]++, nshrdnds[i][3]++, nshrdnds[i][6]++;
      }

      p2 = &(mesh->node[oldtri.vertex[0] - 1]);
      q2 = &(mesh->node[oldtri.vertex[1] - 1]);
      r2 = &(mesh->node[oldtri.vertex[2] - 1]);

      /* Passar nshrdnds para a funçao de interseçao */
      if (GEOM_Tri_tri_overlap_test3D(p1, q1, r1, p2, q2, r2,
        nshrdnds[i]) == TRUE)
        return (FALSE);
      else {
        if (nshrdnds[i][0] > 2)
          break;
      }
    }
  }
  return (TRUE);
}


void SAFT3D_Remove_front_face(unsigned long faceid, unsigned long nshrdnds[3][8],
                Adjacency *fasupt, Heap *front_heap, Heap *reject_heap,
                Surface_mesh *front)
{
  unsigned long lastid, v1id, v2id, v3id;
  Boolean result;

  assert(faceid <= front->nels);
  if (faceid < front->nels) {
    lastid = front->nels;
    front->element[faceid - 1] = front->element[lastid - 1];

    v1id = front->element[faceid - 1].vertex[0];
    ADJAC_Change_neighbor_id(lastid, faceid, v1id, fasupt);

    v2id = front->element[faceid - 1].vertex[1];
    ADJAC_Change_neighbor_id(lastid, faceid, v2id, fasupt);

    v3id = front->element[faceid - 1].vertex[2];
    ADJAC_Change_neighbor_id(lastid, faceid, v3id, fasupt);

    if (lastid == nshrdnds[0][7])
      nshrdnds[0][7] = faceid;

    if (lastid == nshrdnds[1][7])
      nshrdnds[1][7] = faceid;

    if (lastid == nshrdnds[2][7])
      nshrdnds[2][7] = faceid;

    result = HEAP_Change_node_id(lastid, faceid, front_heap);
    if (result != TRUE) {
      result = HEAP_Change_node_id(lastid, faceid, reject_heap);
    }
    assert(result == TRUE);
  }

  front->nels--;
}

void SAFT3D_Print_front(Adjacency *fasupt, Surface_mesh *front, Volume_mesh *mesh)
{
  unsigned long i;
  FILE      *fptr;
  /* Output data file path          */
  #ifdef NDEBUG
    char fpath[150] = "../data/output/";
  #else
    char fpath[150] = "./../../data/output/";
  #endif

  if (front->nels < 1)
    return;

    /* Generate complete file path name */
  strcat(fpath, projname);
  strcat(fpath, "_front-");
  strcat(fpath, status);
  strcat(fpath, ".msh");

  /* Opening files */
  fptr = fopen(fpath, "w");

  fprintf(fptr, "MESH    dimension 3 ElemType Triangle  Nnode 3\n");
  fprintf(fptr, "Coordinates\n");
  for (i = 0; i < mesh->nnds; i++) {
    if (fasupt->node[i].nnbs > 0) {
    fprintf(fptr, "%17lu%17.8lf%17.8lf%17.8lf\n", i+1, mesh->node[i].x,
      mesh->node[i].y, mesh->node[i].z);
    }
  }

  fprintf(fptr, "end coordinates\nElements\n");
  for (i = 0; i < front->nels; i++) {
    fprintf(fptr, "%17lu%17lu%17lu%17lu\n",
       i+1, front->element[i].vertex[0],  front->element[i].vertex[1],
       front->element[i].vertex[2]);
  }
  fprintf(fptr, "end elements\n");
  fclose(fptr);
}

void SAFT3D_Update_data_structures(unsigned long idealptid, unsigned long tid,
                   Triangle *t, Tetrahedron   *newtet,
                   unsigned long nshrdnds[3][8],
                   Heap *front_heap)
{
  int   result;
  double  heapkey;

  /* ATUALIZANDO A FRONTEIRA */
  /* Remover a face base do HEAP */
  HEAP_Remove_root(front_heap);

  /* First, update adjacency list */
  ADJAC_Remove_neighbor(tid, t->vertex[0], fasupt);
  ADJAC_Remove_neighbor(tid, t->vertex[1], fasupt);
  ADJAC_Remove_neighbor(tid, t->vertex[2], fasupt);

  /* Now, update front list */
  SAFT3D_Remove_front_face(tid, nshrdnds, fasupt, front_heap, reject_heap, front);

  /* Caso usei PONTO DA FRONTEIRA */
  if (newtet->vertex[3] != idealptid) {

    if (nshrdnds[0][0] > 2) { /* Face pré-existente */
      result = HEAP_Remove_node(nshrdnds[0][7], front_heap);
      if (result != TRUE) {
        result = HEAP_Remove_node(nshrdnds[0][7], reject_heap);
      }
      assert(result == TRUE);

      ADJAC_Remove_neighbor(nshrdnds[0][7], front->element[nshrdnds[0][7] - 1].vertex[0], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[0][7], front->element[nshrdnds[0][7] - 1].vertex[1], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[0][7], front->element[nshrdnds[0][7] - 1].vertex[2], fasupt);
      SAFT3D_Remove_front_face(nshrdnds[0][7], nshrdnds, fasupt, front_heap, reject_heap, front);
    }
    else { /* Face nova */
      /* Insira no front */
      front->element[front->nels].vertex[0] = newtet->vertex[1];
      front->element[front->nels].vertex[1] = newtet->vertex[3];
      front->element[front->nels].vertex[2] = newtet->vertex[2];
      front->nels++;
      /* Insira no heap */
      heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
      HEAP_Insert(front->nels, heapkey, front_heap);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
    }

    if (nshrdnds[1][0] > 2) { /* Face pré-existente */
      result = HEAP_Remove_node(nshrdnds[1][7], front_heap);
      if (result != TRUE) {
        result = HEAP_Remove_node(nshrdnds[1][7], reject_heap);
      }
      assert(result == TRUE);

      ADJAC_Remove_neighbor(nshrdnds[1][7], front->element[nshrdnds[1][7] - 1].vertex[0], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[1][7], front->element[nshrdnds[1][7] - 1].vertex[1], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[1][7], front->element[nshrdnds[1][7] - 1].vertex[2], fasupt);
      SAFT3D_Remove_front_face(nshrdnds[1][7], nshrdnds, fasupt, front_heap, reject_heap, front);
    }
    else { /* Face nova */
      /* Insira no front */
      front->element[front->nels].vertex[0] = newtet->vertex[0];
      front->element[front->nels].vertex[1] = newtet->vertex[2];
      front->element[front->nels].vertex[2] = newtet->vertex[3];
      front->nels++;
      /* Insira no heap */
      heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
      HEAP_Insert(front->nels, heapkey, front_heap);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
    }

    if (nshrdnds[2][0] > 2) { /* Face pré-existente */
      result = HEAP_Remove_node(nshrdnds[2][7], front_heap);
      if (result != TRUE) {
        result = HEAP_Remove_node(nshrdnds[2][7], reject_heap);
      }
      assert(result == TRUE);

      ADJAC_Remove_neighbor(nshrdnds[2][7], front->element[nshrdnds[2][7] - 1].vertex[0], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[2][7], front->element[nshrdnds[2][7] - 1].vertex[1], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[2][7], front->element[nshrdnds[2][7] - 1].vertex[2], fasupt);
      SAFT3D_Remove_front_face(nshrdnds[2][7], nshrdnds, fasupt, front_heap, reject_heap, front);
    }
    else { /* Face nova */
      /* Insira no front */
      front->element[front->nels].vertex[0] = newtet->vertex[0];
      front->element[front->nels].vertex[1] = newtet->vertex[3];
      front->element[front->nels].vertex[2] = newtet->vertex[1];
      front->nels++;
      /* Insira no heap */
      heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
      HEAP_Insert(front->nels, heapkey, front_heap);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
    }

    if (fasupt->node[newtet->vertex[0] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[0], front_octree, mesh);
    }
    if (fasupt->node[newtet->vertex[1] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[1], front_octree, mesh);
    }
    if (fasupt->node[newtet->vertex[2] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[2], front_octree, mesh);
    }

    if (fasupt->node[newtet->vertex[3] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[3], front_octree, mesh);
    }
  }
  else { /* ESTOU CRIANDO UM NOVO PONTO */
    mesh->nnds++;
    OCTREE_Insert(newtet->vertex[3], front_octree, mesh);

    /* Insira as três novas faces no front */
    front->element[front->nels].vertex[0] = newtet->vertex[1];
    front->element[front->nels].vertex[1] = newtet->vertex[3];
    front->element[front->nels].vertex[2] = newtet->vertex[2];
    front->nels++;
    /* Insira no heap */
    heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
    HEAP_Insert(front->nels, heapkey, front_heap);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);

    front->element[front->nels].vertex[0] = newtet->vertex[0];
    front->element[front->nels].vertex[1] = newtet->vertex[2];
    front->element[front->nels].vertex[2] = newtet->vertex[3];
    front->nels++;
    /* Insira no heap */
    heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
    HEAP_Insert(front->nels, heapkey, front_heap);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);

    front->element[front->nels].vertex[0] = newtet->vertex[0];
    front->element[front->nels].vertex[1] = newtet->vertex[3];
    front->element[front->nels].vertex[2] = newtet->vertex[1];
    front->nels++;
    /* Insira no heap */
    heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
    HEAP_Insert(front->nels, heapkey, front_heap);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
  }

  /* Armazene o elemento, as faces e o novo ponto */
  mesh->element[mesh->nels] = *newtet;
  mesh->nels++;
}



void SAFT3D_Update_data_structuresII(unsigned long idealptid, unsigned long tid,
                   Triangle *t, Tetrahedron   *newtet,
                   unsigned long nshrdnds[3][8],
                   Heap *front_heap, Heap *reject_heap)
{
  int   result;
  double  heapkey;

  /* ATUALIZANDO A FRONTEIRA */
  /* Remover a face base do HEAP */
  HEAP_Remove_root(front_heap);

  /* First, update adjacency list */
  ADJAC_Remove_neighbor(tid, t->vertex[0], fasupt);
  ADJAC_Remove_neighbor(tid, t->vertex[1], fasupt);
  ADJAC_Remove_neighbor(tid, t->vertex[2], fasupt);

  /* Now, update front list */
  SAFT3D_Remove_front_face(tid, nshrdnds, fasupt, front_heap, reject_heap, front);

  /* Caso usei PONTO DA FRONTEIRA */
  if (newtet->vertex[3] != idealptid) {

    if (nshrdnds[0][0] > 2) { /* Face pré-existente */
      result = HEAP_Remove_node(nshrdnds[0][7], front_heap);
      if (result != TRUE) {
        result = HEAP_Remove_node(nshrdnds[0][7], reject_heap);
      }
      assert(result == TRUE);

      ADJAC_Remove_neighbor(nshrdnds[0][7], front->element[nshrdnds[0][7] - 1].vertex[0], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[0][7], front->element[nshrdnds[0][7] - 1].vertex[1], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[0][7], front->element[nshrdnds[0][7] - 1].vertex[2], fasupt);
      SAFT3D_Remove_front_face(nshrdnds[0][7], nshrdnds, fasupt, front_heap, reject_heap, front);
    }
    else { /* Face nova */
      /* Insira no front */
      front->element[front->nels].vertex[0] = newtet->vertex[1];
      front->element[front->nels].vertex[1] = newtet->vertex[3];
      front->element[front->nels].vertex[2] = newtet->vertex[2];
      front->nels++;
      /* Insira no heap */
      heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
      HEAP_Insert(front->nels, heapkey, reject_heap);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
    }

    if (nshrdnds[1][0] > 2) { /* Face pré-existente */
      result = HEAP_Remove_node(nshrdnds[1][7], front_heap);
      if (result != TRUE) {
        result = HEAP_Remove_node(nshrdnds[1][7], reject_heap);
      }
      assert(result == TRUE);

      ADJAC_Remove_neighbor(nshrdnds[1][7], front->element[nshrdnds[1][7] - 1].vertex[0], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[1][7], front->element[nshrdnds[1][7] - 1].vertex[1], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[1][7], front->element[nshrdnds[1][7] - 1].vertex[2], fasupt);
      SAFT3D_Remove_front_face(nshrdnds[1][7], nshrdnds, fasupt, front_heap, reject_heap, front);
    }
    else { /* Face nova */
      /* Insira no front */
      front->element[front->nels].vertex[0] = newtet->vertex[0];
      front->element[front->nels].vertex[1] = newtet->vertex[2];
      front->element[front->nels].vertex[2] = newtet->vertex[3];
      front->nels++;
      /* Insira no heap */
      heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
      HEAP_Insert(front->nels, heapkey, reject_heap);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
    }

    if (nshrdnds[2][0] > 2) { /* Face pré-existente */
      result = HEAP_Remove_node(nshrdnds[2][7], front_heap);
      if (result != TRUE) {
        result = HEAP_Remove_node(nshrdnds[2][7], reject_heap);
      }
      assert(result == TRUE);

      ADJAC_Remove_neighbor(nshrdnds[2][7], front->element[nshrdnds[2][7] - 1].vertex[0], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[2][7], front->element[nshrdnds[2][7] - 1].vertex[1], fasupt);
      ADJAC_Remove_neighbor(nshrdnds[2][7], front->element[nshrdnds[2][7] - 1].vertex[2], fasupt);
      SAFT3D_Remove_front_face(nshrdnds[2][7], nshrdnds, fasupt, front_heap, reject_heap, front);
    }
    else { /* Face nova */
      /* Insira no front */
      front->element[front->nels].vertex[0] = newtet->vertex[0];
      front->element[front->nels].vertex[1] = newtet->vertex[3];
      front->element[front->nels].vertex[2] = newtet->vertex[1];
      front->nels++;
      /* Insira no heap */
      heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
      HEAP_Insert(front->nels, heapkey, reject_heap);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
      ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
    }

    if (fasupt->node[newtet->vertex[0] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[0], front_octree, mesh);
    }
    if (fasupt->node[newtet->vertex[1] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[1], front_octree, mesh);
    }
    if (fasupt->node[newtet->vertex[2] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[2], front_octree, mesh);
    }

    if (fasupt->node[newtet->vertex[3] - 1].nnbs == 0) {
      OCTREE_Remove(newtet->vertex[3], front_octree, mesh);
    }
  }
  else { /* ESTOU CRIANDO UM NOVO PONTO */
    mesh->nnds++;
    OCTREE_Insert(newtet->vertex[3], front_octree, mesh);

    /* Insira as três novas faces no front */
    front->element[front->nels].vertex[0] = newtet->vertex[1];
    front->element[front->nels].vertex[1] = newtet->vertex[3];
    front->element[front->nels].vertex[2] = newtet->vertex[2];
    front->nels++;
    /* Insira no heap */
    heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
    HEAP_Insert(front->nels, heapkey, reject_heap);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);

    front->element[front->nels].vertex[0] = newtet->vertex[0];
    front->element[front->nels].vertex[1] = newtet->vertex[2];
    front->element[front->nels].vertex[2] = newtet->vertex[3];
    front->nels++;
    /* Insira no heap */
    heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
    HEAP_Insert(front->nels, heapkey, reject_heap);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);

    front->element[front->nels].vertex[0] = newtet->vertex[0];
    front->element[front->nels].vertex[1] = newtet->vertex[3];
    front->element[front->nels].vertex[2] = newtet->vertex[1];
    front->nels++;
    /* Insira no heap */
    heapkey = HEAP_Key(&(front->element[front->nels-1]), mesh);
    HEAP_Insert(front->nels, heapkey, reject_heap);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[0], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[1], front->nels, fasupt);
    ADJAC_Insert_neighbor(front->element[front->nels-1].vertex[2], front->nels, fasupt);
  }

  /* Armazene o elemento, as faces e o novo ponto */
  mesh->element[mesh->nels] = *newtet;
  mesh->nels++;
}

Boolean SAFT3D_Dist_cdtnode_to_triangles(double minsqrdist,
                     Point3D *pt,
                     Candidate_triangles *cdttri)
{
  Point3D *nd1, *nd2, *nd3;
  double sqrdist;
  unsigned long i, tid, nd1id, nd2id, nd3id;

  for (i = 0; i < cdttri->ntri; i++) {
    tid = cdttri->triangle[i];
    
    nd1id = front->element[tid-1].vertex[0];
    nd2id = front->element[tid-1].vertex[1];
    nd3id = front->element[tid-1].vertex[2];

    nd1 = &(mesh->node[nd1id - 1]);
    nd2 = &(mesh->node[nd2id - 1]);
    nd3 = &(mesh->node[nd3id - 1]);

    sqrdist = GEOM_Point_tri_sqrdistance(pt, nd1,
      nd2, nd3);

    if (sqrdist < minsqrdist)
      return (FALSE);
  }
  return (TRUE);
}


void SAFT3D_Get_element_size(double  *elemsize,
               Point3D *nd1, Point3D *nd2, Point3D *nd3)
{
  Point3D vec12, vec13, vec23;

  GEOM_SUB3D(vec12,(*nd2),(*nd1))
  GEOM_SUB3D(vec13,(*nd3),(*nd1))
  GEOM_SUB3D(vec23,(*nd3),(*nd2))

  (*elemsize) = 0.0e+000;
  (*elemsize) += sqrt(GEOM_DOT3D(vec12,vec12));
  (*elemsize) += sqrt(GEOM_DOT3D(vec13,vec13));
  (*elemsize) += sqrt(GEOM_DOT3D(vec23,vec23));
  (*elemsize) *= ONETHIRD;
}


Boolean SAFT3D_Dist_faces_to_cdtpts(double minsqrdist,
                  unsigned long idealcdtpos,
                  unsigned long nshrdnds[][8],
                  Point3D *nd1, Point3D *nd2,
                  Point3D *nd3, Point3D *nd4,
                  Candidate_points *cdtpts)
{
  Point3D     n, v1, v2;
  double      dpt, sqrdist;
  unsigned long i;
  Point3D     *pt;

  minsqrdist *= (SQUARED055);

  /* Face 1: 2-4-3 */
  if (nshrdnds[0][0] < 3) {
    for (i = 0; i < idealcdtpos; i++) {
      pt = &(mesh->node[cdtpts->point[i] - 1]);
      /* Verificar se o candidato i está no lado positivo para formar um novo elemento */
      /* Calcular o volume do tetraedro 2-3-4-cdtpt */
      GEOM_SUB3D(v1,(*nd2),(*nd4))
      GEOM_SUB3D(v2,(*nd3),(*nd4))
      GEOM_CROSS3D(n,v1,v2)

      GEOM_SUB3D(v1,(*pt),(*nd4))
      dpt = GEOM_DOT3D(n,v1);
      if (dpt > EPSILON) {
        sqrdist = GEOM_Point_tri_sqrdistance(pt, nd2,
          nd3, nd4);
        if (sqrdist < minsqrdist)
          return (FALSE);
      }
    }
  }
  
  /* Face 2: 1-3-4 */
  if (nshrdnds[1][0] < 3) {
    for (i = 0; i < idealcdtpos; i++) {
      pt = &(mesh->node[cdtpts->point[i] - 1]);
      /* Verificar se o candidato i está no lado positivo para formar um novo elemento */
      /* Calcular o volume do tetraedro 2-3-4-cdtpt */
      GEOM_SUB3D(v1,(*nd3),(*nd4))
      GEOM_SUB3D(v2,(*nd1),(*nd4))
      GEOM_CROSS3D(n,v1,v2)

      GEOM_SUB3D(v1,(*pt),(*nd4))
      dpt = GEOM_DOT3D(n,v1);
      if (dpt > EPSILON) {
        sqrdist = GEOM_Point_tri_sqrdistance(pt, nd1,
          nd4, nd3);
        if (sqrdist < minsqrdist)
          return (FALSE);
      }
    }
  }

  /* Face 3: 1-4-2 */
  if (nshrdnds[2][0] < 3) {
    for (i = 0; i < idealcdtpos; i++) {
      pt = &(mesh->node[cdtpts->point[i] - 1]);
      /* Verificar se o candidato i está no lado positivo para formar um novo elemento */
      /* Calcular o volume do tetraedro 2-3-4-cdtpt */
      GEOM_SUB3D(v1,(*nd1),(*nd4))
      GEOM_SUB3D(v2,(*nd2),(*nd4))
      GEOM_CROSS3D(n,v1,v2)

      GEOM_SUB3D(v1,(*pt),(*nd4))
      dpt = GEOM_DOT3D(n,v1);
      if (dpt > EPSILON) {
        sqrdist = GEOM_Point_tri_sqrdistance(pt, nd1,
          nd2, nd4);
        if (sqrdist < minsqrdist)
          return (FALSE);
      }
    }
  }

  return (TRUE);
}





double TRI_EDGE_ANGLE(Point3D *P1,Point3D *Q1,Point3D *R1,
            Point3D *P2,Point3D *Q2,Point3D *R2,
            double sqtol)
{
  Point3D m1, m2, edge, n, v1, v2, n1, n2, pt1, pt2;
  double dotn, dot1, dot2, ang1, ang2, invnmagn, nmagn;
  
  /* Initialize angles with maximum angle */
  ang1 = HALFPI;
  ang2 = HALFPI;

  /* Calcula o vetor normal do plano contendo o triângulo 01 */
  GEOM_SUB3D(v1,(*R1),(*P1))
  GEOM_SUB3D(v2,(*Q1),(*P1))  
  GEOM_CROSS3D(n,v1,v2)

  /* Calcula os vetores correspondentes às medianas do triângulo 1 */
  GEOM_ADD3D(m1,(*Q1),(*P1))
  GEOM_SCALAR3D((0.5e+000),m1)

  GEOM_ADD3D(m2,(*R1),(*P1))
  GEOM_SCALAR3D((0.5e+000),m2)    

  /* normaliza a normal ao triângulo 1 */
  nmagn = sqrt(GEOM_DOT3D(n,n));
  invnmagn = (1.0e+000)/nmagn;
  GEOM_SCALAR3D(invnmagn,n);

  /* for Q1 */
  pt1.x = nmagn*(n.x) + m1.x;
  pt1.y = nmagn*(n.y) + m1.y;
  pt1.z = nmagn*(n.z) + m1.z;

  /* for R1 */
  pt2.x = nmagn*(n.x) + m2.x;
  pt2.y = nmagn*(n.y) + m2.y;
  pt2.z = nmagn*(n.z) + m2.z;

  /* normal to edges of triangle 1 */
  GEOM_SUB3D(v1,(*Q1),(*P1))
  GEOM_SUB3D(v2,pt1,(*P1))
  GEOM_CROSS3D(n1,v1,v2)

  GEOM_SUB3D(v1,pt2,(*P1))
  GEOM_SUB3D(v2,(*R1),(*P1))    
  GEOM_CROSS3D(n2,v1,v2)

  /* Calcula o vetor correspondente a aresta P2-Q2 */
  GEOM_SUB3D(edge,(*Q2),(*P1))

  /* Checa se o nó terminal da aresta 1 está no lado correto do plano do triângulo 01 */
  dotn = GEOM_DOT3D(n,edge);
  if (dotn >= 0.0e+000) {
    /* check orientation */
    dot1 = GEOM_DOT3D(n1,edge);
    dot2 = GEOM_DOT3D(n2,edge);

    if (dot1 > (0.0e+000) && dot2 > (0.0e+000)) {
      /* calcula ângulo 1 */
      ang1 = HALFPI - GEOM_Angle3D(&n,&edge,sqtol);
    }
  }

  /* Calcula o vetor correspondente a aresta P2-R2 */
  GEOM_SUB3D(edge,(*R2),(*P1))

  /* Checa se o nó terminal da aresta 1 está no lado correto do plano do triângulo 01 */
  dotn = GEOM_DOT3D(n,edge);
  if (dotn >= 0.0e+000) {
    /* check orientation */
    dot1 = GEOM_DOT3D(n1,edge);
    dot2 = GEOM_DOT3D(n2,edge);

    if (dot1 > (0.0e+000) && dot2 > (0.0e+000)) {
      /* calcula ângulo 1 */
      ang2 = HALFPI - GEOM_Angle3D(&n,&edge,sqtol);
    }
  }
  return (GEOM_MIN2(ang1,ang2));
}

double SAFT3D_TRI_EDGE_ANGLE_INITIALPERMUTATE(Point3D *P1, Point3D *Q1, Point3D *R1,
                        Point3D *P2, Point3D *Q2, Point3D *R2,
                        unsigned long *nshrdnds, double sqtol)
{
  if (nshrdnds[4] > 0)
    return TRI_EDGE_ANGLE(P1,Q1,R1,P2,Q2,R2,sqtol);
  else if (nshrdnds[5] > 0)
    return TRI_EDGE_ANGLE(P1,Q1,R1,Q2,R2,P2,sqtol);
  else
    return TRI_EDGE_ANGLE(P1,Q1,R1,R2,P2,Q2,sqtol);
}

double SAFT3D_MINIMUM_TRI_EDGE_ANGLE(Point3D *p1,Point3D *q1,Point3D *r1,
                   Point3D *p2,Point3D *q2,Point3D *r2,
                   unsigned long *nshrdnds, double sqtol)
{
  if (nshrdnds[1] > 0)
    return SAFT3D_TRI_EDGE_ANGLE_INITIALPERMUTATE(p1,q1,r1,p2,q2,r2,nshrdnds,sqtol);
  else if (nshrdnds[2] > 0)
    return SAFT3D_TRI_EDGE_ANGLE_INITIALPERMUTATE(q1,r1,p1,p2,q2,r2,nshrdnds,sqtol);
  else
    return SAFT3D_TRI_EDGE_ANGLE_INITIALPERMUTATE(r1,p1,q1,p2,q2,r2,nshrdnds,sqtol);
}


double SAFT3D_MINIMUM_TRI_TRI_ANGLE(Point3D *P1,Point3D *Q1,Point3D *R1,
                  Point3D *P2,Point3D *Q2,Point3D *R2,double sqtol)
{
  Point3D n1, n2, v1, v2;

  GEOM_SUB3D(v1,(*Q1),(*P1))
  GEOM_SUB3D(v2,(*R1),(*P1))
  GEOM_CROSS3D(n1,v2,v1)

  GEOM_SUB3D(v1,(*Q2),(*P2))
  GEOM_SUB3D(v2,(*R2),(*P2))
  GEOM_CROSS3D(n2,v2,v1)

  return (PI - GEOM_Angle3D(&n1, &n2, sqtol));
}

Boolean SAFT3D_Minimum_angle_tests(double minang_fa2ed, double minang_fa2fa,
                   unsigned long shrdinfo[3][8],
                   Tetrahedron *newtet,
                   Candidate_triangles *cdttri,
                   double *minang)
{
  Point3D     *p1, *q1, *r1;
  Point3D     *p2, *q2, *r2;
  Triangle    newtri[3], oldtri;
  unsigned long nshrdnds[3][8], i, j;
  double      ang1, ang2, sqtol;

  /* Face 01: */
  newtri[0].vertex[0] = newtet->vertex[1];
  newtri[0].vertex[1] = newtet->vertex[3];
  newtri[0].vertex[2] = newtet->vertex[2];

  /* Face 02: */
  newtri[1].vertex[0] = newtet->vertex[0];
  newtri[1].vertex[1] = newtet->vertex[2];
  newtri[1].vertex[2] = newtet->vertex[3];

  /* Face 03: */
  newtri[2].vertex[0] = newtet->vertex[0];
  newtri[2].vertex[1] = newtet->vertex[3];
  newtri[2].vertex[2] = newtet->vertex[1];

  /* Sorting: put smallest index at the first node to preserve orientation */
  for (i = 0; i < 3; i++) {
    while (newtri[i].vertex[0] > newtri[i].vertex[1] ||
      newtri[i].vertex[0] > newtri[i].vertex[2]) {
      j         = newtri[i].vertex[0];
      newtri[i].vertex[0] = newtri[i].vertex[1];
      newtri[i].vertex[1] = newtri[i].vertex[2];
      newtri[i].vertex[2] = j;
    }
  }

  for (i = 0; i < 3; i++) {
    if (shrdinfo[i][0] == 3)
      continue;

    p1 = &(mesh->node[newtri[i].vertex[0] - 1]);
    q1 = &(mesh->node[newtri[i].vertex[1] - 1]);
    r1 = &(mesh->node[newtri[i].vertex[2] - 1]);

    /* Compute the minimum square root tolerance */
    sqtol = epsilon*GEOM_MAX3(GEOM_MAX3(p1->x,p1->y,p1->z),
      GEOM_MAX3(q1->x,q1->y,q1->z),GEOM_MAX3(r1->x,r1->y,r1->z));
    sqtol = sqtol*sqtol;

    for (j = 0; j < cdttri->ntri; j++) {
      nshrdnds[i][0] = 0, nshrdnds[i][1] = 0, nshrdnds[i][2] = 0;
      nshrdnds[i][3] = 0, nshrdnds[i][4] = 0, nshrdnds[i][5] = 0;
      nshrdnds[i][6] = 0, nshrdnds[i][7] = 0;

      nshrdnds[i][7] = cdttri->triangle[j];
      oldtri   = front->element[cdttri->triangle[j] - 1];

      if (newtri[i].vertex[0] == oldtri.vertex[0]) {
        nshrdnds[i][0]++, nshrdnds[i][1]++, nshrdnds[i][4]++;
      }
      else if (newtri[i].vertex[0] == oldtri.vertex[1]) {
        nshrdnds[i][0]++, nshrdnds[i][1]++, nshrdnds[i][5]++;
      }
      else {
        if (newtri[i].vertex[0] == oldtri.vertex[2])
          nshrdnds[i][0]++, nshrdnds[i][1]++, nshrdnds[i][6]++;
      }

      if (newtri[i].vertex[1] == oldtri.vertex[0]) {
        nshrdnds[i][0]++, nshrdnds[i][2]++, nshrdnds[i][4]++;
      }
      else if (newtri[i].vertex[1] == oldtri.vertex[1]) {
        nshrdnds[i][0]++, nshrdnds[i][2]++, nshrdnds[i][5]++;
      }
      else {
        if (newtri[i].vertex[1] == oldtri.vertex[2])
          nshrdnds[i][0]++, nshrdnds[i][2]++, nshrdnds[i][6]++;
      }

      if (newtri[i].vertex[2] == oldtri.vertex[0]) {
        nshrdnds[i][0]++, nshrdnds[i][3]++, nshrdnds[i][4]++;
      }
      else if (newtri[i].vertex[2] == oldtri.vertex[1]) {
        nshrdnds[i][0]++, nshrdnds[i][3]++, nshrdnds[i][5]++;
      }
      else {
        if (newtri[i].vertex[2] == oldtri.vertex[2])
          nshrdnds[i][0]++, nshrdnds[i][3]++, nshrdnds[i][6]++;
      }

      p2 = &(mesh->node[oldtri.vertex[0] - 1]);
      q2 = &(mesh->node[oldtri.vertex[1] - 1]);
      r2 = &(mesh->node[oldtri.vertex[2] - 1]);
      
      if (nshrdnds[i][0] == 0) {
        continue;
      }
      else if (nshrdnds[i][0] == 2) { /* Testar o ângulo entre as duas faces */
        ang1 = SAFT3D_MINIMUM_TRI_TRI_ANGLE(p1,q1,r1,p2,q2,r2,sqtol);
        *minang = GEOM_MIN2(ang1,(*minang));
        if (ang1 < minang_fa2fa) {
          return (FALSE);
        }
        else
          continue;
      }
      else {  /* Testar o ângulo entre a face nova e a aresta */
        ang2 = SAFT3D_MINIMUM_TRI_EDGE_ANGLE(p1,q1,r1,p2,q2,r2,nshrdnds[i],sqtol);
        *minang = GEOM_MIN2(ang2,(*minang));
        if (ang2 < minang_fa2ed) {
          return (FALSE);
        }
        else
          continue;       
      }
    }
  }
  return (TRUE);
}

Boolean SAFT3D_Dist_check_facestonodes(double minsqrdist,
                     unsigned long nd4pos,
                     unsigned long idealptpos,
                     unsigned long nshrdnds[][8],
                     Point3D *nd1, Point3D *nd2,
                     Point3D *nd3, Point3D *nd4,
                     Candidate_points *cdtpts,
                     double *mindist1)
{
  Point3D     centroid, n, v1, v2, vec, *pt;
  unsigned long i, result;
  double  dist;

  result = TRUE;
  if (cdtpts->npts > 1) {
    /*minsqrdist *= (SQUARED055);*/

    if (nd4pos != idealptpos) {
      GEOM_SWAP(cdtpts->point[nd4pos],
        cdtpts->point[idealptpos - 1],i)
    }

    if (nshrdnds[0][0] < 3) {
      /* Calculate triangle centroid */
      GEOM_TRIANGLE_CENTROID3D(centroid,(*nd2),(*nd4),(*nd3));

      GEOM_SUB3D(v1,(*nd3),(*nd2));
      GEOM_SUB3D(v2,(*nd4),(*nd2));
      GEOM_CROSS3D(n,v1,v2);

      for (i = 0; i < idealptpos - 1; i++){
        pt = &(mesh->node[cdtpts->point[i] - 1]);

        /* checar se o ponto está do lado correto */
        GEOM_SUB3D(v1,(*pt),(*nd2));
        if (GEOM_DOT3D(n,v1) < (0.0e+000))
          continue;

        GEOM_SUB3D(vec,(*pt),centroid);

        dist = GEOM_DOT3D(vec,vec);
        *mindist1 = GEOM_MIN2(dist,(*mindist1));
        if (dist < minsqrdist) {
          result = FALSE;
          break;
        }
      }
    }

    if (nshrdnds[1][0] < 3) {
      /* Calculate triangle centroid */
      GEOM_TRIANGLE_CENTROID3D(centroid,(*nd1),(*nd3),(*nd4));

      GEOM_SUB3D(v1,(*nd4),(*nd1));
      GEOM_SUB3D(v2,(*nd3),(*nd1));
      GEOM_CROSS3D(n,v1,v2);

      for (i = 0; i < idealptpos - 1; i++){
        pt = &(mesh->node[cdtpts->point[i] - 1]);

        /* checar se o ponto está do lado correto */
        GEOM_SUB3D(v1,(*pt),(*nd1));
        if (GEOM_DOT3D(n,v1) < (0.0e+000))
          continue;

        GEOM_SUB3D(vec,(*pt),centroid);

        dist = GEOM_DOT3D(vec,vec);
        *mindist1 = GEOM_MIN2(dist,(*mindist1));
        if (dist < minsqrdist) {
          result = FALSE;
          break;
        }
      }
    }

    if (nshrdnds[2][0] < 3) {
      /* Calculate triangle centroid */
      GEOM_TRIANGLE_CENTROID3D(centroid,(*nd1),(*nd4),(*nd2));

      GEOM_SUB3D(v1,(*nd2),(*nd1));
      GEOM_SUB3D(v2,(*nd4),(*nd1));
      GEOM_CROSS3D(n,v1,v2);

      for (i = 0; i < idealptpos - 1; i++){
        pt = &(mesh->node[cdtpts->point[i] - 1]);

        /* checar se o ponto está do lado correto */
        GEOM_SUB3D(v1,(*pt),(*nd1));
        if (GEOM_DOT3D(n,v1) < (0.0e+000))
          continue;

        GEOM_SUB3D(vec,(*pt),centroid);

        dist = GEOM_DOT3D(vec,vec);
        *mindist1 = GEOM_MIN2(dist,(*mindist1));
        if (dist < minsqrdist) {
          result = FALSE;
          break;
        }
      }
    }

    if (nd4pos != idealptpos) {
      GEOM_SWAP(cdtpts->point[nd4pos],
        cdtpts->point[idealptpos - 1],i)
    }
  }

  return (result);
}



double SAFT3D_PointToSegment_SqrDistance(const Point3D *pt,
                     const Point3D *p0,
                     const Point3D *p1)
{
    Point3D diff, direct;
    double  t, sqrlen;

  GEOM_SUB3D(diff,(*pt),(*p0));

  GEOM_SUB3D(direct,(*p1),(*p0));

  t = GEOM_DOT3D(diff,direct);

    if (t > (0.0e+000)) {
        sqrlen = GEOM_DOT3D(direct,direct);
        if (t >= sqrlen) {
      GEOM_SUB3D(diff,diff,direct);
        }
        else {
            t /= sqrlen;
      GEOM_SCALAR3D(t,direct);
      GEOM_SUB3D(diff,diff,direct);
        } 
    }
  return (GEOM_DOT3D(diff,diff));
}

Boolean SAFT3D_Dist_check_edgestonodes(double minsqrdist,
                     unsigned long nd4pos,
                     unsigned long idealptpos,
                     unsigned long nshrdnds[][8],
                     Point3D *nd1, Point3D *nd2,
                     Point3D *nd3, Point3D *nd4,
                     Candidate_points *cdtpts,
                     double *mindist2)
{
  unsigned long i, result;
  Point3D     *pt;
  double      dist;

  result = TRUE;
  if (cdtpts->npts > 1) {
    /*minsqrdist *= (SQUARED055);*/

    if (nd4pos != idealptpos) {
      GEOM_SWAP(cdtpts->point[nd4pos],
        cdtpts->point[idealptpos - 1],i)
    }
    
    if (nshrdnds[0][0] < 3) {
      if (nshrdnds[1][0] < 3) {
        if (nshrdnds[2][0] < 3) {
          /* distâncias entre as três arestas novas e os pontos vizinhos */
          for (i = 0; i < idealptpos - 1; i++) {
            pt = &(mesh->node[cdtpts->point[i] - 1]);
            dist = SAFT3D_PointToSegment_SqrDistance(pt,nd1,nd4);
            *mindist2 = GEOM_MIN2(dist,(*mindist2));
            if (dist < minsqrdist) {
              result = FALSE;
              break;
            }
            dist = SAFT3D_PointToSegment_SqrDistance(pt,nd2,nd4);
            *mindist2 = GEOM_MIN2(dist,(*mindist2));
            if (dist < minsqrdist) {
              result = FALSE;
              break;
            }
            dist = SAFT3D_PointToSegment_SqrDistance(pt,nd3,nd4);
            *mindist2 = GEOM_MIN2(dist,(*mindist2));
            if (dist < minsqrdist) {
              result = FALSE;
              break;
            }
          }
        }
        else {
          /* distâncias entre a aresta nd3-nd4 e os pontos vizinhos */
          for (i = 0; i < idealptpos - 1; i++) {
            pt = &(mesh->node[cdtpts->point[i] - 1]);
            dist = SAFT3D_PointToSegment_SqrDistance(pt,nd3,nd4);
            *mindist2 = GEOM_MIN2(dist,(*mindist2));            
            if (dist < minsqrdist) {
              result = FALSE;
              break;
            }
          }
        }
      }
      else {
        if (nshrdnds[2][0] < 3) {
          /* distâncias entre a aresta nd2-nd4 e os pontos vizinhos */
          for (i = 0; i < idealptpos - 1; i++) {
            pt = &(mesh->node[cdtpts->point[i] - 1]);
            dist = SAFT3D_PointToSegment_SqrDistance(pt,nd2,nd4);
            *mindist2 = GEOM_MIN2(dist,(*mindist2));            
            if (dist < minsqrdist) {
              result = FALSE;
              break;
            }
          }
        }
        /* else {
          não há arestas novas 
        }*/
      }
    }
    else {
      if (nshrdnds[1][0] < 3) {
        if (nshrdnds[2][0] < 3) {
          /* distâncias entre a aresta nd1-nd4 e os pontos vizinhos */
          for (i = 0; i < idealptpos - 1; i++) {
            pt = &(mesh->node[cdtpts->point[i] - 1]);
            dist = SAFT3D_PointToSegment_SqrDistance(pt,nd1,nd4);
            *mindist2 = GEOM_MIN2(dist,(*mindist2));
            if (dist < minsqrdist) {
              result = FALSE;
              break;
            }
          }
        }
        /* else {
          não há arestas novas
        }*/
      }
      /* else {
        if (nshrdnds[2][0] < 3) {
          não há arestas novas
        }
        else {
          não há arestas novas
        }
      }*/
    }

    if (nd4pos != idealptpos) {
      GEOM_SWAP(cdtpts->point[nd4pos],
        cdtpts->point[idealptpos - 1],i)
    }
  }

  return (result);
}


double SAFT3D_SqrDistance_SegmToSegm(const Point3D *p0, const Point3D *q0,
                   const Point3D *p1, const Point3D *q1)
{
    Point3D kDiff, kSeg0Dir, kSeg1Dir;
    double fA00, fA01, fA11;
    double fB0, fC, fDet;
    double fB1, fS, fT, fSqrDist, fTmp;
  double fInvDet;

  GEOM_SUB3D(kDiff,(*p0),(*p1));
  GEOM_SUB3D(kSeg0Dir,(*q0),(*p0));
  GEOM_SUB3D(kSeg1Dir,(*q1),(*p1));
  fA00  = GEOM_DOT3D(kSeg0Dir,kSeg0Dir);
  fA01  = -GEOM_DOT3D(kSeg0Dir,kSeg1Dir);
    fA11  = GEOM_DOT3D(kSeg1Dir,kSeg1Dir);
  fB0   = GEOM_DOT3D(kDiff,kSeg0Dir);
    fC    = GEOM_DOT3D(kDiff,kDiff);
    fDet  = GEOM_ABSDOUBLE(fA00*fA11-fA01*fA01);

    if (fDet >= EPSILON) {
    /* line segments are not parallel */
        fB1 = - GEOM_DOT3D(kDiff,kSeg1Dir);
        fS  = fA01*fB1-fA11*fB0;
        fT  = fA01*fB0-fA00*fB1;

        if (fS >= (0.0e+000)) {
            if (fS <= fDet) {
                if (fT >= (0.0e+000)) {
                    if (fT <= fDet) { /* region 0 (interior) */
                        /* minimum at two interior points of 3D lines */
                        fInvDet = (1.0e+000)/fDet;
                        fS *= fInvDet;
                        fT *= fInvDet;
                        fSqrDist = fS*(fA00*fS + fA01*fT + (2.0e+000)*fB0) +
              fT*(fA01*fS + fA11*fT + (2.0e+000)*fB1) + fC;
                    }
                    else { /* region 3 (side) */
                        fT   = 1.0e+000;
                        fTmp = fA01 + fB0;
                        if (fTmp >= (0.0e+000)) {
              fS = 0.0e+000;
                            fSqrDist = fA11 + (2.0e+000)*fB1 + fC;
                        }
                        else if ((-fTmp) >= fA00) {
                            fS = 1.0e+000;
                            fSqrDist = fA00 + fA11 + fC + (2.0e+000)*(fB1 + fTmp);
                        }
                        else {
                            fS = -fTmp/fA00;
                            fSqrDist = fTmp*fS + fA11 + (2.0e+000)*fB1 + fC;
                        }
                    }
                }
                else { /* region 7 (side) */
                    fT = 0.0e+000;
                    if (fB0 >= (0.0e+000)) {
                        fS = 0.0e+000;
                        fSqrDist = fC;
                    }
                    else if ((-fB0) >= fA00) {
                        fS = 1.0e+000;
                        fSqrDist = fA00 + (2.0e+000)*fB0 + fC;
                    }
                    else {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS + fC;
                    }
                }
            }
            else {
                if (fT >= (0.0e+000)) {
                    if ( fT <= fDet ) { /* region 1 (side) */
                        fS = 1.0e+000;
                        fTmp = fA01 + fB1;
                        if (fTmp >= (0.0e+000)) {
                            fT = 0.0e+000;
                            fSqrDist = fA00 + (2.0e+000)*fB0 + fC;
                        }
                        else if ((-fTmp) >= fA11) {
                            fT = 1.0e+000;
                            fSqrDist = fA00 + fA11 + fC + (2.0e+000)*(fB0 + fTmp);
                        }
                        else {
                            fT = -fTmp/fA11;
                            fSqrDist = fTmp*fT + fA00 + (2.0e+000)*fB0 + fC;
                        }
                    }
                    else { /* region 2 (corner) */
                        fTmp = fA01 + fB0;
                        if ((-fTmp) <= fA00) {
                            fT = 1.0e+000;
                            if (fTmp >= (0.0e+000)) {
                                fS = 0.0e+000;
                                fSqrDist = fA11 + (2.0e+000)*fB1 + fC;
                            }
                            else {
                                 fS = -fTmp/fA00;
                                 fSqrDist = fTmp*fS + fA11 + (2.0e+000)*fB1 + fC;
                            }
                        }
                        else {
                            fS = 1.0e+000;
                            fTmp = fA01 + fB1;
                            if (fTmp >= 0.0e+000) {
                                fT = 0.0e+000;
                                fSqrDist = fA00 + (2.0e+000)*fB0 + fC;
                            }
                            else if ((-fTmp) >= fA11) {
                                fT = 1.0e+000;
                                fSqrDist = fA00 + fA11 + fC +
                  (2.0e+000)*(fB0 + fTmp);
                            }
                            else {
                                fT = -fTmp/fA11;
                                fSqrDist = fTmp*fT + fA00 + (2.0e+000)*fB0 + fC;
                            }
                        }
                    }
                }
                else { /* region 8 (corner) */
                    if ((-fB0) < fA00) {
                        fT = 0.0e+000;
                        if (fB0 >= (0.0e+000)) {
                            fS = 0.0e+000;
                            fSqrDist = fC;
                        }
                        else {
                            fS = -fB0/fA00;
                            fSqrDist = fB0*fS + fC;
                        }
                    }
                    else {
                        fS = 1.0e+000;
                        fTmp = fA01 + fB1;
                        if (fTmp >= (0.0e+000)) {
                            fT = 0.0e+000;
                            fSqrDist = fA00 + (2.0e+000)*fB0 + fC;
                        }
                        else if ((-fTmp) >= fA11) {
                            fT = 1.0e+000;
                            fSqrDist = fA00 + fA11 + fC +
                (2.0e+000)*(fB0 + fTmp);
                        }
                        else {
                            fT = -fTmp/fA11;
                            fSqrDist = fTmp*fT + fA00 + (2.0e+000)*fB0 + fC;
                        }
                    }
                }
            }
        }
        else {
            if (fT >= (0.0e+000)) {
                if (fT <= fDet) { /* region 5 (side) */
                    fS = 0.0e+000;
                    if (fB1 >= (0.0e+000)) {
                        fT = 0.0e+000;
                        fSqrDist = fC;
                    }
                    else if ((-fB1) >= fA11) {
                        fT = 1.0e+000;
                        fSqrDist = fA11 + (2.0e+000)*fB1 + fC;
                    }
                    else {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT + fC;
                    }
                }
                else { /* region 4 (corner) */
                    fTmp = fA01 + fB0;
                    if (fTmp < 0.0e+000) {
                        fT = 1.0e+000;
                        if ((-fTmp) >= fA00) {
                            fS = 1.0e+000;
                            fSqrDist = fA00 + fA11 + fC +
                (2.0e+000)*(fB1 + fTmp);
                        }
                        else {
                            fS = -fTmp/fA00;
                            fSqrDist = fTmp*fS + fA11 + (2.0e+000)*fB1 + fC;
            }
                    }
                    else {
                        fS = 0.0e+000;
                        if (fB1 >= (0.0e+000)) {
                            fT = 0.0e+000;
                            fSqrDist = fC;
                        }
                        else if ((-fB1) >= fA11) {
                            fT = 1.0e+000;
                            fSqrDist = fA11 + (2.0e+000)*fB1 + fC;
                        }
                        else {
                            fT = -fB1/fA11;
                            fSqrDist = fB1*fT + fC;
                        }
                    }
                }
            }
            else { /* region 6 (corner) */
                if (fB0 < (0.0e+000)) {
          fT = 0.0e+000;
                    if ((-fB0) >= fA00) {
                        fS = 1.0e+000;
                        fSqrDist = fA00 + (2.0e+000)*fB0 + fC;
                    }
                    else {
                        fS = -fB0/fA00;
                        fSqrDist = fB0*fS + fC;
                    }
                }
                else {
                    fS = 0.0e+000;
                    if (fB1 >= (0.0e+000)) {
            fT = 0.0e+000;
                        fSqrDist = fC;
                    }
                    else if ((-fB1) >= fA11) {
                        fT = 1.0e+000;
                        fSqrDist = fA11 + (2.0e+000)*fB1 + fC;
                    }
                    else {
                        fT = -fB1/fA11;
                        fSqrDist = fB1*fT + fC;
                    }
                }
            }
        }
    }
    else {
        /* line segments are parallel */
        if (fA01 > (0.0e+000)) {
            /* direction vectors form an obtuse angle */
            if (fB0 >= (0.0e+000)) {
                fS = 0.0e+000;
                fT = 0.0e+000;
                fSqrDist = fC;
            }
            else if ((-fB0) <= fA00) {
        fS = -fB0/fA00;
                fT = 0.0e+000;
                fSqrDist = fB0*fS + fC;
            }
            else {
        fB1 = -GEOM_DOT3D(kDiff,kSeg1Dir);
                fS = 1.0e+000;
                fTmp = fA00 + fB0;
                if ((-fTmp) >= fA01) {
                    fT = 1.0e+000;
                    fSqrDist = fA00 + fA11 + fC +
            (2.0e+000)*(fA01 + fB0 + fB1);
                }
                else {
                    fT = -fTmp/fA01;
                    fSqrDist = fA00 + (2.0e+000)*fB0 + fC +
            fT*(fA11*fT + (2.0e+000)*(fA01 + fB1));
                }
            }
        }
        else {
            /* direction vectors form an acute angle */
            if ((-fB0) >= fA00) {
                fS = 1.0e+000;
                fT = 0.0e+000;
                fSqrDist = fA00 + (2.0e+000)*fB0 + fC;
            }
            else if (fB0 <= (0.0e+000)) {
                fS = -fB0/fA00;
                fT = 0.0e+000;
                fSqrDist = fB0*fS + fC;
            }
            else {
        fB1 = -GEOM_DOT3D(kDiff,kSeg1Dir);
                fS = 0.0e+000;
                if (fB0 >= (-fA01)) {
                    fT = 1.0e+000;
                    fSqrDist = fA11 + (2.0e+000)*fB1 + fC;
                }
                else {
                    fT = -fB0/fA01;
                    fSqrDist = fC + fT*((2.0e+000)*fB1 + fA11*fT);
                }
            }
        }
    }

    return (GEOM_ABSDOUBLE(fSqrDist));
}


#define SAFT3D_EDGE_SHARING_TEST(COUNTER,P1,Q1,P2,Q2) { \
  if (P1 == P2 || P1 == Q2) \
    COUNTER++; \
  else \
    if (Q1 == P2 || Q1 == Q2) \
      COUNTER++; \
}

Boolean SAFT3D_Dist_check_edgestoedges(double minsqrdist,
                     unsigned long nshrdnds[][8],
                     Tetrahedron *newtet,
                     Point3D *nd1, Point3D *nd2,
                     Point3D *nd3, Point3D *nd4,
                     Candidate_triangles *cdttri,
                     double *mindist3)
{
  Triangle    *oldtri;
  unsigned long i, counter;
  Point3D     *pt1, *pt2, *pt3;
  double      dist;

/*  minsqrdist *= (SQUARED030);*/

  if (nshrdnds[0][0] < 3) {
    if (nshrdnds[1][0] < 3) {
      if (nshrdnds[2][0] < 3) {
        /* distâncias entre as três arestas novas e as faces vizinhas */
        for (i = 0; i < cdttri->ntri; i++) {
          oldtri = &(front->element[cdttri->triangle[i] - 1]);

          pt1 = &(mesh->node[oldtri->vertex[0] - 1]);
          pt2 = &(mesh->node[oldtri->vertex[1] - 1]);
          pt3 = &(mesh->node[oldtri->vertex[2] - 1]);
          
          /* NOVA ARESTA NO. 1 */
          counter = 0;
          /* Nova aresta nd1-nd4 com aresta pt1-pt2 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[0],
            newtet->vertex[3],oldtri->vertex[0],oldtri->vertex[1]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd1, nd4, pt1, pt2);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }
          
          /* Nova aresta nd1-nd4 com aresta pt2-pt3 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[0],
            newtet->vertex[3],oldtri->vertex[1],oldtri->vertex[2]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd1, nd4, pt2, pt3);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }
          
          /* Nova aresta nd1-nd4 com aresta pt3-pt1 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[0],
            newtet->vertex[3],oldtri->vertex[2],oldtri->vertex[0]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd1, nd4, pt3, pt1);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* NOVA ARESTA NO. 2 */
          /* Nova aresta nd2-nd4 com aresta pt1-pt2 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[1],
            newtet->vertex[3],oldtri->vertex[0],oldtri->vertex[1]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd2, nd4, pt1, pt2);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd2-nd4 com aresta pt2-pt3 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[1],
            newtet->vertex[3],oldtri->vertex[1],oldtri->vertex[2]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd2, nd4, pt2, pt3);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd2-nd4 com aresta pt3-pt1 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[1],
            newtet->vertex[3],oldtri->vertex[2],oldtri->vertex[0]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd2, nd4, pt3, pt1);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* NOVA ARESTA NO. 3 */
          /* Nova aresta nd3-nd4 com aresta pt1-pt2 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[2],
            newtet->vertex[3],oldtri->vertex[0],oldtri->vertex[1]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd3, nd4, pt1, pt2);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd3-nd4 com aresta pt2-pt3 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[2],
            newtet->vertex[3],oldtri->vertex[1],oldtri->vertex[2]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd3, nd4, pt2, pt3);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd3-nd4 com aresta pt3-pt1 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[2],
            newtet->vertex[3],oldtri->vertex[2],oldtri->vertex[0]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd3, nd4, pt3, pt1);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
        }
      }
      else {
        /* distâncias entre a aresta nd3-nd4 e as faces vizinhas */
        /* NOVA ARESTA NO. 3 */
        for (i = 0; i < cdttri->ntri; i++) {
          oldtri = &(front->element[cdttri->triangle[i] - 1]);

          pt1 = &(mesh->node[oldtri->vertex[0] - 1]);
          pt2 = &(mesh->node[oldtri->vertex[1] - 1]);
          pt3 = &(mesh->node[oldtri->vertex[2] - 1]);

          counter = 0;
          /* Nova aresta nd3-nd4 com aresta pt1-pt2 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[2],
            newtet->vertex[3],oldtri->vertex[0],oldtri->vertex[1]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd3, nd4, pt1, pt2);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd3-nd4 com aresta pt2-pt3 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[2],
            newtet->vertex[3],oldtri->vertex[1],oldtri->vertex[2]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd3, nd4, pt2, pt3);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd3-nd4 com aresta pt3-pt1 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[2],
            newtet->vertex[3],oldtri->vertex[2],oldtri->vertex[0]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd3, nd4, pt3, pt1);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
        }
      }
    }
    else {
      if (nshrdnds[2][0] < 3) {
        /* distâncias entre a aresta nd2-nd4 e as faces vizinhas */
        for (i = 0; i < cdttri->ntri; i++) {
          oldtri = &(front->element[cdttri->triangle[i] - 1]);

          pt1 = &(mesh->node[oldtri->vertex[0] - 1]);
          pt2 = &(mesh->node[oldtri->vertex[1] - 1]);
          pt3 = &(mesh->node[oldtri->vertex[2] - 1]);

          /* NOVA ARESTA NO. 2 */
          counter = 0;
          /* Nova aresta nd2-nd4 com aresta pt1-pt2 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[1],
            newtet->vertex[3],oldtri->vertex[0],oldtri->vertex[1]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd2, nd4, pt1, pt2);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd2-nd4 com aresta pt2-pt3 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[1],
            newtet->vertex[3],oldtri->vertex[1],oldtri->vertex[2]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd2, nd4, pt2, pt3);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd2-nd4 com aresta pt3-pt1 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[1],
            newtet->vertex[3],oldtri->vertex[2],oldtri->vertex[0]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd2, nd4, pt3, pt1);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
        }
      }
      else {
        /* não há arestas novas */
        return (TRUE);
      }
    }
  }
  else {
    if (nshrdnds[1][0] < 3) {
      if (nshrdnds[2][0] < 3) {
        /* distâncias entre a aresta nd1-nd4 e as faces vizinhas */
        for (i = 0; i < cdttri->ntri; i++) {
          oldtri = &(front->element[cdttri->triangle[i] - 1]);

          pt1 = &(mesh->node[oldtri->vertex[0] - 1]);
          pt2 = &(mesh->node[oldtri->vertex[1] - 1]);
          pt3 = &(mesh->node[oldtri->vertex[2] - 1]);
          
          /* NOVA ARESTA NO. 1 */
          counter = 0;
          /* Nova aresta nd1-nd4 com aresta pt1-pt2 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[0],
            newtet->vertex[3],oldtri->vertex[0],oldtri->vertex[1]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd1, nd4, pt1, pt2);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd1-nd4 com aresta pt2-pt3 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[0],
            newtet->vertex[3],oldtri->vertex[1],oldtri->vertex[2]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd1, nd4, pt2, pt3);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
          else {
            counter = 0;
          }

          /* Nova aresta nd1-nd4 com aresta pt3-pt1 da face vizinha */
          SAFT3D_EDGE_SHARING_TEST(counter,newtet->vertex[0],
            newtet->vertex[3],oldtri->vertex[2],oldtri->vertex[0]);
          if (counter == 0) { /* se não compartilham nó */
            dist = SAFT3D_SqrDistance_SegmToSegm(nd1, nd4, pt3, pt1);
            *mindist3 = GEOM_MIN2(dist,(*mindist3));
            if (dist < minsqrdist)
              return (FALSE);
          }
        }
      }
      else {
        /* não há arestas novas */
        return (TRUE);
      }
    }
    else {
      if (nshrdnds[2][0] < 3) {
        /* não há arestas novas */
        return (TRUE);
      }
      else {
        /* não há arestas novas */
        return (TRUE);
      }
    }
  }

  return (TRUE);
}

void SAFT3D_Front_optimum(double elemsize,
              double minidealdist,
              double minvolrate,
              double maxvolrate,
              double minquality,
              double mindist_fa2nd,
              double mindist_ed2nd,
              double mindist_ed2ed,
              double minang_fa2ed,
              double minang_fa2fa,
              Boolean *tets_created)
{
  Candidate_points  cdtpts;
  Candidate_triangles cdttri;
  Tetrahedron     newtet;
  Sphere        range;
  Point3D       idealpt, n, v1, v2;
  unsigned long   i, tpos, bestptid, idealptid, idealcdtpos, NCDTPTS, TEMP;
  unsigned long   nshrdnds[3][8], bp_nshrdnds[3][8];
  double        oldqty, newqty, bestqty;
  double        vol, minvol, maxvol;
  double        mindist1, mindist2, mindist3, minang;
  Triangle      *t;
  Point3D       *nd1,   *nd2,   *nd3,   *nd4;

  while (front_heap->nnds > 0) {
    tpos = HEAP_Get_root(front_heap) - 1;
    t  = &(front->element[tpos]);

    newtet.vertex[0] = t->vertex[0];
    newtet.vertex[1] = t->vertex[2];
    newtet.vertex[2] = t->vertex[1];
    nd1 = &(mesh->node[t->vertex[0] - 1]);
    nd2 = &(mesh->node[t->vertex[2] - 1]);
    nd3 = &(mesh->node[t->vertex[1] - 1]);

/*OK*/  get_ideal_point_3(elemsize, nd1, nd2, nd3, &idealpt);

    cdtpts.npts = 0;
/*OK*/  SAFT3D_Get_range(&idealpt, nd1, nd2, nd3, &range);

    while (cdtpts.npts < 4) {
      OCTREE_Range_search(&range, &cdtpts, front_octree, mesh);
      range.radius *= 1.5;
    }

/*OK*/  SAFT3D_Get_neighbor_faces(&cdtpts, &cdttri, fasupt);
    
/*OK*/  SAFT3D_Remove_selected_face(tpos + 1, &cdttri);

/*OK*/  SAFT3D_Remove_face_vertices(t,&cdtpts);

    if (cdtpts.npts == 0) {
      HEAP_Insert(front_heap->node[1].id, (0.0e+000), reject_heap);
      HEAP_Remove_root(front_heap);
      continue;
    }

    idealptid   = mesh->nnds + 1;
    idealcdtpos = cdtpts.npts;
    if (SAFT3D_Dist_cdtnode_to_triangles(minidealdist,
      &idealpt, &cdttri) == TRUE) {
      cdtpts.point[cdtpts.npts] = idealptid;
      cdtpts.npts++;
      mesh->node[mesh->nnds] = idealpt;
    }

    /* Calculate ideal volume */
    GEOM_SUB3D(v1,(*nd2),(*nd1))
    GEOM_SUB3D(v2,(*nd3),(*nd1))
    GEOM_CROSS3D(n,v1,v2)
    GEOM_SUB3D(v1,idealpt,(*nd1))
    /* Check if the volume is greater than minimum */
    vol = GEOM_ABSDOUBLE(GEOM_DOT3D(n,v1));
    maxvol = maxvolrate * vol;
    minvol = minvolrate * vol;
    assert(minvol > DBL_EPSILON);

    bestptid = 0;
    oldqty   = 0.0e+000;

    /* PONTO IDEAL NÃO FOI CRIADO */
    if (idealcdtpos == cdtpts.npts) {
      NCDTPTS = cdtpts.npts;
    }
    else { /* PONTO IDEAL CRIADO */
      NCDTPTS = cdtpts.npts - 1;
    }
    
    
    for (i = 0; i < NCDTPTS; i++) {
      mindist1 = DBL_MAX;
      mindist2 = DBL_MAX;
      mindist3 = DBL_MAX;
      minang   = DBL_MAX;

      newtet.vertex[3] = cdtpts.point[i];
      nd4 = &(mesh->node[newtet.vertex[3] - 1]);
      
      /* Volume deve ser positivo */
      GEOM_SUB3D(v1,(*nd2),(*nd1))
      GEOM_SUB3D(v2,(*nd3),(*nd1))
      GEOM_CROSS3D(n,v1,v2)
      GEOM_SUB3D(v1,(*nd4),(*nd1))
      vol = GEOM_DOT3D(n,v1);
      if (vol < minvol || vol > maxvol)
        continue;
      assert(vol > DBL_EPSILON);

      /* Teste topológico local */
/*OK*/    newqty = QUALI_Tetrahedron_solid_angle_ratio(nd1, nd2,  nd3, nd4);
      if (newqty < minquality)
        continue;

      /* Teste geométrico Local */
/*OK*/    if (SAFT3D_Empty_tetrahedron_test(i, idealcdtpos, nd1, nd2,
      nd3, nd4, &cdtpts) != TRUE)
        continue;

      /* Teste geométrico de vizinhança */
      clockbefore_intersec = clock();
/*OK*/    if (SAFT3D_Intersection_tests(&newtet, &cdttri,
      nshrdnds) != TRUE) {
        clockafter_intersec = clock();
        cpu_time_intersec += (clockafter_intersec - clockbefore_intersec);
        continue;
      }
      
      /* Teste topológico de vizinhança por distâncias */
      if (SAFT3D_Dist_check_facestonodes(mindist_fa2nd, i, idealcdtpos,
        nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist1) != TRUE)
        continue;

      if (SAFT3D_Dist_check_edgestonodes(mindist_ed2nd, i, idealcdtpos,
        nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist2) != TRUE)
        continue;

      if (SAFT3D_Dist_check_edgestoedges(mindist_ed2ed, nshrdnds,
        &newtet, nd1, nd2, nd3, nd4, &cdttri, &mindist3) != TRUE)
        continue;

      /* Teste topológico de vizinhança por ângulos */
      if (SAFT3D_Minimum_angle_tests(minang_fa2ed, minang_fa2fa,
        nshrdnds, &newtet, &cdttri, &minang) != TRUE)
        continue;

      newqty *= (GEOM_MIN3(mindist1,mindist2,mindist3) + minang);

      if (newqty > oldqty) {
        memcpy(bp_nshrdnds[0], nshrdnds[0], 8*sizeof(unsigned long));
        memcpy(bp_nshrdnds[1], nshrdnds[1], 8*sizeof(unsigned long));
        memcpy(bp_nshrdnds[2], nshrdnds[2], 8*sizeof(unsigned long));
        bestqty  = newqty;
        bestptid = newtet.vertex[3];
        /***
        SAFT3D_Update_data_structures(idealptid, tpos+1, t, &newtet,
        nshrdnds, front_heap);  atencao, antes era bp_nshrdnds 
        break;
         ***/
      }
    }
    newtet.vertex[3] = bestptid;

    if (bestptid > 0) {
      SAFT3D_Update_data_structuresII(idealptid, tpos+1, t, &newtet,
      bp_nshrdnds, front_heap, reject_heap);
      *tets_created = TRUE;
    }
    else {
      if (idealcdtpos != cdtpts.npts) {/* ponto ideal foi criado */
        mindist1 = DBL_MAX;
        mindist2 = DBL_MAX;
        mindist3 = DBL_MAX;
        minang   = DBL_MAX;
        newtet.vertex[3] = cdtpts.point[idealcdtpos];
        nd4 = &(mesh->node[newtet.vertex[3] - 1]);
        
        /* Volume deve ser positivo */
        GEOM_SUB3D(v1,(*nd2),(*nd1))
        GEOM_SUB3D(v2,(*nd3),(*nd1))
        GEOM_CROSS3D(n,v1,v2)
        GEOM_SUB3D(v1,(*nd4),(*nd1))
        vol = GEOM_DOT3D(n,v1);
        assert(vol > DBL_EPSILON);
        TEMP = TRUE;
        if (vol < minvol)
          TEMP = FALSE;

        /* Teste topológico local */
        newqty = QUALI_Tetrahedron_mean_ratio(nd1, nd2, nd3, nd4);
        if (newqty < minquality)
          TEMP = FALSE;

        /* Teste geométrico Local */
        if (SAFT3D_Empty_tetrahedron_test(i, idealcdtpos, nd1, nd2,
        nd3, nd4, &cdtpts) != TRUE)
          TEMP = FALSE;

        /* Teste geométrico de vizinhança */
        clockbefore_intersec = clock();
        if (SAFT3D_Intersection_tests(&newtet, &cdttri,
          nshrdnds) != TRUE) {
          clockafter_intersec = clock();
          cpu_time_intersec += (clockafter_intersec - clockbefore_intersec);
          TEMP = FALSE;
        }

        /* Teste topológico de vizinhança por distâncias */
        if (SAFT3D_Dist_check_facestonodes(mindist_fa2nd, i,
          idealcdtpos, nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist1) != TRUE)
          TEMP = FALSE;

        if (SAFT3D_Dist_check_edgestonodes(mindist_ed2nd, i,
          idealcdtpos, nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist2) != TRUE)
          TEMP = FALSE;

        if (SAFT3D_Dist_check_edgestoedges(mindist_ed2ed,
          nshrdnds, &newtet, nd1, nd2, nd3, nd4, &cdttri, &mindist3) != TRUE)
          TEMP = FALSE;

        /* Teste topológico de vizinhança por ângulos */
        if (SAFT3D_Minimum_angle_tests(minang_fa2ed, minang_fa2fa,
          nshrdnds, &newtet, &cdttri, &minang) != TRUE)
          TEMP = FALSE;

        if (TEMP == TRUE) {
          SAFT3D_Update_data_structuresII(idealptid, tpos+1, t, &newtet,
          nshrdnds, front_heap, reject_heap);
          *tets_created = TRUE;
        }
        else {
          /* Deficient convergence */
          HEAP_Insert(front_heap->node[1].id, (0.0e+000), reject_heap);
          HEAP_Remove_root(front_heap);
        }
      }
      else {
        /* Deficient convergence */
        HEAP_Insert(front_heap->node[1].id, (0.0e+000), reject_heap);
        HEAP_Remove_root(front_heap);
      }
    }   
    free(cdtpts.point);
    free(cdttri.triangle);
  }
}

void SAFT3D_Front_passable(double elemsize,
              double minidealdist,
              double minvolrate,
              double maxvolrate,
              double minquality,
              double mindist_fa2nd,
              double mindist_ed2nd,
              double mindist_ed2ed,
              double minang_fa2ed,
              double minang_fa2fa,
              Boolean *tets_created)
{
  Candidate_points  cdtpts;
  Candidate_triangles cdttri;
  Tetrahedron     newtet;
  Sphere        range;
  Point3D       idealpt, n, v1, v2, vec, centroid;
  unsigned long   i, tpos, bestptid, idealptid, idealcdtpos;
  unsigned long   nshrdnds[3][8], bp_nshrdnds[3][8];
  int         j;
  double        oldqty, newqty, bestqty;
  double        vol, idealvol, minvol, maxvol;
  double        mindist1, mindist2, mindist3, minang;
  Triangle      *t;
  Point3D       *nd1,   *nd2,   *nd3,   *nd4;
  Boolean TEMP;

  while (front_heap->nnds > 0) {
    tpos = HEAP_Get_root(front_heap) - 1;
    t  = &(front->element[tpos]);

    newtet.vertex[0] = t->vertex[0];
    newtet.vertex[1] = t->vertex[2];
    newtet.vertex[2] = t->vertex[1];
    nd1 = &(mesh->node[t->vertex[0] - 1]);
    nd2 = &(mesh->node[t->vertex[2] - 1]);
    nd3 = &(mesh->node[t->vertex[1] - 1]);

    get_ideal_point_3(elemsize, nd1, nd2, nd3, &idealpt);

    cdtpts.npts = 0;
    SAFT3D_Get_range(&idealpt, nd1, nd2, nd3, &range);

    while (cdtpts.npts < 4) {
      OCTREE_Range_search(&range, &cdtpts, front_octree, mesh);
      range.radius *= 1.5;
    }

    SAFT3D_Get_neighbor_faces(&cdtpts, &cdttri, fasupt);
    SAFT3D_Remove_selected_face(tpos + 1, &cdttri);
    SAFT3D_Remove_face_vertices(t,&cdtpts);

    if (cdtpts.npts == 0) {
      HEAP_Insert(front_heap->node[1].id, (0.0e+000), reject_heap);
      HEAP_Remove_root(front_heap);
      continue;
    }

    idealptid   = mesh->nnds + 1;
    idealcdtpos = cdtpts.npts;
    cdtpts.point[cdtpts.npts] = idealptid;
    cdtpts.npts++;
    mesh->node[mesh->nnds] = idealpt;

    /* Calculate ideal volume */
    GEOM_SUB3D(v1,(*nd2),(*nd1))
    GEOM_SUB3D(v2,(*nd3),(*nd1))
    GEOM_CROSS3D(n,v1,v2)
    GEOM_SUB3D(v1,idealpt,(*nd1))
    /* Check if the volume is greater than minimum */
    idealvol = GEOM_ABSDOUBLE(GEOM_DOT3D(n,v1));
    maxvol = maxvolrate * idealvol;
    minvol = minvolrate * idealvol;
    assert(minvol > DBL_EPSILON);
    
    /* Calculate triangle centroid */
    GEOM_TRIANGLE_CENTROID3D(centroid,(*nd1),(*nd2),(*nd3));

    bestptid = 0;
    oldqty   = 0.0e+000;
    
    for (i = 0; i < cdtpts.npts - 1; i++) {
      mindist1 = DBL_MAX;
      mindist2 = DBL_MAX;
      mindist3 = DBL_MAX;
      minang   = DBL_MAX;

      newtet.vertex[3] = cdtpts.point[i];
      nd4 = &(mesh->node[newtet.vertex[3] - 1]);
      
      /* Volume deve ser positivo */
      GEOM_SUB3D(v1,(*nd2),(*nd1))
      GEOM_SUB3D(v2,(*nd3),(*nd1))
      GEOM_CROSS3D(n,v1,v2)
      GEOM_SUB3D(v1,(*nd4),(*nd1))
      vol = GEOM_DOT3D(n,v1);
      if (vol < minvol)
        continue;
      assert(vol > DBL_EPSILON);

      /* Teste topológico local */
      newqty = QUALI_Tetrahedron_mean_ratio(nd1, nd2, nd3, nd4);
      if (newqty < minquality)
        continue;

      /* Teste geométrico Local */
      if (SAFT3D_Empty_tetrahedron_test(i, idealcdtpos, nd1, nd2,
            nd3, nd4, &cdtpts) != TRUE)
        continue;

      /* Teste geométrico de vizinhança
      if (SAFT3D_Intersection_tests(&newtet, &cdttri,
            nshrdnds) != TRUE)
        continue;*/

      clockbefore_intersec = clock();
      if (SAFT3D_Intersection_tests(&newtet, &cdttri, nshrdnds) != TRUE) {
        clockafter_intersec = clock();
        cpu_time_intersec += (clockafter_intersec - clockbefore_intersec);
        continue;
      }

      /* Teste topológico de vizinhança por distâncias */
      if (SAFT3D_Dist_check_facestonodes(mindist_fa2nd, i,
        idealcdtpos, nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist1) != TRUE)
        continue;

      if (SAFT3D_Dist_check_edgestonodes(mindist_ed2nd, i,
        idealcdtpos, nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist2) != TRUE)
        continue;

      if (SAFT3D_Dist_check_edgestoedges(mindist_ed2ed,
        nshrdnds, &newtet, nd1, nd2, nd3, nd4, &cdttri, &mindist3) != TRUE)
        continue;

      /* Teste topológico de vizinhança por ângulos */
      if (SAFT3D_Minimum_angle_tests(minang_fa2ed, minang_fa2fa,
        nshrdnds, &newtet, &cdttri, &minang) != TRUE)
        continue;

      newqty *= (GEOM_MIN3(mindist1,mindist2,mindist3) + minang);

      if (newqty > oldqty) {
        memcpy(bp_nshrdnds[0], nshrdnds[0], 8*sizeof(unsigned long));
        memcpy(bp_nshrdnds[1], nshrdnds[1], 8*sizeof(unsigned long));
        memcpy(bp_nshrdnds[2], nshrdnds[2], 8*sizeof(unsigned long));
        bestqty  = newqty;
        bestptid = newtet.vertex[3];
        /***
        SAFT3D_Update_data_structures(idealptid, tpos+1, t, &newtet,
        nshrdnds, front_heap); atencao, antes era bp_nshrdnds 
        break;
         ***/
      }
    }
    newtet.vertex[3] = bestptid;

    if (bestptid > 0) {
      SAFT3D_Update_data_structuresII(idealptid, tpos+1, t, &newtet,
      bp_nshrdnds, front_heap, reject_heap);
      *tets_created = TRUE;
    }
    else {
      newtet.vertex[3] = cdtpts.point[idealcdtpos];
      nd4 = &(mesh->node[newtet.vertex[3] - 1]);
      GEOM_SUB3D(vec,idealpt,centroid);
      for (j = 4; j >= 0; j--) {
        mindist1 = DBL_MAX;
        mindist2 = DBL_MAX;
        mindist3 = DBL_MAX;
        minang   = DBL_MAX;
        
        /* Volume deve ser positivo */
        GEOM_SUB3D(v1,(*nd2),(*nd1))
        GEOM_SUB3D(v2,(*nd3),(*nd1))
        GEOM_CROSS3D(n,v1,v2)
        GEOM_SUB3D(v1,(*nd4),(*nd1))
        vol = GEOM_DOT3D(n,v1);
        assert(vol > DBL_EPSILON);
        TEMP = TRUE;
        if (vol < minvol)
          TEMP = FALSE;

        /* Teste topológico local */
        newqty = QUALI_Tetrahedron_mean_ratio(nd1, nd2, nd3, nd4);
        if (newqty < minquality)
          TEMP = FALSE;

        /* Teste geométrico Local */
        if (SAFT3D_Empty_tetrahedron_test(i, idealcdtpos, nd1, nd2,
        nd3, nd4, &cdtpts) != TRUE)
          TEMP = FALSE;

        /* Teste geométrico de vizinhança
        if (SAFT3D_Intersection_tests(&newtet, &cdttri,
        nshrdnds) != TRUE)
          TEMP = FALSE;*/

        clockbefore_intersec = clock();
        if (SAFT3D_Intersection_tests(&newtet, &cdttri,
          nshrdnds) != TRUE) {
          clockafter_intersec = clock();
          cpu_time_intersec += (clockafter_intersec - clockbefore_intersec);
          TEMP = FALSE;
        }

        if (SAFT3D_Dist_cdtnode_to_triangles(minidealdist,
          &idealpt, &cdttri) != TRUE)
          TEMP = FALSE;

        /* Teste topológico de vizinhança por distâncias */
        if (SAFT3D_Dist_check_facestonodes(mindist_fa2nd, i,
          idealcdtpos, nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist1) != TRUE)
          TEMP = FALSE;
  
        if (SAFT3D_Dist_check_edgestonodes(mindist_ed2nd, i,
          idealcdtpos, nshrdnds, nd1, nd2, nd3, nd4, &cdtpts, &mindist2) != TRUE)
          TEMP = FALSE;

        if (SAFT3D_Dist_check_edgestoedges(mindist_ed2ed,
          nshrdnds, &newtet, nd1, nd2, nd3, nd4, &cdttri, &mindist3) != TRUE)
          TEMP = FALSE;

        /* Teste topológico de vizinhança por ângulos */
        if (SAFT3D_Minimum_angle_tests(minang_fa2ed, minang_fa2fa,
          nshrdnds, &newtet, &cdttri, &minang) != TRUE)
          TEMP = FALSE;

        if (TEMP == TRUE) {
          SAFT3D_Update_data_structuresII(idealptid, tpos+1, t, &newtet,
          nshrdnds, front_heap, reject_heap);
          bestptid = idealptid;
          *tets_created = TRUE;
          break;
        }
        else {
          /* Decremente a distancia do ponto ideal a face */
          idealpt.x = centroid.x + ((0.2e+000)*j*vec.x);
          idealpt.y = centroid.y + ((0.2e+000)*j*vec.y);
          idealpt.z = centroid.z + ((0.2e+000)*j*vec.z);
          mesh->node[mesh->nnds] = idealpt;           
        }
      }

      if (bestptid < 1) {
        /* Deficient convergence */
        HEAP_Insert(front_heap->node[1].id, (0.0e+000), reject_heap);
        HEAP_Remove_root(front_heap);
      }
    }   
    free(cdtpts.point);
    free(cdttri.triangle);
  }
}

int generate_mesh()
{
  double elemsize, sqrdelemsize;
  double minidealdist, minvolrate, maxvolrate, minquality;
  double mindist_fa2nd, mindist_ed2nd;
  double mindist_ed2ed;
  double minang_fa2ed, minang_fa2fa;
  unsigned long oldntet, oldnnds;
  Boolean tets_created;
  clock_t startime, endtime;
  clock_t t0_fase1, t1_fase1;
  clock_t t0_fase2, t1_fase2;
  Heap *swap_aux;

  /* Grid the boundary surfaces */
  SAFT2D_Grid_boundary();
  SAFT2D_Normalize_coordinates(boundary->nnds, boundary->node);
  SAFT3D_Initialize_mesh();
  SAFT3D_Initialize_front();
  ADJAC_Initialize_fasupt(boundary);
  OCTREE_Initialize_front_octree(mesh);
  HEAP_Initialize_front_heap(front, mesh);

  /* Allocate a block of memory to the main heap structure */
  reject_heap = (Heap *) malloc(sizeof(Heap)); 
  assert(reject_heap != NULL);
  reject_heap->nnds = 0;
  reject_heap->maxnnds = front_heap->maxnnds;
  /* Allocate the heap array (nodes + SENTINEL) */
  reject_heap->node = (Heap_node *)
    malloc((reject_heap->maxnnds + 1)*sizeof(Heap_node));
  assert(reject_heap->node != NULL);
  /* Initialize the sentinel node */
  reject_heap->node[0].id  = 0;
  reject_heap->node[0].key = (double) SENTINEL;

  strcpy(status, "initial");
  SAFT3D_Print_front(fasupt, front, mesh);

  OCTREE_Print(front_octree);

  printf("Element size: ");
  scanf("%lf", &elemsize);

  /* Calculate parameters */
  sqrdelemsize  = elemsize*elemsize;
  minidealdist  = SQUARED067*sqrdelemsize;
  minvolrate    = 9.0e-001;
  maxvolrate    = 1.1e+000;
  minquality    = 6.0e-001;
  mindist_fa2nd = SQUARED067*sqrdelemsize;
  mindist_ed2nd = SQUARED070*sqrdelemsize;
  mindist_ed2ed = SQUARED060*sqrdelemsize;
  minang_fa2ed  = ANGLE45_INRAD;
  minang_fa2fa  = ANGLE60_INRAD;

  /* Timing */
  startime = clock();

  t0_fase1 = clock();

  fase1_nels0 = 0;

  printf("\n****************************************\n");
  printf("*************    FRONT 1    ************\n");
  printf("****************************************\n");
  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_optimum(elemsize, minidealdist, minvolrate,
        maxvolrate, minquality, mindist_fa2nd, mindist_ed2nd,
        mindist_ed2ed, minang_fa2ed, minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");

  } while (tets_created == TRUE);

  /* Calculate parameters */
  minidealdist  = SQUARED067*sqrdelemsize;
  minvolrate    = 9.0e-001;
  maxvolrate    = 1.1e+000;
  minquality    = 6.0e-001;
  mindist_fa2nd = SQUARED055*sqrdelemsize;
  mindist_ed2nd = SQUARED060*sqrdelemsize;
  mindist_ed2ed = SQUARED050*sqrdelemsize;
  minang_fa2ed  = ANGLE45_INRAD;
  minang_fa2fa  = ANGLE60_INRAD;

  printf("\n****************************************\n");
  printf("*************    FRONT 2    ************\n");
  printf("****************************************\n");
  

  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_optimum(elemsize, minidealdist, minvolrate,
        maxvolrate, minquality, mindist_fa2nd, mindist_ed2nd,
        mindist_ed2ed, minang_fa2ed, minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");

  } while (tets_created == TRUE);

  printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
  printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
  printf("   - No. of left front faces: %lu\n", front->nels);

  /* Calculate parameters */
  minidealdist  = SQUARED067*sqrdelemsize;
  minvolrate    = 8.0e-001;
  maxvolrate    = 1.2e+000;
  minquality    = 4.0e-001;
  mindist_fa2nd = SQUARED045*sqrdelemsize;
  mindist_ed2nd = SQUARED045*sqrdelemsize;
  mindist_ed2ed = SQUARED045*sqrdelemsize;
  minang_fa2ed  = ANGLE25_INRAD;
  minang_fa2fa  = ANGLE25_INRAD;

  printf("\n****************************************\n");
  printf("*************    FRONT 3    ************\n");
  printf("****************************************\n");
  
  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_optimum(elemsize, minidealdist, minvolrate,
        maxvolrate, minquality, mindist_fa2nd, mindist_ed2nd,
        mindist_ed2ed, minang_fa2ed, minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");

  } while (tets_created == TRUE);

  /* Calculate parameters */
  minidealdist  = SQUARED065*sqrdelemsize;
  minvolrate    = 7.0e-001;
  maxvolrate    = 1.3e+000;
  minquality    = 4.0e-001;
  mindist_fa2nd = SQUARED035*sqrdelemsize;
  mindist_ed2nd = SQUARED035*sqrdelemsize;
  mindist_ed2ed = SQUARED035*sqrdelemsize;
  minang_fa2ed  = ANGLE17_INRAD;
  minang_fa2fa  = ANGLE17_INRAD;

  printf("\n****************************************\n");
  printf("*************    FRONT 4    ************\n");
  printf("****************************************\n");
  
  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_optimum(elemsize, minidealdist, minvolrate,
        maxvolrate, minquality, mindist_fa2nd, mindist_ed2nd,
        mindist_ed2ed, minang_fa2ed, minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");

  } while (tets_created == TRUE);

  fase1_nels1 = mesh->nels;

  t1_fase1 = clock();

  t0_fase2 = clock();

  fase2_nels0 = fase1_nels1;

  /* Calculate parameters */
  minidealdist  = SQUARED067*sqrdelemsize;
  minvolrate    = 3.0e-001;
  maxvolrate    = 1.7e+000;
  minquality    = 5.0e-001;
  mindist_fa2nd = SQUARED025*sqrdelemsize;
  mindist_ed2nd = SQUARED025*sqrdelemsize;
  mindist_ed2ed = SQUARED025*sqrdelemsize;
  minang_fa2ed  = ANGLE17_INRAD;
  minang_fa2fa  = ANGLE17_INRAD;

  printf("\n****************************************\n");
  printf("*************    FRONT 5    ************\n");
  printf("****************************************\n");
  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_passable(elemsize, minidealdist, minvolrate,
        maxvolrate, minquality, mindist_fa2nd, mindist_ed2nd,
        mindist_ed2ed, minang_fa2ed, minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");
  } while (tets_created == TRUE);

  /*SAFT3D_Print_mesh(mesh);
  SAFT3D_Print_front(fasupt, front, mesh);
  fflush(stdin);
  getchar();*/

  /* Calculate parameters */
  minidealdist  = SQUARED035*sqrdelemsize;
  minvolrate    = 3.0e-001;
  maxvolrate    = 1.7e+000;
  minquality    = 3.0e-001;
  mindist_fa2nd = SQUARED025*sqrdelemsize;
  mindist_ed2nd = SQUARED025*sqrdelemsize;
  mindist_ed2ed = SQUARED025*sqrdelemsize;
  minang_fa2ed  = ANGLE14_INRAD;
  minang_fa2fa  = ANGLE14_INRAD;

  printf("\n****************************************\n");
  printf("*************    FRONT 6    ************\n");
  printf("****************************************\n");
  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_passable(elemsize, minidealdist, minvolrate, maxvolrate, minquality,
      mindist_fa2nd, mindist_ed2nd, mindist_ed2ed, minang_fa2ed,
      minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");
  } while (tets_created == TRUE);

  /* Calculate parameters */
  minidealdist  = SQUARED035*sqrdelemsize;
  minvolrate    = 5.0e-001;
  maxvolrate    = 1.5e+000;
  minquality    = 2.0e-001;
  mindist_fa2nd = SQUARED025*sqrdelemsize;
  mindist_ed2nd = SQUARED025*sqrdelemsize;
  mindist_ed2ed = SQUARED025*sqrdelemsize;
  minang_fa2ed  = ANGLE10_INRAD;
  minang_fa2fa  = ANGLE10_INRAD;

  printf("\n****************************************\n");
  printf("*************    FRONT 7    ************\n");
  printf("****************************************\n");

  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_passable(elemsize, minidealdist, minvolrate, maxvolrate, minquality,
      mindist_fa2nd, mindist_ed2nd, mindist_ed2ed, minang_fa2ed,
      minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");
  } while (tets_created == TRUE);
  /*SAFT3D_Print_mesh(mesh);
  SAFT3D_Print_front(fasupt, front, mesh);
  fflush(stdin);
  getchar();*/

  /* Calculate parameters */
  minidealdist  = 0.0000001;
  minvolrate    = 0.0e-001;
  maxvolrate    = 5.0e+000;
  minquality    = 0.0e+000;
  mindist_fa2nd = 0.000000;
  mindist_ed2nd = 0.000000;
  mindist_ed2ed = 0.000000;
  minang_fa2ed  = 0.000000;
  minang_fa2fa  = 0.000000;

  printf("\n****************************************\n");
  printf("*************    FRONT 8    ************\n");
  printf("****************************************\n");

  do {
    oldntet  = mesh->nels;
    oldnnds  = mesh->nnds;
    tets_created = FALSE;

    SAFT3D_Front_passable(elemsize, minidealdist, minvolrate, maxvolrate, minquality,
      mindist_fa2nd, mindist_ed2nd, mindist_ed2ed, minang_fa2ed,
      minang_fa2fa, &tets_created);

    swap_aux    = front_heap;
    front_heap  = reject_heap;
    reject_heap = swap_aux;

    printf("   - No. of generated tetrahedra: %lu\n", (mesh->nels)-oldntet);
    printf("   - No. of generated points: %lu\n", (mesh->nnds)-oldnnds);
    printf("   - No. of left front faces: %lu\n", front->nels);
    printf("\n");
  } while (tets_created == TRUE);
  
  t1_fase2 = clock();

  fase2_nels1 = mesh->nels;

  cpu_time_fase1 = (double) (t1_fase1 - t0_fase1)/CLOCKS_PER_SEC;
  cpu_time_fase2 = (double) (t1_fase2 - t0_fase2)/CLOCKS_PER_SEC;

  endtime = clock();
  cpu_time = (double) (endtime - startime)/CLOCKS_PER_SEC;


  if (front->nels > 0) {
    strcpy(status, "final");
    SAFT3D_Print_front(fasupt, front, mesh);
    printf("        Convergence deficiency!\n");
  }
  else {
    SAFT2D_Restore_coordinates(mesh->nnds, mesh->node);
    printf("     Writing tetrahedral mesh file...");
/*    SAFT3D_Print_mesh(mesh);*/
    printf(" completed.\n\n");
    printf("     Creating mesh quality report ...");
    QUALI_Print_quality_report(elemsize*deltamax,mesh);
    printf(" completed.\n");
    printf("\n****************************************\n");
    printf("************  FINAL SUMMARY  ***********\n");
    printf("****************************************\n");
    printf("   - Domain:\n");
    printf("        No. of tetrahedra: %lu\n", mesh->nels);
    printf("                in Fase 1: %lu\n", fase1_nels1 - fase1_nels0);
    printf("                in Fase 2: %lu\n", fase2_nels1 - fase2_nels0);
    printf("        No. of faces:      %lu\n", mesh->nfaces);
    printf("        No. of edges:      %lu\n", mesh->nedges);
    printf("        No. of nodes:      %lu\n", mesh->nnds);
    printf("   - Elapsed time:\n");
    printf("       Stage 1 %9.6lf sec\n", cpu_time_fase1);
    printf("       Stage 2 %9.6lf sec\n", cpu_time_fase2);
    printf("       Intersection test %9.6lf sec\n",
      cpu_time_intersec/CLOCKS_PER_SEC);
    printf("       Total   %9.6lf sec\n\n", cpu_time);
  }

  /* Finalize the heap structure  */
  HEAP_Finalize(front_heap);

  /* Finalize the front array */
  SAFT2D_Finalize();

  return 0;
}

/*
 * SAFT3D_Finalize()
 *
 * Prototype:
 *    void SAFT3D_Finalize(Domain *);
 *
 * Required header:
 *    saft3d.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to an allocated mesh structure.
 *
 * Comments:
 *    It deallocates the main mesh structure and respective sub-structures.
 *
 * Last modified: 21/06/2004.
 */
void SAFT3D_Finalize_mesh()
{
  free(mesh->node);
  mesh->node  = NULL;
  assert(mesh->node == NULL);
  
  free(mesh->element);
  mesh->element = NULL;
  assert(mesh->element == NULL);

  free(mesh);
  mesh = NULL;
  assert(mesh == NULL);
}
