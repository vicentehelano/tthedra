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
#include <saft2d.h>

/* Global Variables */
Surface_mesh *boundary;
char         projname[150];
double     deltamax;
Point3D      bbmin, bbmax, bbdelta;

void SAFT2D_Normalize_coordinates(unsigned long nnds, Point3D *node)
{
  unsigned long i;

  bbmin = node[0];
  bbmax = node[0];
  for (i = 1; i < nnds; i++) {
    if (node[i].x < bbmin.x)
      bbmin.x = node[i].x;
    if (node[i].y < bbmin.y)
      bbmin.y = node[i].y;
    if (node[i].z < bbmin.z)
      bbmin.z = node[i].z;

    if (node[i].x > bbmax.x)
      bbmax.x = node[i].x;
    if (node[i].y > bbmax.y)
      bbmax.y = node[i].y;
    if (node[i].z > bbmax.z)
      bbmax.z = node[i].z;
  }

  GEOM_SUB3D(bbdelta,bbmax,bbmin);

  deltamax = GEOM_MAX3(bbdelta.x,bbdelta.y,bbdelta.z);

  for (i = 0; i < nnds; i++) {
    node[i].x =  (node[i].x - bbmin.x) / deltamax;
    node[i].y =  (node[i].y - bbmin.y) / deltamax;
    node[i].z =  (node[i].z - bbmin.z) / deltamax;
  }
}

void SAFT2D_Restore_coordinates(unsigned long nnds, Point3D *node)
{
  unsigned long i;

  for (i = 0; i < nnds; i++) {
    node[i].x =  (node[i].x)*deltamax + bbmin.x;
    node[i].y =  (node[i].y)*deltamax + bbmin.y;
    node[i].z =  (node[i].z)*deltamax + bbmin.z;
  }
}

/*
 * SAFT2D_Initialize()
 *
 * Prototype:
 *    void SAFT2D_Initialize(Boundary **);
 *
 * Required header:
 *    saft2d.h and assert.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives the address of a pointer to a boundary structure yet not allocated.
 *
 * Comments:
 *    It allocates the main boundary structure and respective sub-structures.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Initialize_boundary()
{
  /* Allocate a block of memory to the main boundary structure */
  boundary = (Surface_mesh *) malloc(sizeof(Surface_mesh));
  assert(boundary != NULL);         /* Check if memory was allocated    */

  /* Initialize the number of nodes */
  boundary->nnds = 0;
  /* Initial maximum number of nodes ((2^15) - 1) + 1 */
  /*boundary->maxnnds = 32768;*/
  boundary->maxnnds = 150000;

  /* Allocate the array for the nodes */
  boundary->node = (Point3D *)
    malloc(boundary->maxnnds*sizeof(Point3D));
  assert(boundary->node != NULL);     /* Check if memory was allocated    */

  /* Initialize the number of elements */
  boundary->nels    = 0;
  /* Initial maximum number of elements ((2^15) - 1) + 1 */
  /*boundary->maxnels = 32768;*/
  boundary->maxnels = 150000;

  /* Allocate the array for the elements */
  boundary->element = (Triangle *)
    malloc(boundary->maxnels*sizeof(Triangle));
  assert(boundary->element != NULL);      /* Check if memory was allocated    */
}

/*
 * SAFT2D_Get_boundary_mesh()
 *
 * Prototype:
 *    void SAFT2D_Grid_boundary(Boundary *);
 *
 * Required header:
 *    saft2d.h and assert.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a boundary structure.
 *
 * Comments:
 *    Read the triangulation of the surfaces (internal and external)
 *    of some domain from a file and print the global surface mesh.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Grid_boundary()
{
  SAFT2D_Initialize_boundary();

  /* Read external surface mesh */
  SAFT2D_Read_external_surfs();

  /* Read internal surface mesh */
  SAFT2D_Read_internal_surfs();
}

/*
 * SAFT2D_Read_external_surfs()
 *
 * Prototype:
 *    void SAFT2D_Read_external_surfs(Boundary *);
 *
 * Required header:
 *    saft2d.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a boundary structure.
 *
 * Comments:
 *    It calls the SAFT2D_Read_surface_mesh() with the parameter ".esf"
 *    that indicates a valid extension of a external surface mesh file.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Read_external_surfs()
{
  SAFT2D_Read_surface_mesh(".esf");
}

/*
 * SAFT2D_Read_internal_surfs()
 *
 * Prototype:
 *    void SAFT2D_Read_internal_surfs(Boundary *);
 *
 * Required header:
 *    saft2d.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a boundary structure.
 *
 * Comments:
 *    It calls the SAFT2D_Read_surface_mesh() with the parameter ".isf"
 *    that indicates a valid extension of a internal surface mesh file.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Read_internal_surfs()
{
  SAFT2D_Read_surface_mesh(".isf");
}

/*
 * SAFT2D_Read_surface_mesh()
 *
 * Prototype:
 *    void SAFT2D_Read_surface_mesh(char*, Boundary *);
 *
 * Required header:
 *    saft2d.h, string.h, stdlib.h, stdio.h and assert.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a boundary structure and a pointer to a
 *    file extension string.
 *
 * Comments:
 *    It reads the surface meshes (internal and external).
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Read_surface_mesh(char *extension)
{
  /* Definitions of local variables */
  unsigned long i, tag;           /* Counter and a tag number       */
  unsigned long oldnnds, oldnels;     /* Old last node and element      */
  unsigned long v1, v2, v3;         /* Auxiliary vertices id        */
  char      ftype[]  = "external";    /* File type              */
  char      fname[150], fsectr[20];   /* File name and file sector      */
  FILE      *fptr;            /* File pointer             */
  /* Input data file path         */
  #ifdef NDEBUG
    char fpath[150] = "../data/input/";
  #else
    char fpath[150] = "./../../data/input/";
  #endif

  /* Check file type and change it if necessary */
  if (strcmp(extension, ".esf") != 0) {
    ftype[0] = 'i';
    ftype[1] = 'n';
  }

  /* Initialize file pointer and name */
  fptr = NULL;
  fname[0] = '\0';

  /* Read name of the file with the surface mesh */
  while (fptr == NULL) {
    printf("File with %s surface mesh: ", ftype);
    fflush(stdin);
    scanf("%s", fname);

    /* Check for 'none' option */
    if (strcmp(extension, ".isf") == 0) {
      if (strcmp(fname, "none") == 0)
        return;
      else if (strcmp((fname + (strlen(fname)-4)), ".isf") != 0) {
        printf("    Not a extension of internal surface (.isf)\n");
        continue;
      }
    }
    else {
      strcpy(projname,fname);
      projname[strlen(projname) - 4] = '\0';
    }

    /* Generate complete file path name */
    strcat(fpath, fname);

    /* Open file with the surface mesh */
    fptr = fopen(fpath, "r");

    if (fptr == NULL) {
      printf("File not found!\n");
      fpath[strlen(fpath)-strlen(fname)] = '\0';
      fname[0]  = '\0';
    }
  }
  /* Check if the file is valid */
  assert(fptr != NULL);

  /* Read old number of nodes and elements */
  oldnnds = boundary->nnds;
  oldnels = boundary->nels;

  /* Read begin of coordinate sector */
  fscanf(fptr, "%s %lu\n", fsectr, &(boundary->nnds));
  assert(strcmp(fsectr, "coor") == 0);      /* Check if the file sector is valid*/

  /* Check if the available memory block is sufficient */
  boundary->nnds += oldnnds;            /* Update the number of nodes     */
  SAFT2D_Check_memory();          /* Check memory             */

  /* Read vertex coordinates */
  for (i = oldnnds; i < boundary->nnds; i++) {
    fscanf(fptr, "\t%lu\t%lf\t%lf\t%lf\n", &tag, &(boundary->node[i].x),
      &(boundary->node[i].y), &(boundary->node[i].z));
  }
  /* End of coordinate sector */

  /* Read begin of element sector */
  fscanf(fptr, "%s  %lu\n", fsectr, &(boundary->nels));
  assert(strcmp(fsectr, "tria3") == 0);     /* Check if the file sector is valid*/

  /* Check if the available memory block is sufficient */
  boundary->nels += oldnels;            /* Update the number of elements    */
  SAFT2D_Check_memory();          /* Check memory             */

  /* Read the table of connectivities */
  for (i = oldnels; i < boundary->nels; i++) {
    fscanf(fptr, "\t%lu\t%lu\t%lu\t%lu\n", &tag, &v1, &v2, &v3);
    /* Put the smallest index in the first node preserving orientation */
    boundary->element[i].vertex[0] = v1 + oldnnds;
    boundary->element[i].vertex[1] = v2 + oldnnds;
    boundary->element[i].vertex[2] = v3 + oldnnds;
  }
  /* End of element sector */

  /* Close file with the surface mesh */
  fclose(fptr);
}

/*
 * SAFT2D_Print_global_surfs()
 *
 * Prototype:
 *    void SAFT2D_Print_global_surfs(Boundary *);
 *
 * Required header:
 *    saft2d.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a boundary structure.
 *
 * Comments:
 *    It calls SAFT2D_Print_surface_mesh() to print the global surface
 *    meshes (internal and external).
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Print_boundary()
{
  /* Definitions of local variables */
  unsigned long i;              /* Counter and a tag number       */
  char      fname[150];         /* File name and file sector      */
  FILE      *fptr;            /* File pointer             */
  /* Output data file path          */
  #ifdef NDEBUG
    char fpath[150] = "../data/output/";
  #else
    char fpath[150] = "./../../data/output/";
  #endif

  /* Initialize file pointer and name */
  fptr = NULL;
  fname[0] = '\0';

  /* Generate complete file path name */
  strcat(fpath, projname);
  strcat(fpath, "_boundary.msh");

  /* Open file with the surface mesh */
  fptr = fopen(fpath, "w");

  if (fptr == NULL)
    printf("Boundary mesh file could not be opened!\n");

  /* Check if the file is valid */
  assert(fptr != NULL);

  /* begin GiD standard */
  /* Print the begin of coordinate sector */
  fprintf(fptr, "MESH    dimension 3 ElemType Triangle  Nnode 3 \nCoordinates");

  /* Print vertex coordinates */
  for (i = 0; i < boundary->nnds; i++) {
    fprintf(fptr, "\t%lu\t%14.8lf\t%14.8lf\t%14.8lf\n", (i + 1),
      boundary->node[i].x, boundary->node[i].y, boundary->node[i].z);
  }

  /* Print the begin of coordinate sector */
  fprintf(fptr, "end coordinates\nElements");

  /* Print the table of connectivities */
  for (i = 0; i < boundary->nels; i++) {
    fprintf(fptr, "\t%lu\t%lu\t%lu\t%lu\n", (i + 1),
      boundary->element[i].vertex[0], boundary->element[i].vertex[1],
      boundary->element[i].vertex[2]);
  }
  fprintf(fptr, "end elements\n");
  /* end GiD standard */

  /* Close file with the global surface mesh */
  fclose(fptr);
}

/*
 * HEAP_Check_memory()
 *
 * Prototype:
 *    void SAFT2D_Check_memory(Boundary *);
 *
 * Required Header:
 *    saft2d.h, stdlib.h and assert.h.
 *
 * Return Value:
 *    It returns void.
 *
 * Actual Parameters:
 *    It receives a pointer to a boundary structure.
 *
 * Comments:
 *    It checks if the currently available memory block is enough to
 *    insert a new node or element in the boundary structure. If necessary,
 *    it allocates a new memory block.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Check_memory()
{
  /* Check node array */
  if (boundary->nnds > boundary->maxnnds) {
    /* Update the maximum number of nodes */
    boundary->maxnnds = 2*(boundary->nnds) - 1;
    /* Reallocate the node array */
    boundary->node = (Point3D *)
      realloc(boundary->node, (boundary->maxnnds)*sizeof(Point3D));
    assert(boundary->node != NULL);     /* Check if memory was allocated    */
  }

  /* Check element array */
  if (boundary->nels > boundary->maxnels) {
    /* Update the maximum number of elements */
    boundary->maxnels = 2*(boundary->nels) - 1;
    /* Reallocate the element array */
    boundary->element = (Triangle *)
      realloc(boundary->element, (boundary->maxnels)*sizeof(Triangle));
    assert(boundary->element != NULL);      /* Check if memory was allocated    */
  }
}

/*
 * SAFT2D_Finalize()
 *
 * Prototype:
 *    void SAFT2D_Finalize(Boundary *);
 *
 * Required header:
 *    saft2d.h and stdlib.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to an allocated boundary structure.
 *
 * Comments:
 *    It deallocates the main boundary structure and respective sub-structures.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Finalize()
{
  /* Deallocate the array of nodes and set it to NULL */
  free(boundary->node);
  boundary->node = NULL;
  assert(boundary->node == NULL);       /* Check if memory was deallocated    */

  /* Deallocate the array of elements and set it to NULL */
  free(boundary->element);
  boundary->element = NULL;
  assert(boundary->element == NULL);        /* Check if memory was deallocated    */

  /* Deallocate the main boundary structure and set it to NULL */
  free(boundary);
  boundary = NULL;
  assert(boundary == NULL);           /* Check if memory was deallocated    */
}
