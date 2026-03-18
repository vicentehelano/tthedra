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
#ifndef SAFT2D_H
#define SAFT2D_H

#include <geom.h>

extern Surface_mesh *boundary;
extern char projname[150];
extern double deltamax;

/*
 * Initialize()
 *
 * Prototype:
 *    void Initialize(Boundary **);
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
void SAFT2D_Initialize_boundary(void);

void SAFT2D_Normalize_coordinates(unsigned long, Point3D *);

void SAFT2D_Restore_coordinates(unsigned long, Point3D *);

/*
 * Get_boundary_mesh()
 *
 * Prototype:
 *    void Grid_boundary(Boundary *);
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
 *	  Read the triangulation of the surfaces (internal and external)
 *	  of some domain from a file and print the global surface mesh.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Grid_boundary(void);

/*
 * Read_external_surfs()
 *
 * Prototype:
 *    void Read_external_surfs(Boundary *);
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
 *	  It calls the Read_surface_mesh() with the parameter ".esf"
 *	  that indicates a valid extension of a external surface mesh file.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Read_external_surfs(void);

/*
 * Read_internal_surfs()
 *
 * Prototype:
 *    void Read_internal_surfs(Boundary *);
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
 *	  It calls the Read_surface_mesh() with the parameter ".isf"
 *	  that indicates a valid extension of a internal surface mesh file.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Read_internal_surfs(void);

/*
 * Read_surface_mesh()
 *
 * Prototype:
 *    void Read_surface_mesh(char*, Boundary *);
 *
 * Required header:
 *    saft2d.h, string.h, stdlib.h, stdio.h and assert.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a boundary structure and a pointer to a
 *	  file extension string.
 *
 * Comments:
 *	  It reads the surface meshes (internal and external).
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Read_surface_mesh(char *);

/*
 * Print_global_surfs()
 *
 * Prototype:
 *    void Print_global_surfs(Boundary *);
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
 *	  It calls Print_surface_mesh() to print the global surface
 *	  meshes (internal and external).
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Print_boundary(void);

/*
 * Print_surface_mesh()
 *
 * Prototype:
 *    void Print_surface_mesh(char *, Boundary *);
 *
 * Required header:
 *    saft2d.h, stdio.h, string.h and assert.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a boundary structure.
 *
 * Comments:
 *	  It prints a surface mesh.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Print_surface_mesh(char *, Surface_mesh *);

/*
 * HEAP_Check_memory()
 *
 * Prototype:
 *    void Check_memory(Boundary *);
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
 *	  insert a new node or element in the boundary structure. If necessary,
 *	  it allocates a new memory block.
 *
 * Last modified: 06/07/2004.
 */
void SAFT2D_Check_memory(void);

/*
 * Finalize()
 *
 * Prototype:
 *    void Finalize(Boundary *);
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
void SAFT2D_Finalize(void);

#endif /* SAFT2D_H */
