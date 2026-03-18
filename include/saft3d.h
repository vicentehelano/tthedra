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
#ifndef SAFT3D_H
#define SAFT3D_H

#include <geom.h>

extern double cpu_time;

/*
 * SAFT3D_Initialize()
 *
 * Prototype:
 *    void SAFT3D_Initialize(Domain **);
 *
 * Required header:
 *    saft3d.h and assert.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives the address of a pointer to a domain structure yet not allocated.
 *
 * Comments:
 *    It allocates the main domain structure and respective sub-structures.
 *
 * Last modified: 21/06/2004.
 */
void SAFT3D_Initialize_mesh(void);

void SAFT3D_Initialize_front(void);

/*
 * SAFT3D_Grid_domain()
 *
 * Prototype:
 *    void SAFT3D_Grid_domain(Domain *);
 *
 * Required header:
 *    saft3d.h, saft2d.h, heap.h and assert.h.
 *
 * Return value:
 *    It returns void.
 *
 * Actual parameters:
 *    It receives a pointer to a domain structure.
 *
 * Comments:
 *	  It applies the advancing-front technique to generate tetrahedra
 *	  from the given boundary triangulation.
 *
 * Last modified: 21/06/2004.
 */
int generate_mesh(void);

void SAFT3D_Get_ideal_point(double, Point3D *, Point3D *,
							Point3D *, Point3D *);
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
 *    It receives a pointer to an allocated domain structure.
 *
 * Comments:
 *    It deallocates the main domain structure and respective sub-structures.
 *
 * Last modified: 21/06/2004.
 */
void SAFT3D_Finalize_mesh(void);

#endif /* SAFT3D_H */
