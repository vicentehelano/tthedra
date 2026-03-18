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
#ifndef QUALI_H
#define QUALI_H

#include <assert.h>
#include <float.h>
#include <geom.h>

/* Minimum geometric tolerance (quali.c) */
extern const double	epsilon;

/* Inverse of the minimum dihedral angle of a regular tetrahedron	*/
#define INVMINDIHEANG_REGTET (8.12374466544386e-001)

/* Inverse of sine of the minimum solid angle of a regular tetrahedron	*/
#define INVOFSIGMA_REGTET (3.67423461417476e+000)

/* minimum solid angle of a regular tetrahedron
	31.586338096527925892262001980113 degrees
	 0.551285598432530807942144151464 radians
*/

/*
 * Inline functions
 *
 * Last modified: 22/10/2004.
 */

/*
 * Functions prototypes
 *
 * Last modified: 22/10/2004.
 */

/*
 * GEOM_Angle3D()
 *
 * Prototype:
 *    double GEOM_Angle3D(Point3D *, Point3D *, double);
 *
 * Return value:
 *    It returns a double precision floating-point number equal to the
 *	  minimum dihedral angle of the tetrahedron in radians.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Last modified: 22/10/2004.
 */
double GEOM_Angle3D(Point3D *, Point3D *, double);


/*
 * QUALI_Tetrahedron_dihedral_angle_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_dihedral_angle_ratio(Point3D *, Point3D *,
 *			Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point number equal to the
 *	  ratio between the minimum dihedral angle of the input tetrahedron
 *	  and the minimum dihedral angle of a regular tetrahedron.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_dihedral_angle_ratio(Point3D *, Point3D *,
											  Point3D *, Point3D *);

/*
 * QUALI_Tetrahedron_solid_angle_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_solid_angle_ratio(Point3D *, Point3D *,
 *			Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point number equal to the
 *    ratio between the sin(minimum solid angle/2) of the input
 *	  tetrahedron and a regular tetrahedron.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_solid_angle_ratio(Point3D *, Point3D *,
										   Point3D *, Point3D *);

/*
 * QUALI_Tetrahedron_radius_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_radius_ratio(Point3D *, Point3D *,
 *			Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point equal number to the
 *	  radius (or aspect) ratio of the tetrahedron.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_radius_ratio(Point3D *, Point3D *,
									  Point3D *, Point3D *);

/*
 * QUALI_Tetrahedron_mean_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_mean_ratio(Point3D *, Point3D *,
 *			Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point equal number to the
 *	  mean ratio.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_mean_ratio(Point3D *, Point3D *,
									Point3D *, Point3D *);

/*
 * QUALI_Tetrahedron_gamma_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_gamma_ratio(Point3D *, Point3D *,
 *			Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point equal number to the
 *	  gamma ratio.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Last modified: 01/11/2004.
 */
double QUALI_Tetrahedron_gamma_ratio(Point3D *, Point3D *,
									 Point3D *, Point3D *);

void QUALI_Print_quality_report(double, Volume_mesh *);

#endif /* QUALI_H */
