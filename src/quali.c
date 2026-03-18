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
#include <quali.h>
#include <saft3d.h>
#include <saft2d.h>

/*
 * Global variables
 *
 * Last modified: 22/10/2004.
 */

/* Edges */
Edge *edg;
Triangle faces[40000];
unsigned long nedge, *point_help, nface;
double *edglens, *elemvols, *solids, *radius, *means;
long neighbours[4][50000];

/*
 * Angle3D()
 *
 * Prototype:
 *    double Angle3D(Point3D *, Point3D *, double);
 *
 * Return value:
 *    It returns a double precision floating-point number equal to the
 *    minimum dihedral angle of the tetrahedron in radians.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Last modified: 22/10/2004.
 */
double Angle3D(Point3D *vec1, Point3D *vec2, double sqtol)
{
  /* Definitions of local variables */
  double  dot;            /* Store the dot product result */
  double  cosine;           /* Cosine between vec1 and vec2 */
  double  sqrdlen1, sqrdlen2;     /* Squared vector lengths */

  dot    = GEOM_DOT3D((*vec1),(*vec2));
  sqrdlen1 = GEOM_DOT3D((*vec1),(*vec1));
  sqrdlen2 = GEOM_DOT3D((*vec2),(*vec2));

  /* Check for minimum geometric tolerance */
  if (sqrdlen1 > sqtol && sqrdlen2 > sqtol) {
    cosine = dot/sqrt(sqrdlen1*sqrdlen2);
    if (fabs(cosine) > (1.0e+000 - epsilon))
      cosine = GEOM_SIGNDOUBLE((1.0e+000),cosine);
    return (acos(cosine));
  }
  else
    return (PI);
}

/*
 * QUALI_Tetrahedron_dihedral_angle_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_dihedral_angle_ratio(Point3D *, Point3D *,
 *      Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point number equal to the
 *    ratio between the minimum dihedral angle of the input tetrahedron
 *    and the minimum dihedral angle of a regular tetrahedron.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Comments:
 *    It computes the minimum dihedral angle of a tetrahedron applying the
 *    formula 2.2 of the article:
 *    @article{LIU94,
 *      author={A. Liu and B. Joe},
 *      title={Relationship between Tetrahedron Shape Mesaures},
 *      journal={BIT},
 *      year={1994},
 *      month={January},
 *      volume={34},
 *      pages={268--287},
 *      key={tetrahedron shape measure; mesh generation; finite
 *        element analysis},
 *      annote={Printed copy (UFRJ/NCE)}
 *    }
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_dihedral_angle_ratio(Point3D *nd1, Point3D *nd2,
                        Point3D *nd3, Point3D *nd4)
{
  /* Definitions of local variables */
  Point3D v12, v13, v14, v23, v24;  /* Tetrahedron edge vectors   */
  double  sqtol;            /* Minimum square root tolerance*/
  Point3D n123, n124, n134, n234;   /* Normal vectors of faces    */
  double  diang[6];         /* Dihedral angle       */
  double  diangmin;         /* Minimum dihedral angle   */

  /* Compute edge vectors */
  GEOM_SUB3D(v12,(*nd2),(*nd1))
  GEOM_SUB3D(v13,(*nd3),(*nd1))
  GEOM_SUB3D(v14,(*nd4),(*nd1))
  GEOM_SUB3D(v23,(*nd3),(*nd2))
  GEOM_SUB3D(v24,(*nd4),(*nd2))

  /* Compute the minimum square root tolerance */
  sqtol = epsilon*GEOM_MAX12((fabs(nd1->x)),
    (fabs(nd1->y)),(fabs(nd1->z)),
    (fabs(nd2->x)),(fabs(nd2->y)),
    (fabs(nd2->z)),(fabs(nd3->x)),
    (fabs(nd3->y)),(fabs(nd3->z)),
    (fabs(nd4->x)),(fabs(nd4->y)),
    (fabs(nd4->z)));
  sqtol = sqtol*sqtol;

  /* Compute normal vectors of faces */
  GEOM_CROSS3D(n123,v13,v12)
  GEOM_CROSS3D(n124,v12,v14)
  GEOM_CROSS3D(n134,v14,v13)
  GEOM_CROSS3D(n234,v23,v24)

  /* Compute the six dihedral angles */
  diang[0] = PI - Angle3D(&n123,&n124,sqtol);
  diang[1] = PI - Angle3D(&n123,&n134,sqtol);
  diang[2] = PI - Angle3D(&n124,&n134,sqtol);
  diang[3] = PI - Angle3D(&n123,&n234,sqtol);
  diang[4] = PI - Angle3D(&n124,&n234,sqtol);
  diang[5] = PI - Angle3D(&n134,&n234,sqtol);

  diangmin = GEOM_MIN6((diang[0]),(diang[1]),(diang[2]),
    (diang[3]),(diang[4]),(diang[5]));

  return (diangmin*INVMINDIHEANG_REGTET);
}

/*
 * QUALI_Tetrahedron_solid_angle_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_solid_angle_ratio(Point3D *, Point3D *,
 *      Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point number equal to the
 *    ratio between the sin(minimum solid angle/2) of the input
 *    tetrahedron and a regular tetrahedron.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Comments:
 *    It computes the sin(minimum solid angle/2) of a tetrahedron
 *    applying the formula 2.2 of the article:
 *    @article{LIU94,
 *      author={A. Liu and B. Joe},
 *      title={Relationship between Tetrahedron Shape Mesaures},
 *      journal={BIT},
 *      year={1994},
 *      month={January},
 *      volume={34},
 *      pages={268--287},
 *      key={tetrahedron shape measure; mesh generation; finite
 *        element analysis},
 *      annote={Printed copy (UFRJ/NCE)}
 *    }
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_solid_angle_ratio(Point3D *nd1, Point3D *nd2,
                       Point3D *nd3, Point3D *nd4)
{
  /* Definitions of local variables */
  Point3D n, v1, v2;          /* Auxiliary normal vectors   */
  double  leng12, leng13, leng14;   /* Tetrahedron edge lengths   */
  double  leng23, leng24, leng34;   /* Tetrahedron edge lengths   */
  double  parc1, parc2, parc3;    /* Partial sums         */
  double  numer, denom;       /* Numerator and denominator  */
  double  sigma[4];         /* Store sin(solid angle/2.0) */
  double  sigmamin;         /* Store minimum sigma      */

  /* Compute edge lengths */
  GEOM_SUB3D(v2,(*nd3),(*nd2))
  leng23 = sqrt(GEOM_DOT3D(v2,v2));

  GEOM_SUB3D(v1,(*nd4),(*nd2))
  GEOM_SUB3D(v2,(*nd4),(*nd3))
  leng24 = sqrt(GEOM_DOT3D(v1,v1));
  leng34 = sqrt(GEOM_DOT3D(v2,v2));

  GEOM_SUB3D(v1,(*nd2),(*nd1))
  GEOM_SUB3D(v2,(*nd3),(*nd1))
  leng12 = sqrt(GEOM_DOT3D(v1,v1));
  leng13 = sqrt(GEOM_DOT3D(v2,v2));

  GEOM_CROSS3D(n,v1,v2);       /* Used to compute the volume */
  GEOM_SUB3D(v1,(*nd4),(*nd1))
  leng14 = sqrt(GEOM_DOT3D(v1,v1));

  /* Compute (12.0 x volume of the tetrahedron) */
  numer = fabs((2.0e+000*GEOM_DOT3D(v1,n)));

  /* Compute partial sums */
  parc1 = leng12 + leng13;
  parc2 = leng12 + leng14;
  parc3 = leng13 + leng14;
  denom = (parc1+leng23)*(parc1-leng23)*(parc2+leng24)*
    (parc2-leng24)*(parc3+leng34)*(parc3-leng34);

  /* Check for precision degeneracies */
  if (denom <= DBL_EPSILON) {
    sigma[0] = 0.0e+000;
  }
  else {
    sigma[0] = numer/sqrt(denom);
  }

  /* Compute partial sums */
  parc1 = leng12 + leng23;
  parc2 = leng12 + leng24;
  parc3 = leng23 + leng24;
  denom = (parc1+leng13)*(parc1-leng13)*(parc2+leng14)*
    (parc2-leng14)*(parc3+leng34)*(parc3-leng34);

  /* Check for precision degeneracies */
  if (denom <= DBL_EPSILON) {
    sigma[1] = 0.0e+000;
  }
  else {
    sigma[1] = numer/sqrt(denom);
  }

  /* Compute partial sums */
  parc1 = leng13 + leng23;
  parc2 = leng13 + leng34;
  parc3 = leng23 + leng34;
  denom = (parc1+leng12)*(parc1-leng12)*(parc2+leng14)*
    (parc2-leng14)*(parc3+leng24)*(parc3-leng24);

  /* Check for precision degeneracies */
  if (denom <= DBL_EPSILON) {
    sigma[2] = 0.0e+000;
  }
  else {
    sigma[2] = numer/sqrt(denom);
  }

  /* Compute partial sums */
  parc1 = leng14 + leng24;
  parc2 = leng14 + leng34;
  parc3 = leng24 + leng34;
  denom = (parc1+leng12)*(parc1-leng12)*(parc2+leng13)*
    (parc2-leng13)*(parc3+leng23)*(parc3-leng23);

  /* Check for precision degeneracies */
  if (denom <= DBL_EPSILON) {
    sigma[3] = 0.0e+000;
  }
  else {
    sigma[3] = numer/sqrt(denom);
  }

  /* Get the minimum sine of half of the solid angles */
  sigmamin = GEOM_MIN4(sigma[0],sigma[1],sigma[2],sigma[3]);

  return (sigmamin*INVOFSIGMA_REGTET);
}

/*
 * QUALI_Tetrahedron_radius_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_radius_ratio(Point3D *, Point3D *,
 *      Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point equal number to the
 *    radius (or aspect) ratio of the tetrahedron.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Comments:
 *    It computes the radius (or aspect) ratio of a tetrahedron applying
 *    the formula 3.1 of the article:
 *    @article{LIU94,
 *      author={A. Liu and B. Joe},
 *      title={Relationship between Tetrahedron Shape Mesaures},
 *      journal={BIT},
 *      year={1994},
 *      month={January},
 *      volume={34},
 *      pages={268--287},
 *      key={tetrahedron shape measure; mesh generation; finite
 *        element analysis},
 *      annote={Printed copy (UFRJ/NCE)}
 *    }
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_radius_ratio(Point3D *nd1, Point3D *nd2,
                    Point3D *nd3, Point3D *nd4)
{
  /* Definitions of local variables */
  Point3D vec12, vec13, vec14;    /* Tetrahedron edge vectors   */
  Point3D vec23, vec24, vec34;    /* Tetrahedron edge vectors   */
  Point3D n;              /* Auxiliary normal vector    */
  double  sarea;            /* Total superficial area   */
  double  sqrdlen12, sqrdlen13;   /* Squared edge lengths     */
  double  sqrdlen14, sqrdlen23;   /* Squared edge lengths     */
  double  sqrdlen24, sqrdlen34;   /* Squared edge lengths     */
  double  prod1, prod2, prod3;    /* Product of squared lengths */
  double  parc1, parc2;       /* Partial sums         */
  double  numer, denom;       /* Numerator and denominator  */

  /* Compute edge vectors */
  GEOM_SUB3D(vec12,(*nd2),(*nd1))
  GEOM_SUB3D(vec13,(*nd3),(*nd1))
  GEOM_SUB3D(vec14,(*nd4),(*nd1))
  GEOM_SUB3D(vec23,(*nd3),(*nd2))
  GEOM_SUB3D(vec24,(*nd4),(*nd2))
  GEOM_SUB3D(vec34,(*nd4),(*nd3))

  /* Compute total superficial area */
  GEOM_CROSS3D(n,vec12,vec14)
  sarea = sqrt(GEOM_DOT3D(n,n));

  GEOM_CROSS3D(n,vec12,vec13)
  sarea += sqrt(GEOM_DOT3D(n,n));

  GEOM_CROSS3D(n,vec23,vec24)
  sarea += sqrt(GEOM_DOT3D(n,n));

  GEOM_CROSS3D(n,vec13,vec14)
  sarea += sqrt(GEOM_DOT3D(n,n));

  /* Compute squared edge lengths */
  sqrdlen12 = GEOM_DOT3D(vec12,vec12);
  sqrdlen13 = GEOM_DOT3D(vec13,vec13);
  sqrdlen14 = GEOM_DOT3D(vec14,vec14);
  sqrdlen23 = GEOM_DOT3D(vec23,vec23);
  sqrdlen24 = GEOM_DOT3D(vec24,vec24);
  sqrdlen34 = GEOM_DOT3D(vec34,vec34);

  /* Compute product of squared lengths */
  prod1 = sqrt(sqrdlen12*sqrdlen34);
  prod2 = sqrt(sqrdlen13*sqrdlen24);
  prod3 = sqrt(sqrdlen14*sqrdlen23);

  /* Compute partial sums and denominator */
  parc1 = prod1 + prod2;
  parc2 = prod1 - prod2;
  denom = (sarea)*sqrt(fabs((parc1+prod3)*(parc1-prod3)*
    (prod3+parc2)*(prod3-parc2)));

  /* Check for precision degeneracies */
  if (denom <= DBL_EPSILON)
    return (0.0e+000);
  else {
    numer = GEOM_DOT3D(vec12,n);
    return ((12.0e+000)*numer*numer/denom);
  }
}

/*
 * QUALI_Tetrahedron_mean_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_mean_ratio(Point3D *, Point3D *,
 *      Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point equal number to the
 *    mean ratio.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Comments:
 *    It computes the mean ratio of a tetrahedron applying the formula
 *    3.4 of the article:
 *    @article{LIU94,
 *      author={A. Liu and B. Joe},
 *      title={Relationship between Tetrahedron Shape Mesaures},
 *      journal={BIT},
 *      year={1994},
 *      month={January},
 *      volume={34},
 *      pages={268--287},
 *      key={tetrahedron shape measure; mesh generation; finite
 *        element analysis},
 *      annote={Printed copy (UFRJ/NCE)}
 *    }
 *
 * Last modified: 22/10/2004.
 */
double QUALI_Tetrahedron_mean_ratio(Point3D *nd1, Point3D *nd2,
                  Point3D *nd3, Point3D *nd4)
{
  /* Definitions of local variables */
  Point3D n, v1, v2;          /* Auxiliary normal vectors   */
  double  denom, numer;       /* Numerator and denominator  */

  /* Compute the denominator of the formula 3.4 in LIU94 */
  GEOM_SUB3D(n,(*nd3),(*nd2))
  denom = GEOM_DOT3D(n,n);

  GEOM_SUB3D(n,(*nd4),(*nd2))
  denom += GEOM_DOT3D(n,n);
  GEOM_SUB3D(n,(*nd4),(*nd3))
  denom += GEOM_DOT3D(n,n);
  
  GEOM_SUB3D(v1,(*nd2),(*nd1))
  denom += GEOM_DOT3D(v1,v1);
  GEOM_SUB3D(v2,(*nd3),(*nd1))
  denom += GEOM_DOT3D(v2,v2);

  /* Compute the denominator of the formula 3.4 in LIU94 */
  GEOM_CROSS3D(n,v1,v2);
  GEOM_SUB3D(v1,(*nd4),(*nd1))
  numer = fabs((GEOM_DOT3D(v1,n)));

  denom += GEOM_DOT3D(v1,v1);

  /* Check for precision degeneracies */
  if (denom <= DBL_EPSILON)
    return (0.0e+000);
  else
    return (((7.55952629936924e+000)/denom)*
      pow(numer,(6.66666666666667e-001)));
}

/*
 * QUALI_Tetrahedron_gamma_ratio()
 *
 * Prototype:
 *    double QUALI_Tetrahedron_gamma_ratio(Point3D *, Point3D *,
 *      Point3D *, Point3D *);
 *
 * Required header:
 *    quali.h.
 *
 * Return value:
 *    It returns a double precision floating-point equal number to the
 *    gamma ratio.
 *
 * Actual parameters:
 *    It receives pointers to the vertices of a tetrahedron.
 *
 * Comments:
 *    It computes the gamma ratio of a tetrahedron applying the formula
 *    in the article:
 *    @article{PAR93,
 *      author={V. N. Parthasarathy and C. M. Graichen and
 *        A. F. Hathaway},
 *      title={A Comparison of Tetrahedral Quality Mesures},
 *      journal={Finite Elements in Analysis and Design},
 *      year={1993},
 *      month={Jun.},
 *      volume={15},
 *      pages={255--261},
 *      key={tetrahedral quality measures},
 *      annote={Printed and digital copy}
 *    }
 *
 * Last modified: 01/11/2004.
 */
double QUALI_Tetrahedron_gamma_ratio(Point3D *nd1, Point3D *nd2,
                   Point3D *nd3, Point3D *nd4)
{
  /* Definitions of local variables */
  Point3D n, v1, v2;          /* Auxiliary normal vectors   */
  double  denom, numer;       /* Numerator and denominator  */

  /* Compute the numerator of the formula in PAR93 */
  GEOM_SUB3D(n,(*nd3),(*nd2))
  numer = GEOM_DOT3D(n,n);

  GEOM_SUB3D(n,(*nd4),(*nd2))
  numer += GEOM_DOT3D(n,n);
  GEOM_SUB3D(n,(*nd4),(*nd3))
  numer += GEOM_DOT3D(n,n);
  
  GEOM_SUB3D(v1,(*nd2),(*nd1))
  numer += GEOM_DOT3D(v1,v1);
  GEOM_SUB3D(v2,(*nd3),(*nd1))
  numer += GEOM_DOT3D(v2,v2);

  /* Compute the denominator of the formula in PAR93 */
  GEOM_CROSS3D(n,v1,v2);
  GEOM_SUB3D(v1,(*nd4),(*nd1))
  denom = fabs((GEOM_DOT3D(v1,n)));  /* 6*Volume */

  numer += GEOM_DOT3D(v1,v1);

  /* Check for precision degeneracies */
  if (denom <= DBL_EPSILON)
    return (0.0e+000);
  else
    return ((SQRT_INV_SIX*numer*sqrt(numer))/denom);
}

void QUALI_Get_elems_surround_elems(Volume_mesh *domain)
{
  unsigned long i, j, k;
  Triangle face[4];

  point_help = (unsigned long *)
    malloc(domain->nnds*sizeof(unsigned long));
  assert(point_help != NULL); /* Check if memory was */
  memset(point_help, 0, domain->nnds*sizeof(unsigned long));

  for (j = 0; j < domain->nels; j++) {
    for (i = 0; i < 4; i++) {
      neighbours[i][j] = 0;
    }
  }
  
  for (i = 0; i < domain->nels; i++) {
    face[0].vertex[0] = domain->element[i].vertex[1];
    face[0].vertex[1] = domain->element[i].vertex[2];
    face[0].vertex[2] = domain->element[i].vertex[3];

    face[1].vertex[0] = domain->element[i].vertex[0];
    face[1].vertex[1] = domain->element[i].vertex[3];
    face[1].vertex[2] = domain->element[i].vertex[2];

    face[2].vertex[0] = domain->element[i].vertex[0];
    face[2].vertex[1] = domain->element[i].vertex[1];
    face[2].vertex[2] = domain->element[i].vertex[3];

    face[3].vertex[0] = domain->element[i].vertex[0];
    face[3].vertex[1] = domain->element[i].vertex[2];
    face[3].vertex[2] = domain->element[i].vertex[1];
    for (j = 0; j < 4; j++) {

      point_help[face[j].vertex[0]-1] = 1;
      point_help[face[j].vertex[1]-1] = 1;
      point_help[face[j].vertex[2]-1] = 1;
      for (k = 0; k < domain->nels; k++) {
        
        if ((k != i) && ((domain->element[k].vertex[0] == face[j].vertex[0]) ||
        (domain->element[k].vertex[1] == face[j].vertex[0]) ||
        (domain->element[k].vertex[2] == face[j].vertex[0]) ||
        (domain->element[k].vertex[3] == face[j].vertex[0]))) {
          
          if ((k != i) && ((domain->element[k].vertex[0] == face[j].vertex[1]) ||
          (domain->element[k].vertex[1] == face[j].vertex[1]) ||
          (domain->element[k].vertex[2] == face[j].vertex[1]) ||
          (domain->element[k].vertex[3] == face[j].vertex[1]))) {
            
            if ((k != i) && ((domain->element[k].vertex[0] == face[j].vertex[2]) ||
            (domain->element[k].vertex[1] == face[j].vertex[2]) ||
            (domain->element[k].vertex[2] == face[j].vertex[2]) ||
            (domain->element[k].vertex[3] == face[j].vertex[2]))) {
              neighbours[j][i] = k + 1;
              
            }
            
          }
          
        }
        
      }
      point_help[face[j].vertex[0]-1] = 0;
      point_help[face[j].vertex[1]-1] = 0;
      point_help[face[j].vertex[2]-1] = 0;
    }
  }
}

void QUALI_Get_faces(Volume_mesh *domain)
{
  unsigned long i, j;
  Triangle f[4];

  QUALI_Get_elems_surround_elems(domain);
  free(point_help);
  
  nface = 0;
  for (i = 0; i < domain->nels; i++) {
    f[0].vertex[0] = domain->element[i].vertex[1];
    f[0].vertex[1] = domain->element[i].vertex[2];
    f[0].vertex[2] = domain->element[i].vertex[3];

    f[1].vertex[0] = domain->element[i].vertex[0];
    f[1].vertex[1] = domain->element[i].vertex[3];
    f[1].vertex[2] = domain->element[i].vertex[2];

    f[2].vertex[0] = domain->element[i].vertex[0];
    f[2].vertex[1] = domain->element[i].vertex[1];
    f[2].vertex[2] = domain->element[i].vertex[3];

    f[3].vertex[0] = domain->element[i].vertex[0];
    f[3].vertex[1] = domain->element[i].vertex[2];
    f[3].vertex[2] = domain->element[i].vertex[1];
    for (j = 0; j < 4; j++) {
      if (neighbours[j][i] == 0) {
        faces[nface] = f[j];
        nface++;
      }
      else if (neighbours[j][i] > 0) {
        faces[nface] = f[j];
        nface++;
        if (neighbours[0][neighbours[j][i]-1] == (i+1))
          neighbours[0][neighbours[j][i]-1] = -1;
        if (neighbours[1][neighbours[j][i]-1] == (i+1))
          neighbours[1][neighbours[j][i]-1] = -1;
        if (neighbours[2][neighbours[j][i]-1] == (i+1))
          neighbours[2][neighbours[j][i]-1] = -1;
        if (neighbours[3][neighbours[j][i]-1] == (i+1))
          neighbours[3][neighbours[j][i]-1] = -1;
      }
      else
        continue;
    }
  }
  domain->nfaces = nface;
}

void QUALI_Get_edges(Volume_mesh *domain)
{
  unsigned long i, j, k, iedge, nedge0;

  edg = (Edge *) malloc(12*domain->nnds*sizeof(Edge));
  assert(edg != NULL);        /* Check if memory was */

  point_help = (unsigned long *) malloc(domain->nnds*sizeof(unsigned long));
  assert(point_help != NULL);       /* Check if memory was */
  memset(point_help, 0, domain->nnds*sizeof(unsigned long));

  nedge = 0;
  for (i = 0; i < domain->nnds; i++) {
    nedge0 = nedge;

    for (j = 0; j < domain->nels; j++) {
      if ((i + 1) == domain->element[j].vertex[0] || (i + 1) == domain->element[j].vertex[1] ||
        (i + 1) == domain->element[j].vertex[2] || (i + 1) == domain->element[j].vertex[3]) {
        for (k = 0; k < 4; k++) {
          if (domain->element[j].vertex[k] > (i + 1) &&
            point_help[domain->element[j].vertex[k] - 1] == 0) {
              nedge++;
              edg[nedge-1].vertex[0] = (i + 1);
              edg[nedge-1].vertex[1] = domain->element[j].vertex[k];
              point_help[domain->element[j].vertex[k] - 1] = 1;
          }
        }
      }
    }

    for (iedge = nedge0; iedge < nedge; iedge++) {
      point_help[edg[iedge].vertex[1] - 1] = 0;
    }
  }
  domain->nedges = nedge;
}

void QUALI_Print_quality_report(double elemsize, Volume_mesh *mesh_ptr)
{
  unsigned long i, v1, v2, v3, v4;
  double      minedglen, avgedglen, maxedglen, edglendev;
  double      minelemvol, avgelemvol, maxelemvol, elemvoldev;
  double      minsolid, avgsolid, maxsolid, soliddev;
  double      minradius, avgradius, maxradius, radiusdev;
  double      minmean, avgmean, maxmean, meandev;
  Point3D     *nd1, *nd2, *nd3, *nd4, delta, n, vec1, vec2;
  FILE      *fptr;
  time_t      t0;

  /* Output data file path          */
  #ifdef NDEBUG
    char fpath[150] = "../data/output/";
  #else
    char fpath[150] = "./../../data/output/";
  #endif


  QUALI_Get_faces(mesh_ptr);


  /* Generate complete file path name */
  strcat(fpath, projname);
  strcat(fpath, ".qty");

  fptr = NULL;
  while (fptr == NULL) {
    fptr = fopen(fpath, "w");
  }

  /* Edge Lengths */
  QUALI_Get_edges(mesh_ptr);
  free(point_help);
  edglens = (double *) malloc(nedge*sizeof(double));
  minedglen = DBL_MAX;
  maxedglen = -1.0e+000;
  avgedglen = 0.0e+000;
  for (i = 0; i < nedge; i++) {
    GEOM_SQRDDIST3D(edglens[i],mesh_ptr->node[edg[i].vertex[0]-1],
      mesh_ptr->node[edg[i].vertex[1]-1]);
    edglens[i] = sqrt(edglens[i]);
    if (edglens[i] < minedglen)
      minedglen = edglens[i];
    if (edglens[i] > maxedglen)
      maxedglen = edglens[i];
    avgedglen += edglens[i];
  }
  avgedglen /= nedge;
  /* Standard Deviation */
  edglendev = 0.0e+000;
  for (i = 0; i < nedge; i++) {
    edglendev += (edglens[i] - avgedglen)*(edglens[i] - avgedglen);
  }
  edglendev = sqrt(edglendev/(nedge - 1));

  /* Element Volumes */
  elemvols = (double *) malloc(mesh_ptr->nels*sizeof(double));
  minelemvol = DBL_MAX;
  maxelemvol = -1.0e+000;
  avgelemvol = 0.0e+000;
  for (i = 0; i < mesh_ptr->nels; i++) {
    v1 = mesh_ptr->element[i].vertex[0] - 1;
    v2 = mesh_ptr->element[i].vertex[1] - 1;
    v3 = mesh_ptr->element[i].vertex[2] - 1;
    v4 = mesh_ptr->element[i].vertex[3] - 1;

    nd1 = &(mesh_ptr->node[v1]);
    nd2 = &(mesh_ptr->node[v2]);
    nd3 = &(mesh_ptr->node[v3]);
    nd4 = &(mesh_ptr->node[v4]);

    GEOM_SUB3D(vec1,(*nd2),(*nd1))
    GEOM_SUB3D(vec2,(*nd3),(*nd1))
    GEOM_CROSS3D(n,vec1,vec2)
    GEOM_SUB3D(vec1,(*nd4),(*nd1))

    elemvols[i] = 1.66666666666667e-001*GEOM_DOT3D(n,vec1);

    if (minelemvol > elemvols[i])
      minelemvol = elemvols[i];
    if (maxelemvol < elemvols[i])
      maxelemvol = elemvols[i];
    avgelemvol += elemvols[i];
  }
  avgelemvol /= mesh_ptr->nels;
  /* Standard Deviation */
  elemvoldev = 0.0e+000;
  for (i = 0; i < mesh_ptr->nels; i++) {
    elemvoldev += (elemvols[i] - avgelemvol)*(elemvols[i] - avgelemvol);
  }
  elemvoldev = sqrt(elemvoldev/(mesh_ptr->nels - 1));


  /* Element Quality Shape Metrics */
  solids = (double *) malloc(mesh_ptr->nels*sizeof(double));
  minsolid = DBL_MAX;
  maxsolid = -1.0e+000;
  avgsolid = 0.0e+000;
  radius = (double *) malloc(mesh_ptr->nels*sizeof(double));
  minradius = DBL_MAX;
  maxradius = -1.0e+000;
  avgradius = 0.0e+000;
  means = (double *) malloc(mesh_ptr->nels*sizeof(double));
  minmean = DBL_MAX;
  maxmean = -1.0e+000;
  avgmean = 0.0e+000;

  for (i = 0; i < mesh_ptr->nels; i++) {
    v1 = mesh_ptr->element[i].vertex[0] - 1;
    v2 = mesh_ptr->element[i].vertex[1] - 1;
    v3 = mesh_ptr->element[i].vertex[2] - 1;
    v4 = mesh_ptr->element[i].vertex[3] - 1;

    nd1 = &(mesh_ptr->node[v1]);
    nd2 = &(mesh_ptr->node[v2]);
    nd3 = &(mesh_ptr->node[v3]);
    nd4 = &(mesh_ptr->node[v4]);

    means[i] = QUALI_Tetrahedron_mean_ratio(nd1, nd2, nd3, nd4);
    radius[i] = QUALI_Tetrahedron_radius_ratio(nd1, nd2, nd3, nd4);
    solids[i] = QUALI_Tetrahedron_solid_angle_ratio(nd1, nd2, nd3, nd4);

    if (minsolid > solids[i])
      minsolid = solids[i];
    if (maxsolid < solids[i])
      maxsolid = solids[i];
    avgsolid += solids[i];

    if (minradius > radius[i])
      minradius = radius[i];
    if (maxradius < radius[i])
      maxradius = radius[i];
    avgradius += radius[i];

    if (minmean > means[i])
      minmean = means[i];
    if (maxmean < means[i])
      maxmean = means[i];
    avgmean += means[i];
  }
  avgsolid /= mesh_ptr->nels;
  avgradius /= mesh_ptr->nels;
  avgmean /= mesh_ptr->nels;
  /* Standard Deviation */
  soliddev = 0.0e+000;
  radiusdev = 0.0e+000;
  meandev = 0.0e+000;
  for (i = 0; i < mesh_ptr->nels; i++) {
    soliddev += (solids[i] - avgsolid)*(solids[i] - avgsolid);
    radiusdev += (radius[i] - avgradius)*(radius[i] - avgradius);
    meandev += (means[i] - avgmean)*(means[i] - avgmean);
  }
  soliddev = sqrt(soliddev/(mesh_ptr->nels - 1));
  radiusdev = sqrt(radiusdev/(mesh_ptr->nels - 1));
  meandev = sqrt(meandev/(mesh_ptr->nels - 1));

  fprintf(fptr,"Saft3D Mesh Quality Report\n\n");
  time(&t0);
  fprintf(fptr,"Mesh file: %s.msh\n", projname);
  fprintf(fptr,"\n");
  fprintf(fptr,"1. Mesh Entities\n");
  fprintf(fptr,"        No. of nodes: %lu\n", mesh_ptr->nnds);
  fprintf(fptr,"        No. of edges: %lu \n", nedge);
  fprintf(fptr,"        No. of elems: %lu \n", mesh_ptr->nels);
  fprintf(fptr,"\n");

  fprintf(fptr,"2. CPU Time\n");
  fprintf(fptr,"        %-8.4lf secs\n", cpu_time);
  fprintf(fptr,"\n");

  fprintf(fptr,"3. Statistical Outputs\n");
  fprintf(fptr,"        Target element size: %4.3e\n\n", elemsize);
  fprintf(fptr,"                      &      min       &      avrg      &");
  fprintf(fptr,"       max      &     stddev     &      |min|     &");
  fprintf(fptr,"      |max|     &     |stddev|   \n");
  fprintf(fptr,"      ----------------&----------------&----------------&");
  fprintf(fptr,"----------------&----------------&----------------&");
  fprintf(fptr,"----------------&----------------\n");

  fprintf(fptr,"        {\\small $\\ell$}         &  $%4.3e$  &  $%4.3e$  &", minedglen,
    avgedglen);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  & & $%4.3e$  &", maxedglen, edglendev,
    minedglen/avgedglen);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  \\\\\n", maxedglen/avgedglen, edglendev/avgedglen);

  fprintf(fptr,"        {\\small $V$}            &  $%4.3e$  &  $%4.3e$  &", minelemvol,
    avgelemvol);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  & & $%4.3e$  &", maxelemvol, elemvoldev,
    minelemvol/avgelemvol);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  \\\\\n", maxelemvol/avgelemvol, elemvoldev/avgelemvol);

  fprintf(fptr,"        {\\small $\\sigma_{min}$} &  $%4.3e$  &  $%4.3e$  &", minsolid,
    avgsolid);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  & & $%4.3e$  &", maxsolid, soliddev,
    minsolid/avgsolid);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  \\\\\n", maxsolid/avgsolid, soliddev/avgsolid);

  fprintf(fptr,"        {\\small $\\rho$}         &  $%4.3e$  &  $%4.3e$  &", minradius,
    avgradius);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  & & $%4.3e$  &", maxradius, radiusdev,
    minradius/avgradius);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  \\\\\n", maxradius/avgradius, radiusdev/avgradius);

  fprintf(fptr,"        {\\small $\\eta$}         &  $%4.3e$  &  $%4.3e$  &", minmean,
    avgmean);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  & & $%4.3e$  &", maxmean, meandev,
    minmean/avgmean);
  fprintf(fptr,"  $%4.3e$  &  $%4.3e$  \\\\\n", maxmean/avgmean, meandev/avgmean);
  fprintf(fptr,"\n");

  fprintf(fptr,"4. Free Form Data Sources\n");
  fprintf(fptr,"             Id          Elem. Volumes    Solid Angles   ");
  fprintf(fptr,"  Radius Ratio      Mean Ratio\n");
  for (i = 0; i < mesh_ptr->nels; i++) {
    fprintf(fptr,"             %-10ld   %6.5e    %6.5e     %6.5e     %6.5e\n",
      i+1, elemvols[i], solids[i], radius[i], means[i]);
  }
  fprintf(fptr,"\n");

  fprintf(fptr,"             Id          Edge Lengths     \n");
  for (i = 0; i < nedge; i++) {
    fprintf(fptr,"             %-10ld  %6.5e    \n",
      i+1, edglens[i]);
  }

  fprintf(fptr,"\n5. Normalized Data Sources\n");
  fprintf(fptr,"             Id         |Elem. Volumes|\n");
  for (i = 0; i < mesh_ptr->nels; i++) {
    fprintf(fptr,"             %-10ld   %6.5e\n",
      i+1, elemvols[i]/avgelemvol);
  }
  fprintf(fptr,"\n");

  fprintf(fptr,"             Id         |Edge Lengths|    \n");
  for (i = 0; i < nedge; i++) {
    fprintf(fptr,"             %-10ld  %6.5e    \n",
      i+1, edglens[i]/avgedglen);
  }

  fprintf(fptr,"\n\n");
  fprintf(fptr,"Creation date: %s\n", asctime(localtime(&t0)));
  free(edglens);
  free(elemvols);
  free(solids);
  free(radius);
  free(means);

  fclose(fptr);
}
