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
#ifndef GEOM_H
#define GEOM_H

#include <stdlib.h>
#include <math.h>

#define EPSILON 1.0E-16

extern const double epsilon;

/*
 * 3D Point structure
 *
 */
typedef enum boolean {
  TRUE    = 1,
  FALSE   = 0
} Boolean;

typedef struct  point2D {
  double x, y;
} Point2D;

typedef struct  point3D {
  double x, y, z;
} Point3D;

/*
 * Edge structure
 *
 */
typedef struct  edge  {
  unsigned long vertex[2];          /* Id number of the two vertex       */
} Edge;

/*
 * Linear triangle structure
 *
 */
typedef struct  triangle  {
  unsigned long vertex[3];          /* Id number of the three vertex       */
} Triangle;

/*
 * Linear tetrahedron structure
 *
 */
typedef struct  tetrahedron {
  unsigned long vertex[4];          /* Id number of the four vertex        */
} Tetrahedron;

/*
 * Sphere structure
 *
 */
typedef struct  sphere {
  Point3D center;
  double  radius;
} Sphere;

/*
 * Bounding box structure
 *
 */
typedef struct  bounding_box {
  Point3D min, max;
} Bounding_box;

/*
 * Candidate points structure
 *
 */
typedef struct  candidate_triangles {
  unsigned long ntri;
  unsigned long maxntri;
  unsigned long *triangle;
} Candidate_triangles;


/*
 * Advancing front structure
 *
 */
typedef struct  advancing_front {
  unsigned long nels;           /* Number of elements                  */
  unsigned long maxnels;          /* Maximum number of elements          */
  Triangle    *element;         /* Elements array                      */
} Advancing_front;

/*
 * Surface_mesh structure
 *
 */
typedef struct  surface_mesh  {
  unsigned long nnds;           /* Number of coordinate nodes          */
  unsigned long nels;           /* Number of elements                  */
  unsigned long maxnnds;          /* Maximum number of coordinate nodes  */
  unsigned long maxnels;          /* Maximum number of elements          */
  Point3D     *node;            /* Coordinate nodes array              */
  Triangle    *element;         /* Elements array                      */
} Surface_mesh;

/*
 * Volume_mesh structure
 *
 */
typedef struct  volume_mesh   {
  unsigned long nnds;           /* Number of coordinate nodes          */
  unsigned long nels;           /* Number of elements                  */
  unsigned long nfaces;
  unsigned long nedges;
  unsigned long maxnnds;          /* Maximum number of coordinate nodes  */
  unsigned long maxnels;          /* Maximum number of elements          */
  Point3D     *node;            /* Coordinate nodes array              */
  Tetrahedron   *element;         /* Elements array                      */
} Volume_mesh;


#define SQR(a) (a*a)

#define GEOM_SQRDDIST3D(dist,pt1,pt2) {\
  GEOM_SUB3D(delta,pt1,pt2) \
  dist = SQR(delta.x) + SQR(delta.y) + SQR(delta.z);\
}


enum contants {
  EMPTY = 0,
  VERTEX = 1,
  EDGE = 2
} Contants;


#define GEOM_TRIANGLE_CENTROID3D(CENTROID,V1,V2,V3) \
  CENTROID.x=((V1.x+V2.x+V3.x)*ONETHIRD);\
  CENTROID.y=((V1.y+V2.y+V3.y)*ONETHIRD);\
  CENTROID.z=((V1.z+V2.z+V3.z)*ONETHIRD);

/* CONSTANTS */
#define PI (3.14159265358979e+000)

#define HALFPI (1.57079632679489e+000)

#define SQUARED010 (0.100e-001)

#define SQUARED015 (0.225e-001)

#define SQUARED020 (0.400e-001)

#define SQUARED025 (0.625e-001)

#define SQUARED030 (0.900e-001)

#define SQUARED035 (1.225e-001)

#define SQUARED037 (1.369e-001)

#define SQUARED040 (1.600e-001)

#define SQUARED045 (2.025e-001)

#define SQUARED046 (2.116e-001)

#define SQUARED047 (2.209e-001)

#define SQUARED048 (2.304e-001)

#define SQUARED049 (2.401e-001)

#define SQUARED050 (2.500e-001)

#define SQUARED055 (3.025e-001)

#define SQUARED060 (3.600e-001)

#define SQUARED065 (4.225e-001)

#define SQUARED066 (4.356e-001)

#define SQUARED067 (4.489e-001)

#define SQUARED068 (4.624e-001)

#define SQUARED069 (4.761e-001)

#define SQUARED070 (4.900e-001)

#define SQUARED075 (5.652e-001)

#define SQUARED080 (6.400e-001)

#define SQUARED085 (7.225e-001)

#define SQUARED090 (8.100e-001)

#define SQUARED095 (9.025e-001)

#define SQUARED0100 (1.000e+000)

#define ONETHIRD (3.33333333333333e-001)

#define ANGLE70_INRAD (1.22173047639603e+000)

#define ANGLE60_INRAD (1.04719755119659e+000)

#define ANGLE54_INRAD (9.42477796076937e-001)

#define ANGLE45_INRAD (7.85398163397448e-001)

#define ANGLE35_INRAD (6.10865238198015e-001)

#define ANGLE30_INRAD (5.23598775598298e-001)

#define ANGLE29_INRAD (5.06145483078355e-001)

#define ANGLE27_INRAD (4.71238898038469e-001)

#define ANGLE25_INRAD (4.36332312998582e-001)

#define ANGLE20_INRAD (3.49065850398865e-001)

#define ANGLE17_INRAD (2.96705972839036e-001)

#define ANGLE14_INRAD (2.44346095279206e-001)

#define ANGLE10_INRAD (1.74532925199432e-001)

/* sqrt(1/6) */
#define SQRT_INV_SIX (4.08248290463863e-001)

/* ABSOLUTE */
#define GEOM_ABSINTEGER(A)    ((A)<(0)?(-(A)):(A))

#define GEOM_ABSFLOAT(A)    ((A)<(0.0f)?(-(A)):(A))

#define GEOM_ABSDOUBLE(A)   ((A)<(0.0e+000)?(-(A)):(A))

/* SIGN MACRO FUNCTIONS (like in Fortran) */
#define GEOM_SIGNDOUBLE(A,B) ((B)>=(0.0e+000))?\
  (GEOM_ABSFLOAT(A)):(-GEOM_ABSFLOAT(A))

/* MINIMUM */
#define GEOM_MIN2(A,B)      ((A)<(B)?(A):(B))

#define GEOM_MIN3(A,B,C)    GEOM_MIN2(GEOM_MIN2(A,B),C)

#define GEOM_MIN4(A,B,C,D)    GEOM_MIN2(GEOM_MIN2(A,B),GEOM_MIN2(C,D))

#define GEOM_MIN6(A,B,C,D,E,F)  GEOM_MIN2(GEOM_MIN2(GEOM_MIN2(A,B),\
  GEOM_MIN2(C,D)),GEOM_MIN2(E,F))

/* MAXIMUM */
#define GEOM_MAX2(A,B)    ((A>B)?(A):(B))

#define GEOM_MAX3(A,B,C)  GEOM_MAX2(GEOM_MAX2(A,B),C)

#define GEOM_MAX4(A,B,C,D)  GEOM_MAX2(GEOM_MAX3(A,B,C),D)

#define GEOM_MAX12(A,B,C,D,E,F,G,H,I,J,K,L) \
  GEOM_MAX4((GEOM_MAX3(A,B,C)),(GEOM_MAX3(D,E,F)),\
  (GEOM_MAX3(G,H,I)),(GEOM_MAX3(J,K,L)))

/* VECTOR OPERATIONS */
#define GEOM_SCALAR3D(ALPHA,V)      /* 3D Scalar multiplication   */\
  V.x*=ALPHA;        \
  V.y*=ALPHA;        \
  V.z*=ALPHA;

#define GEOM_ADD3D(DEST,V1,V2)      /* 3D Vector addition     */\
  DEST.x=V1.x+V2.x;        \
  DEST.y=V1.y+V2.y;        \
  DEST.z=V1.z+V2.z;

#define GEOM_SUB3D(DEST,V1,V2)      /* 3D Vector subtraction    */\
  DEST.x=V1.x-V2.x;        \
  DEST.y=V1.y-V2.y;        \
  DEST.z=V1.z-V2.z;

#define GEOM_DOT3D(V1,V2)       /* 3D Dot product       */\
  (V1.x*V2.x+V1.y*V2.y+V1.z*V2.z)

#define GEOM_CROSS3D(DEST,V1,V2)    /* 3D Cross product       */\
  DEST.x=V1.y*V2.z-V1.z*V2.y;  \
  DEST.y=V1.z*V2.x-V1.x*V2.z;  \
  DEST.z=V1.x*V2.y-V1.y*V2.x;

#define GEOM_SWAP(A,B,AUX) {\
  AUX = B; \
  B   = A; \
  A = AUX; \
}

#define GEOM_CHECK_MIN_MAX3D(P1,Q1,R1,P2,Q2,R2) { \
  GEOM_SUB3D(V1,P2,Q1) \
  GEOM_SUB3D(V2,P1,Q1) \
  GEOM_CROSS3D(N,V1,V2) \
  GEOM_SUB3D(V1,Q2,Q1) \
\
  if (GEOM_DOT3D(V1,N) > 0.0f) \
    return (FALSE); \
\
  GEOM_SUB3D(V1,P2,P1) \
  GEOM_SUB3D(V2,R1,P1) \
  GEOM_CROSS3D(N,V1,V2) \
  GEOM_SUB3D(V1,R2,P1) \
\
  if (GEOM_DOT3D(V1,N) > 0.0f) \
    return (FALSE); \
  else \
    return (TRUE); \
}

#define GEOM_TRI_TRI_GENERAL3D(P1,Q1,R1,P2,Q2,R2,DP2,DQ2,DR2) { \
  /* Permutation in a canonical form of T2's vertices */ \
  if (DP2 > 0.0f) { \
    if (DQ2 > 0.0f) \
      GEOM_CHECK_MIN_MAX3D(P1,R1,Q1,R2,P2,Q2) \
    else if (DR2 > 0.0f) \
      GEOM_CHECK_MIN_MAX3D(P1,R1,Q1,Q2,R2,P2) \
    else \
      GEOM_CHECK_MIN_MAX3D(P1,Q1,R1,P2,Q2,R2) \
  } \
  else if (DP2 < 0.0f) { \
    if (DQ2 < 0.0f) \
      GEOM_CHECK_MIN_MAX3D(P1,Q1,R1,R2,P2,Q2) \
    else if (DR2 < 0.0f) \
      GEOM_CHECK_MIN_MAX3D(P1,Q1,R1,Q2,R2,P2) \
    else \
      GEOM_CHECK_MIN_MAX3D(P1,R1,Q1,P2,Q2,R2) \
  } \
  else { \
    if (DQ2 < 0.0f) { \
      if (DR2 >= 0.0f) \
        GEOM_CHECK_MIN_MAX3D(P1,R1,Q1,Q2,R2,P2) \
      else \
        GEOM_CHECK_MIN_MAX3D(P1,Q1,R1,P2,Q2,R2) \
    } \
    else if (DQ2 > 0.0f) { \
      if (DR2 > 0.0f) \
        GEOM_CHECK_MIN_MAX3D(P1,R1,Q1,P2,Q2,R2) \
      else \
        GEOM_CHECK_MIN_MAX3D(P1,Q1,R1,Q2,R2,P2) \
    } \
    else { \
      if (DR2 > 0.0f) \
        GEOM_CHECK_MIN_MAX3D(P1,Q1,R1,R2,P2,Q2) \
      else \
        GEOM_CHECK_MIN_MAX3D(P1,R1,Q1,R2,P2,Q2) \
    } \
  } \
}

#define  GEOM_SEGM_TRI_OVERLAP3D(A,B,P,Q,R) \
{ \
  GEOM_SUB3D(V3,A,B) \
\
  GEOM_SUB3D(V1,P,B) \
  GEOM_SUB3D(V2,Q,B) \
  GEOM_CROSS3D(N,V1,V2) \
  if (GEOM_DOT3D(V3,N) >= 0.0f) {\
    GEOM_SUB3D(V1,Q,B) \
    GEOM_SUB3D(V2,R,B) \
    GEOM_CROSS3D(N,V1,V2) \
    if (GEOM_DOT3D(V3,N) >= 0.0f) {\
      GEOM_SUB3D(V1,R,B) \
      GEOM_SUB3D(V2,P,B) \
      GEOM_CROSS3D(N,V1,V2) \
      if (GEOM_DOT3D(V3,N) >= 0.0f) \
        return (TRUE); \
    } \
  } \
}

#define GEOM_TRI_TRI_SHARING_VERTEX3D(P1,Q1,R1,P2,Q2,R2,DQ1,DR1,DQ2,DR2) \
{ \
  /* Test line segment q1-r1 (or r1-q1) against T2 */ \
  if (DQ1 > 0.0f) { \
    if (DR1 > 0.0f) \
      return (FALSE); \
    else \
      GEOM_SEGM_TRI_OVERLAP3D(Q1, R1, P2, Q2, R2) \
  } \
  else if (DQ1 < 0.0f) { \
    if (DR1 < 0.0f) \
      return (FALSE); \
    else \
      GEOM_SEGM_TRI_OVERLAP3D(R1, Q1, P2, Q2, R2) \
  } \
  else { \
    if (DR1 < 0.0f) \
      GEOM_SEGM_TRI_OVERLAP3D(Q1, R1, P2, Q2, R2) \
    else \
      GEOM_SEGM_TRI_OVERLAP3D(R1, Q1, P2, Q2, R2) \
  } \
\
  /* Test line segment q2-r2 (or r2-q2) against T2 */ \
  if (DQ2 > 0.0f) { \
    if (DR2 > 0.0f) \
      return (FALSE); \
    else \
      GEOM_SEGM_TRI_OVERLAP3D(Q2, R2, P1, Q1, R1) \
  } \
  else if (DQ2 < 0.0f) { \
    if (DR2 < 0.0f) \
      return (FALSE); \
    else \
      GEOM_SEGM_TRI_OVERLAP3D(R2, Q2, P1, Q1, R1) \
  } \
  else { \
    if (DR2 < 0.0f) \
      GEOM_SEGM_TRI_OVERLAP3D(Q2, R2, P1, Q1, R1) \
    else \
      GEOM_SEGM_TRI_OVERLAP3D(R2, Q2, P1, Q1, R1) \
  } \
\
  return (FALSE); \
}


/* some 2D macros */

#define ORIENT_2D(a,b,c)  ((a->x-c->x)*(b->y-c->y)-(a->y-c->y)*(b->x-c->x))

#define GEOM_TRI_TRI_SHARING_EDGE2D(P1,Q1,R1,P2,Q2,R2,NSHRDNDS) \
{ \
  if (NSHRDNDS[4] > 0) { \
    if (NSHRDNDS[5] > 0) \
      if ((ORIENT_2D(P1,Q1,R1)*ORIENT_2D(P1,Q1,R2)) < 0.0f) \
        return (FALSE); \
      else \
        return (TRUE); \
    else \
      if ((ORIENT_2D(P1,Q1,R1)*ORIENT_2D(P1,Q1,Q2)) < 0.0f) \
        return (FALSE); \
      else \
        return (TRUE); \
  } \
  else { \
    if ((ORIENT_2D(P1,Q1,R1)*ORIENT_2D(P1,Q1,P2)) < 0.0f) \
      return (FALSE); \
    else \
      return (TRUE); \
  } \
}

#define GEOM_TRI_TRI_CHECK_EDGES2D(P1,Q1,R1,P2,Q2,R2) \
{ \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) { \
    if (ORIENT_2D(P2,Q2,Q1) >= 0.0f) \
      return (TRUE); \
    else \
      if (ORIENT_2D(P2,Q2,R1) >= 0.0f) \
        return (TRUE); \
      else \
        return (FALSE); \
  } \
  else { \
    if (ORIENT_2D(P2,Q2,Q1) >= 0.0f) \
      return (FALSE); \
    else \
      if (ORIENT_2D(P2,Q2,R1) >= 0.0f) \
        return (TRUE); \
      else \
        return (FALSE); \
  } \
}

#define GEOM_TRI_TRI_SHARING_VERTEX2D(P1,Q1,R1,P2,Q2,R2,NSHRDNDS) \
{ \
  if (NSHRDNDS[1] > 0) { \
    if (NSHRDNDS[4] > 0) \
      GEOM_TRI_TRI_CHECK_EDGES2D(P1,Q1,R1,P2,Q2,R2) \
    else if (NSHRDNDS[5] > 0) \
      GEOM_TRI_TRI_CHECK_EDGES2D(P1,Q1,R1,Q2,R2,P2) \
    else \
      GEOM_TRI_TRI_CHECK_EDGES2D(P1,Q1,R1,R2,P2,Q2) \
  } \
  else if (NSHRDNDS[2] > 0) { \
    if (NSHRDNDS[4] > 0) \
      GEOM_TRI_TRI_CHECK_EDGES2D(Q1,R1,P1,P2,Q2,R2) \
    else if (NSHRDNDS[5] > 0) \
      GEOM_TRI_TRI_CHECK_EDGES2D(Q1,R1,P1,Q2,R2,P2) \
    else \
      GEOM_TRI_TRI_CHECK_EDGES2D(Q1,R1,P1,R2,P2,Q2) \
  } \
  else { \
    if (NSHRDNDS[4] > 0) \
      GEOM_TRI_TRI_CHECK_EDGES2D(R1,P1,Q1,P2,Q2,R2) \
    else if (NSHRDNDS[5] > 0) \
      GEOM_TRI_TRI_CHECK_EDGES2D(R1,P1,Q1,Q2,R2,P2) \
    else \
      GEOM_TRI_TRI_CHECK_EDGES2D(Q1,P1,Q1,R2,P2,Q2) \
  } \
}

#define GEOM_TRI_TRI_GENERAL_VERTEX2D(P1,Q1,R1,P2,Q2,R2) \
{ \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) { \
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f) { \
      if (ORIENT_2D(P1,P2,Q1) > 0.0f) { \
        if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) \
          return (TRUE); \
        else \
          return (FALSE); \
      } \
      else { \
        if (ORIENT_2D(P1,P2,R1) >= 0.0f) { \
          if (ORIENT_2D(Q1,R1,P2) >= 0.0f) \
            return (TRUE); \
          else \
            return (FALSE); \
        } \
        else \
          return (FALSE); \
      } \
    } \
    else { \
      if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) { \
        if (ORIENT_2D(R2,Q2,R1) <= 0.0f) { \
          if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) \
            return (TRUE); \
          else \
            return (FALSE); \
        } \
        else \
          return (FALSE); \
      } \
      else \
        return (FALSE); \
    } \
  } \
  else { \
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) { \
      if (ORIENT_2D(Q1,R1,R2) >= 0.0f) { \
        if (ORIENT_2D(P1,P2,R1) >= 0.0f) \
          return (TRUE); \
        else \
          return (FALSE); \
      } \
      else  { \
        if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) { \
          if (ORIENT_2D(R2,R1,Q2) >= 0.0f) \
            return (TRUE); \
          else \
            return (FALSE); \
        } \
        else \
          return (FALSE); \
      } \
    } \
    else \
      return (FALSE); \
  } \
}

#define GEOM_TRI_TRI_GENERAL_EDGE2D(P1,Q1,R1,P2,Q2,R2) \
{ \
  if (ORIENT_2D(R2,P2,Q1) >= 0.0f) { \
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
      if (ORIENT_2D(P1,Q1,R2) >= 0.0f) \
        return (TRUE); \
      else \
        return (FALSE); \
    } \
    else { \
      if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
        if (ORIENT_2D(R1,P1,P2) >= 0.0f) \
          return (TRUE); \
        else \
          return (FALSE); \
      } \
      else \
        return (FALSE); \
    } \
  } \
  else { \
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) { \
      if (ORIENT_2D(P1,P2,R1) >= 0.0f) { \
        if (ORIENT_2D(P1,R1,R2) >= 0.0f) \
          return (TRUE); \
        else { \
          if (ORIENT_2D(Q1,R1,R2) >= 0.0f) \
            return (TRUE); \
          else \
            return (FALSE); \
        } \
      } \
      else \
        return (FALSE); \
    } \
    else \
      return (FALSE); \
  } \
}

#define GEOM_TRI_TRI_CCW2D(P1,Q1,R1,P2,Q2,R2) \
{ \
  if (ORIENT_2D(P2,Q2,P1) >= 0.0f) { \
    if (ORIENT_2D(Q2,R2,P1) >= 0.0f) { \
      if (ORIENT_2D(R2,P2,P1) >= 0.0f) \
        return (TRUE); \
      else \
        GEOM_TRI_TRI_GENERAL_EDGE2D(P1,Q1,R1,P2,Q2,R2) \
    } \
    else { \
      if (ORIENT_2D(R2,P2,P1) >= 0.0f) \
        GEOM_TRI_TRI_GENERAL_EDGE2D(P1,Q1,R1,R2,P2,Q2) \
      else \
        GEOM_TRI_TRI_GENERAL_VERTEX2D(P1,Q1,R1,P2,Q2,R2) \
    } \
  } \
  else { \
    if (ORIENT_2D(Q2,R2,P1) >= 0.0f) { \
      if (ORIENT_2D(R2,P2,P1) >= 0.0f) \
        GEOM_TRI_TRI_GENERAL_EDGE2D(P1,Q1,R1,Q2,R2,P2) \
      else \
        GEOM_TRI_TRI_GENERAL_VERTEX2D(P1,Q1,R1,Q2,R2,P2) \
    } \
    else \
      GEOM_TRI_TRI_GENERAL_VERTEX2D(P1,Q1,R1,R2,P2,Q2) \
  } \
}

/* Function Prototypes */
Boolean GEOM_In_tetrahedron(Point3D *p, Point3D *nd1, Point3D *nd2,
              Point3D *nd3, Point3D *nd4);




Boolean GEOM_Tri_tri_overlap_test3D(Point3D *, Point3D *, Point3D *,
                  Point3D *, Point3D *, Point3D *,
                  unsigned long *);

Boolean GEOM_Tri_tri_noncoplanar3D(Point3D *, Point3D *, Point3D *,
                   Point3D *, Point3D *, Point3D *,
                   double,    double,    double,
                   double,    double,    double,
                   unsigned long *);

Boolean GEOM_Tri_tri_coplanar3D(Point3D *, Point3D *, Point3D *,
                Point3D *, Point3D *, Point3D *,
                Point3D *,
                unsigned long *);

Boolean GEOM_Tri_tri_overlap_test2D(Point2D *, Point2D *, Point2D *,
                  Point2D *, Point2D *, Point2D *,
                  unsigned long *);

double GEOM_Point_tri_sqrdistance(const Point3D *, const Point3D *,
                  const Point3D *, const Point3D *);

#endif /* GEOM_H */
