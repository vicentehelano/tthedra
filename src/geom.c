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
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <geom.h>
#include <float.h>

const double epsilon = GEOM_MAX2((0.0e+000),((100.0e+000)*DBL_EPSILON));

double GEOM_Angle3D(Point3D *vec1, Point3D *vec2, double sqtol)
{
  double  dot;
  double  cosine;
  double  sqrdlen1, sqrdlen2;

  dot    = GEOM_DOT3D((*vec1),(*vec2));
  sqrdlen1 = GEOM_DOT3D((*vec1),(*vec1));
  sqrdlen2 = GEOM_DOT3D((*vec2),(*vec2));

  /* Check for minimum geometric tolerance */
  if (sqrdlen1 > sqtol && sqrdlen2 > sqtol) {
    cosine = dot/sqrt(sqrdlen1*sqrdlen2);
    if (GEOM_ABSFLOAT(cosine) > (1.0e+000 - epsilon))
      cosine = GEOM_SIGNDOUBLE((1.0e+000),cosine);
    return (acos(cosine));
  }
  else
    return (PI);
}

Boolean GEOM_In_tetrahedron(Point3D *p, Point3D *nd1, Point3D *nd2,
              Point3D *nd3, Point3D *nd4)
{
  double  dp1, dp2, dp3, dp4;
  Point3D v1, v2, v3, n;

  /* Face 04: triangle (nd1, nd3, nd2) */
  GEOM_SUB3D(v1,(*nd1),(*nd2))
  GEOM_SUB3D(v2,(*nd3),(*nd2))
  GEOM_CROSS3D(n,v1,v2)

  GEOM_SUB3D(v3,(*p),(*nd2))
  dp1 = GEOM_DOT3D(v3,n);

  if (GEOM_ABSDOUBLE(dp1) < EPSILON)
    dp1 = 0.0e+000;

  /* Common vector */
  GEOM_SUB3D(v3,(*p),(*nd4))

  /* Face 01: triangle (nd1, nd2, nd4) */
  GEOM_SUB3D(v1,(*nd1),(*nd4))
  GEOM_SUB3D(v2,(*nd2),(*nd4))
  GEOM_CROSS3D(n,v1,v2)
  dp2 = GEOM_DOT3D(v3,n);

  /* Robustness check */
  if (GEOM_ABSDOUBLE(dp2) < EPSILON)
    dp2 = 0.0e+000;

  /* Face 02: triangle (nd2, nd3, nd4) */
  GEOM_SUB3D(v1,(*nd2),(*nd4))
  GEOM_SUB3D(v2,(*nd3),(*nd4))
  GEOM_CROSS3D(n,v1,v2)
  dp3 = GEOM_DOT3D(v3,n);

  /* Robustness check */
  if (GEOM_ABSDOUBLE(dp3) < EPSILON)
    dp3 = 0.0e+000;

  /* Face 03: triangle (nd3, nd1, nd4) */
  GEOM_SUB3D(v1,(*nd3),(*nd4))
  GEOM_SUB3D(v2,(*nd1),(*nd4))
  GEOM_CROSS3D(n,v1,v2)
  dp4 = GEOM_DOT3D(v3,n);

  /* Robustness check */
  if (GEOM_ABSDOUBLE(dp4) < EPSILON)
    dp4 = 0.0e+000;

  /* Is 'p' in tetrahedron? */
  if (dp1 >= 0.0e+000 && dp2 >= 0.0e+000 &&
    dp3 >= 0.0e+000 && dp4 >= 0.0e+000)
    return (TRUE);
  else
    return (FALSE);
}

Boolean GEOM_Tri_tri_overlap_test3D(Point3D *p1, Point3D *q1, Point3D *r1,
                  Point3D *p2, Point3D *q2, Point3D *r2,
                  unsigned long *nshrdnds)
{
  /* Definitions of local variables */
  double dp1, dq1, dr1, dp2, dq2, dr2;
  Point3D v1, v2;
  Point3D n1, n2; 

  /* Compute the normal vector of T2 */
  GEOM_SUB3D(v1,(*p2),(*r2))
  GEOM_SUB3D(v2,(*q2),(*r2))
  GEOM_CROSS3D(n2,v1,v2)
  /* Compute distance signs of p1, q1 and r1 to the plane of T2 */
  GEOM_SUB3D(v1,(*p1),(*r2))
  dp1 = GEOM_DOT3D(v1,n2);
  GEOM_SUB3D(v1,(*q1),(*r2))
  dq1 = GEOM_DOT3D(v1,n2);
  GEOM_SUB3D(v1,(*r1),(*r2))
  dr1 = GEOM_DOT3D(v1,n2);

  #ifdef EPSILON /* then apply epsilon test */
    if(GEOM_ABSDOUBLE(dp1) < EPSILON)
      dp1 = 0.0e+000;
    if(GEOM_ABSDOUBLE(dq1) < EPSILON)
      dq1 = 0.0e+000;
    if(GEOM_ABSDOUBLE(dr1) < EPSILON)
      dr1 = 0.0e+000;
  #endif

  /* Check for coplanarity */
  if (dp1 == 0.0e+000 && dq1 == 0.0e+000 && dr1 == 0.0e+000)
    return (GEOM_Tri_tri_coplanar3D(p1,q1,r1,p2,q2,r2,&n2,nshrdnds));

  /* Check if p1, q1 and r1 are on the same side of the plane of T2 */
  if (((dp1 * dq1) > 0.0e+000) && ((dp1 * dr1) > 0.0e+000)) 
    return (FALSE);
  
  /* Compute the normal vector of T1 */
  GEOM_SUB3D(v1,(*q1),(*p1))
  GEOM_SUB3D(v2,(*r1),(*p1))
  GEOM_CROSS3D(n1,v1,v2)
  /* Compute distance signs of p2, q2 and r2 to the plane of T1 */
  GEOM_SUB3D(v1,(*p2),(*r1))
  dp2 = GEOM_DOT3D(v1,n1);
  GEOM_SUB3D(v1,(*q2),(*r1))
  dq2 = GEOM_DOT3D(v1,n1);
  GEOM_SUB3D(v1,(*r2),(*r1))
  dr2 = GEOM_DOT3D(v1,n1);

  #ifdef EPSILON /* then apply epsilon test */
    if(GEOM_ABSDOUBLE(dp2) < EPSILON)
      dp2 = 0.0e+000;
    if(GEOM_ABSDOUBLE(dq2) < EPSILON)
      dq2 = 0.0e+000;
    if(GEOM_ABSDOUBLE(dr2) < EPSILON)
      dr2 = 0.0e+000;
  #endif

  /* Check for coplanarity */
  if (dp2 == 0.0e+000 && dq2 == 0.0e+000 && dr2 == 0.0e+000)
    return (GEOM_Tri_tri_coplanar3D(p1,q1,r1,p2,q2,r2,&n1,nshrdnds));

  /* Check if p2, q2 and r2 are on the same side of the plane of T1 */
  if (((dp2 * dq2) > 0.0e+000) && ((dp2 * dr2) > 0.0e+000))
    return (FALSE);

  /* Non-coplanar triangles */
  return (GEOM_Tri_tri_noncoplanar3D(p1,q1,r1,p2,q2,r2,dp1,dq1,dr1,
      dp2,dq2,dr2,nshrdnds));
}

Boolean GEOM_Tri_tri_noncoplanar3D(Point3D *p1, Point3D *q1, Point3D *r1,
                   Point3D *p2, Point3D *q2, Point3D *r2,
                   double  dp1, double  dq1, double  dr1,
                   double  dp2, double  dq2, double  dr2,
                   unsigned long *nshrdnds)
{
  Point3D V1, V2, V3, N;

  if (nshrdnds[0] == EMPTY) { /* Not sharing vertex or edge */
    /* Permutation in a canonical form of T1's vertices */
    if (dp1 > 0.0e+000) {
      if (dq1 > 0.0e+000)
        GEOM_TRI_TRI_GENERAL3D((*r1),(*p1),(*q1),(*p2),(*r2),(*q2),dp2,dr2,dq2)
      else if (dr1 > 0.0e+000)
        GEOM_TRI_TRI_GENERAL3D((*q1),(*r1),(*p1),(*p2),(*r2),(*q2),dp2,dr2,dq2)
      else
        GEOM_TRI_TRI_GENERAL3D((*p1),(*q1),(*r1),(*p2),(*q2),(*r2),dp2,dq2,dr2)
    }
    else if (dp1 < 0.0e+000) {
      if (dq1 < 0.0e+000)
        GEOM_TRI_TRI_GENERAL3D((*r1),(*p1),(*q1),(*p2),(*q2),(*r2),dp2,dq2,dr2)
      else if (dr1 < 0.0e+000)
        GEOM_TRI_TRI_GENERAL3D((*q1),(*r1),(*p1),(*p2),(*q2),(*r2),dp2,dq2,dr2)
      else
        GEOM_TRI_TRI_GENERAL3D((*p1),(*q1),(*r1),(*p2),(*r2),(*q2),dp2,dr2,dq2)
    }
    else {
      if (dq1 < 0.0e+000) {
        if (dr1 >= 0.0e+000)
          GEOM_TRI_TRI_GENERAL3D((*q1),(*r1),(*p1),(*p2),(*r2),(*q2),dp2,dr2,dq2)
        else
          GEOM_TRI_TRI_GENERAL3D((*p1),(*q1),(*r1),(*p2),(*q2),(*r2),dp2,dq2,dr2)
      }
      else if (dq1 > 0.0e+000) {
        if (dr1 > 0.0e+000)
          GEOM_TRI_TRI_GENERAL3D((*p1),(*q1),(*r1),(*p2),(*r2),(*q2),dp2,dr2,dq2)
        else
          GEOM_TRI_TRI_GENERAL3D((*q1),(*r1),(*p1),(*p2),(*q2),(*r2),dp2,dq2,dr2)
      }
      else {
        if (dr1 > 0.0e+000)
          GEOM_TRI_TRI_GENERAL3D((*r1),(*p1),(*q1),(*p2),(*q2),(*r2),dp2,dq2,dr2)
        else
          GEOM_TRI_TRI_GENERAL3D((*r1),(*p1),(*q1),(*p2),(*r2),(*q2),dp2,dr2,dq2)
      }
    }
  }
  else if (nshrdnds[0] == VERTEX) { /* Sharing vertex */
    /* Permutation such that the shared nodes are moved to p1 and p2 */
    if (nshrdnds[1] > 0) {
      if (nshrdnds[4] > 0)
        GEOM_TRI_TRI_SHARING_VERTEX3D((*p1),(*q1),(*r1),(*p2),(*q2),(*r2),dq1,dr1,dq2,dr2)
      else if (nshrdnds[5] > 0)
        GEOM_TRI_TRI_SHARING_VERTEX3D((*p1),(*q1),(*r1),(*q2),(*r2),(*p2),dq1,dr1,dr2,dp2)
      else 
        GEOM_TRI_TRI_SHARING_VERTEX3D((*p1),(*q1),(*r1),(*r2),(*p2),(*q2),dq1,dr1,dp2,dq2)
    }
    else if (nshrdnds[2] > 0) {
      if (nshrdnds[4] > 0)
        GEOM_TRI_TRI_SHARING_VERTEX3D((*q1),(*r1),(*p1),(*p2),(*q2),(*r2),dr1,dp1,dq2,dr2)
      else if (nshrdnds[5] > 0)
        GEOM_TRI_TRI_SHARING_VERTEX3D((*q1),(*r1),(*p1),(*q2),(*r2),(*p2),dr1,dp1,dr2,dp2)
      else 
        GEOM_TRI_TRI_SHARING_VERTEX3D((*q1),(*r1),(*p1),(*r2),(*p2),(*q2),dr1,dp1,dp2,dq2)
    }
    else {
      if (nshrdnds[4] > 0)
        GEOM_TRI_TRI_SHARING_VERTEX3D((*r1),(*p1),(*q1),(*p2),(*q2),(*r2),dp1,dq1,dq2,dr2)
      else if (nshrdnds[5] > 0)
        GEOM_TRI_TRI_SHARING_VERTEX3D((*r1),(*p1),(*q1),(*q2),(*r2),(*p2),dp1,dq1,dr2,dp2)
      else 
        GEOM_TRI_TRI_SHARING_VERTEX3D((*r1),(*p1),(*q1),(*r2),(*p2),(*q2),dp1,dq1,dp2,dq2)
    }
  }
  else { /* Sharing an edge */
    return (FALSE);
  }
}

Boolean GEOM_Tri_tri_coplanar3D(Point3D *p1, Point3D *q1, Point3D *r1,
                Point3D *p2, Point3D *q2, Point3D *r2,
                Point3D *n,
                unsigned long *nshrdnds)
{
  Point2D P1,Q1,R1;
  Point2D P2,Q2,R2;
  double nx, ny, nz;
  unsigned long NSHRDNDS[8];
  
  nx = GEOM_ABSDOUBLE(n->x);
  ny = GEOM_ABSDOUBLE(n->y);
  nz = GEOM_ABSDOUBLE(n->z);
  
  /* Projection of the triangles in 3D onto 2D such that the area of
  the projection is maximized. */
  NSHRDNDS[0] = nshrdnds[0];
  NSHRDNDS[7] = nshrdnds[7];

  if ((nx > nz) && (nx >= ny)) {
    /* Project onto plane YZ */
    P1.x = q1->z; P1.y = q1->y;
    Q1.x = p1->z; Q1.y = p1->y;
    R1.x = r1->z; R1.y = r1->y;
    NSHRDNDS[1] = nshrdnds[2];
    NSHRDNDS[2] = nshrdnds[1];
    NSHRDNDS[3] = nshrdnds[3];
    
    P2.x = q2->z; P2.y = q2->y;
    Q2.x = p2->z; Q2.y = p2->y;
    R2.x = r2->z; R2.y = r2->y;
    NSHRDNDS[4] = nshrdnds[5];
    NSHRDNDS[5] = nshrdnds[4];
    NSHRDNDS[6] = nshrdnds[6];
  }
  else if ((ny > nz) && (ny >= nx)) {
    /* Project onto plane XZ */
    P1.x = q1->x; P1.y = q1->z;
    Q1.x = p1->x; Q1.y = p1->z;
    R1.x = r1->x; R1.y = r1->z;
    NSHRDNDS[1] = nshrdnds[2];
    NSHRDNDS[2] = nshrdnds[1];
    NSHRDNDS[3] = nshrdnds[3];
    
    P2.x = q2->x; P2.y = q2->z;
    Q2.x = p2->x; Q2.y = p2->z;
    R2.x = r2->x; R2.y = r2->z;
    NSHRDNDS[4] = nshrdnds[5];
    NSHRDNDS[5] = nshrdnds[4];
    NSHRDNDS[6] = nshrdnds[6];
  }
  else {
    /* Project onto plane XY */
    P1.x = p1->x; P1.y = p1->y;
    Q1.x = q1->x; Q1.y = q1->y;
    R1.x = r1->x; R1.y = r1->y;
    NSHRDNDS[1] = nshrdnds[1];
    NSHRDNDS[2] = nshrdnds[2];
    NSHRDNDS[3] = nshrdnds[3];
    
    P2.x = p2->x; P2.y = p2->y;
    Q2.x = q2->x; Q2.y = q2->y;
    R2.x = r2->x; R2.y = r2->y;
    NSHRDNDS[4] = nshrdnds[4];
    NSHRDNDS[5] = nshrdnds[5];
    NSHRDNDS[6] = nshrdnds[6];
  }
  return (GEOM_Tri_tri_overlap_test2D(&P1,&Q1,&R1,&P2,&Q2,&R2,NSHRDNDS));
}

Boolean GEOM_Tri_tri_overlap_test2D(Point2D *p1, Point2D *q1, Point2D *r1,
                  Point2D *p2, Point2D *q2, Point2D *r2,
                  unsigned long *nshrdnds)
{
  unsigned long NSHRDNDS[8];

  NSHRDNDS[0] = nshrdnds[0];
  NSHRDNDS[7] = nshrdnds[7];

  if (nshrdnds[0] == 0) {
    if (ORIENT_2D(p1,q1,r1) < 0.0e+000) {
      if (ORIENT_2D(p2,q2,r2) < 0.0e+000)
        GEOM_TRI_TRI_CCW2D(p1,r1,q1,p2,r2,q2)
      else
        GEOM_TRI_TRI_CCW2D(p1,r1,q1,p2,q2,r2)
    }
    else {
      if (ORIENT_2D(p2,q2,r2) < 0.0e+000)
        GEOM_TRI_TRI_CCW2D(p1,q1,r1,p2,r2,q2)
      else
        GEOM_TRI_TRI_CCW2D(p1,q1,r1,p2,q2,r2)
    }
  }
  else if (nshrdnds[0] == 1) {
    if (ORIENT_2D(p1,q1,r1) < 0.0e+000) {
      if (ORIENT_2D(p2,q2,r2) < 0.0e+000) {
        NSHRDNDS[1] = nshrdnds[1], NSHRDNDS[2] = nshrdnds[3], NSHRDNDS[3] = nshrdnds[2];
        NSHRDNDS[4] = nshrdnds[4], NSHRDNDS[5] = nshrdnds[6], NSHRDNDS[6] = nshrdnds[5];
        GEOM_TRI_TRI_SHARING_VERTEX2D(p1,r1,q1,p2,r2,q2,NSHRDNDS)
      }
      else {
        NSHRDNDS[1] = nshrdnds[1], NSHRDNDS[2] = nshrdnds[3], NSHRDNDS[3] = nshrdnds[2];
        NSHRDNDS[4] = nshrdnds[4], NSHRDNDS[5] = nshrdnds[5], NSHRDNDS[6] = nshrdnds[6];
        GEOM_TRI_TRI_SHARING_VERTEX2D(p1,r1,q1,p2,q2,r2,NSHRDNDS)
      }
    }
    else {
      if (ORIENT_2D(p2,q2,r2) < 0.0e+000) {
        NSHRDNDS[1] = nshrdnds[1], NSHRDNDS[2] = nshrdnds[2], NSHRDNDS[3] = nshrdnds[3];
        NSHRDNDS[4] = nshrdnds[4], NSHRDNDS[5] = nshrdnds[6], NSHRDNDS[6] = nshrdnds[5];
        GEOM_TRI_TRI_SHARING_VERTEX2D(p1,q1,r1,p2,r2,q2,NSHRDNDS)
      }
      else {
        NSHRDNDS[1] = nshrdnds[1], NSHRDNDS[2] = nshrdnds[2], NSHRDNDS[3] = nshrdnds[3];
        NSHRDNDS[4] = nshrdnds[4], NSHRDNDS[5] = nshrdnds[5], NSHRDNDS[6] = nshrdnds[6];
        GEOM_TRI_TRI_SHARING_VERTEX2D(p1,q1,r1,p2,q2,r2,NSHRDNDS)
      }
    }
  }
  else if (nshrdnds[0] == 2) {
    if (nshrdnds[1] > 0) {
      if (nshrdnds[2] > 0)
        GEOM_TRI_TRI_SHARING_EDGE2D(p1,q1,r1,p2,q2,r2,nshrdnds)
      else
        GEOM_TRI_TRI_SHARING_EDGE2D(r1,p1,q1,p2,q2,r2,nshrdnds)
    }
    else
      GEOM_TRI_TRI_SHARING_EDGE2D(q1,r1,p1,p2,q2,r2,nshrdnds)
  }
  else { /* Triangulos coincidentes */
    return (FALSE);
  }
}

double GEOM_Point_tri_sqrdistance(const Point3D *pt,  const Point3D *nd1,
                  const Point3D *nd2, const Point3D *nd3)
{
  Point3D diff, v12, v13;
  double  a11, a12, a22, b1, b2, b3;
  double  s, t, det, sqrddist;
  double  tmp1, tmp2, numer, denom;
  double  invdet;

  GEOM_SUB3D(diff,(*nd1),(*pt))
  GEOM_SUB3D(v12,(*nd2),(*nd1))
  GEOM_SUB3D(v13,(*nd3),(*nd1))

  a11 = GEOM_DOT3D(v12,v12);
  a12 = GEOM_DOT3D(v12,v13);
  a22 = GEOM_DOT3D(v13,v13);

  b1  = GEOM_DOT3D(diff,v12);
  b2  = GEOM_DOT3D(diff,v13);
  b3   = GEOM_DOT3D(diff,diff);

  s   = a12*b2 - a22*b1;
    t   = a12*b1 - a11*b2;
  det = GEOM_ABSDOUBLE((a11*a22 - a12*a12));

    if ((s + t) <= det) {
        if (s < 0.0e+000) {
            if (t < 0.0e+000) { /* region 4 */
                if (b1 < 0.0e+000) {
                    if (-b1 >= a11)
                        sqrddist = a11 + (2.0e+000)*b1 + b3;
                    else {
                        s = -b1/a11;
                        sqrddist = b1*s + b3;
                    }
                }
                else {
                    if (b2 >= 0.0e+000)
                        sqrddist = b3;
          else if (-b2 >= a22)
                        sqrddist = a22 + (2.0e+000)*b2 + b3;
                    else {
                        t = -b2/a22;
                        sqrddist = b2*t + b3;
                    }
                }
            }
    else { /* region 3 */
                if (b2 >= 0.0e+000)
                    sqrddist = b3;
                else if (-b2 >= a22)
                    sqrddist = a22 + (2.0e+000)*b2 + b3;
                else {
                    t = -b2/a22;
                    sqrddist = b2*t + b3;
                }
            }
        }
        else if (t < 0.0e+000) {  /* region 5 */
            if (b1 >= 0.0e+000)
                sqrddist = b3;
            else if (-b1 >= a11)
                sqrddist = a11 + (2.0e+000)*b1 + b3;
            else {
                s = -b1/a11;
                sqrddist = b1*s + b3;
            }
        }
        else {  /* region 0 (minimum at interior point)*/
            invdet = (1.0e+000)/det;
            s *= invdet;
            t *= invdet;
            sqrddist = s*(a11*s + a12*t + (2.0e+000)*b1) +
                t*(a12*s + a22*t + (2.0e+000)*b2) + b3;
        }
    }
    else {
        if (s < 0.0e+000) { /* region 2 */
            tmp1 = a12 + b1;
            tmp2 = a22 + b2;
            if (tmp2 > tmp1) {
                numer = tmp2 - tmp1;
                denom = a11 - (2.0e+000)*a12 + a22;
                if (numer >= denom)
                    sqrddist = a11 + (2.0e+000)*b1 + b3;
                else {
                    s = numer/denom;
                    t = 1.0e+000 - s;
                    sqrddist = s*(a11*s + a12*t + (2.0e+000)*b1) +
            t*(a12*s + a22*t + (2.0e+000)*b2) + b3;
                }
            }
            else {
                if (tmp2 <= 0.0e+000)
                    sqrddist = a22 + (2.0e+000)*b2 + b3;
                else if (b2 >= 0.0e+000)
                    sqrddist = b3;
                else {
                    t = -b2/a22;
                    sqrddist = b2*t + b3;
                }
            }
        }
        else if (t < 0.0e+000) {  /* region 6 */
            tmp1 = a12 + b2;
            tmp2 = a11 + b1;
            if (tmp2 > tmp1) {
                numer = tmp2 - tmp1;
                denom = a11 - (2.0e+000)*a12 + a22;
                if (numer >= denom)
                    sqrddist = a22 + (2.0e+000)*b2 + b3;
                else {
                    t = numer/denom;
                    s = 1.0e+000 - t;
                    sqrddist = s*(a11*s + a12*t + (2.0e+000)*b1) +
            t*(a12*s + a22*t + (2.0e+000)*b2) + b3;
                }
            }
            else {
                if (tmp2 <= 0.0e+000)
                    sqrddist = a11 + (2.0e+000)*b1 + b3;
                else if (b1 >= 0.0e+000)
                    sqrddist = b3;
                else {
                    s = -b1/a11;
                    sqrddist = b1*s + b3;
                }
            }
        }
        else {  /* region 1 */
            numer = a22 + b2 - a12 - b1;
            if (numer <= 0.0e+000)
                sqrddist = a22 + (2.0e+000)*b2 + b3;
            else {
                denom = a11 - (2.0e+000)*a12 + a22;
                if (numer >= denom)
                    sqrddist = a11 + (2.0e+000)*b1 + b3;
                else {
                    s = numer/denom;
                    t = 1.0e+000 - s;
                    sqrddist = s*(a11*s + a12*t + (2.0e+000)*b1) +
            t*(a12*s + a22*t + (2.0e+000)*b2) + b3;
                }
            }
        }
    }
    return (GEOM_ABSDOUBLE(sqrddist));
}
