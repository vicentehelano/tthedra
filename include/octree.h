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
#ifndef OCTREE_H
#define OCTREE_H

#include <geom.h>

typedef enum octree_const {
	BUCKETSIZE		= 8,
	INIT_NLEVELS	= 6,
	MAX_NLEVELS		= 11,
	INIT_MAXCDTNDS	= 500
} Octree_const;

typedef enum octant {
	OCT1 = 1, OCT2 = 2, OCT3 = 3, OCT4 = 4,
	OCT5 = 5, OCT6 = 6, OCT7 = 7, OCT8 = 8
} Octant;


/*
 * Octree nodes structure
 *
 */
typedef struct octree_node {
	int				status;
	unsigned long	*point;
} Octree_node;


/*
 * Candidate points structure
 *
 */
typedef struct	candidate_points {
	unsigned long	npts;
	unsigned long	maxnpts;
	unsigned long	*point;
} Candidate_points;

/*
 * Candidate nodes structure
 *
 */
typedef struct	candidate_nodes {
	unsigned long	nnds;
	unsigned long	maxnnds;
	unsigned long	*node;
} Candidate_nodes;



/*
 * Main octree structure
 *
 */
typedef struct octree {
	int				nlevels;
	unsigned long	maxnnds;
	Octree_node		*node;
	Bounding_box	bounds;
	Point3D			center;
} Octree;


extern Octree *front_octree;			/* Pointer to a main octree structure	*/

/*
 * Functions Prototypes
 */
void			OCTREE_Initialize_front_octree(Volume_mesh *);

void			OCTREE_Get_bounds(Octree *, Volume_mesh *);

void			OCTREE_Insert(unsigned long, Octree *, Volume_mesh *);

unsigned long	OCTREE_Get_leaf(Point3D *, Octree *);

void			OCTREE_Split_node(unsigned long, Octree *);

void			OCTREE_Check_memory(unsigned long, Octree *);

Boolean			OCTREE_Remove(unsigned long, Octree *, Volume_mesh *);

void			OCTREE_Range_search(Sphere *, Candidate_points *,
									Octree	*, Volume_mesh		*);

void			OCTREE_Get_candidate_nodes(Sphere		 *,	Bounding_box	*,
										   unsigned long  ,	Candidate_nodes	*,
										   Octree		 *);

void			OCTREE_Get_candidate_points(Sphere			*, Candidate_points	*,
											Candidate_nodes *, Octree			*,
											Volume_mesh		*);


Boolean			OCTREE_Box_sphere_intersect(Bounding_box *, Sphere *);

void			OCTREE_Print(Octree *);

void			OCTREE_Print_temps(Bounding_box	*, unsigned long  ,
								   unsigned long *, unsigned long *,
								   FILE *,  FILE *, Octree *);

Octant			OCTREE_Compare(Point3D *, Point3D *);

#endif /* OCTREE_H */
