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
#ifndef	ADJAC_H
#define	ADJAC_H

typedef enum adjac_consts {
	INIT_NNEIGHBORS_FACE = 8
} Adjac_consts;

typedef struct adjac_node {
	int				nnbs;
	int				maxnnbs;
	unsigned long	*neighbor;
} Adjac_node;

typedef struct adjacency {
	unsigned long	nnds;
	unsigned long	maxnnds;
	Adjac_node		*node;
} Adjacency;

extern Adjacency *fasupt;

void ADJAC_Initialize_fasupt(Surface_mesh *);

void ADJAC_Insert_neighbor(unsigned long, unsigned long, Adjacency *);

void ADJAC_Remove_neighbor(unsigned long, unsigned long, Adjacency *);


void ADJAC_Change_neighbor_id(unsigned long, unsigned long,
							  unsigned long, Adjacency *);

#endif /* ADJAC_H */
