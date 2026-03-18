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
#ifndef MAIN_H
#define MAIN_H

#include <stdbool.h>
#include <limits.h>

#include <config.h>

#include <getopt.h>
#include <gettext.h>

#ifdef ENABLE_NLS
  #define _(String) gettext(String)
  #define gettext_noop(String) String
  #define N_(String) gettext_noop(String)
#else
  #define _(String) (String)
  #define N_(String) String
#endif

/* The official name of this program */
#define PROGRAM_NAME "tthedra"
#define AUTHORS "Vicente H. F. Batista"

/* The name of the executable file */
extern char *program_name;

/* Some helpful enumerations */
enum {
  COPYRIGHT_YEAR = 2010,
  GETOPT_HELP_CHAR = (CHAR_MIN - 2),
  GETOPT_VERSION_CHAR = (CHAR_MIN - 3)
};

struct user_options {
  bool print_mesh, generate_mesh, print_report, force, be_verbose, read_mesh;
};

extern struct user_options it_must;

#endif /* MAIN_H */
