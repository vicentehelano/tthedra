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
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>

#include <error.h>
#include <errno.h>

#include <main.h>

/* The name of the executable file */
char *program_name;

struct user_options it_must;

static struct option const long_options[] = {
  {"file", required_argument, NULL, 'f'},
  {"print-report", no_argument, NULL, 'q'},
  {"force", no_argument, NULL, 'f'},
  {"verbose", no_argument, NULL, 'v'},
  {"help", no_argument, NULL, GETOPT_HELP_CHAR},
  {"version", no_argument, NULL, GETOPT_VERSION_CHAR},
  {NULL, 0, NULL, 0}
};

static void
usage(int status)
{
  if (status != EXIT_SUCCESS)
    fprintf(stderr, _("Try `%s --help' for more information.\n"),
        program_name);
  else {
    printf(_("Usage: %s [OPTION] FILE\n"), program_name);
    fputs(_("Construct tetrahedral meshes from triangulated surfaces.\n\n"),
        stdout);
    fputs (_("\
Mandatory arguments to long options are mandatory for short options too.\n"),
  stdout);
    fputs(_("      --help     display this help and exit\n"), stdout);
    fputs(_("      --version  output version information and exit\n"), stdout);
    printf(_("\nReport bugs to <%s>\n"), PACKAGE_BUGREPORT);
  }
  exit(status);
}

static void
version()
{
  fprintf(stdout, "%s %s\n", PROGRAM_NAME, VERSION);
  fprintf(stdout, "Copyright (C) %d %s", COPYRIGHT_YEAR, AUTHORS);
  fputs(_("\n\
This is free software. You may redistribute it and/or modify it under terms\n\
the GNU General Public License. There is NO warranty; not even for\n\
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\
See COPYING for more details.\n\n"), stdout);
  fprintf(stdout, _("Written by %s.\n"), AUTHORS);
}

/* Set global variables here */

int
main(int argc, char *argv[])
{
  bool ok;
  int  opt;
  char *input_file = NULL, *output_file = NULL;
/*  SurfaceMesh_2 *front;*/

  program_name = argv[0];

#ifdef HAVE_SETLOCALE
  setlocale (LC_ALL, "");
#endif
#ifdef ENABLE_NLS
  bindtextdomain (PACKAGE, LOCALEDIR);
  textdomain (PACKAGE);
#endif

  ok = false;
  it_must.generate_mesh = false;
  it_must.print_report = false;
  it_must.print_mesh = false;
  it_must.be_verbose = false;
  it_must.read_mesh = false;

  while ((opt = getopt_long(argc, argv, "f:qv",
                            long_options, (int *)0)) != -1) {
    switch (opt) {
      case 'f':
        it_must.generate_mesh = true;
        it_must.read_mesh = true;
        input_file = argv[optind];
        output_file = optarg;
        it_must.print_mesh = true;
        break;
      case 'q':
        it_must.print_report = true;
        break;
      case 'v':
        it_must.be_verbose = true;
        break;
      case GETOPT_HELP_CHAR:
        usage(EXIT_SUCCESS);
        break;
      case GETOPT_VERSION_CHAR:
        version();
        exit(EXIT_SUCCESS);
        break;
      default:
        usage(EXIT_FAILURE);
    }
  }

  if (it_must.read_mesh)
    ;
/*    front = read_mesh_2(input_file);*/
/*
//  if (it_must.generate_mesh) {
//    ok = generate_mesh();
//  }

//  if (it_must.print_report)
//    ok = write_quality_report(mesh[0], argv[argc - 3]);

//  if (it_must.print_mesh)
//    ok = write_mesh(output_file, mesh[0]);
*/

  if (argc <= optind) {
    error(0, 0, _("too few arguments"));
    usage(EXIT_FAILURE);
  }

  exit(ok ? EXIT_FAILURE:EXIT_SUCCESS);
}
