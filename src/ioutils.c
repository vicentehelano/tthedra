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
#include <stdlib.h>
#include <errno.h>

#include <error.h>
#include <unistd.h>
#include <regex.h>

#include <main.h>

bool prompt_for_overwriting(char *file_name)
{
  char response[4];

  fprintf(stderr, _("%s: overwrite `%s'? "), program_name, file_name);
  flockfile(stdin);
  fgets(response, 4, stdin);
  funlockfile(stdin);
#if ENABLE_NLS
  int bad_regex, matched;
  static regex_t compiled_yes;

  bad_regex = regcomp(&compiled_yes, _("^[yY]"), REG_EXTENDED);
  if (bad_regex) {
    error(0, 0, _("cannot compile regular expression `%s'\n"), _("^[yY]"));
    exit(EXIT_FAILURE);
  }
  matched = regexec(&compiled_yes, response, 0, NULL, 0);
  if (!matched) return true;
    else return false;
#else
  return (*response == 'y' || *response == 'Y' ? true : false);
#endif
}

FILE *get_file(char *file_name, char *mode_str)
{
  FILE *fstream;

  if (mode_str[0] == 'r') {
    fstream = fopen(file_name, mode_str);
  } else {
    errno = 0;
    fstream = fopen(file_name, "wx");
    if (errno == EEXIST) {
      if (it_must.force) {
        if (remove(file_name) != 0) {
          error(0, errno, _("cannot remove file `%s'"), file_name);
          exit(EXIT_FAILURE);
        }
      } else {
        error(0, errno, _("cannot create regular file `%s'"), file_name);
        it_must.force = prompt_for_overwriting(file_name);
        if (it_must.force) {
          if (remove(file_name) != 0) {
            error(0, errno, _("cannot remove file `%s'"), file_name);
            exit(EXIT_FAILURE);
          }
          it_must.force = false;
        } else {
          exit(EXIT_SUCCESS);
        }
      }
      fstream = fopen(file_name, mode_str);
    }
  }

  return fstream;
}
