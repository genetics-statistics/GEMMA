/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright (C) 2011-2017, Xiang Zhou

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "gemma.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

int main(int argc, char *argv[]) {
  GEMMA cGemma;
  PARAM cPar;

  gsl_set_error_handler (&gemma_gsl_error_handler);

  if (argc <= 1) {
    cGemma.PrintHeader();
    cGemma.PrintHelp(0);
    return EXIT_SUCCESS;
  }
  if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
    cGemma.PrintHeader();
    cGemma.PrintHelp(0);
    return EXIT_SUCCESS;
  }
  if (argc == 3 && argv[1][0] == '-' && argv[1][1] == 'h') {
    string str;
    str.assign(argv[2]);
    cGemma.PrintHeader();
    cGemma.PrintHelp(atoi(str.c_str()));
    return EXIT_SUCCESS;
  }
  if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'l') {
    cGemma.PrintHeader();
    cGemma.PrintLicense();
    return EXIT_SUCCESS;
  }

  cGemma.Assign(argc, argv, cPar);

  ifstream check_dir((cPar.path_out).c_str());
  if (!check_dir) {
    #ifdef WINDOWS
      mkdir((cPar.path_out).c_str());
    #else
      mkdir((cPar.path_out).c_str(), S_IRWXU | S_IRGRP | S_IROTH);
    #endif
  }

  if (!is_quiet_mode())
    cGemma.PrintHeader();

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  if (is_quiet_mode()) {
    stringstream ss;
    cout.rdbuf(ss.rdbuf());
  }

  cPar.CheckParam();

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  cGemma.BatchRun(cPar);

  if (cPar.error == true) {
    return EXIT_FAILURE;
  }

  cGemma.WriteLog(argc, argv, cPar);

  return EXIT_SUCCESS;
}
