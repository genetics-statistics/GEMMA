/*
    Genome-wide Efficient Mixed Model Association API (GEMMA API)
    Copyright © 2018, Pjotr Prins

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

    MODULE INFO

    This file includes the api functions that can be called from R,
    Python etc. Functions starting with 'internal_' are an internal
    complement that typically handle GEMMA state. These functions
    should not be used outside GEMMA and do not expose a "C" interface
    for general use.
*/

#include "int_api.h"
#include "debug.h"
#include "faster_lmm_d.h"

char *api_faster_lmm_d_version(char *buf) {
  #ifdef FASTER_LMM_D
  if (use_faster_lmm_d())
    return flmmd_api_version(buf);
  #endif
  return NULL;
}

void api_compute_and_write_K(const char* target, const char* file_geno, int is_centered) {
  flmmd_compute_and_write_K(target, file_geno, is_centered);
}

void api_write_K(string filen, const gsl_matrix *K) {
  // if (use_fast_lmm_d())
  //   flmmd_write_K(inds,G,is_centered);
  // else
    fail_msg("Unsupported function without faster-lmm-d");
}

// Handles internal state
void int_api_write_K(const PARAM *cPar, const gsl_matrix *K, bool is_centered) {
  // if (use_faster_lmm_d())
  //   api_write_K(filename, cPar.setKSnps(), K, is_centered);
  // else {
    cPar->WriteMatrix(K, (is_centered ? "cXX" : "sXX"));
  // }
}
