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
*/

#ifndef __GEMMA_API_H
#define __GEMMA_API_H

// Exposed C API

extern "C" {
  char *api_faster_lmm_d_version(char *buf);
  void api_compute_and_write_K(const char* target, const char* file_geno, const char * file_anno, int is_loco, int is_centered, double maf_level);
}


#endif
