/*
    Genome-wide Efficient Mixed Model Association (GEMMA)
    Copyright © 2011-2017, Xiang Zhou
    Copyright © 2017, Peter Carbonetto
    Copyright © 2017, Pjotr Prins

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

#ifndef __GEMMA_H__
#define __GEMMA_H__

#include "param.h"

using namespace std;

// OPTIONS
// -------
// gk:      21-22
// gs:      25-26
// gq:      27-28
// eigen:   31-32
// lmm:     1-5
// bslmm:   11-15
// predict: 41-43
// lm:      51
// vc:      61
// ci:      66-67
// calccor: 71
// gw:      72

enum M_MODE { M_LMM1=1, M_LMM2=2, M_LMM3=3, M_LMM4=4, M_LMM5=5,
  M_LMM9=9,  // GeneNetwork mode
  M_BSLMM5=15,
  M_KIN=21, M_KIN2=22, M_EIGEN=31
};

class GEMMA {

public:
  // Parameters.
  string version;
  string date;
  string year;

  // Constructor.
  GEMMA(void);

  // Functions.
  void PrintHeader(void);
  void PrintHelp(size_t option);
  void PrintLicense(void);
  void Assign(int argc, char **argv, PARAM &cPar);
  void BatchRun(PARAM &cPar);
  void WriteLog(int argc, char **argv, PARAM &cPar);
};

#endif
