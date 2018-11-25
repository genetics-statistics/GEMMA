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

    Internal API mainly splits logic between calling GEMMA legacy and
    faster-lmm-d code. Do not add logic - just what is needed for
    passing around data.

    Functions starting with 'internal_' are an internal complement
    that typically handle GEMMA state. These functions should not be
    used outside GEMMA and do not expose a "C" interface for general
    use.
*/

#include "int_api.h"
#include "debug.h"
#include "faster_lmm_d.h"
