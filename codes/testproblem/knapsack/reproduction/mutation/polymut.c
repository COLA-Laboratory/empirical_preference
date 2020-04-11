/*
 * polymut.c:
 *  This file contains the functions to perform polynomial mutation operation.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Ke Li, Renzhi Chen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

# include "../../header/global.h"
# include "../../header/rand.h"
# include "../../header/reproduction.h"

/* Routine for real polynomial mutation of an individual */
void polymut_ind (individual_real *ind)
{
    int i;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    for (i = 0; i < number_variable; i++)
    {
        if (randomperc() <= pmut_real) {
            ind->xreal[i] = 1 - ind->xreal[i];

        }
    }

    return;
}
