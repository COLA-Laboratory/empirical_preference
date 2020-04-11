/*
 * sbx.c:
 *  This file contains the functions to perform Simulated Binary Crossover (SBX) operation.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Renzhi Chen, Ke Li
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
# include <time.h>
# include "../../header/reproduction.h"
double *weights_obj;

/* Routine for real variable SBX crossover */
void sbx_crossover (individual_real *parent1, individual_real *parent2, individual_real *child1,individual_real *child2)
{
    int i,duandian;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;

    if (randomperc() <= pcross_real) {

        duandian = rnd(0,number_variable-1);

        for(i=0;i < duandian;i++)
        {
            child1->xreal[i] = parent1->xreal[i];
            child2->xreal[i] = parent2->xreal[i];
        }
        for(i=duandian;i < number_variable;i++)
        {
            child1->xreal[i] = parent2->xreal[i];
            child2->xreal[i] = parent1->xreal[i];
        }


    }
}

void sbx_mcrossover (individual_real *parent1, individual_real *parent2, individual_real *child1,int sub_problem_id)
{
    int i;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
    individual_real *child2;
    int duandian;
    child2 = child1;
    if (randomperc() <= pcross_real) {

        duandian = rnd(0,number_variable-1);

        for(i=0;i < duandian;i++)
        {
            child1->xreal[i] = parent1->xreal[i];
            child2->xreal[i] = parent2->xreal[i];
        }
        for(i=duandian;i < number_variable;i++)
        {
            child1->xreal[i] = parent2->xreal[i];
            child2->xreal[i] = parent1->xreal[i];
        }


    }

    return;
}