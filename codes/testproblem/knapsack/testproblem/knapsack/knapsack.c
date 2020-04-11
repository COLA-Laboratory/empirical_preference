/*
 * DTLZ2.c
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
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

#include "../../header/problems.h"
# include "../../header/rand.h"

int **p_nap,**w_nap,*capa,**qtemp;
void knapsack (individual_real *ind)
{
    int i, j, k;

    double *xreal, *obj;
    obj   = ind->obj;
    xreal = ind->xreal;



    for(i = 0; i < number_objective; i++)
    {
        obj[i] = 0.0;
        double temp_sum=0;
        for(j = 0; j < number_variable;j++)
        {
            obj[i] += w_nap[i][j] *  p_nap[i][j] * xreal[j];
            if(xreal[j]==1)
               temp_sum += w_nap[i][j];
        }
        if(temp_sum>capa[i])
            obj[i]=0;




    }
    for (i=2;i<number_objective;i++)
    {
        obj[i] = 0.9 * obj[i] + 0.1 *obj[0];
        i = i + 1;
    }
    for (i=3;i < number_objective; i++)
    {
        obj[i] = 0.9 * obj[i] + 0.1 *obj[1];
        i = i + 1;
    }
    for (i=0;i<number_objective;i++)
    {
        obj[i] = -1 * obj[i];

    }


}
