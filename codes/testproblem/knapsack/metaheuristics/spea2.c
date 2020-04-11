/*
 * spea2.c:
 *  This file contains the main procedures of the standard SPEA2.
 *
 * Authors:
 *  Qi Xu <qixu.student@gmail.com>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Qi Xu, Renzhi Chen, Ke Li
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

# include "../header/metaheuristics.h"

void SPEA2(population_real *parent_pop, population_real *archive, population_real *mixed_pop)
{
    // Here, the offspring is actually used as archive, so the name is changed for convenience.
    int i,j;
    int p[8][500];
    for(i = 0; i <8; i++)
    {
        for(j = 0; j <500;j++)
        {
            p[i][j]=rnd(10,100);
        }
    }
    FILE *fp1;
    if ((fp1 = fopen("w10.dat", "w")) == NULL)
    {
        printf("file cannot open \n");
        exit(0);
    }

    for(j = 0; j <500;j++)
    {
        fprintf(fp1, "%d\n",p[0][j]);
    }

    fclose(fp1);
}















