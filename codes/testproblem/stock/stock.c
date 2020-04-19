/*
 * stock.c
 *
 * Authors:
 *  Minhui Liao <minhui.liao1@gmail.com>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  COLA-Laboratory @ University of Exeter | http://cola-laboratory.github.io
 *
 * Copyright (c) 2020 Minhui Liao, Ke Li
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
#include "../../header/selection.h"

void kron (double a[],double b[],int dimension, double kron_mat[])
{
    int i, j;

    for (i = 0; i < dimension; i++)
    {
        for (j = 0; j < dimension; j++)
            kron_mat[i * dimension + j] = a[i] * b[j];
    }
}

void stock (individual_real *ind)
{
    int i, j;
    double sum;
    double turnover[number_variable], mean[number_variable], weights[number_variable], kron_value[number_variable * number_variable], cov[number_variable][number_variable];
    double temp1[number_variable];
    double co_sk[number_variable][number_variable * number_variable];

    double *D_kron_value;
    double *xreal, *obj;
    double **co_ku;

    FILE *fp1, *fp2, *fp3, *fp4, *fp5;

    D_kron_value = (double *) malloc(number_variable * number_variable * number_variable * sizeof(double));
    co_ku        = (double **) malloc (number_variable * sizeof(double *));
    for (i = 0; i < number_variable; i++)
        co_ku[i] = (double *) malloc(number_variable * number_variable * number_variable * sizeof(double));

    obj   = ind->obj;
    xreal = ind->xreal;
    for (i = 0; i < number_objective; i++)
        obj[i] = 0.0;

    sum = 0.0;
    for (i = 0; i< number_variable; i++)
        sum += xreal[i];
    for (i = 0; i < number_variable; i++)
    {
         xreal[i] = xreal[i] / sum;
         weights[i] = xreal[i];
    }

    if ((fp1 = fopen ("data/mean.dat","r")) == NULL)
    {
        printf ("error\n");
        exit (1);
    }
    for (i = 0; i < number_variable; i++)
        fscanf (fp1, "%lf", &(mean[i]));
    fclose (fp1);
    obj[0] =0.0;
    for (i = 0; i < number_variable; i++)
        obj[0] += xreal[i] * mean[i];

    if ((fp2 = fopen ("data/cov.dat","r")) == NULL)
    {
        printf ("error\n");
        exit (1);
    }
    for (i = 0; i < number_variable; i++)
    {
        for (j = 0; j < number_variable; j++)
        fscanf (fp2, "%lf", &(cov[i][j]));
    }
    fclose (fp2);

    for (i = 0; i < number_variable; i++)
    {
        temp1[i] = 0.0;
        for (j = 0; j < number_variable; j++)
           temp1[i] += xreal[j] * cov[i][j];
    }
    obj[1] = 0.0;
    for (i = 0; i < number_variable; i++)
        obj[1] += xreal[i] * temp1[i];

    kron (weights, weights, number_variable, kron_value);
    if ((fp3 = fopen ("data/skewness.dat","r")) == NULL)
    {
        printf ("error\n");
        exit (1);
    }
    for (i = 0; i < number_variable; i++)
    {
        for (j =0; j < number_variable*number_variable;j++)
            fscanf (fp3, "%lf", &(co_sk[i][j]));
    }
    for (i = 0; i < number_variable; i++)
    {
        temp1[i] = 0.0;
        for (j = 0; j < number_variable * number_variable; j++)
            temp1[i] += co_sk[i][j] * kron_value[j];
    }
    obj[2] = 0.0;
    for (i = 0; i < number_variable; i++)
    {
        obj[2] += xreal[i] * temp1[i];
    }

    if ((fp4 = fopen ("data/ku.dat","r")) == NULL)
    {
        printf ("error\n");
        exit (1);
    }
    for (i = 0; i < number_variable; i++)
    {
        for (j = 0; j < number_variable * number_variable * number_variable; j++)
            fscanf (fp4, "%lf", &(co_ku[i][j]));

    }

    for (i = 0; i < number_variable; i++)
    {
        for (j = 0; j < number_variable * number_variable; j++)
            D_kron_value[i * number_variable * number_variable + j] = kron_value[j] * weights[i];
    }

    for (i = 0; i < number_variable; i++)
        temp1[i]=0.0;

    for (i = 0; i < number_variable; i++)
    {
        for (j = 0; j < (number_variable * number_variable * number_variable); j++)
            temp1[i] += co_ku[i][j] * D_kron_value[j];
    }
    obj[3] = 0.0;
    for (i = 0; i < number_variable; i++)
        obj[3] += xreal[i] * temp1[i];

    if ((fp5 = fopen ("data/turnover.dat","r")) == NULL)
    {
        printf ("error\n");
        exit (1);
    }
    for (i = 0; i < number_variable; i++)
        fscanf (fp5, "%lf", &(turnover[i]));
    fclose (fp5);
    obj[4] =0.0;
    for (i = 0; i < number_variable; i++)
        obj[4] += xreal[i] * turnover[i];

    obj[0] = -obj[0];
    obj[2] = -obj[2];
    obj[4] = -obj[4];

    for (i = 0 ; i < number_variable ; i++)
        free (co_ku[i]);
    free (co_ku);
    free (D_kron_value);

    return;
}
