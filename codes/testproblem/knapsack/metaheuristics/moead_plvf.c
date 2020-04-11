/*
 * moead.c:
 *  This file contains the main procedures of the standard MOEAD.
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
 * MERCHANTABILITY or FITNESS FOR central_vector PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received central_vector copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

# include "../header/metaheuristics.h"

double *weights_obj,*reference_point;
//double *reference_point;
void PLVF (population_real *pop, population_real *offspring_pop, population_real *mixed_pop)
{
    int i, j ,k = 0,m;
    int generation;
    int subproblem_id, neighbor_type;
    double *rbf;
    int *rbf_index;
    int best_index=0;
    int r[2];
/* common paramters */
    int miu,landa;                                //number of incumbent candidates.
    int tau;                                // number of generations between two consecutive consultation.//  number of generations in first consecutive consultation.
    double eta;
    miu = 2 * number_objective + 1;//step size.// the width of the Gaussian function.
    int best_number=5;
    int pairs=20;
    int flag;

    generation       = 1;
    evaluation_count = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", run_index);

   //nums_weight();
    initialize_uniform_weight ();


    print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_neighborhood ();
    initialize_population_real (pop);
    evaluate_population (pop);
    initialize_idealpoint (pop);
    initialize_nadirpoint(pop);

    track_evolution (pop, generation, 0);
    double *weights_rbf;
    weights_rbf=(double *)malloc(40 * sizeof(double));




    // initialize parameter settings
    permutation = malloc (number_weight * sizeof(int));

    tau = 5;
    landa = 5;
    eta = 0.25;
    individual_real* offspring = &(offspring_pop->ind[0]);
    rbf = (double *) malloc (popsize * sizeof(double));

    double temp;

    rbf_index = (int *) malloc (popsize * sizeof(int));
    double **seed;
    seed = (double **) malloc (miu * sizeof(double *));
    for (i = 0; i < miu; i++)
        seed[i] = (double *) malloc(number_objective * sizeof(double));
    while (evaluation_count < max_evaluation)
    {
        print_progress();
        //if ((generation) % tau ==0)
        if ((generation % tau == 0))
        {



                FILE *fp1;
                if ((fp1 = fopen("plvf_train.dat", "w")) == NULL)
                {
                    printf("file cannot open \n");
                    exit(0);
                }
            for (i = 0; i < popsize; i++)
            {
                //rbf[i] = cos_angle_inverse(pop->ind[i].obj, reference_point, number_objective);
                 rbf[i] = Rbf_fitnessFunction(&(pop->ind[i]),weights_obj);
            }

            int a, b;
            a = rnd(0,popsize-1);
            for (i = 0; i < pairs; i++)
            {

                b = rnd(0, popsize - 1);
                while (b == a)
                    b = rnd(0, popsize - 1);
                if (rbf[a] < rbf[b])
                {
                    fprintf(fp1, "%d\t%d\n", a+1, b+1);
                    // a  = rnd(0,popsize-1);
                }

                if (rbf[a] > rbf[b])
                {
                    fprintf(fp1, "%d\t%d\n", b+1, a+1);
                    a = b;
                }
            }
            fclose(fp1);



            FILE *fp3;
            if ((fp3 = fopen("plvf_samples.dat", "w")) == NULL) {
                printf("file cannot open \n");
                exit(0);
            }
            for (i = 0; i < popsize; i++)
            {
                for (j = 0; j < number_objective; j++)
                    fprintf(fp3, "%lf\t", pop->ind[i].obj[j]);
                fprintf(fp3, "\n");
            }
            fclose(fp3);
            //printf("test\n");
            system("./gpref plvf_train.dat plvf_samples.dat");
            FILE *fid;
            while ((fid = fopen("test.dat", "r")) == NULL )
            {
                if ((fp1 = fopen("plvf_train.dat", "w")) == NULL)
                {
                    printf("file cannot open \n");
                    exit(0);
                }

                a = rnd(0,popsize-1);
                for (i = 0; i < pairs; i++)
                {

                    b = rnd(0,popsize-1);
                    while(b==a)
                        b = rnd(0,popsize-1);
                    if (rbf[a] < rbf[b])
                    {
                        fprintf(fp1, "%d\t%d\n", a+1, b+1);
                       // a  = rnd(0,popsize-1);
                    }

                    if (rbf[a] > rbf[b])
                    {
                        fprintf(fp1, "%d\t%d\n", b+1, a+1);
                        a = b;
                    }


                }
                fclose(fp1);
                system("./gpref plvf_train.dat plvf_samples.dat");

            }
            system("rm test.dat");
            FILE *fp;
            if ((fp = fopen("plvf_train.dat.func", "r")) == NULL) {
                printf("error\n");
                exit(1);
            }
            for (i = 0; i < popsize; i++)
                fscanf(fp, "%lf", &(rbf[i]));
            fclose(fp);
            for (i = 0; i < popsize; i++)
                rbf[i] = -1 * rbf[i];

            // Find the best incumbent candidates.
            for (i = 0; i < best_number; i++) {
                temp = INF;
                for (j = 0; j < popsize; j++) {
                    if (rbf[j] < temp) {
                        temp = rbf[j];
                        rbf_index[i] = j;
                    }
                }
                rbf[rbf_index[i]] = INF;
            }

            int best_index1[best_number];
            for (i = 0; i < best_number; i++)
                best_index1[i] =-1;
            int flag1 ;

            for (i = 0; i < best_number; i++)
            {
                best_index = rbf_index[i];
                double temple_angle=INF;
                double angle[number_weight];
                for (j = 0;j < number_weight;j++)
                {
                    flag1 = 0;
                    angle[j] = cos_angle_inverse(lambda[j],pop->ind[best_index].obj,number_objective);
                    for (m=0;m < best_number;m++)
                    {
                        if (j == best_index1[m])
                            flag1 =1;
                    }
                    if ((angle[j] < temple_angle)&&(flag1 == 0))
                    {
                        temple_angle = angle[j];
                        k = j;
                    }
                }
              best_index1[i] = k;
               // printf("k=%d\trbf_index[i]=%d\n",k,rbf_index[i]);

            }




            //Tune the other weights to the best weights.


            for (i = 0; i < best_number; i++)
            {
               // rbf_index[i] = best_index1[i];
                best_index = rbf_index[i];

                //for (j = 0; j < number_objective; j++)
                     lambda[best_index][0] += 10000;

            }


            //Tune the other weights to the best weights.

            for (i = 0; i < best_number; i++)
            {
                j = rbf_index[i];
               // for (k = 0; k < number_objective; k++)
               lambda[j][0] -= 10000;

                update_neighborhood(lambda[j], best_number, i, eta);
            }
            neighbor_size = 20;
            for (i = 0; i < number_weight; i++) {
                for (j = 0; j < number_objective; j++) {
                    lambda[i][j] -= 10000;
                }
            }

        }

        if ((generation % landa == 0)&& (generation % tau != 0)&& (generation > tau)&&(generation > 100))
        {
            FILE *fp4;
            if ((fp4 = fopen("samples.dat", "w")) == NULL)
            {
                printf("file cannot open \n");
                exit(0);
            }

            for (i = 0; i < popsize; i++) {
                for (j = 0; j < number_objective; j++)
                    fprintf(fp4, "%lf\t", pop->ind[i].obj[j]);
                fprintf(fp4, "\n");
            }
            fclose(fp4);
            system("./gpref1 plvf_train.dat plvf_samples.dat");
            FILE *fp;
            if ((fp = fopen("plvf_train.dat.func", "r")) == NULL) {
                printf("error\n");
                exit(1);
            }
            for (i = 0; i < popsize; i++)
                fscanf(fp, "%lf", &(rbf[i]));
            fclose(fp);
            for (i = 0; i < popsize; i++)
                rbf[i] = -1 * rbf[i];
            // Find the best incumbent candidates.
            for (i = 0; i < best_number; i++)
            {
                temp = INF;
                for (j = 0; j < popsize; j++)
                {
                    if (rbf[j] < temp) {
                        temp = rbf[j];
                        rbf_index[i] = j;
                    }
                }
                rbf[rbf_index[i]] = INF;
            }



            for (i = 0; i < best_number; i++)
            {
                best_index = rbf_index[i];
                lambda[best_index][0] += 10000;
            }


            //Tune the other weights to the best weights.

            for (i = 0; i < best_number; i++) {
                j = rbf_index[i];
                lambda[j][0] -= 10000;
                update_neighborhood(lambda[j], best_number, i, eta);
            }
            neighbor_size = 20;
            for (i = 0; i < number_weight; i++) {
                for (j = 0; j < number_objective; j++) {
                    lambda[i][j] -= 10000;
                }

            }


        }

        random_permutation (permutation, number_weight);
        for (i = 0; i < number_weight; i++)
        {
            subproblem_id = permutation[i];

            // crossover and mutation
            crossover_moead_real (pop, offspring, subproblem_id, &neighbor_type);
            mutation_ind (offspring);
            evaluate_individual (offspring);

            // update ideal point
            update_ideal_point (offspring);
            update_nadir_point (offspring);
            for (j = 0; j < number_objective; j++)
            {
                offspring->obj[j] = (offspring->obj[j] - ideal_point[j]) / (nadir_point[j] - ideal_point[j] + 0.01);
            }
            // update subproblem
            update_subproblem (pop, offspring, subproblem_id, neighbor_type);
        }
        generation++;

        track_evolution (pop, generation, evaluation_count >= max_evaluation);

    }


    free (permutation);
    moead_free ();
    free (rbf_index);
    free (rbf);
    free (weights_rbf);
    for (i = 0; i < miu; i++)
        free (seed[i]);
    free (seed);
    return;
}
