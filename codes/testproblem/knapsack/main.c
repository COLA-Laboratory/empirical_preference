/*
 * main.c:
 *  This is the main procedures of a general EMO algorithm (generational evolution model).
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

# include "header/rand.h"
# include "header/metaheuristics.h"
# include "header/reproduction.h"
#include <unistd.h>
#include <time.h>

/* common paramters */
int run_index;
int run_index_begin;
int run_index_end;
int max_evaluation;              // maximum number of evaluations (stopping criterion)
int evaluation_count;            // evaluation counter
int popsize;                     // population size
int number_variable;             // number of variables
int number_objective;            // number of objectives
double* ideal_point;             // ideal point
double* nadir_point;             // nadir point
double* variable_lowerbound;     // variable lower bound
double* variable_upperbound;     // variable upper bound
char dummy[BUFSIZE_S];
char problem_name[BUFSIZE_S];
char algorithm_name[BUFSIZE_S];
char analyse_stream[BUFSIZE_L];
char problem_param_stream[BUFSIZE_L];
/* crossover and mutation */
double eta_c;                    // eta_c in SBX
double eta_m;                    // eta_m in polynomial mutation
double pcross_real;              // crossover rate for real encoded
double pmut_real;                // mutation rate for real encoded
double CR;                       // CR in DE
double F;                        // F in DE
double K;

/* performance metrics */
int PF_size;                 // size of the true Pareto-optimal Front
double **PF_data;            // true Pareto-optimal front data
double *ref_point;           // reference point for Hypervolume calculation

/* MOEA/D variants */
int neighbor_size;                           // neighborhood length
int number_weight;                           // number of weight vectors
char *weight_file;
int function_type;                           // type of the aggregation function
int maximumNumberOfReplacedSolutions;        // the maximum replacement number of a superior offspring
double neighborhood_selection_probability;   // probability to replace in the neighborhood
double **lambda;                             // weight vectors
int **neighborhood;                          // neighborhood structure
int *permutation;                            // subproblem index permutation
int *frequency;                              // subproblem usages counter arrary
double *utility;                             // subproblem utility array
struct int_vector *selected;
struct int_vector *candidate;
double **rbf_lambda;
/* analysis platform */
int runtime_output;
int output_interval;
int analyse_list[BUFSIZE_S];
FILE *pythonplot;
pthread_t *plot_thread;
double *reference_point, *weights_obj;
int   directory_index;

int **p_nap,**w_nap,*capa,**qtemp;
double **q_nap;
int main(int argc, char *argv[])
{
    int i,j;
    reference_point = (double *) malloc (sizeof(double) * 10);
    directory_index=3;

    //reference_point[0]=0.511,reference_point[1]=0.383,reference_point[2]=0.383;
    initialization_real (argc,argv);
    p_nap = (int **) malloc (number_objective * sizeof(int *));
    for (i = 0; i < number_objective; i++)
        p_nap[i] = (int *) malloc(number_variable * sizeof(int));
    qtemp = (int **) malloc (number_objective * sizeof(int *));
    for (i = 0; i < number_objective; i++)
        qtemp[i] = (int *) malloc(number_variable * sizeof(int));
    q_nap = (double **) malloc (number_objective * sizeof(double *));
    for (i = 0; i < number_objective; i++)
        q_nap[i] = (double *) malloc(number_variable * sizeof(double));
    capa=(int *) malloc(number_objective * sizeof(int));
    w_nap = (int **) malloc (number_objective * sizeof(int *));
    for (i = 0; i < number_objective; i++)
        w_nap[i] = (int *) malloc(number_variable * sizeof(int));
    FILE *fp1;
    if ((fp1 = fopen("data/w1.dat", "r")) == NULL) {
        printf("file cannot open \n");
        exit(0);
    }

    for (i = 0; i < number_variable; i++)
        fscanf(fp1, "%d", &(w_nap[0][i]));
    fclose(fp1);


    FILE *fp2;
    if ((fp2 = fopen("data/w2.dat", "r")) == NULL) {
        printf("file cannot open \n");
        exit(0);
    }

    for (i = 0; i < number_variable; i++)
        fscanf(fp2, "%d", &(w_nap[1][i]));
    fclose(fp2);

    FILE *fp23;
    if ((fp23 = fopen("data/w3.dat", "r")) == NULL) {
        printf("file cannot open \n");
        exit(0);
    }

    for (i = 0; i < number_variable; i++)
        fscanf(fp23, "%d", &(w_nap[2][i]));
    fclose(fp23);

//    FILE *fp24;
//    if ((fp24 = fopen("w4.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp24, "%d", &(w_nap[3][i]));
//    fclose(fp24);
//
//    FILE *fp25;
//    if ((fp25 = fopen("w5.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp25, "%d", &(w_nap[4][i]));
//    fclose(fp25);
//    FILE *fp26;
//    if ((fp26 = fopen("w6.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp26, "%d", &(w_nap[5][i]));
//    fclose(fp26);
//
//    FILE *fp27;
//    if ((fp27 = fopen("w7.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp27, "%d", &(w_nap[6][i]));
//    fclose(fp27);
//    FILE *fp28;
//    if ((fp28 = fopen("w8.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp28, "%d", &(w_nap[7][i]));
//    fclose(fp28);

//===========================
    FILE *fp3;
    if ((fp3 = fopen("data/p1.dat", "r")) == NULL) {
        printf("file cannot open \n");
        exit(0);
    }

    for (i = 0; i < number_variable; i++)
        fscanf(fp3, "%d", &(p_nap[0][i]));
    fclose(fp3);
    FILE *fp4;
    if ((fp4 = fopen("data/p2.dat", "r")) == NULL) {
        printf("file cannot open \n");
        exit(0);
    }

    for (i = 0; i < number_variable; i++)
        fscanf(fp4, "%d", &(p_nap[1][i]));
    fclose(fp4);


    FILE *fp43;
    if ((fp43 = fopen("data/p3.dat", "r")) == NULL) {
        printf("file cannot open \n");
        exit(0);
    }

    for (i = 0; i < number_variable; i++)
        fscanf(fp43, "%d", &(p_nap[2][i]));
    fclose(fp43);
//    FILE *fp44;
//    if ((fp44 = fopen("data/p4.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp44, "%d", &(p_nap[3][i]));
//    fclose(fp44);
//    FILE *fp45;
//    if ((fp45 = fopen("p5.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp45, "%d", &(p_nap[4][i]));
//    fclose(fp45);
//
//    FILE *fp46;
//    if ((fp46 = fopen("p6.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp46, "%d", &(p_nap[5][i]));
//    fclose(fp46);
//
//    FILE *fp47;
//    if ((fp47 = fopen("p7.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp47, "%d", &(p_nap[6][i]));
//    fclose(fp47);
//    FILE *fp48;
//    if ((fp48 = fopen("p8.dat", "r")) == NULL) {
//        printf("file cannot open \n");
//        exit(0);
//    }
//
//    for (i = 0; i < number_variable; i++)
//        fscanf(fp48, "%d", &(p_nap[7][i]));
//    fclose(fp48);


    for (i = 0; i < number_objective; i++)
    {
        capa[i]=0;
        for (j = 0; j < number_variable; j++)
            capa[i] += w_nap[i][j];
        capa[i] =capa[i]/2;

    }


    // initialize parameter settings




    population_real *parent_pop;
    population_real *offspring_pop;
    population_real *mixed_pop;
    parent_pop    = (population_real *) malloc (sizeof(population_real));
    offspring_pop = (population_real *) malloc (sizeof(population_real));
    mixed_pop     = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (offspring_pop, popsize);
    allocate_memory_pop (mixed_pop, 2 * popsize);
    randomize ();
    int m;
    /* parameter settings for preference-based methods */

   // printf ("Please input a reference point:\n"); //reference point.
//    for (i = 0; i < number_objective; i++)
//       scanf ( "%lf", &(reference_point[i]));
//        reference_point[0]=-900000,reference_point[1]=-1000000;
    //DTLZ1-3-M
//    reference_point[0]=-945000,reference_point[1]=-1050000,reference_point[2]=-945000;
    //DTLZ2-3-M

//DTLZ2-5-M
       reference_point[0]=-900000,reference_point[1]=-1000000,reference_point[2]=-900000,reference_point[3]=-1000000,reference_point[4]=-910000,reference_point[5]=-1015000,reference_point[6]=-920000,reference_point[7]=-1030000;
    //DTLZ1-5-M
    //reference_point[0]=0.1,reference_point[1]=0.09,reference_point[2]=0.12,reference_point[3]=0.09,reference_point[4]=0.1;
    //DTLZ1-8-M
//   reference_point[0]=0.065,reference_point[1]=0.065,reference_point[2]=0.06,reference_point[3]=0.065,reference_point[4]=0.07,
//    reference_point[5]=0.06,reference_point[6]=0.055,reference_point[7]=0.055;
   //DTLZ1-8-e
//    reference_point[0]=0.250,reference_point[1]=0.035,reference_point[2]=0.04,reference_point[3]=0.035,reference_point[4]=0.035,
//   reference_point[5]=0.030,reference_point[6]=0.040,reference_point[7]=0.035;
//DTLZ2 -8-m
//    reference_point[0]=0.370,reference_point[1]=0.370,reference_point[2]=0.342,reference_point[3]=0.370,reference_point[4]=0.399,
//    reference_point[5]=0.342,reference_point[6]=0.313,reference_point[7]=0.313;
//DTLZ2 -8-e
//    reference_point[0]=0.131,reference_point[1]=0.112,reference_point[2]=0.150,reference_point[3]=0.131,reference_point[4]=0.935,
//    reference_point[5]=0.131,reference_point[6]=0.150,reference_point[7]=0.131;
//DTLZ1-10-M
//    reference_point[0]=0.055,reference_point[1]=0.06,reference_point[2]=0.055,reference_point[3]=0.04,reference_point[4]=0.05,
//   reference_point[5]=0.045,reference_point[6]=0.050,reference_point[7]=0.045,reference_point[8]=0.055,reference_point[9]=0.045;
//DTLZ1-10-e
//    reference_point[0]=0.025,reference_point[1]=0.03,reference_point[2]=0.025,reference_point[3]=0.03,reference_point[4]=0.03,
//    reference_point[5]=0.250,reference_point[6]=0.025,reference_point[7]=0.03,reference_point[8]=0.025,reference_point[9]=0.03;
//DTLZ2-10-M
//    reference_point[0]=0.345,reference_point[1]=0.377,reference_point[2]=0.345,reference_point[3]=0.251,reference_point[4]=0.314,
//    reference_point[5]=0.283,reference_point[6]=0.314,reference_point[7]=0.283,reference_point[8]=0.345,reference_point[9]=0.283;
//DTLZ2-10-e
//    reference_point[0]=0.114,reference_point[1]=0.114,reference_point[2]=0.948,reference_point[3]=0.095,reference_point[4]=0.114,
//    reference_point[5]=0.095,reference_point[6]=0.114,reference_point[7]=0.095,reference_point[8]=0.114,reference_point[9]=0.095;
//DTLZ2-5-e
//reference_point[0]=0.119,reference_point[1]=0.965,reference_point[2]=0.149,reference_point[3]=0.134,reference_point[4]=0.119;
//DTLZ5-5-M
//reference_point[0]=0.3,reference_point[1]=0.3,reference_point[2]=0.4242,reference_point[3]=0.6,reference_point[4]=0.5293;
//DTLZ5 -8-m
//    reference_point[0]=0.1145,reference_point[1]=0.1145,reference_point[2]=0.1620,reference_point[3]=0.2291,reference_point[4]=0.3240,
//    reference_point[5]=0.4582,reference_point[6]=0.6479,reference_point[7]=0.4004;
    //DTLZ5-10
//        reference_point[0]=-86000,reference_point[1]=-88000,reference_point[2]=-81000,reference_point[3]=-80000,reference_point[4]=-78000,
//    reference_point[5]=-75000,reference_point[6]=-82000,reference_point[7]=-77000,reference_point[8]=-84000,reference_point[9]=-76000;
    weights_obj = (double *) malloc (sizeof(double) * number_objective);
  // weights_obj[0]=0.2,weights_obj[1]=0.18,weights_obj[2]=0.24,weights_obj[3]=0.18,weights_obj[4]=0.2;//DTLZ 5-M
 //   weights_obj[0]=0.13,weights_obj[1]=0.13,weights_obj[2]=0.12,weights_obj[3]=0.13,weights_obj[4]=0.14,weights_obj[5]=0.12,weights_obj[6]=0.11,weights_obj[7]=0.11;//DTLZ-8-M
  //  weights_obj[0]=0.07,weights_obj[1]=0.08,weights_obj[2]=0.08,weights_obj[3]=0.07,weights_obj[4]=0.5,weights_obj[5]=0.07,weights_obj[6]=0.08,weights_obj[7]=0.07;//DTLZ2-8-e
    //weights_obj[0]=0.1,weights_obj[1]=0.09,weights_obj[2]=0.08,weights_obj[3]=0.08,weights_obj[4]=0.65;//DTLZ1-5-E
   // weights_obj[0]=0.08,weights_obj[1]=0.65,weights_obj[2]=0.1,weights_obj[3]=0.09,weights_obj[4]=0.1; //DTLZ2-5-E
   // weights_obj[0]=0.4,weights_obj[1]=0.3,weights_obj[2]=0.3;
  for (i = 0; i < number_objective; i++)
      //scanf ( "%lf", &(weights_obj[i]));
        weights_obj[i] = 1 / (double) number_objective; // weights_obj for different objectives
   //     reference_point[i]=0;

    double specificity = 0.001;     // the parameter of PBEA
    double sigma       = 0.2;       // the parameter of r-NSGA-II
    double epsilon     = 0.000001;    // the parameter of R-NSGA-II
    double radius   =   2;
    // run experiments
    for (run_index = run_index_begin; run_index <= run_index_end; run_index++)
    {
        printf ("-----------------------------\n");
        printf ("|\tThe %d run\t|\t", run_index);
        if (!strcmp (algorithm_name, "NSGA2"))
            NSGA2 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD"))
            MOEAD (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD_DRA"))
            MOEAD_DRA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD_STM"))
            MOEAD_STM (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "MOEAD_STM_DRA"))
            MOEAD_STM_DRA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "SMSEMOA"))
            SMSEMOA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "IBEA"))
            IBEA (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "HYPE"))
            HypE (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "NSGA3"))
            NSGA3 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "SPEA2"))
            SPEA2 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "rNSGA2"))
            r2NSGA2 (parent_pop, offspring_pop, mixed_pop, sigma);
        else if(!strcmp (algorithm_name, "gNSGA2"))
            gNSGA2 (parent_pop, offspring_pop, mixed_pop);
        else if(!strcmp (algorithm_name, "RNSGA2"))
            RNSGA2 (parent_pop, offspring_pop, mixed_pop, reference_point, weights_obj, epsilon);
        else if(!strcmp (algorithm_name, "PBEA"))
            PBEA (parent_pop, offspring_pop, mixed_pop, reference_point, weights_obj, specificity);
        else if(!strcmp (algorithm_name, "RMEAD2"))
            RMEAD2(parent_pop,offspring_pop,mixed_pop,reference_point,radius);
        else if(!strcmp (algorithm_name, "PLVF"))
            PLVF (parent_pop, offspring_pop, mixed_pop);
        else
            print_error (1, 2, "UNKNOWN algorithm:", algorithm_name);

        printf ("\n");
    }
    printf ("-----------------------------\n");

    // free memory
    if (number_variable != 0)
    {
        free (variable_lowerbound);
        free (variable_upperbound);
    }
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (offspring_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2 * popsize);
    free (parent_pop);
    free (offspring_pop);
    free (mixed_pop);

    for (i = 0; i < PF_size; i++)
        free (PF_data[i]);
    free (PF_data);
    free(reference_point);
    free(weights_obj);

    return 0;
}
