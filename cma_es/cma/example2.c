#include <stdio.h>
#include <string.h> /* strncmp */
#include <math.h>
#include <stdlib.h>
#include "cmaes_interface.h"
/*___________________________________________________________________________
 *
 * Function Declarations
 *___________________________________________________________________________
*/
double **OrthogonalBasis(int DIM);
double f_rosenbrock( double const *x);
double f_rand( double const *x);
double f_constant( double const *x);
double f_kugelmin1( double const *x);
double f_sphere( double const *x);
double f_stepsphere( double const *x);
double f_cigar( double const *x);
double f_cigtab( double const *x);
double f_tablet( double const *x);
double f_elli( double const *x);
double f_ellirot( double const *x);
double f_elli100( double const *x);
double f_ellinumtest( double const *x);
double f_parabR( double const *x);
double f_sharpR( double const *x);
double f_diffpow( double const *x);
double f_diffpowrot( double const *x);
double f_gleichsys5( double const *x);

double * optimize(double(*pFun)(double const *), int number_of_restarts,
                  double increment_factor_for_population_size,
                  char *input_parameter_filename);

extern void   random_init( random_t *, long unsigned seed /*=0=clock*/);
extern void   random_exit( random_t *);
extern double random_Gauss( random_t *); /* (0,1)-normally distributed */

/*___________________________________________________________________________
//___________________________________________________________________________
//
// reads from file "initials.par" here and in cmaes_init()
//___________________________________________________________________________
*/
int main(int argn, char **args)
{
  typedef double (*pfun_t)(double const *);
  pfun_t rgpFun[99];  /* array (range) of pointer to objective function */
  char *filename = "initials.par"; /* input parameter file */
  FILE *fp = NULL;
  int nb = 0, nbrestarts = 0;
  double incpopsize = 2;
  int maxnb, ret=1;
  char c;
  double *x;

  /* Put together objective functions */
  rgpFun[0] = f_sphere;
  rgpFun[1] = f_elli;
  rgpFun[2] = f_cigar;
  rgpFun[3] = f_cigtab;
  rgpFun[4] = f_tablet;
  rgpFun[5] = f_rosenbrock;
  rgpFun[6] = f_parabR;
  rgpFun[7] = f_sharpR;
  rgpFun[8] = f_diffpow;
  rgpFun[9] = f_kugelmin1;
  rgpFun[10] = f_ellinumtest;
  rgpFun[11] = f_elli100;
  rgpFun[18] = f_gleichsys5;
  rgpFun[19] = f_rand;
  rgpFun[20] = f_constant;
  rgpFun[21] = f_stepsphere;
  rgpFun[22] = f_ellirot;
  rgpFun[23] = f_diffpowrot;
  maxnb = 23;

  /* Read objective function number and number of restarts from file */
  fp = fopen(filename, "r");
  if (fp) {
    fscanf(fp, " function number %d ", &nb);
    /* go to next line, a bit sloppy */
    for (c = ' ', ret = 1; c != '\n' && c != '\0' && c != EOF && ret && ret != EOF;
         ret=fscanf(fp, "%c", &c))
      ;
    fscanf(fp, " restarts %d %lf", &nbrestarts, &incpopsize);
    fclose(fp);
    if (nb < 0 || nb > maxnb)
      nb = 0;
    if (nbrestarts < 0)
      nbrestarts = 0;
  } else
    printf("main(): could not open %s to read function number", filename);

  /* Optimize function */

  x = optimize(rgpFun[nb], nbrestarts, incpopsize, filename);

  /* here we could utilize the solution x, and finally free memory */

  free(x);

  return 0;

} /* main() */

/*___________________________________________________________________________
//___________________________________________________________________________
//
// Somewhat extended interface for optimizing pFun with cmaes_t
// implementing a restart procedure with increasing population size
//___________________________________________________________________________
*/

double * optimize(double(*pFun)(double const *), int nrestarts, double incpopsize, char * filename)
{
  cmaes_t evo;       /* the optimizer */
  double *const*pop; /* sampled population */
  double *fitvals;   /* objective function values of sampled population */
  double fbestever=0, *xbestever=NULL; /* store best solution */
  double fmean;
  int i, irun,
    lambda = 0,      /* offspring population size, 0 invokes default */
    countevals = 0;  /* used to set for restarts */
  char const * stop; /* stop message */

  for (irun = 0; irun < nrestarts+1; ++irun) /* restarts */
    {
      /* Parameters can be set in three ways. Here as input parameter
       * to cmaes_init, as value read from initials.par in readpara_init
       * during initialization, and as value read from signals.par by
       * calling cmaes_ReadSignals explicitely.
       */
      fitvals = cmaes_init(&evo, 0, NULL, NULL, 0, lambda, filename); /* allocs fitvals */
      printf("%s\n", cmaes_SayHello(&evo));
      evo.countevals = countevals; /* a hack, effects the output and termination */
      cmaes_ReadSignals(&evo, "signals.par"); /* write initial values, headers in case */

      while(!(stop=cmaes_TestForTermination(&evo)))
        {
          /* Generate population of new candidate solutions */
          pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

          /* Here optionally handle constraints etc. on pop. You may
           * call cmaes_ReSampleSingle(&evo, i) to resample the i-th
           * vector pop[i], see below.  Do not change pop in any other
           * way. You may also copy and modify (repair) pop[i] only
           * for the evaluation of the fitness function and consider
           * adding a penalty depending on the size of the
           * modification.
           */

          /* Compute fitness value for each candidate solution */
          for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i) {
            /* You may resample the solution i until it lies within the
               feasible domain here, e.g. until it satisfies given
               box constraints (variable boundaries). The function
               is_feasible() needs to be user-defined.
               Assumptions: the feasible domain is convex, the optimum
               is not on (or very close to) the domain boundary,
               initialX is feasible (or in case typicalX +- 2*initialStandardDeviations
               is feasible) and initialStandardDeviations is (are)
               sufficiently small to prevent quasi-infinite looping.
            */
            /* while (!is_feasible(pop[i]))
                 cmaes_ReSampleSingle(&evo, i);
            */
            fitvals[i] = (*pFun)(pop[i]);
          }

          /* update search distribution */
          cmaes_UpdateDistribution(&evo, fitvals);

          /* read control signals for output and termination */
          cmaes_ReadSignals(&evo, "signals.par"); /* from file signals.par */

          fflush(stdout);
        } /* while !cmaes_TestForTermination(&evo) */

      lambda = incpopsize * cmaes_Get(&evo, "lambda");   /* needed for the restart */
      countevals = cmaes_Get(&evo, "eval");     /* ditto */

      /* print some "final" output */
      printf("%.0f generations, %.0f fevals (%.1f sec): f(x)=%g\n",
             cmaes_Get(&evo, "gen"), cmaes_Get(&evo, "eval"),
             evo.eigenTimings.totaltime,
             cmaes_Get(&evo, "funval"));
      printf("  (axis-ratio=%.2e, max/min-stddev=%.2e/%.2e)\n",
             cmaes_Get(&evo, "maxaxislen") / cmaes_Get(&evo, "minaxislen"),
             cmaes_Get(&evo, "maxstddev"), cmaes_Get(&evo, "minstddev")
             );
      printf("Stop (run %d):\n%s\n",  irun+1, cmaes_TestForTermination(&evo));

      /* write some data */
      cmaes_WriteToFile(&evo, "all", "allcmaes.dat");

      /* keep best ever solution */
      if (irun == 0 || cmaes_Get(&evo, "fbestever") < fbestever) {
        fbestever = cmaes_Get(&evo, "fbestever");
        xbestever = cmaes_GetInto(&evo, "xbestever", xbestever); /* alloc mem if needed */
      }
      /* best estimator for the optimum is xmean, therefore check */
      if ((fmean = (*pFun)(cmaes_GetPtr(&evo, "xmean"))) < fbestever) {
        fbestever = fmean;
        xbestever = cmaes_GetInto(&evo, "xmean", xbestever);
      }

      cmaes_exit(&evo); /* does not effect the content of stop string and xbestever */

      /* abandon restarts if target fitness value was achieved or MaxFunEvals reached */
      if (stop) /* as it can be NULL */ {
        if (strncmp(stop, "Fitness", 7) == 0 || strncmp(stop, "MaxFunEvals", 11) == 0)
          break;
      }
      if (strncmp(stop, "Manual", 6) == 0) {
        printf("Press RETURN to start next run\n"); fflush(stdout);
        getchar();
      }
    } /* for restarts */

  return xbestever; /* was dynamically allocated, should be freed in the end */
}

#if 1
/*___________________________________________________________________________
//___________________________________________________________________________
*/
double f_rand( double const *x)
{
  double d = (double)rand() / RAND_MAX;
  while (d == 0.)
    d = (double)rand() / RAND_MAX;
  return d;
}
double f_constant( double const *x)
{
  return 1;
}
#endif

static double SQR(double d)
{
  return (d*d);
}

/* ----------------------------------------------------------------------- */
double f_stepsphere( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);
  for (i = 0; i < DIM; ++i)
    sum += floor(x[i]*x[i]);
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_sphere( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);
  for (i = 0; i < DIM; ++i)
    sum += x[i]*x[i];
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_cigar( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);

  for (i = 1; i < DIM; ++i)
    sum += x[i]*x[i];
  sum *= 1e6;
  sum += x[0]*x[0];
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_cigtab( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);

  sum = x[0]*x[0] + 1e8*x[DIM-1]*x[DIM-1];
  for (i = 1; i < DIM-1; ++i)
    sum += 1e4*x[i]*x[i];
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_tablet( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);

  sum = 1e6*x[0]*x[0];
  for (i = 1; i < DIM; ++i)
    sum += x[i]*x[i];
  return sum;
}

/* ----------------------------------------------------------------------- */
/* a hack, memory is never released */
double **OrthogonalBasis(int DIM) {
  static int b_dim;
  static double **b;
  double sp;
  int i,j,k;
  random_t R;

  if(b_dim != 0) { /* Initialization was done */

    if (b_dim != DIM) {
      printf("function OrthogonalBasis cannot change dimensionality in file example2.c");
      exit(0);
    }

    return b;
  }

  /* Otherwise initialize basis b */
  random_init(&R, 2); /* TODO: choose not always the same basis? */

  /* allocate b */
  b = (double **) calloc((unsigned) DIM, sizeof(double*));
  if (!b) {
    printf("calloc failed in function OrthogonalBasis in file example2.c");
    exit(0);
  }
  for (i = 0; i < DIM; ++i) {
    b[i] = (double *) calloc((unsigned) DIM, sizeof(double));
    if (!b[i]) {
      printf("calloc failed in function Orthogonalbasis in file example2.c");
      exit(0);
    }
  }
  b_dim = DIM;

  /* generate orthogonal basis */
  for (i = 0; i < DIM; ++i) {
    /* sample components gaussian */
    for (j = 0; j < DIM; ++j)
      b[i][j] = random_Gauss(&R);
    /* substract projection of previous vectors */
    for (j = i-1; j >= 0; --j) {
      for (sp = 0., k = 0; k < DIM; ++k)
        sp += b[i][k]*b[j][k]; /* scalar product */
      for (k = 0; k < DIM; ++k)
        b[i][k] -= sp * b[j][k]; /* substract */
    }
    /* normalize */
    for (sp = 0., k = 0; k < DIM; ++k)
      sp += b[i][k]*b[i][k]; /* squared norm */
    for (k = 0; k < DIM; ++k)
      b[i][k] /= sqrt(sp);
  }
  random_exit(&R);

  return b;

} /* OrthogonalBasis(int DIM) */

/* ----------------------------------------------------------------------- */
double f_ellirot( double const *x)
{
  int i, k;
  double sum = 0., y;
  int DIM = (int)(x[-1]);
  double **B = OrthogonalBasis(DIM);

  if (DIM == 1)
    return x[0] * x[0];
  for (i = 0; i < DIM; ++i) {
    for (y = 0., k = 0; k < DIM; ++k)
      y += B[i][k] * x[k];
    sum += exp(log(1e6) * 2. * (double)(i)/(DIM-1)) * y*y;
  }
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_elli( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);

  if (DIM == 1)
    return x[0] * x[0];
  for (i = 0; i < DIM; ++i)
    sum += exp(log(1000.) * 2. * (double)(i)/(DIM-1)) * x[i]*x[i];
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_elli100( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);

  if (DIM == 1)
    return x[0] * x[0];
  for (i = 0; i < DIM; ++i)
    sum += exp(log(100.) * 2. * (double)(i)/(DIM-1)) * x[i]*x[i];
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_diffpow( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);

  if (DIM == 1)
    return x[0] * x[0];
  for (i = 0; i < DIM; ++i)
    sum += pow(fabs(x[i]), 2.+10*(double)(i)/(DIM-1));
  return sum;
}
/* ----------------------------------------------------------------------- */
double f_diffpowrot( double const *x)
{
  int i, k;
  double sum = 0., y;
  int DIM = (int)(x[-1]);
  double **B = OrthogonalBasis(DIM);

  if (DIM == 1)
    return x[0] * x[0];
  for (i = 0; i < DIM; ++i) {
    for (y = 0., k = 0; k < DIM; ++k)
      y += B[i][k] * x[k];
    sum += pow(fabs(y), 2.+10*(double)(i)/(DIM-1));
  }
  return sum;
}
/* ----------------------------------------------------------------------- */
double f_kugelmin1( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);

  for (i = 1; i < DIM; ++i)
    sum += x[i]*x[i];
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_rosenbrock( double const *x)
/*
        Rosenbrock's Function, generalized.
*/
{
  double qualitaet;
  int i;
  int DIM = (int)(x[-1]);
        qualitaet = 0.0;

        for( i = DIM-2; i >= 0; --i)
          qualitaet += 100.*SQR(SQR(x[i])-x[i+1]) + SQR(1.-x[i]);
        return ( qualitaet);
} /* f_rosenbrock() */

/* ----------------------------------------------------------------------- */
double f_parabR( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);
  for (i = 1; i < DIM; ++i)
    sum += x[i]*x[i];
  return -x[0] + 100.*sum;
}

/* ----------------------------------------------------------------------- */
double f_sharpR( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);
  for (i = 1; i < DIM; ++i)
    sum += x[i]*x[i];
  return -x[0] + 100*sqrt(sum);
}

/* ----------------------------------------------------------------------- */
double f_ellinumtest( double const *x)
{
  int i;
  double sum = 0.;
  int DIM = (int)(x[-1]);
  static double maxVerhaeltnis = 0.;
  if (maxVerhaeltnis == 0.)
    {
      for (maxVerhaeltnis = 1.;
           maxVerhaeltnis < 1e99 && maxVerhaeltnis < 2. * maxVerhaeltnis;
           maxVerhaeltnis *= 2.)
        if (maxVerhaeltnis == maxVerhaeltnis + 1.)
          break;
      maxVerhaeltnis *= 10.;
      maxVerhaeltnis = sqrt (maxVerhaeltnis);
    }
  if (DIM < 3)
    return x[0] * x[0];
  for (i = 1; i < DIM; ++i)
    sum += exp(log(maxVerhaeltnis) * 2. * (double)(i-1)/(DIM-2)) * x[i]*x[i];
  return sum;
}

/* ----------------------------------------------------------------------- */
double f_gleichsys5( double const *x)
/*
        Gleichungssystem 5-dimensional von Juergen Bremer
        Fuer jede Zeile soll gelten:
         c_1*x[1] + c_2*x[2] + c_3*x[3] + c_4*x[4] + c_5*x[5] + c_0 = 0
         Deswegen geht das Quadrat der linken Seite in die
         Qualitaetsfunktion ein.
*/
{
  double qualitaet = 0.0;

#if 1
  static double koeff[5][6] =
    {/* c_1,   c_2,  c_3,   c_4,  c_5,   c_0 */
      {   4,   191,   27,   199,   21,   172},
      { 191, 10883, 1413,  5402,  684, -8622},
      {  27,  1413,  191,  1032,  118,   -94},
      { 199,  5402, 1032, 29203, 2331, 78172},
      {  21,   684,  118,  2331,  199,  5648}
    };
  int i, j;
  double sum;

  for( i = 0; i < 5; ++i)
    {
      sum = koeff[i][5];
      for ( j = 0; j < 5; ++j)
        {
          sum += koeff[i][j] * x[j];
        }
      qualitaet += sum * sum;
    }
#endif
  return qualitaet;
} /* f_gleichsys5() */


/*
  05/10/05: revised buggy comment on handling constraints by resampling
*/
