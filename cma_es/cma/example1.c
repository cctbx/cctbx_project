/* --------------------------------------------------------- */
/* --------------- A Very Short Example -------------------- */
/* --------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h> /* free() */
#include "cmaes_interface.h"

double fitfun(double const *x, int dim);

/* the objective (fitness) function to be minized */
double fitfun(double const *x, int N) { /* function "cigtab" */
  int i;
  double sum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];
  for(i = 2; i < N; ++i)
    sum += x[i]*x[i];
  return sum;
}

/* the optimization loop */
int main(int argn, char **args) {
  cmaes_t evo; /* an CMA-ES type struct or "object" */
  double *arFunvals, *const*pop, *xfinal;
  int i;

  /* Initialize everything into the struct evo, 0 means default */
  arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "initials.par");
  printf("%s\n", cmaes_SayHello(&evo));
  cmaes_ReadSignals(&evo, "signals.par");  /* write header and initial values */

  /* Iterate until stop criterion holds */
  while(!cmaes_TestForTermination(&evo))
    {
      /* generate lambda new search points, sample population */
      pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */

      /* Here you may resample each solution point pop[i] until it
         becomes feasible, e.g. for box constraints (variable
         boundaries). function is_feasible(...) needs to be
         user-defined.
         Assumptions: the feasible domain is convex, the optimum is
         not on (or very close to) the domain boundary, initialX is
         feasible and initialStandardDeviations are sufficiently small
         to prevent quasi-infinite looping.
      */
      /* for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i)
           while (!is_feasible(pop[i]))
             cmaes_ReSampleSingle(&evo, i);
      */

      /* evaluate the new search points using fitfun from above */
      for (i = 0; i < cmaes_Get(&evo, "lambda"); ++i) {
        arFunvals[i] = fitfun(pop[i], (int) cmaes_Get(&evo, "dim"));
      }

      /* update the search distribution used for cmaes_SampleDistribution() */
      cmaes_UpdateDistribution(&evo, arFunvals);

      /* read instructions for printing output or changing termination conditions */
      cmaes_ReadSignals(&evo, "signals.par");
      fflush(stdout); /* useful in MinGW */
    }
  printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */
  cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

  /* get best estimator for the optimum, xmean */
  xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
  cmaes_exit(&evo); /* release memory */

  /* do something with final solution and finally release memory */
  free(xfinal);

  return 0;
}
