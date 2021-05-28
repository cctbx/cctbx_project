/*                    suitenutil.c                         */
#include "suitename.h"
#define SUITENUTIL
#include "suitenutil.h"
#undef  SUITENUTIL

#include "suiteninit.h"

/****vector7ab()**************************************************************/
void vector7ab(float* atob, float* a, float* b)
{
   int   j=0,ncnt=1,ndim=7; /*could pass these values ... */

   if(ncnt==1) {atob[0] = 0;} /*count from 1*/
   for(j=ncnt; j<=ndim; j++)
   {
      atob[j] = b[j] - a[j];
   }
}
/*___vector7ab()_____________________________________________________________*/

/****hyperellipsoiddist()*****************************************************/
float hyperellipsoiddist(int i, int j, int nang, float* warray)
{
   /*Piet Hein superellipse, hyperellipsoids */
   /*library power function: */
   /*http://www.cs.utah.edu/dept/old/texinfo/glibc-manual-0.02/library_17.html*/
   /*double pow (double base, double power)*/

   int k=0,kmin=0,kmax=0,n=0;
   double del=0,delpower=0,dpower=0;

   /*power is global, can be specified on input: -power #.# */

   if(nang == 4){kmin=2; kmax=5;}
   else {kmin=1; kmax=7;}

   dpower = 0;
   for(k=kmin; k<=kmax; k++) /*evaluate distance to this i,j cluster*/
   {
      /*    1      2      3      4      5      6      7   */
      /* deltam epsilon zeta   alpha  beta   gamma  delta */
      /*    X                                  X      X   */
      /*    X not used in 4 angle distance calc*/
      /*del =  suiteptr->ang[k] - clusterav[i][j][k];*/ /*globals*/
      del =  suiteptr->ang[k] - bin[i].clst[j].ang[k]; /*globals*/
      if(del < 0){del = -del;}
      del = del/warray[k]; /*normalize, del < 1 inside the ellipsoid*/
      delpower = 1;
      /*for(n=1; n<=power; n++)*/
         /*{delpower = del*delpower;}*/ /*this component*/
      delpower = pow(del,power);
      dpower = dpower + delpower;   /*accummulate components*/
   }
   dpower = pow(dpower,1/power);
/*
fprintf(stderr,"suite: %s nangle==%d  power== %f, distance== %f\n",suiteptr->ptID,nang,power,dpower);
*/
   return((float)dpower);
}
/*___hyperellipsoiddist()____________________________________________________*/

/****dotproduct()*************************************************************/
float dotproduct(float* first, float* secnd, int nang)
{
   int k=0,kmin=0,kmax=0;
   float answer=0;

   //if(nang = 4){kmin=2; kmax=5;}
   if(nang == 4){kmin=2; kmax=5;}
   else {kmin=1; kmax=7;}

   answer = 0;
   for(k=kmin; k<=kmax; k++) /*evaluate distance to this i,j cluster*/
   {
      /*    1      2      3      4      5      6      7   */
      /* deltam epsilon zeta   alpha  beta   gamma  delta */
      /*    X                                  X      X   */
      /*    X not used in 4 angle calc*/

      answer = answer + (first[k] * secnd[k]);
   }
   return(answer);
}
/*___dotproduct()____________________________________________________________*/

/****resetcoordw()************************************************************/
void resetcoordw(float* widptr, int ith, int jth, int Lspecial)
{
         widptr[1] = deltamw  ;
         widptr[6] = gammaw   ;
         widptr[7] = deltaw   ;
   if(Lspecial && satellites[ith][jth] > 0)
   {
         widptr[2] = epsilonsatw ;
         widptr[3] = zetasatw    ;
         widptr[4] = alphasatw   ;
         widptr[5] = betasatw    ;
   }
   else
   {
         widptr[2] = epsilonw ;
         widptr[3] = zetaw    ;
         widptr[4] = alphaw   ;
         widptr[5] = betaw    ;
   }
}
/*___resetcoordw()___________________________________________________________*/

/****confirmsuite()***********************************************************/
int  confirmsuite(void)
{
   if(   suiteptr->deltam  >= 0 && suiteptr->deltam  <= 360
      && suiteptr->epsilon >= 0 && suiteptr->epsilon <= 360
      && suiteptr->zeta    >= 0 && suiteptr->zeta    <= 360
      && suiteptr->alpha   >= 0 && suiteptr->alpha   <= 360
      && suiteptr->beta    >= 0 && suiteptr->beta    <= 360
      && suiteptr->gamma   >= 0 && suiteptr->gamma   <= 360
      && suiteptr->delta   >= 0 && suiteptr->delta   <= 360
     )
        {return(1);}
   else {return(0);}
}
/*___confirmsuite()__________________________________________________________*/

/****CompArgStr()***********jmw/utility***************************************/
int CompArgStr(char *str, char *arg, int min)
{
   /*requires exact (caseless) match up through min number of characters */
   /* past that, must break or continue exactly matching to be valid*/

   int i, max;
   char s, a;

   if (!str || !arg) return 0;

   max = strlen(arg);

   for(i=0; i<max; i++)
   {
      s = toupper(str[i]);
      a = toupper(arg[i]);

      if (i >= min && (   s == '\0' || s == ',' || s == '(' || s == '.'
                       || s == '+'  || s == '-' || isdigit(s)))
      {
          break; /* good ending point */
      }
      else if (s != a)
      {
          i = 0; /* failed to match */
          break;
      }
   }
   return i;
}
/*___CompArgStr()____________________________________________________________*/

