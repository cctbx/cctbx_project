/*                       suitename.h                           */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define MAXBINS 14  /*1--12 named bins, bin 0 triaged, bin 13 incomplete*/
#define MAXCLST 16  /*practical, <observed, limit of clusters in a bin*/
         /*clst indexed from 1: bin 33p had 10, index to 11, as of 070428 */

/*  main defines EXTERN as nothing ( "" ), so really does the declarations */
/* if EXTERN not defined, then it is defined as "extern"  so only referenced*/

#ifdef  EXTERN
#undef  EXTERN
#define EXTERN   /*no prefix, so do definition */
#else
#undef  EXTERN
#define EXTERN extern  /*extern prefix so just declaration*/
#endif


#define EOLO "\n" /*for UNIX_X11 */
#define CRLF "\n" /*for UNIX_X11 */

EXTERN int Ltest;

EXTERN FILE  *fpin,*fpout; /*070210 input stdin, output stdout*/
EXTERN char  version[256];
EXTERN int   LatEOF;
EXTERN int   Lhitend;
EXTERN int   itext;
EXTERN char  texts[256];

typedef struct suitestruct {
   char  ptID[256];
   char  basechr[2];  /*070412*/
   float chim;
   float deltam;
   float epsilon;
   float zeta;
   float alpha;
   float beta;
   float gamma;
   float delta;
   float chi;
   float ang[9]; /*8 to allow count from 1, */
                 /*9 to save chi-1, chi OR theta, eta angles 070524*/
   /*stores as names for readability, and indexed for computations*/
}suitestruct;
EXTERN struct suitestruct* suiteptr;

typedef struct residuestruct { 
   char  ptID[10][32];
   char  basechr[2];  /*070412*/
   float alpha;
   float beta;
   float gamma;
   float delta;
   float epsilon;
   float zeta;
}residuestruct;
EXTERN residuestruct* residueptr;
EXTERN residuestruct* newresidueptr;
EXTERN residuestruct* oldresidueptr;

EXTERN float matches[MAXBINS][MAXCLST];
    /*incl 1--12 named bins, # defined clsters070429*/

EXTERN float coordw[8];     /*general weightings*/
EXTERN float dominantw[8];  /*weights of dominant next to a satellite*/
EXTERN float satellitew[8]; /*weights of satellite next to a dominant*/

/*prototypes*/
/*main(int argc, char** argv)*/
int   evaluatesuite(void);
void  membership(int);
void  clearmatches(void);

