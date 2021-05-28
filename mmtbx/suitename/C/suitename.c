/*                           suitename.c                         */
/*****************************************************************/
/* NOTICE: This is free software and the source code is freely   */
/* available. You are free to redistribute or modify under the   */
/* conditions that (1) this notice is not removed or modified    */
/* in any way and (2) any modified versions of the program are   */
/* also available for free.                                      */
/*               ** Absolutely no Warranty **                    */
/* Copyright (C) 2007 David C. Richardson                        */
/*****************************************************************/

#define EXTERN
    /* that is, EXTERN is defined, so suitename.h will do main declarations*/
#include "suitename.h"
#undef  EXTERN

#include "suitenscrt.h"
#include "suitenutil.h"
#include "suiteninit.h"
#include "suiteninpt.h"
#include "suitenout.h"

/*0.2.070524 preserve chi-1 and chi, so could preserve eta, theta */
/*0.3.070525 general read dangle record for, e.g.,  eta, theta */
/*0.3.070628 triage reports zeta-1, epsilon-1, delta-1,... Ltriage codes */
/*0.3.070803 notes: rearranged suitenhead.h/janesviews ...  */
/*                  put something in to say what veiws mean */
/*  what are masters e and d ???? */
/*0.3.070919 3g wannabe (tRNA TpseudoUC loop) */
/*0.3.110606 range of delta updated by S.J. */
/* 01/07/2014 S.J. updated so that it can take input with alternate conformations, *nd by default will calculate the suite for altA*/
/* 09/18/2014 S.J. updated so that suitename will ignore DNA residues*/
/****main()*******************************************************************/
int main(int argc, char** argv)
{
   int  LOK=1,ibin=0,jclst=0;
   char sour[32];
   float suiteness=0, distance=0;
   char ptmaster[4];
   char ptcolor[16];

   int i=0,j=0,k=0,m=0;

   ptcolor[0]  = '\0';  /*default is no point color*/
   ptmaster[0] = '\0';  /*default is no point master*/


   sprintf(version,"suitename.0.4.130509 ");  /*  VERSION  */

   /* default values before parsecommandline */
   Ltest = 0; /*compile in... */
   Ltestout = 0; /*commandline -test */
   fpin = stdin;
   fpout = stdout;

     /*NOTE: for consistency, edit defaults text in usageout() */
   Lreportout = 1; /* suite by suite suiteness report, & summary */
   Lchart = 0;     /* summary-less report for MolProbity multichart 070521*/
   Ldangle = 0;  /* read straight dangle records  070525*/
   Lsourout = 0; /* extra info in kinemage ptIDs. optional as of 070524*/
   Letatheta = 0; /* theta,eta instead of chi-1,chi kinemage angles 070524*/
   Lkinemageout=0; /* kinemage of the clusters */
   Lstringout = 0; /* 3 char string instead of cluster kinemage */
   Lsuitesin=0;
   Lresiduesin=1;
   NptIDfields=6;  /*for dangle residue input*/
   altIDfield=0; /* by default, no altID field - S.J. 01/07/2014*/ 
   altID[0]='\0'; altID[1]='\0';
   Nanglefields=9; /*for 9D kinemage edited suite input*/
   Lnewfile=0;
   Lhelpout=0;
   LNptIDfields=0; /* initializing the four variables to check correct usage of altIDfield - S.J. 01/07/2014*/
   LaltIDfield=0;
   Luseincorrect=0; 
   LaltID=0;
   Lchangesout=0;
   NameStr[0] = '\0';
   Lgeneralsatw = 0; /*flag for special general case satellite widths 070328*/
   Lwannabe = 1;/* -wannabe input flag 070429, default 070525, else -nowannabe*/

   Lsequence = 1; /*output 1 letter Base code sequence part of string 070409*/
   Loverlap = 0;  /*overlap string lines, e.g. 1-20, 11-30, 21-40, ... 070409*/
   Loneline = 0;  /*string output all as oneline 070409*/

   if(argc > 1)
   {
     LOK = parsecommandline(&argc, argv);
   }

   if(LOK)
   {
      /*command OK*/
     inittextblock(&mainscratch); /*__open scratch "tapes" */
     /*--rewind scratch "tapes", inplicit with inittextblock() */
     /*--allocation of space automatic when attempting write*/

     newresidueptr = (residuestruct*)malloc(sizeof(struct residuestruct));
     oldresidueptr = (residuestruct*)malloc(sizeof(struct residuestruct));
     suiteptr      = (suitestruct*)malloc(sizeof(struct suitestruct));
     if(   newresidueptr != NULL
         && oldresidueptr != NULL
         && suiteptr != NULL
        )  {LOK=1;}
     else {LOK=0;}

     if(LOK)
     {/*storage allocated*/
       LOK = initializations(); /*070325*/
       if(LOK)
       {/*LOK to work*/

         clearnewresidue(); /*only at beginning so first-1 handled correctly*/
/*sudo pre-Algorithm: accummulate suites or suites from residues */
         while(!LatEOF) /*allow EOF on end of last line, see Lhitend */
         {/*loop over all residues in the file*/
/*
if(Lresiduesin)
{
fprintf(stderr,"LOOP:newresidue:");
for(j=0;j<10;j++) fprintf(stderr,"%s,",newresidueptr->ptID[j]);
fprintf(stderr,";%s;",newresidueptr->basechr);
fprintf(stderr,"%8.3f ",newresidueptr->alpha);
fprintf(stderr,"%8.3f ",newresidueptr->beta);
fprintf(stderr,"%8.3f ",newresidueptr->gamma);
fprintf(stderr,"%8.3f ",newresidueptr->delta);
fprintf(stderr,"%8.3f ",newresidueptr->epsilon);
fprintf(stderr,"%8.3f ",newresidueptr->zeta);
fprintf(stderr,"\n");
}
*/
            if(Lsuitesin) {LOK = getsuite();}
            else {LOK = getresidue(); Lresiduesin = 1;}
            if(LOK)
            {
               if(Lresiduesin) {LOK = loadsuite();}
               if(LOK)
               {
                  LOK = confirmsuite();
                  if(!LOK)
                  {
                     ibin = 13; /*angles not fully specified*/
                     jclst = 0; /*cluster not assigned in dummy bin*/
                     sprintf(sour," tangled "); /*suite incomplete angles*/
                     writesuite(ibin,jclst,sour,distance,suiteness,ptmaster,ptcolor);
                       /*binname "inc "  clustername "__" 070414 */
                  }
                  else     /*now have to do some work*/
                  {
/*sudo Algorithm  */
/*sudo FOR all RNA suites si in S */
  /*sudo program body: evaluate suite, cluster membership */

                     ibin = evaluatesuite();

                     if(ibin > 0){membership(ibin);}

                  }
               }/*loaded suite*/
            }/*got a residue (or a suite from a kinemage)*/
         }/*loop over all residues in the file*/

         writeoutput();  /*has logicals for what output possible*/

       }/*LOK to work*/
     }/*storage allocated*/
   }/*command OK*/
   else
   {
      if(Lkinemageout) /*dump of actual cluster information   070421*/
      {
         initializations();
         L33out = L32out = L23out = L22out = 1;
         for(i=1; i<=12; i++)
         {
            binout[i]=1;
            /*for(j=1; j<ddgmax; j++)*/
            j=0; /*zeroth cluster of a bin holds its outliers*/
            while(bin[i].clst[j].LOK || j==0)  /*070429*/
            {
               clusterout[i][j] = 1;
               j++;
            }
         }
         writeoutput();
      }
      else if(Ltestout) /*test of alternate storage of cluster info  070421*/
      {  /*brute force show what is in the suiteninit.h defined arrays */
         for(i=0; i<14; i++)
         {
            for(j=0; j<MAXCLST; j++) /*stupidly loop over whole array*/
            {
               if(bin[i].clst[j].LOK) /*however, only take defined clusters*/
               {
                  fprintf(stderr,"%s %s :%d: %s %s %s: "
                  "%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n"
                  ,bin[i].binname
                  ,bin[i].clst[j].clustername
                  ,bin[i].clst[j].LOK
                  ,bin[i].clst[j].status
                  ,bin[i].clst[j].clustercolor
                  ,bin[i].clst[j].domsatness
                  ,bin[i].clst[j].ang[1]
                  ,bin[i].clst[j].ang[2]
                  ,bin[i].clst[j].ang[3]
                  ,bin[i].clst[j].ang[4]
                  ,bin[i].clst[j].ang[5]
                  ,bin[i].clst[j].ang[6]
                  ,bin[i].clst[j].ang[7]
                 );
               }
            }
         }
     /*format spacings for output layout columns: */
     fprintf(stderr,"triage axes limits version: %s\n",axeslimitsversion);
     fprintf(stderr,"epsilonmin %.0f, epsilonmax %.0f\n",epsilonmin,epsilonmax);
     fprintf(stderr,"delta3min   %.0f, delta3max  %.0f\n",delta3min,delta3max);
     fprintf(stderr,"delta2min  %.0f, delta2max  %.0f\n",delta2min,delta2max);
     fprintf(stderr,"gammapmin   %.0f, gammapmax   %.0f\n",gammapmin,gammapmax);
     fprintf(stderr,"gammatmin  %.0f, gammatmax  %.0f\n",gammatmin,gammatmax);
     fprintf(stderr,"gammammin  %.0f, gammammax  %.0f\n",gammammin,gammammax);
     fprintf(stderr,"alphamin    %.0f, alphamax   %.0f\n",alphamin,alphamax);
     fprintf(stderr,"betamin     %.0f, betamax    %.0f\n",betamin,betamax);
     fprintf(stderr,"zetamin     %.0f, zetamax    %.0f\n",zetamin,zetamax);
     fprintf(stderr,"cluster half-widths, & special satellite widths");
     fprintf(stderr,"    version: %s\n",clusterhalfwidthsversion);
     fprintf(stderr,"deltamw  %.0f\n",deltamw);
     fprintf(stderr,"epsilonw %.0f, epsilonsatw %.0f\n",epsilonw,epsilonsatw);
     fprintf(stderr,"zetaw    %.0f, zetasatw    %.0f\n",zetaw,zetasatw);
     fprintf(stderr,"alphaw   %.0f, alphasatw   %.0f\n",alphaw,alphasatw);
     fprintf(stderr,"betaw    %.0f, betasatw    %.0f\n",betaw,betasatw);
     fprintf(stderr,"gammaw   %.0f\n",gammaw);
     fprintf(stderr,"deltaw   %.0f\n",deltaw);

         for(k=0; k<MAXSAT; k++)
         {
            if(satinfo[k].name[0] != '\0')
            {
               fprintf(stderr,"satellite: %s, satw: ",satinfo[k].name);
               for(m=1; m<8; m++)
               {
                  fprintf(stderr,"%3.0f ",satinfo[k].satw[m]);
               }
               fprintf(stderr,", domw: ");
               for(m=1; m<8; m++)
               {
                  fprintf(stderr,"%3.0f ",satinfo[k].domw[m]);
               }
               fprintf(stderr,": %s\n",satinfo[k].doma);
            }
         }
      }
      else if (Luseincorrect) /* if both the point field and alt field are not specified, then exit - S.J. 01/07/2014*/
      {
	fprintf(stderr,"Please specify both -altIDfield and -pointIDfields\n");
        usageout();
      }
      else /*Lhelpout*/
      {
         usageout();
      }
   }
   return(0);
}
/*___main()__________________________________________________________________*/

/****evaluatesuite()**********************************************************/
int  evaluatesuite(void)
{
   int  ibin=0;  /* 0 for no bin assignment*/
   int  jclst=0; /* 0 for no cluster assignment, i.e. for triaged suite*/
   char sour[32];      /*flag type of non-suiteness */
   int  LOK=1; /*presume OK */
   int  puckerdm=0,puckerd=0;
   char gammaname=' ';
   float suiteness=0, distance=0;
   int  ddg=0; /*deltam,delta,gamma   12 group numbers */
   char ptmaster[4];
   char ptcolor[16];

   /*enters with LOK = 1 */
   Ltriage=0; Liswannabe=0; Lcomment=0; /*reset flags 070521,070525,070628*/
   ptcolor[0]  = '\0';  /*default is no point color*/
   ptmaster[0] = '\0';  /*default is no point master*/

   sour[0] = '\0'; /* flagless */

  /*sudo program body: evaluate suite is a decision cascade*/
   /*sudo order important, epsilon and delta triage faulty pucker */

   if(LOK) /*------ filter on epsilon -----------------------------*/
   {
    /*sudo IF si(e) < e_min OR si(e) > e_max */
       /*sudo THEN s_bini = e_triage AND !OK */

      if(suiteptr->epsilon < epsilonmin || suiteptr->epsilon > epsilonmax)
      { sprintf(sour," e out  "); LOK = 0;
        Ltriage = EPSILONM; sprintf(ptmaster,"'E'");}
   }
   if(LOK) /*------ filter on delta minus 1 -----------------------*/
   {
    /*sudo IF      OK AND si(dm) >= d3_min AND si(dm) <= d3_max */
       /*sudo THEN s_dm = 3 AND OK */
    /*sudo ELSE IF OK AND si(dm) >= d2_min AND si(dm) <= d2_max */
       /*sudo THEN s_dm = 2 AND OK */
    /*sudo ELSE s_bini = dm_triage AND !OK */

      if(suiteptr->deltam >=  delta3min && suiteptr->deltam <= delta3max)
      { puckerdm = 3; LOK = 1; }
      else if(suiteptr->deltam >= delta2min && suiteptr->deltam <= delta2max)
      { puckerdm = 2; LOK = 1; }
      else
      { puckerdm = 0; sprintf(sour," bad deltam "); LOK = 0;
        Ltriage = DELTAM; sprintf(ptmaster,"'D'");}
   }
   if(LOK) /*------ filter on delta        -----------------------*/
   {
    /*sudo IF      OK AND si(d) >= d3_min AND si(d) <= d3_max */
       /*sudo THEN s_d = 3 AND OK */
    /*sudo ELSE IF OK AND si(d) >= d2_min AND si(d) <= d2_max */
       /*sudo THEN s_d = 2 AND OK */
    /*sudo ELSE s_bini = d_triage AND !OK */

      if(suiteptr->delta >=  delta3min && suiteptr->delta <= delta3max)
      { puckerd = 3; LOK = 1; }
      else if(suiteptr->delta  >= delta2min && suiteptr->delta  <= delta2max)
      { puckerd = 2; LOK = 1; }
      else
      { puckerd = 0; sprintf(sour," bad delta "); LOK = 0;
        Ltriage = DELTA; sprintf(ptmaster,"'D'");}
      if(LOK)
      {/*both deltas in range*/
         sprintf(sour," %1d%1d delta ",puckerdm,puckerd);
      }
   }
   if(LOK) /*------ filter on gamma        -----------------------*/
   {
    /*sudo IF      OK AND si(g) >= gp_min AND si(g) <= gp_max */
       /*sudo THEN s_g = p AND OK */
    /*sudo ELSE IF OK AND si(g) >= gt_min AND si(g) <= gt_max */
       /*sudo THEN s_g = t AND OK */
    /*sudo ELSE IF OK AND si(g) >= gm_min AND si(g) <= gm_max */
       /*sudo THEN s_g = m AND OK */
    /*sudo ELSE s_bini = g_triage AND !OK */

      if(suiteptr->gamma >= gammapmin && suiteptr->gamma <= gammapmax)
      { gammaname = 'p'; LOK = 1; }
      else if(suiteptr->gamma >= gammatmin && suiteptr->gamma <= gammatmax)
      { gammaname = 't'; LOK = 1; }
      else if(suiteptr->gamma >= gammammin && suiteptr->gamma <= gammammax)
      { gammaname = 'm'; LOK = 1; }
      else
      { gammaname = 'o'; LOK = 0; sprintf(sour," g out  ");
        Ltriage = GAMMA; sprintf(ptmaster,"'T'");}
   }
   if(LOK) /*------ filter on alpha        -----------------------*/
   {
    /*sudo IF OK AND si(a) >= a_min AND si(a) <= a_max */
       /*sudo THEN OK */
    /*sudo ELSE s_bini = a_triage AND !OK */

      if(suiteptr->alpha >= alphamin && suiteptr->alpha <= alphamax)
      {LOK = 1; }
      else
      {LOK = 0; sprintf(sour," a out  ");
        Ltriage = ALPHA; sprintf(ptmaster,"'T'");}
   }
   if(LOK) /*------ filter on beta         -----------------------*/
   {
    /*sudo IF OK AND si(b) >= b_min AND si(b) <= b_max */
       /*sudo THEN OK */
    /*sudo ELSE s_bini = b_triage AND !OK */

      if(suiteptr->beta  >= betamin && suiteptr->beta  <= betamax)
      {LOK = 1; }
      else
      {LOK = 0; sprintf(sour," b out  ");
        Ltriage = BETA; sprintf(ptmaster,"'T'");}
   }
   if(LOK) /*------ filter on zeta         -----------------------*/
   {
    /*sudo IF OK AND si(z) >= z_min AND si(z) <= z_max */
       /*sudo THEN OK */
    /*sudo ELSE s_bini = z_triage AND !OK */

      if(suiteptr->zeta  >= zetamin && suiteptr->zeta  <= zetamax)
      {LOK = 1; }
      else
      {LOK = 0; sprintf(sour," z out  ");
        Ltriage = ZETAM; sprintf(ptmaster,"'T'");}
   }

   ddg = 0; /*NO bin assigned yet...*/
   /*sudo IF OK ASSIGN s_bini to one of 12 delta,delta,gamma bins */
   if(LOK) /*------ categorize as one of the 12 ddg groups ------*/
   {/*OK to place in one of the 12 bins*/
      /*    ddg     */
      /*    0   out */
      /* 33 123 ptm */
      /* 32 456 ptm */
      /* 23 789 ptm */
      /* 22 ABC ptm */


      /*sudo IF s_dm == 3 */
      if     (puckerdm == 3)
      {
         /*sudo IF s_d == 3 */
         if     (puckerd == 3)
         {
            L33out = 1;
            /*sudo IF      s_g == p THEN s_bini = 33p */
            if     (gammaname == 'p') {ddg =  1; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == t THEN s_bini = 33t */
            else if(gammaname == 't') {ddg =  2; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == m THEN s_bini = 33m */
            else if(gammaname == 'm') {ddg =  3; binout[ddg] = 1;}
         }
         /*sudo ELSE IF s_d == 2 */
         else if(puckerd == 2)
         {
            L32out = 1;
            /*sudo IF      s_g == p THEN s_bini = 32p */
            if     (gammaname == 'p') {ddg =  4; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == t THEN s_bini = 32t */
            else if(gammaname == 't') {ddg =  5; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == m THEN s_bini = 32m */
            else if(gammaname == 'm') {ddg =  6; binout[ddg] = 1;}
         }
      }
      /*sudo ELSE IF OK AND s_dm == 2 */
      else if(puckerdm == 2)
      {
         /*sudo IF s_d == 3 */
         if     (puckerd == 3)
         {
            L23out = 1;
            /*sudo IF      s_g == p THEN s_bini = 23p */
            if     (gammaname == 'p') {ddg =  7; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == t THEN s_bini = 23t */
            else if(gammaname == 't') {ddg =  8; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == m THEN s_bini = 23m */
            else if(gammaname == 'm') {ddg =  9; binout[ddg] = 1;}
         }
         /*sudo ELSE IF s_d == 2 */
         else if(puckerd == 2)
         {
            L22out = 1;
            /*sudo IF      s_g == p THEN s_bini = 22p */
            if     (gammaname == 'p') {ddg = 10; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == t THEN s_bini = 22t */
            else if(gammaname == 't') {ddg = 11; binout[ddg] = 1;}
            /*sudo ELSE IF s_g == m THEN s_bini = 22m */
            else if(gammaname == 'm') {ddg = 12; binout[ddg] = 1;}

         }
      }
      sprintf(sour," ddg== %2d ",ddg);
   }/*OK to place in one of the 12 bins*/
   else
   {/*culled by single angle triage*/

    /*sudo ELSE {s_bini already set to be some kind of triage} */

      Ltriageout = 1;
      binout[0] = 1;
      ibin = 0;
      jclst = 0;
      suiteness = 0;       /*ibin,jclst,suiteness default to zero anyway*/
      clusterout[0][0] = 1; /*triage bin, cluster number 0 */
      writesuite(ibin,jclst,sour,distance,suiteness,ptmaster,ptcolor);
        /*binname "trig"  clustername "!!" 070414 */
   }
   return(ddg); /* ddg is the ibin number! */
}
/*___evaluatesuite()_________________________________________________________*/

/****membership()*************************************************************/
void membership(int ibin)
{
   int  i=0,j=0,k=0,minj=0,matchcnt=0,domj=0,thej=0;
   float dsq = 0, del = 0, distance = 0, mindistance = 999;
   float disttodom=0, disttosat=0, thedist=0;
   float penmindist = 999;
   float pi = 3.14159;
   float suiteness=0;
   char sour[96],sourness[64]; /*070311 accum report*/
   int  Ldominant=0;
   int  LOK = 0;
   float sattodom[8], domtosat[8],sattopt[8],domtopt[8];
   char ptmaster[4];
   int   closestj = 0, nextclosej=0;
   float closestd = 999, nextclosed = 999;
   char ptcolor[16];
   int  Lassigned = 0; /*070430 status of this point...*/

   ptcolor[0]  = '\0';  /*default is no point color*/
   sour[0] = '\0'; /*extra comment in outputed ptID*/
   sourness[0] = '\0'; /*extra comment in outputed ptID*/
   ptmaster[0] = '\0';  /*default is no point master*/

   clearmatches();

   matchcnt = 0;
   i = ibin;

   /*sudo */
   /*sudo program body: cluster membership in assigned s_bini */

   /*logic: find all clusters in this bin that match this conformation*********/

   /*bin[].clst[].ang  working-arrays count from 1 for known clusters*/

   /* some bins have a major cluster with satellites */
   /*--first clear the flags */
   domj = 0; /*default, no dominant cluster in the bin*/
   Ldominant = 0; /*conformation close to a major (dominant) cluster*/

   mindistance = 999; /*especially note closest cluster*/
   closestd = 999;

   /*sudo FOR each cluster cj in this delta,delta,gamma bin */
   /*for(j=1; j<=ddgcnt[i]; j++)*/
   j = 1; /*start with index 1 for first observed clst in this bin  070428*/

   while(bin[i].clst[j].LOK) /*070428*/
   {/*loop over all observed clusters in this ith bin*/
      if(   !Lwannabe
         && strcmp(bin[i].clst[j].status,"wannabe")==0) {j++;continue;}/*match*/
      /*empty bin dummy cluster at 0,0,0,0,0,0,0 is out of angle range*/
      /*but safer to explicity skip*/
      /*if(i==3 || i==9){continue;}  not needed  070428*/

     /*sudo COMPUTE d=scaled 4D distance of si to cj (epsilon,zeta,alpha,beta)*/

      resetcoordw(coordw,i,j,Lgeneralsatw); /*revert to general defaults*/
         /*insert this resetcoordw() here 070414 */
      distance = hyperellipsoiddist(i,j,4,coordw);
         /*4D scaled, normalized distance*/
      if(distance < closestd)   /*070311*/
      {
         closestd = distance;
         closestj = j;
      }

      if(distance >= 0) /* store all distances, needed for diagnostics */
      {
        matches[i][j] = distance;
        /*sudo IF 0 <= d < 1 THEN */
        if(distance < 1) /* suite could be a member of this cluster */
        {   /*APPEND to list of clusters of which suite could be a member*/
            /*sudo APPEND cj to s_list_ci */
            /*sudo FLAG when possible_cluster is a dominant_cluster */

            matchcnt++; /*only count possible clusters*/
            /*if(j== domj) {Ldominant = 1;}*/
            if(strcmp(bin[i].clst[j].domsatness,"dom")==0) /*match 070429*/
            {
               domj = j; /*this j is a dominant cluster*/
               Ldominant = 1; /*there is a close dominant cluster to consider*/
            }
        }
        /*sudo END IF */
        /*sudo FLAG closest_cluster */
        if(distance < mindistance)
        {
           mindistance = distance;
           minj = j;
        }
      }
      j++; /*increment index for next cycle...*/
   }/*loop over all clusters in this ith bin*/
   /*sudo END FOR each cluster cj in this delta,delta,gamma bin */

   /*now find the next closest cluster (for debugging purposes)*/
   nextclosed = 999; /*next closest cluster*/

   j = 1; /*start with index 1 for first observed clst in this bin  070428*/
   while(bin[i].clst[j].LOK) /*070428*/
   {/*loop over all clusters in this ith bin*/
      if(   !Lwannabe
         && strcmp(bin[i].clst[j].status,"wannabe")==0) {j++;continue;}/*match*/
      if(matches[i][j] < nextclosed && j != closestj)
      {
         nextclosed = matches[i][j];
         nextclosej = j;
      }
      j++; /*increment index for next cycle...*/
   }/*loop over all clusters in this ith bin*/

   /*end logic: find all clusters that match this conformation_______________*/

   /*now assign suite to a particular cluster*/
   Lassigned = 0;

   /*sudo */
   /*sudo IF      |s_list_ci| == 1 THEN */
   if(matchcnt == 1)
   { /*output the only distance match*/

      /*sudo ASSIGN s_cluster = cj in s_list_ci */

      thej = minj; /*cluster number in this bin*/
      Lassigned = 1;
      sprintf(sourness,"%d-only-one",matchcnt);
      thedist = matches[ibin][thej];
   }
   /*sudo ELSE IF |s_list_ci| > 1 THEN */
   else if(matchcnt > 1 )
   {/*suite close to more than one cluster*/

      /*sudo IF all cj in s_list_ci are not dominant THEN */
      if(!Ldominant) /*dominant_cluster is not a member of possible_clusters */
      { /*output the minimum distance match*/

         /*sudo ASSIGN s_cluster = cj in s_list_ci minimizing d */

         thej = minj; /*cluster number in this bin*/
         Lassigned = 1;
         sprintf(sourness,"%d-none-dom",matchcnt);
      }

      /*sudo ELSE IF s_list_ci contains a dominant cj THEN */
      else if(Ldominant)
      {/*match dominant cluster and highest quality other cluster*/

         /*sudo ASSIGN domj = dominant cj in s_list_ci */
         /*sudo ASSIGN thej = cj in s_list_ci minimizing d for all non-dominant cj in s_list_ci */

         penmindist = 999;

         j = 1;/*start with index 1 for first observed clst in this bin 070428*/
         while(bin[i].clst[j].LOK) /*070428*/
         {/*loop over all clusters in this ith bin*/
            if(   !Lwannabe
               && strcmp(bin[i].clst[j].status,"wannabe")==0) {j++;continue;}

            if(   ( strcmp(bin[i].clst[j].domsatness,"dom")!=0 )/*i.e. NOT dom*/
               && (matches[i][j] < penmindist))
            {/*find shortest distance to non-dominant and satellite clusters*/
               penmindist = matches[i][j];
               thej = j; /*j index of closest other cluster*/
            }
            j++; /*increment index for next cycle...*/
         }

         /*distinguish best match cluster as satellite or just another cluster*/
         /*this bin has a dominant cluster and thus satellite clusters*/

         /*sudo IF thej is a satellite of domj THEN */

         if( strcmp(bin[i].clst[thej].domsatness,"sat")==0 )/*i.e. matches sat*/
         {  /*the other cluster is a satellite cluster*/
            /*point close to both a satellite and the dominant cluster*/

            /*NEED TO KNOW IF INBETWEEN DOMINANT AND SATELLITE CLUSTERS*/

            /*if dotproducts both positive, then inbetween*/
            /*     p                  p      */
            /* dom/___sat  and   dom___\sat  */

            vector7ab(domtopt, bin[i].clst[domj].ang, suiteptr->ang);

            vector7ab(sattopt, bin[i].clst[thej].ang, suiteptr->ang);

            vector7ab(domtosat, bin[i].clst[domj].ang, bin[i].clst[thej].ang);

            vector7ab(sattodom, bin[i].clst[thej].ang, bin[i].clst[domj].ang);

            /*sudo IF si is located between the means of thej and domj THEN */
            if(   (dotproduct(domtopt,domtosat,4) > 0)
                &&(dotproduct(sattopt,sattodom,4) > 0) )
            {/*pt between dom and sat*/

               /*sudo CALCULATE rescaled distances d* of si to thej and domj */
               /*sudo ASSIGN s_clusteri = thej or domj minimizing rescaled distance d* */

               /*recalc distances to dominant cluster and to satellite cluster*/
               /*considering different weights between special cluster pairs */
               resetcoordw(dominantw, i,domj,0); /*set to dominant  widths*/
               resetcoordw(satellitew,i,thej,Lgeneralsatw);
                 /*set to satellite widths, using Lgeneralsatw  070414*/
               assignweights(i,thej,dominantw,satellitew);
                  /*special widths for this satellite: weight dom--sat pair*/
/*
fprintf(stderr,"satj: %d satw: ",thej);
for(k=1; k<=7; k++) {fprintf(stderr,"%3.0f ",satellitew[k]);}
fprintf(stderr,", domj: %d domw: ",domj);
for(k=1; k<=7; k++) {fprintf(stderr,"%3.0f ",dominantw[k]);}
fprintf(stderr,"\n");
*/
               disttodom = hyperellipsoiddist(i,domj,4,dominantw);
               disttosat = hyperellipsoiddist(i,thej,4,satellitew);

               if(disttosat <= disttodom) {thej = thej;} /*bookkeeping*/
               else {thej = domj;} /*in effect, the default*/
               Lassigned = 1;
               sprintf(sourness,"%d-BETWEEN-dom-sat(%7.3f|%7.3f)"
                        ,matchcnt,disttodom,disttosat);
            }/*pt between dom and sat*/
            /*sudo ELSE (not between the means of thej and domj) */
            else
            {/*pt NOT inbetween*/
               /*just assign by closest standard distance evaluation*/

               /*sudo ASSIGN s_clusteri = thej or domj minimizing d */

               if(matches[i][domj] < matches[i][thej]) {thej = domj;}
               Lassigned = 1;
               sprintf(sourness,"%d-OUTSIDE-dom-sat",matchcnt);
            }/*pt NOT inbetween*/
            /*sudo END IF (between or not between satellite and dominant) */

         }/*point close to both a satellite and the dominant cluster*/

         /*sudo ELSE (thej is NOT a satellite of domj) */
         else /* near cluster is not a satellite */
         {  /*point just belongs to the closest cluster*/
            /*output the minimum distance match*/

            /*sudo ASSIGN s_clusteri = thej or domj minimizing d */

            thej = minj; /*cluster number in this bin*/
            Lassigned = 1;
            sprintf(sourness,"%d-not-sat",matchcnt);
         }
         /*sudo END IF (satellite or not a satellite) */
      }/*match dominant cluster and highest quality other cluster*/
      /*sudo END IF (a cj in s_list_ci is or is not dominant) */

      thedist = matches[ibin][thej];
   }/*suite close to more than one cluster*/
   else
   /*sudo ELSE IF |s_list_ci| == 0 THEN */
   { /*NO match in this bin*/

      /*sudo ASSIGN s_clusteri = outlier */

      thej = 0; /*outlier named !!*/
      Lassigned = 0;
      Loutlier = 1;
      thedist = closestd;
      sprintf(sourness,"outlier distance %.3f",closestd); /*070311*/
      sprintf(ptmaster,"'O'"); /*outlier within a bin*/
      sprintf(ptcolor," white"); /*leading space*/
      /*thej = closestj; NOT viable to have outlier named by closest cluster*/
      /* some are very far from any cluster, */
      /* and even a singleton defines an extant cluster and makes a ring*/
   }
   /*sudo END IF ( |s_list_ci| == ? ) */

   sprintf(sour," %s:%s, %7.3f:%s, %7.3f: "
           ,sourness
           ,bin[i].clst[closestj].clustername,matches[i][closestj]
           ,bin[i].clst[nextclosej].clustername,matches[i][nextclosej]);

   clusterout[i][thej] = 1; /*this specific cluster has an entry*/

   if(Lwannabe && strcmp(bin[i].clst[thej].status,"wannabe")==0)
   {
      Lwannabeout = 1; /*once set, stays set*/
      Liswannabe = 1;  /*set only for this cluster*/
   }

   /*----------- final suiteness ------------------------------------------*/
   /*sudo COMPUTE dist_si = scaled 7D distance of si to cluster s_clusteri */
   /*sudo possible reassignment based on 7D distance */

   resetcoordw(coordw,i,thej,Lgeneralsatw); /*revert to general defaults*/

   /*recompute distance for all 7 dimensions, calc suiteness*/
   if(Lassigned){distance = hyperellipsoiddist(i,thej,7,coordw);}
   else         {distance = hyperellipsoiddist(i,closestj,7,coordw);}

   /*sudo IF s_clusteri not outlier THEN */
   /*sudo IF dist_si <= 1 THEN */
   if(distance <= 1)  /* less than ...   or equal added 070430 */
   {
      /*sudo COMPUTE "suiteness" suiti = (cos(pi*dist_si) +1)/2 */
      /*sudo limit   "suiteness" suiti >= 0.01                  */
      /*sudo ASSIGN s_clusteri = thej                           */

      suiteness = (cos(pi*distance) +1)/2;
      if(suiteness < 0.01) {suiteness = 0.01;} /*floor  070430*/

      if(!Lassigned)
      {
         thej = closestj; /*assign in this case*/
         if(Lreportout)
         {
            /*"7D distance forces assignment to closest cluster\n");*/
            sprintf(commentstr," by 7Ddist"); /*070628*/
            Lcomment=1;
         }
      }
   }
   /*sudo ELSE "suiteness" suiti = 0 */
      /*sudo ASSIGN s_clusteri = outlier */
   else
   {
      if(Lassigned)
      {
         if(Lreportout)
         {  /*"7D distance forces %s assignment to be outlier\n"*/
            sprintf(commentstr," 7D dist %s",bin[ibin].clst[thej].clustername);
            Lcomment=1; /*070628*/
         }
      }
      thej = 0; /*make not assigned in any case! */
      suiteness = 0;
   }
   /*Decisions are made piecemeal (ddg) or in 4D (abez)                       */
   /*  inclusion is definite to a particular cluster, based on 4D distances   */
   /*  exclusion is general, except for knowledge of closest in 4D            */
   /* final suiteness is based on 7D distance                                 */
   /* but it is very difficult to understand shape in 7D                      */
   /*Use 7D suiteness to do final choice about inclusion:                     */
   /*case: currently assigned membership in a cluster:                        */
   /*   calc 7D dist to assigned cluster:  retain as member or make an outlier*/
   /*case: currently an outlier:                                              */
   /*   calc 7D dist to closest  cluster: assign as member or leave as outlier*/
   /*----------------------------------------------------------------------*/

   writesuite(ibin,thej,sour,distance,suiteness,ptmaster,ptcolor);
        /*binname "## c"  clustername "#c" ; outlier "!!" 070414 */
   if(Ltestout)
   {
      fprintf(fpout," [suite: %s %s power== %4.2f, 4Ddist== %f, 7Ddist== %f, suiteness==%f] \n",bin[ibin].clst[thej].clustername,suiteptr->ptID,power,thedist,distance,suiteness);
   }
}
/*___membership()____________________________________________________________*/

/*sudo END FOR */

/****clearmatches()***********************************************************/
void clearmatches(void)  /*called from membership()*/
{
   int  i=0,j=0;

   for(i=1; i<=12; i++) /*12 actual named bins, this is robust within MAXBINS*/
   {
      /*named clusters start at index 1, no tally of actual # defined clusters*/
      for(j=1; j<MAXCLST; j++)  /*MAXCLST should be > # defined clusters*/
      {
         matches[i][j] = 999; /*going for minimum of match distance*/
      }
   }
}
/*___clearmatches()__________________________________________________________*/

