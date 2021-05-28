/*                           suiteninit.c                                    */

#include "suitename.h"
#define SUITENINIT
#include "suiteninit.h"
#undef  SUITENINIT
#include "suitenscrt.h"
#include "suitenutil.h"
#include "suiteninpt.h"
#include "suitenout.h"

/****initializations()********************************************************/
int  initializations(void)
{
   int   i=0,j=0,LOK=1; /*070426 not using LOK here, was for buildclusterav*/

   Lcomment = 0; /*report 7D distance forces outlier or assignment  070628*/
   Ltriage = 0; Loutlier = 0; /*diagnostic flags*/
   L33out = 0; L32out = 0; L23out = 0; L22out = 0; Ltriageout = 0; /*group*/
   LTdeltam=0;LTdelta=0;LTepsilon=0;LTzeta=0;LTalpha=0;LTbeta=0;LTgamma=0;
      /* LT...flags accummulate for pointmaster definition 070628*/
   Liswannabe=0; Lwannabeout = 0; /*070429*/
   clearbinout();       /*flags each found bin: subgroup */
   clearclusterout();   /*flags each found cluster: list */
   resetcoordw(coordw,0,0,0); /*default general cluster half-widths along axes*/

   return(LOK);
}
/*___initializations()_______________________________________________________*/

/****assignweights()**********************************************************/
void assignweights(int ibin, int jth, float* domWarray, float* satWarray)
{  /*called from membership() for a dom-sat pair: empirical values to edit */
   
   /*observation is only one dominant cluster in any one bin*/
   /* so test on common bin number and number of satellite cluster*/
   /*    1       2       3       4       5       6       7   */
   /* deltam  epsilon  zeta    alpha   beta    gamma   delta */
   /*special pair half-widths*/
   /*dominant Width toward satellite, satellite Width toward dominant*/
   /* dom    sat   w,    sat    dom   w;    ibin domj     satj  angle(s)*/

   /*rewritten to use initialized suiteninit.h arrays  070429*/
   int k=0,m=0;

   /*sanity check: */
   if(strcmp(bin[ibin].clst[jth].domsatness,"sat")==0) /*match*/
   {
      for(k=0; k<MAXSAT; k++)
      {
         if(strcmp(satinfo[k].name,bin[ibin].clst[jth].clustername)==0)/*match*/
         {
/*fprintf(stderr,"assignweights called with bin[%d].clst[%d].domsatness: %s, clstname: %s\n",ibin,jth,bin[ibin].clst[jth].domsatness,satinfo[k].name);*/
            for(m=0; m<9; m++)
            {
               if(satinfo[k].satw[m] > 0) 
{
satWarray[m] = satinfo[k].satw[m];
/*fprintf(stderr,"sat ang %d = %3.0f ",m,satWarray[m]);*/ 
}
               if(satinfo[k].domw[m] > 0) 
{
domWarray[m] = satinfo[k].domw[m];
/*fprintf(stderr,"dom ang %d = %3.0f \n",m,domWarray[m]);*/ 
}
            }
            break; /*just find one*/
         }
      }      
   }
   /*no match just leaves arrays unmodified*/
}
/*___assignweights()_________________________________________________________*/

/****parsecommandline()*******************************************************/
int  parsecommandline(int *argc, char** argv)
{ 
   int LOK=0;
   char *p;
   int  i=0,j=0,k=0;
   char numstr[256];

   for(i=1; i<*argc; i++)
   {/*loop over arguments*/
      p = argv[i];
      if(p[0] == '-')
      {/*flag*/
         if(p[1] == '\0')
         {/*naked flag ignored*/
            ;
         }/*naked flag ignored*/
         else
         {/*interpret flag*/
            if(   CompArgStr(p+1,"residuein",  7)
               || CompArgStr(p+1,"residuesin", 7))
            {
               Lsuitesin = 0;
               Lresiduesin = 1;
            }
            else if(   CompArgStr(p+1,"suitein",  5)
                    || CompArgStr(p+1,"suitesin", 5))
            {
               Lsuitesin = 1;
               Lresiduesin = 0;
            }
            else if(CompArgStr(p+1,"string", 6))
            {
               Lstringout = 1; 
               Lreportout = 0; 
               Lkinemageout = 0; 
               Lchart = 0;
            }
            else if(CompArgStr(p+1,"report", 6))
            {
               Lreportout = 1; 
               Lstringout = 0; 
               Lkinemageout = 0; 
               Lchart = 0;
            }
            else if(CompArgStr(p+1,"chart", 5))
            {
               Lreportout = 1; 
               Lstringout = 0; 
               Lkinemageout = 0; 
               Lchart = 1;  /*no summary of report, MolProbity multichart use*/
            }
            else if(CompArgStr(p+1,"kinemage", 3))
            {
               Lkinemageout = 1; 
               Lreportout = 0; 
               Lstringout = 0; 
            }
            else if(CompArgStr(p+1,"satellites", 9))
            {
               Lgeneralsatw = 1; 
            }
            else if(CompArgStr(p+1,"wannabes", 7)) /*070429*/
            {
               Lwannabe = 1; 
            }
            else if(CompArgStr(p+1,"nowannabes", 9)) /*070525*/
            {
               Lwannabe = 0; 
            }
            else if(CompArgStr(p+1,"nosequence", 5))
            {
               Lsequence = 0; 
            }
            else if(CompArgStr(p+1,"overlaps", 7))
            {
               Loverlap = 1; 
            }
            else if(CompArgStr(p+1,"oneline", 7))
            {
               Loneline = 1; 
            }
            else if(CompArgStr(p+1,"help", 1))
            {
               Lhelpout = 1; 
            }
            else if(CompArgStr(p+1,"changes", 7))
            {
               Lchangesout = 1; 
            }
            else if(CompArgStr(p+1,"test", 4))
            {
               Ltestout = 1; 
            }
            else if(CompArgStr(p+1,"sour", 4)) /*070524*/
            {
               Lsourout = 1; 
            }
            else if(CompArgStr(p+1,"dangle", 6)) /*070525*/
            {
               Ldangle = 1; 
            }
            else if(   CompArgStr(p+1,"thetaeta", 5)  /*070524*/
                    || CompArgStr(p+1,"etatheta", 3)) /*070524*/
            {
               Letatheta = 1; 
            }
            else if(   CompArgStr(p+1,"pointIDfields", 5) 
                    || CompArgStr(p+1,"ptID",  2) )
            {
	       LNptIDfields=1; /* flag to see if this is specified or not - S.J. 01/07/2014*/
               i = i+1; /*number of pointID fields is in the next string*/
               p = argv[i];
               j=0;
               k=0;
               while(j<256)
               {/*strobe out number*/
                 if( isdigit(p[j]) ) { numstr[k] = p[j]; k++; j++; }
                 else if(k==0) {/*no char yet*/ j++; }
                 else
                 {/*presume whole number is in*/
                    numstr[k]='\0'; /*end number string*/
                    sscanf(numstr,"%d",&NptIDfields);
                    j=999; /*end the while loop*/
                 }
               }/*strobe out number*/
            }
            else if(CompArgStr(p+1,"altIDfield", 10)) /* to check for altIDfield and storing it S.J. 01/07/2014*/
            {
	       LaltIDfield=1; 
               i = i+1; /*altID field is in the next string*/
               p = argv[i];
               j=0;
               k=0;
               while(j<256)
               {/*strobe out number*/
                 if( isdigit(p[j]) ) { numstr[k] = p[j]; k++; j++; }
                 else if(k==0) {/*no char yet*/ j++; }
                 else
                 {/*presume whole number is in*/
                    numstr[k]='\0'; /*end number string*/
                    sscanf(numstr,"%d",&altIDfield);
                    j=999; /*end the while loop*/
                 }
               }/*strobe out number*/
            }
	    else if(CompArgStr(p+1,"altIDval",8)) /* to get the value of the altID - S.J. 01/07/2014*/
            {
                  LaltID=1;
                  i = i+1;
                  p=argv[i];
                  altID[0]=p[0];
            }
            else if(   CompArgStr(p+1,"angles",  5) 
                    || CompArgStr(p+1,"anglefields", 5) 
                    || CompArgStr(p+1,"nangles", 6) )
            {
               i = i+1; /*number of angle fields is in the next string*/
               p = argv[i];
               j=0;
               k=0;
               while(j<256)
               {/*strobe out number*/
                 if( isdigit(p[j]) ) { numstr[k] = p[j]; k++; j++; }
                 else if(k==0) {/*no char yet*/ j++; }
                 else
                 {/*presume whole number is in*/
                    numstr[k]='\0'; /*end number string*/
                    sscanf(numstr,"%d",&Nanglefields);
                    j=999; /*end the while loop*/
                 }
               }/*strobe out number*/
            }
            else if(   CompArgStr(p+1,"hyper", 5) 
                    || CompArgStr(p+1,"power", 5) )
            {
               i = i+1; /*hyperellipsoid power is in the next string*/
               p = argv[i];
               j=0;
               k=0;
               while(j<256)
               {/*strobe out number*/
                 if( isdigit(p[j]) || p[j]=='.') { numstr[k] = p[j]; k++; j++; }
                 else if(k==0) {/*no char yet*/ j++; }
                 else
                 {/*presume whole number is in*/
                    numstr[k]='\0'; /*end number string*/
                    sscanf(numstr,"%lf",&power);
                    j=999; /*end the while loop*/
                 }
               }/*strobe out number*/
            }
         }/*interpret flag*/
      }/*flag*/
      else
      {/*presume an input file name*/
        /*THIS IS NOT IMPLEMENTED: NO CODE TO OPEN A FILE !!!! */
         strcpy(NameStr,argv[i]); /*copy name into input file Name*/
         Lnewfile=1;   /*file is present*/
      }/*presume an input file name*/
   }/*loop over arguments*/
   /*070222 much of this not yet implemented*/
   if(LaltIDfield == 1 && LNptIDfields == 0) {Luseincorrect=1;} /* added to check is both altID and pointID are specified - S.J. 01/07/2014*/
   if( Lhelpout || Ltestout || Lchangesout || Lnewfile || Luseincorrect) {LOK = 0;} /* added Luseincorrect check S.J. 01/07/2014*/
   else {LOK = 1;}
   return(LOK);
}
/*___parsecommandline()______________________________________________________*/

