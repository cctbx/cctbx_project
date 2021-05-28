/*                  suitenout.c                  */

#include "suitename.h"
#define SUITENOUT
#include "suitenout.h"
#include "suitenhead.h" /* kinemageframe[] 070414*/
#undef  SUITENOUT
#include "suitenscrt.h"
#include "suiteninit.h"

/****writesuite()*************************************************************/
void writesuite(int ibin,int jclst,char* sour, float distance, float suiteness,char* ptmaster,char* ptcolor)
{
   int  j=0,n=0;
   /*static int nout=1; Lstringout residue/suite counter, suitenout.h 070421*/
   char resstr[32];
   char basestr[2]={'\0','\0'};
   char clststr[4]={'\0','\0','\0','\0'};
   static char lappedstr[256]={'\0'};
   char reason[16]; /*used in report: reason for triage*/
   char stray[16]; /*used in report: stray wannabe*/
   char  sourpuss[1]; /*070524*/
   char* sourptr; /*070524*/
   
   sourpuss[0] = '\0'; /*070524*/

   if(Lkinemageout)
   {/*kinemage showing clusters*/

     if(Lsourout) {sourptr = sour;} /*070524*/
     else         {sourptr = sourpuss;}
     sprintf(temps,"{%s %s %s:D==%5.3f:S==%5.3f: %s} %s%s ,%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%7.2f,%7.2f",
     bin[ibin].binname,
     bin[ibin].clst[jclst].clustername,
     sourptr,
     distance,
     suiteness,
     suiteptr->ptID,
     ptmaster,
     ptcolor,
     suiteptr->chim,
     suiteptr->deltam,
     suiteptr->epsilon,
     suiteptr->zeta,
     suiteptr->alpha,
     suiteptr->beta,
     suiteptr->gamma,
     suiteptr->delta,
     suiteptr->chi   );

     if(ibin > 12 || ibin < 0)
     {
        fprintf(stderr,"%s\n",temps); /*e.g. 13: angles not fully specified*/
     }
     else
     {
        if(Ltest){fprintf(stderr,"PUT %s\n",temps);}

        putonetextblockline(&mainscratch,temps);
     }
   }
   else if(Lstringout)
   {/*3char string in order of input for all entries*/
      if(Lsequence) 
      { 
         basestr[0] = suiteptr->basechr[0]; /*070412 found previously*/ 
      }
      else /*colon instead of single-character Base code 070409*/
      {
         basestr[0] = ':';
      }

      fprintf(fpout,"%s%s",bin[ibin].clst[jclst].clustername,basestr);
      
      if(!Loneline)  /*070409*/
      {
         if(Loverlap && nout >= 11)
         {
            sprintf(clststr,"%s%s",bin[ibin].clst[jclst].clustername,basestr);
            j = n = 0;
            while(lappedstr[j] != '\0') {j++;} /*find end of lappedstr*/
            while((lappedstr[j++] = clststr[n++]) != '\0') ; /*copy on end*/
         }
         if(nout == 20)
         {
            fprintf(fpout," {%s}\n",suiteptr->ptID); 
            if(Loverlap)
            {
               fprintf(fpout,"%s",lappedstr);
               lappedstr[0] = '\0';
               nout = 11;
            }
            else
            {
               nout = 1;
            }
         }
         else {nout++;}
      }
   }
   else if(Lreportout)
   {
      /*report on all entries, even if suite is incomplete*/
      reason[0] = '\0'; /*default to none*/
      stray[0] = '\0'; /*default to none*/
      if(Ltriage==EPSILONM)   {sprintf(reason," epsilon-1");} /*070628*/
      else if(Ltriage==DELTAM){sprintf(reason," delta-1");} /*070628*/
      else if(Ltriage==DELTA) {sprintf(reason," delta");} /*070628*/
      else if(Ltriage==GAMMA) {sprintf(reason," gamma");} /*070521*/
      else if(Ltriage==BETA)  {sprintf(reason," beta");} /*070521*/
      else if(Ltriage==ALPHA) {sprintf(reason," alpha");} /*070521*/
      else if(Ltriage==ZETAM) {sprintf(reason," zeta-1");} /*070628*/
      else if(Lcomment){sprintf(reason,"%s",commentstr);} /*070628*/
      if(Liswannabe){sprintf(stray," wannabe");} /*070525*/
      fprintf(fpout,"%s %s %s %5.3f%s%s\n",
         suiteptr->ptID,
         bin[ibin].binname,
         bin[ibin].clst[jclst].clustername,
         suiteness,
         reason,
         stray   /*070628*/
         );

      reportcountall++; /*everything that suitename brought in */

      /*now accummulate sum and count for complete suites*/

      if(ibin == 0)
      {
         triagecountall++; /*070328*/
      }
      else if(ibin >= 0 && ibin < 13) /*ibin==13 when angles not all defined*/
      {
         /*defineable suites, including triaged and outliers */ 
         suitenesssumall = suitenesssumall + suiteness; 
            /*adding (float)zero won't hurt this sum*/
         binnedsuitecountall++; 
      }
      /*distribution of suiteness values...*/
      /*j==0 pseudo clusters always have a suiteness == 0, so sum not changed*/
      if(jclst == 0) /*augment the k==11 outlier place for this cluster*/ 
      {  /*outliers in bins 1 to 12, triages in bin 0, not defined in bin 13*/
         /*jclst membership 4D distance < 1, suiteness > 0 if  7D distance > 1*/
         /*outliers in jclst 0, suiteness can be == 0 for others if 7D dist>1 */
         suitenesscnt[ibin][jclst][11]++;
      }
      else
      {
         suitenesssum[ibin][jclst] = suitenesssum[ibin][jclst] + suiteness;

         /*3rd index for suiteness intervals: 10 at 10ths + extra at zero*/
      
         if(suiteness == 0                   ) {suitenesscnt[ibin][jclst][0]++;}
         if(suiteness >  0  && suiteness < .1) {suitenesscnt[ibin][jclst][1]++;}
         if(suiteness >= .1 && suiteness < .2) {suitenesscnt[ibin][jclst][2]++;}
         if(suiteness >= .2 && suiteness < .3) {suitenesscnt[ibin][jclst][3]++;}
         if(suiteness >= .3 && suiteness < .4) {suitenesscnt[ibin][jclst][4]++;}
         if(suiteness >= .4 && suiteness < .5) {suitenesscnt[ibin][jclst][5]++;}
         if(suiteness >= .5 && suiteness < .6) {suitenesscnt[ibin][jclst][6]++;}
         if(suiteness >= .6 && suiteness < .7) {suitenesscnt[ibin][jclst][7]++;}
         if(suiteness >= .7 && suiteness < .8) {suitenesscnt[ibin][jclst][8]++;}
         if(suiteness >= .8 && suiteness < .9) {suitenesscnt[ibin][jclst][9]++;}
         if(suiteness >= .9) {suitenesscnt[ibin][jclst][10]++;}
      }
   }
}
/*___writesuite()____________________________________________________________*/

/****suitenessaverage()******************************************************/
void suitenessaverage(int mode) /*mode for all, just Aform, etc. */
{                  
   int  i=0,j=0,k=0;
   int  ibin=0,nbin=0,jclst=0,nclst=0,xbin = -1,xclst = -1;
   double sum = 0;
   int    num[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
   float  average = 0;
   int    number = 0;
   char   comment[256];
   int    Ntriaged = 0; /*070328*/

   if(mode == 0) /*all complete suites, all ways*/
   {  /*bin==13 has the incomplete pseudo-suites with incomplete angles*/
      ibin = 1; nbin = 12; jclst = 0; nclst = 11;
      xbin = -1; xclst = -1;  /*exception when >= 0*/
      /*ibin 0 contains triaged clusters, never use these*/
   
      fprintf(fpout,"Found %d complete suites derived from %d entries\n"
        ,binnedsuitecountall+triagecountall,reportcountall);
      fprintf(fpout,"%d suites were triaged, leaving %d assigned to bins\n"
        ,triagecountall,binnedsuitecountall); /*070328*/
      sprintf(comment,"For all");

   }
   else if(mode == 1) /* A form   1,1*/
   {  /*bin==13 has the incomplete, pseudo-suites with incomplete angles*/
      ibin = 1; nbin = 1; jclst = 1; nclst = 1;
                         /*no outliers (j==0)*/
      xbin = -1; xclst = -1;  /*exception for this clst in this bin when >= 0*/
   
      sprintf(comment," A form (1a)");

   }
   else if(mode == 2) /* everything except A form, i.e. except 1,1*/
   {  /*bin==13 has the incomplete pseudo-suites with incomplete angles*/
      ibin = 1; nbin = 12; jclst = 1; nclst = 11; /*only complete suite pts*/
                         /*no outliers (j==0)*/
      xbin = 1; xclst = 1;  /*exception when >= 0*/
   
      sprintf(comment," non-1a  has");

   }
   for(i=ibin; i<=nbin; i++)
   {
      for(j=jclst; j<=nclst; j++) /*some clusters not defined in some bins*/
      {
            if(i != xbin || (i == xbin && j != xclst))
            {
               sum = sum + suitenesssum[i][j];
               for(k=0; k<12; k++)
               {
                  num[k] = num[k] + suitenesscnt[i][j][k];
               }
            }  
      }
   }

   for(k=0; k<12; k++) { number = number + num[k]; }

   if(number > 1){average = sum/number;}
   else {average = 0;}

   fprintf(fpout,"%s %d suites: average suiteness== %5.3f (power==%4.2f)\n"
     ,comment,number,average,power);
   if(mode == 0) {fprintf(fpout,"%6d suites are  outliers\n",num[11]);} 
   fprintf(fpout,"%6d suites have suiteness == 0    \n",num[0]); 
   fprintf(fpout,"%6d suites have suiteness >  0 <.1\n",num[1]); 
   fprintf(fpout,"%6d suites have suiteness >=.1 <.2\n",num[2]); 
   fprintf(fpout,"%6d suites have suiteness >=.2 <.3\n",num[3]); 
   fprintf(fpout,"%6d suites have suiteness >=.3 <.4\n",num[4]); 
   fprintf(fpout,"%6d suites have suiteness >=.4 <.5\n",num[5]); 
   fprintf(fpout,"%6d suites have suiteness >=.5 <.6\n",num[6]); 
   fprintf(fpout,"%6d suites have suiteness >=.6 <.7\n",num[7]); 
   fprintf(fpout,"%6d suites have suiteness >=.7 <.8\n",num[8]); 
   fprintf(fpout,"%6d suites have suiteness >=.8 <.9\n",num[9]); 
   fprintf(fpout,"%6d suites have suiteness >=.9    \n",num[10]); 
}
/*___suitenessaverage()______________________________________________________*/

/****writeoutput()************************************************************/
void   writeoutput(void)
{
   char commentstr[128];

   if(Lgeneralsatw) 
   {
       sprintf(commentstr,
         " special general case satellite widths, power = %4.2f",power);
   }
   else
   {
       sprintf(commentstr,
         " all general case widths, power = %4.2f",power);
   }
   if(Lreportout && !Lchart)
   {
      fprintf(fpout,"%s\n",commentstr);
      suitenessaverage(0); /* 0 mode reports on all suites */
      if(clusterout[1][1] > 0) /*Aform 1a    070325*/
      {
         suitenessaverage(1); /* 1 mode reports on A form suites */
         suitenessaverage(2); /* 2 mode reports on non-A suites */
      }
   }
   if(Lkinemageout &&(L33out ||  L32out ||  L23out ||  L22out || Ltriageout))
   {/*stuff for a kinemage*/
      /*Ltriageout for intact suites that fail individual angle test(s) */
      /*each populated bin gets a subgroup*/

      kinemageheader(commentstr);
      kinemagestuffer(janesviews);  /*070421*/
      kinemagestuffer(kinemageframe);  /*070414, 070421*/

      /*work through all possible bins and their suites in standard order*/
      if(L33out)
      {
         fprintf(fpout,"@group {33} recessiveon dimension=9 wrap=360 select animate off\n");
         binstuffout(1,3); /*bins in this group*/
      }
      if(L32out)
      {
         fprintf(fpout,"@group {32} recessiveon dimension=9 wrap=360 select animate off\n");
         binstuffout(4,6); /*bins in this group*/
      }
      if(L23out)
      {
         fprintf(fpout,"@group {23} recessiveon dimension=9 wrap=360 select animate off\n");
         binstuffout(7,9); /*bins in this group*/
      }
      if(L22out)
      {
         fprintf(fpout,"@group {22} recessiveon dimension=9 wrap=360 select animate off\n");
         binstuffout(10,12); /*bins in this group*/
      }
      if(Ltriageout)
      {
         fprintf(fpout,"@group {triaged} dominant dimension=9 wrap=360 select off\n");
         binstuffout(0,0); /*one bin (0th) in this group*/
      }
   }/*stuff for a kinemage*/
   if(Lstringout)
   {
      if(nout != 1) {fprintf(fpout," {%s}\n",suiteptr->ptID);}
      else          {fprintf(fpout,"\n");}  /*just final EOL*/
   }
}
/*___writeoutput()___________________________________________________________*/

/****binstuffout()************************************************************/
void binstuffout(int nbin, int mbin)
{
   int i=0,j=0,ncnt=0;
   char  ctrl[10]; /* 7 actual characters in curly braces*/
   char  extras[32];

   for(i=nbin; i<=mbin; i++)
   {/*loop over the three bins in this delta delta group*/
      if(binout[i])
      {/*ith bin has points*/
         fprintf(fpout,"@subgroup {%s} recessiveon \n",bin[i].binname);
         
         j=1; /*zeroth cluster is for bin outliers, if any, do them later */
         while(bin[i].clst[j].LOK)  /*named clusters j>0 ... 070429*/
         {/*loop over named clusters in this ith bin, NOT including outliers*/
            if(clusterout[i][j])
            {
               if(strcmp(bin[i].clst[j].status,"wannabe")==0) /*070429*/
               {
                  sprintf(extras," master= {wannabees}");
                  Lwannabeout = 1; /*redundant, already been set*/
               }
               else {extras[0] = '\0';}

               fprintf(fpout,
 "@balllist {%s %s} color= %s radius= 1 nohilite master= {data}%s\n"
 ,bin[i].binname,bin[i].clst[j].clustername,bin[i].clst[j].clustercolor,extras);
               sprintf(ctrl,"{%s %s}"
 ,bin[i].binname,bin[i].clst[j].clustername);
               ncnt = transferout(ctrl); /*gets just pts starting with ctrl*/

              if(j>0) /*obselete, but is a safety*/
              {/*belong to a named cluster  (av & labels initially on 070414)*/
               fprintf(fpout,
 "@ringlist {%s %s} color= %s radius= 10 width= 1 nobutton master= {avsigma}%s\n"
 ,bin[i].binname,bin[i].clst[j].clustername,bin[i].clst[j].clustercolor,extras);
               fprintf(fpout,
 "{%s %s} 180 %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f 180\n"
 ,bin[i].binname,bin[i].clst[j].clustername
 ,bin[i].clst[j].ang[1]
 ,bin[i].clst[j].ang[2]
 ,bin[i].clst[j].ang[3]
 ,bin[i].clst[j].ang[4]
 ,bin[i].clst[j].ang[5]
 ,bin[i].clst[j].ang[6]
 ,bin[i].clst[j].ang[7]);

               fprintf(fpout,
 "@labellist {%s %s} color= %s nobutton master= {labels}%s\n"
 ,bin[i].binname,bin[i].clst[j].clustername,bin[i].clst[j].clustercolor,extras);
               fprintf(fpout,
 "{%s %s} 180 %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f 180\n"
 ,bin[i].binname,bin[i].clst[j].clustername
 ,bin[i].clst[j].ang[1]
 ,bin[i].clst[j].ang[2]
 ,bin[i].clst[j].ang[3]
 ,bin[i].clst[j].ang[4]
 ,bin[i].clst[j].ang[5]
 ,bin[i].clst[j].ang[6]
 ,bin[i].clst[j].ang[7]);
             }/*belong to a named cluster*/
            }/*members in i,j cluster*/
            j++; /*increment cluster index for the while condition 070429*/
         }/*loop over named clusters in this ith bin*/
         j = 0; /*look for outliers... (set j==0 to reuse code from above)*/
         if(clusterout[i][j]) /*outliers in this bin*/
         {
               extras[0] = '\0';
               fprintf(fpout,
 "@balllist {%s %s} color= %s radius= 1 nohilite master= {data}%s\n"
 ,bin[i].binname,bin[i].clst[j].clustername,bin[i].clst[j].clustercolor,extras);
               sprintf(ctrl,"{%s %s}"
 ,bin[i].binname,bin[i].clst[j].clustername);
               ncnt = transferout(ctrl); /*gets just pts starting with ctrl*/
         }/*outliers in this bin*/
      }/*ith bin has points*/
   }/*loop over bins*/
}
/*___binstuffout()___________________________________________________________*/

/****transferout()************************************************************/
int  transferout(char* ctrl)
{
   /*cribbed from PKINCOUT/transferout() 070210 */
   /* for later flexibility, include some controls not now needed 070210 */

   int  LatEOF=0, LOK=0;
   int  j=0, iscrt=0, iexist=0;
   char scrts[256];

   rewindtextblock(&mainscratch);
   LatEOF = 0;

   while(!LatEOF)
   {/*scan through scratch block*/
      getonetextblockline(&mainscratch,temps);
      if(temps[0] != '\0')
      {/*entries still in scratch block*/
         if(   (temps[0]==ctrl[0])  /*check chars at beginning of ptID*/
             &&(temps[1]==ctrl[1])
             &&(temps[2]==ctrl[2])
             &&(temps[3]==ctrl[3])
             &&(temps[4]==ctrl[4])
             &&(temps[5]==ctrl[5])
             &&(temps[6]==ctrl[6])
             &&(temps[7]==ctrl[7]) )
         {
            LOK = 1; /*this point belongs in this list*/
            iscrt = 0; 
            if(Ltest){fprintf(stderr,"GOT %s\n",temps);}
         }
         else
         {/*unrecognized, can't use this point in this list*/
            LOK = 0;
         }
         if(LOK)
         {/*this point belongs in this list*/
            for( j=iscrt ; j<=255 ; j++) /*extract content from temps[] */
            {
                scrts[j-iscrt] = temps[j];
                if(temps[j] == '\0') break;
            }
            iexist++; /*increase counter*/
            /*suitename writeoutput has already written list header*/
            /* so do not write this here when iexist==1 */

            fprintf(fpout,"%s\n",scrts);

         /*this point belongs in this list*/}
      }/*entries still in scratch block*/
      else
      {/*reached end of scratch block*/
          LatEOF = 1;
      }
   }
   return(iexist);
}
/*___transferout()___________________________________________________________*/

/****kinemageheader()*********************************************************/
void kinemageheader(char* textstr)
{
   fprintf(fpout,"@text\n %s\n %s\n",version,textstr);
   fprintf(fpout,"@kinemage 1\n");
   fprintf(fpout,"@onewidth\n");
   if(Letatheta) /*070524*/
      {fprintf(fpout,"@dimension {theta} {delta-1} {epsilon-1} {zeta-1} {alpha} {beta} {gamma} {delta} {eta}\n");}
   else
      {fprintf(fpout,"@dimension {chi-1} {delta-1} {epsilon-1} {zeta-1} {alpha} {beta} {gamma} {delta} {chi}\n");}
   fprintf(fpout,"@dimminmax 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000 0.000 360.000\n");

   if(LTepsilon) {fprintf(fpout,"@pointmaster 'E' {epsilon bad}\n");}
   else if(LTdelta || LTdeltam) 
                         {fprintf(fpout,"@pointmaster 'D' {delta bad}\n");}
   else if(LTzeta || LTalpha || LTbeta || LTgamma)
                 {fprintf(fpout,"@pointmaster 'T' {various bad}\n");}
   if(Loutlier)  {fprintf(fpout,"@pointmaster 'O' {outliers}\n");}
   if(Lwannabeout)  {fprintf(fpout,"@master {wannabees}\n");}
}
/*___kinemageheader()________________________________________________________*/

/****kinemagestuffer()********************************************************/
void kinemagestuffer(char* kinemagestuff[])
{
   int  more = 1, j=0, Nth=0;
   char temps[256];

   while(more)
   {/*loop over all text for this call*/
      for(j=0; j<255; j++) /*save a space for safety '\0' */
      {/*load transfer string*/
         temps[j] =  kinemagestuff[Nth][j];
         if(temps[j] == '\r')
         {/*ASCII Carraige Return*/
            temps[j] = EOL; /*platform specific End-Of-Line character*/
         }
         if(temps[j] == '\0'){ break;} /*separate end of text piece*/
      }/*load transfer string*/
      /*check for END before writing to stdout*/
      if(temps[0] =='E' && temps[1] =='N' && temps[2] =='D'){more = 0;}
      temps[j+1] = '\0'; /*safety if runs off end of temps str*/
      if(more != 0) {fprintf(fpout,"%s",temps);}
      Nth++; /*local, not static, so will get reset ==0 on new call*/
   }/*loop over all text for this call*/
}
/*___kinemagestuffer()_______________________________________________________*/

/****clearbinout()************************************************************/
void clearbinout(void)
{
   int  i=0,j=0;

   for(i=1; i<MAXBINS; i++)
   {
      binout[i] = 0; /*set 1 if an entry is encountered*/
   }
}
/*___clearbinout()___________________________________________________________*/

/****clearclusterout()********************************************************/
void clearclusterout(void)
{                    /* clear storage arrays, re suitenout.h */
   int  i=0,j=0,k=0;

   for(i=0; i<MAXBINS; i++)
   {
      for(j=0; j<MAXCLST; j++)
      {
         clusterout[i][j] = 0; /*set 1 if an entry is encountered*/
         suitenesssum[i][j] = 0;
         for(k=0; k<12; k++) /*intervals: 10 10ths + extras at zero & 12*/
         {
            suitenesscnt[i][j][k] = 0;
         }
      }
   }
   suitenesssumall = 0;
   binnedsuitecountall = 0;
   triagecountall = 0;
   reportcountall = 0;
}
/*___clearclusterout()_______________________________________________________*/

/****usageout()***************************************************************/
void usageout(void) /* -h  -help */
{
  fprintf(stderr,"%s\n",version);
  fprintf(stderr,"suitename -flags <stdin >stdout\n");
  fprintf(stderr,"output flags: [ -report || -string || -kinemage ]\n");
  fprintf(stderr,"default:  -report -residuein \n");
  fprintf(stderr,"input flags: [ -residuein || -suitein  ]\n");
  fprintf(stderr,"flags: [ -residuein [ -pointIDfields # ] ] default#==%d\n"
                        ,NptIDfields);
  fprintf(stderr," OR ");
  fprintf(stderr,"flags: [ -suitein [ -anglefields # ] ]   default#==%d\n"
                        ,Nanglefields);
  fprintf(stderr,"\n");
  fprintf(stderr,"defaults: -residuein  -pointIDfields %d\n",NptIDfields);
  fprintf(stderr," as made by dangle \n");
  fprintf(stderr,
     "dangle \"alpha, beta, gamma, delta, epsilon, zeta\" in.pdb >out.dngl\n");
  fprintf(stderr,
     "label:model:chain:number:ins:type:alpha:beta:gamma:delta:epsilon:zeta\n");
  fprintf(stderr,"if the file has alternate conformations, then use both -pointIDfields # -altIDfield # to specify the number of pointID fields and which field contains the altID\n"); /* S.J. 01/07/2014*/
  fprintf(stderr,"use -altIDval <altID> to specify which alternate conformation to calculate suite for. By default calculated for alt A\n");
  fprintf(stderr,"-suitein presumes point records from a kinemage\n");
  fprintf(stderr,"{pointID} 7 or 9 anglefields \n");
  fprintf(stderr,
      "{ptID} [chi] deltam epsilon zeta alpha beta gamma delta [chi] \n");
  fprintf(stderr," Note that all other kinemage lines must be stripped off.\n");
  fprintf(stderr,"-thetaeta  kinemage labels theta,eta instead of chi-1,chi\n");
  fprintf(stderr,"Note dangle trick to make theta,...,eta suites directly\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"flag: -report [ -chart ]\n");
  fprintf(stderr," suites in order of input, suiteness summary at end\n");  
  fprintf(stderr,"( -chart : NO summary at end, for MolProbity multichart\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"flag: -string\n");
  fprintf(stderr," 3 character per suite string in order of input\n");  
  fprintf(stderr,"   20 per line, ptID of n*20th at end of line\n");  
  fprintf(stderr,"   flag: -nosequence\n");
  fprintf(stderr,"      only suite names, no Base sequence character\n");  
  fprintf(stderr,"   flag: -oneline\n");
  fprintf(stderr,"      string all one line, no point IDs\n");  
  fprintf(stderr,"   flag: -overlap\n");
  fprintf(stderr,"      20 per line: overlap 10 each line, 10 new per line\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"flag: -kinemage\n");
  fprintf(stderr," kinemage of clusters grouped by pucker,pucker ... \n");  
  fprintf(stderr," group {delta,delta},subgroup {gamma},list {cluster name}\n");
  fprintf(stderr,"\n");
  /*fprintf(stderr,"flag: -satellite\n"); DO NOT ADVERTISE THIS OPTION*/
  /*fprintf(stderr,"  use special general case satellite widths\n");*/
  fprintf(stderr," assigns to designated wannabe clusters, default: wannabe\n");
  /*fprintf(stderr,"flag: -wannabe\n");*/
  fprintf(stderr,"flag: -nowannabe   to not assign them\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"[ -power #.#] default# %4.2f multi-dimension distance calc\n"
                        ,power);
  fprintf(stderr,"[ -test ] dump cluster centers, halfwidths,... to stderr\n");
  fprintf(stderr,"cluster averages version: %s\n",clusteraveragesversion);
  fprintf(stderr,"cluster half-widths version: %s\n",clusterhalfwidthsversion);
  fprintf(stderr,"axes limits version: %s\n",axeslimitsversion);
  fprintf(stderr,"suitename is not paying attention to chains...\n");
  fprintf(stderr,"%s\n",version);
}
/*___usageout()______________________________________________________________*/


