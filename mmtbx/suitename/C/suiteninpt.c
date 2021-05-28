/*                           suiteninpt.c                                    */

#include "suitename.h"
#define SUITENINPT
#include "suiteninpt.h"
#undef  SUITENINPT

#include "suiteninit.h"
#include "suitenscrt.h"
#include "suitenutil.h" 

/****getresidue()*************************************************************/ 
int  getresidue(void)
{
   int LOK=0;
   int j=0;

   while(LOK == 0 && !LatEOF) /*==1 for valid residue line*/
   {
      LOK = Getoneinputline(); /* ==0 for comment or too-short line */
      /*return global texts[] and itext */
      if(LOK && (itext<=1 || texts[0]=='#')) {LOK = 0;} /* # dangle comment*/
   }
   if(LOK)
   {
      LOK = interpretresiduerecord();
      /* ptID: : :alpha:beta:gamma:delta:epsilon:zeta */
      /*loads ptID[NptIDfields], angle[6] */ 
      if(LOK)
      {
         movenewtooldresidue();
         clearnewresidue();
         LOK = loadnewresidue();
      }
   }
   return(LOK);
}
/*___getresidue()____________________________________________________________*/

/****getsuite()***************************************************************/
int  getsuite(void)
{
   int LOK=0;

   while(LOK == 0 && !LatEOF) /*==1 for valid suite line*/
   {
      LOK = Getoneinputline(); /* ==0 for comment or too-short line */
      /*return global texts[] and itext */
      if(LOK && (itext<=1 || texts[0]=='#')) {LOK = 0;} /* # dangle comment*/
   }
   if(LOK) /*body is 7D or 9D kinemage with {ptID}*/
   {
      LOK = interpretsuiterecord();
      /* {ptID} [chi],deltam,epsilon,zeta,alpha,beta,gamma,delta,[chi] */
      /* and loads suiteptr so do not need to loadsuite()  */
   }
   return(LOK);
}
/*___getsuite()______________________________________________________________*/

/****loadsuite()**************************************************************/
int  loadsuite(void)
{
    int  j=0,k=0,n=0;

    n=0;
    for(j=1; j<=NptIDfields; j++) 
    {
       for(k=0; k<32; k++)
       {
          if(newresidueptr->ptID[j][k] == '\0') {break;}
          suiteptr->ptID[n++] = newresidueptr->ptID[j][k];
       }
    }
    suiteptr->ptID[n] = '\0'; /*terminate ptID string*/
    suiteptr->basechr[0] = newresidueptr->basechr[0]; /*070412*/
    suiteptr->basechr[1] = '\0';
    suiteptr->ang[1] = suiteptr->deltam  = oldresidueptr->delta;
    suiteptr->ang[2] = suiteptr->epsilon = oldresidueptr->epsilon;
    suiteptr->ang[3] = suiteptr->zeta    = oldresidueptr->zeta;
    suiteptr->ang[4] = suiteptr->alpha   = newresidueptr->alpha;
    suiteptr->ang[5] = suiteptr->beta    = newresidueptr->beta;
    suiteptr->ang[6] = suiteptr->gamma   = newresidueptr->gamma;
    suiteptr->ang[7] = suiteptr->delta   = newresidueptr->delta;

    return(1);
}
/*___loadsuite()_____________________________________________________________*/

/****movenewtooldresidue()****************************************************/
int  movenewtooldresidue(void)
{
    int  j=0;

    for(j=1; j<=9; j++)
    {
       strcpy(oldresidueptr->ptID[j],newresidueptr->ptID[j]);
    }
    //oldresidueptr->basechr[0],newresidueptr->basechr[0]; /*070412*/
    oldresidueptr->basechr[0] = '\0';
    oldresidueptr->alpha   = newresidueptr->alpha;
    oldresidueptr->beta    = newresidueptr->beta;
    oldresidueptr->gamma   = newresidueptr->gamma;
    oldresidueptr->delta   = newresidueptr->delta;
    oldresidueptr->epsilon = newresidueptr->epsilon;
    oldresidueptr->zeta    = newresidueptr->zeta;
    return(1);
}
/*___movenewtooldresidue()___________________________________________________*/

/****printresidue()***********************************************************/
void printresidue(char* type)
{
    int  j=0;
    struct residuestruct* theresidueptr;

    if(strcmp(type,"old")==0) {theresidueptr = oldresidueptr;}
    else {theresidueptr = newresidueptr;}
    fprintf(stderr,"%s ",type);
    for(j=1; j<=NptIDfields; j++)
    {
       fprintf(stderr,"%s",theresidueptr->ptID[j]);
    }
    fprintf(stderr," %7.2f, %7.2f, %7.2f, %7.2f, %7.2f, %7.2f\n"
       ,theresidueptr->alpha
       ,theresidueptr->beta
       ,theresidueptr->gamma
       ,theresidueptr->delta
       ,theresidueptr->epsilon
       ,theresidueptr->zeta);
}
/*___printresidue()__________________________________________________________*/

/****clearnewresidue()********************************************************/
int  clearnewresidue(void)
{
    int  j=0;

    for(j=0; j<=9; j++)
    {
       newresidueptr->ptID[j][0] = '\0';
    }
    newresidueptr->basechr[0] = '\0'; /*070412*/
    newresidueptr->basechr[1] = '\0';
    newresidueptr->alpha   = 9999.99;
    newresidueptr->beta    = 9999.99;
    newresidueptr->gamma   = 9999.99;
    newresidueptr->delta   = 9999.99;
    newresidueptr->epsilon = 9999.99;
    newresidueptr->zeta    = 9999.99;
    return(1);
}
/*___clearnewresidue()_______________________________________________________*/

/****loadnewresidue()*********************************************************/
int  loadnewresidue(void)
{
    int  j=0;

    for(j=1; j<=NptIDfields; j++)
    {
       strcpy(newresidueptr->ptID[j],ptID[j]);
    }
    newresidueptr->basechr[0] = basechr[0]; /*070412*/
    newresidueptr->basechr[1] = '\0';
    for(j=1; j<=9; j++)
    {  /*work with angles 0 to 360, dang calc angles -180 to +180 */
       if(angle[j] < 0)
       {
          angle[j] = angle[j] + 360;
       }
    }
    newresidueptr->alpha   = angle[1];
    newresidueptr->beta    = angle[2];
    newresidueptr->gamma   = angle[3];
    newresidueptr->delta   = angle[4];
    newresidueptr->epsilon = angle[5];
    newresidueptr->zeta    = angle[6];
    return(1);
}
/*___loadnewresidue()________________________________________________________*/

/****interpretdanglerecord()**************************************************/
int  interpretdanglerecord(void) 
{
   int  i=0,j=0,k=0,n=0,ns=0,nan=0;
   char numstr[256];

   /*ptID: :  : : :...:alpha:beta:gamma:delta:epsilon:zeta*/

   n = 1; /*count fields from 1 */
   ns = 0; /*first j of current field*/
   k = 0; /*counter within a field*/

   for(j=0 ; j<= itext; j++)
   {/*loop through inputed line*/
      if(texts[j] == ':' || j==itext)
      {/* : */
         {/*close nth field*/
            if(n <= NptIDfields)
            {/*ID fields*/
               {
                  /*k not decremented, k is index of the : character*/
                  for(i=0; i<=k; i++)
                  {
                     ptID[n][i] = texts[ns+i];
                  }
                  if(n == NptIDfields)
                  {
                     ptID[n][k] = '\0'; /*overwrite last : character*/
                              /*so last char of last field is a Base name char*/
                      /*and this last field is the full 3char Base name*/
                     if(strlen(ptID[n]) == 3)
                     { /*interpret as Base name --> 1char   070412 */
                        if     (strstr(NAListA, ptID[n])) {basechr[0] = 'A';}
                        else if(strstr(NAListG, ptID[n])) {basechr[0] = 'G';}
                        else if(strstr(NAListC, ptID[n])) {basechr[0] = 'C';}
                        else if(strstr(NAListU, ptID[n])) {basechr[0] = 'U';}
                        else if(strstr(NAListT, ptID[n])) {basechr[0] = 'T';}
                        else {basechr[0] = 'Y';}
                     }
                     else {basechr[0] = 'Z';}
                     basechr[1] = '\0';
                  }
                  else {ptID[n][k+1] = '\0';} /* not overwrite : character*/
               }
               /*reset for another field*/
               {
                  k=0;
                  nan=0;
                  ns = j+1;
                  n++;
               }
            }/*ID fields*/
            else
            {/*angle coord fields*/
               if(nan) {angle[n-NptIDfields] = 9999.99;}
               else {angle[n-NptIDfields] = floatfromstr(numstr);}
               k=0;
               nan=0;
               ns = j+1;
               n++;
            }
         }/*close nth field*/
      }/* : */
      else
      {/* input character , only actually use numstr for angle fields*/
         numstr[k++] = texts[j];
         numstr[k] = '\0'; /*keep string terminated*/
         if(texts[j] == '?') {nan = 1;} /*NOT a number*/
      }
   }/*loop through inputed line*/
   return(1);
}
/*___interpretdanglerecord()________________________________________________*/

/****interpretresiduerecord()*************************************************/
int  interpretresiduerecord(void)
{
   int  i=0,j=0,k=0,n=0,ns=0,nan=0;
   char numstr[256];

   /*ptID: :  : : :...:alpha:beta:gamma:delta:epsilon:zeta*/

   n = 1; /*count fields from 1 */
   ns = 0; /*first j of current field*/
   k = 0; /*counter within a field*/

   for(j=0 ; j<= itext; j++)
   {/*loop through inputed line*/
      if(texts[j] == ':' || j==itext)
      {/* : */
         {/*close nth field*/
            if(n <= NptIDfields)
            {/*ID fields*/
               {
                  /*k not decremented, k is index of the : character*/
                  for(i=0; i<=k; i++)
                  {
                     ptID[n][i] = texts[ns+i];
                  }
                  if(n == NptIDfields)
                  {
                     ptID[n][k] = '\0'; /*overwrite last : character*/
                              /*so last char of last field is a Base name char*/
                      /*and this last field is the full 3char Base name*/
                      if(!strcmp(ptID[n]," DA") || !strcmp(ptID[n]," DG") || !strcmp(ptID[n]," DC") || !strcmp(ptID[n]," DT")) // SJ - 09/17/2014 to ignore DNA residues
                      {
                          return(0);
                      }
                     if(strlen(ptID[n]) == 3)
                     { /*interpret as Base name --> 1char   070412 */
                        if     (strstr(NAListA, ptID[n])) {basechr[0] = 'A';}
                        else if(strstr(NAListG, ptID[n])) {basechr[0] = 'G';}
                        else if(strstr(NAListC, ptID[n])) {basechr[0] = 'C';}
                        else if(strstr(NAListU, ptID[n])) {basechr[0] = 'U';}
                        else if(strstr(NAListT, ptID[n])) {basechr[0] = 'T';}
                        else {basechr[0] = 'Y';}
                     }
                     else {basechr[0] = 'Z';}
                     basechr[1] = '\0';
                  }
                  else {ptID[n][k+1] = '\0';} /* not overwrite : character*/
               }
               /*reset for another field*/
               {
                  k=0;
                  nan=0;
                  ns = j+1;
                  n++;
               }
            }/*ID fields*/
            else
            {/*angle coord fields*/
               if(nan) {angle[n-NptIDfields] = 9999.99;}
               else {angle[n-NptIDfields] = floatfromstr(numstr);}
               k=0;
               nan=0;
               ns = j+1;
               n++;
            }
         }/*close nth field*/
      }/* : */
      else
      {/* input character , only actually use numstr for angle fields*/
         numstr[k++] = texts[j];
         numstr[k] = '\0'; /*keep string terminated*/
         if(texts[j] == '?') {nan = 1;} /*NOT a number*/
      }
   }/*loop through inputed line*/
   if(LaltIDfield) /* if the altIDfield specified - S.J. 01/07/2014*/
   { 
        if(!LaltID) /* if the user has not specified which alt to calculate for, calculate for altA*/
        	altID[0]='A';
	if(!CompArgStr(ptID[altIDfield],altID,1) && !CompArgStr(ptID[altIDfield]," ",1)) 
         		return(0);
   }
   return(1);
}
/*___interpretresiduerecord()________________________________________________*/

/****interpretsuiterecord()*************************************************/
int  interpretsuiterecord(void)
{  /*hack unprotected way to read well-formated kinemage pt records*/

   int  i=0,j=0,k=0,n=0,nan=0; /*nan 070525*/
   char numstr[256];
   int LptID = 0, iptID = 0;
   char ptIDstr[256];
   char basechr[2];  /*070412*/

   for(j=0; j<256; j++){ptIDstr[j] = '\0';}
   basechr[0] = 'X'; /*070412*/
   basechr[1] = '\0';

   /* {ptID} deltam,epsilon,zeta,alpha,beta,gamma,delta */

   for(j=0 ; j<= itext; j++)
   {/*loop through inputed line*/
      if(texts[j] == '{') /* }  balance curly braces */
      {/*pointID */
         LptID = 1;
         iptID = 0;
      }
/*{*/ else if(texts[j] == '}') /* balance curly braces */
      {
         ptIDstr[iptID] = '\0'; /*end ptIDstr*/
         LptID = 0;
         basechr[0] = ptIDstr[iptID -1]; /*last char of kin ptID is base char*/
      }
      else if(LptID)
      {
         ptIDstr[iptID++] = texts[j];
      }
      else if( (k > 0 && texts[j] == ' ') || texts[j] == ',' || j==itext)
      {/*white space ends a field*/
         numstr[k] = '\0';
         if(nan) {angle[n] = 9999.99;} /*070525*/
         else {angle[n] = floatfromstr(numstr);} /*n starts at 0*/
         if(angle[n] < 0) {angle[n] = 360 + angle[n];} /*scope 0 to 360 */
         if(texts[j] == ' ' || texts[j] == ',') /*interior record white space*/
         {/*expect another angle*/
            k=0;
            n++;
            nan = 0; /*reset not-a-number flag  070525*/
         }
      }
      else
      {/* input character */
         if(texts[j] != ' ') /*ignore leading blanks*/
         if(texts[j] == '_' || texts[j] == '?' ){nan = 1;} /*070525*/
         {numstr[k++] = texts[j];}
      }
   }/*loop through inputed line*/
   /*presume suite kinemage record was intact and good...*/

    strcpy(suiteptr->ptID,ptIDstr);

    if(Nanglefields == 9) /*rearranged 070524 to preserve chim,chi (eta,theta)*/
    {
       i=1;
       suiteptr->ang[0] = suiteptr->chim = angle[0];
       suiteptr->ang[8] = suiteptr->chi  = angle[7+i];
    }
    else
    {
       i=0;
       suiteptr->ang[0] = suiteptr->chim = 180.0;
       suiteptr->ang[8] = suiteptr->chi  = 180.0;
    }
    suiteptr->ang[1] = suiteptr->deltam  = angle[0+i];
    suiteptr->ang[2] = suiteptr->epsilon = angle[1+i];
    suiteptr->ang[3] = suiteptr->zeta    = angle[2+i];
    suiteptr->ang[4] = suiteptr->alpha   = angle[3+i];
    suiteptr->ang[5] = suiteptr->beta    = angle[4+i];
    suiteptr->ang[6] = suiteptr->gamma   = angle[5+i];
    suiteptr->ang[7] = suiteptr->delta   = angle[6+i];
    /*now DO NOT need to loadsuite() */
    return(1);
}
/*___interpretsuiterecord()________________________________________________*/

/****Getoneinputline()********************************************************/
int  Getoneinputline(void) /*cribbed from PKININPT.c*/
{
   int   c=0,ireturn=0; /*NOTE: normal return is from first part of while loop*/

   for(itext=0; itext<256; itext++) {texts[itext] = '\0';}
   /*this will avoid left over garbage at ends of short records*/
   itext = 0; /*itext is a global, numbers inputed characters on each line*/
   texts[0] = '\0';/*initially a NULL string*/

   if(!Lhitend)
   {/*bypass for EOF on previous line*/
      while ((c = fgetc(fpin)) != EOF)
      {
         if (c == '\r' || c == '\n')
         {
            texts[itext] = '\0';

            LatEOF = 0;    /* not at EOF */
            strcat(texts,EOLO);
            if(itext < 2) {ireturn = 0;} /*ignore empty or very short lines*/
            else ireturn = 1;
            return(ireturn);  /* break out of this subroutine */
         }
         else
         {
            texts[itext] = c;
            if (itext > 254)
            {
               texts[itext+1] = '\0';
               /*make array of characters into properly ended C string*/
               itext = 0;

               /* Warn the reader about line too long */
               fprintf(stderr,"%s\n",texts);
               fprintf(stderr,"above line is too long");
            }
            else itext++;
         }
      }/* loop while file not at EOF */
      /* EOF : finish out anything in text buffer */
      /*if line is long enough to be coords*/
      if(itext > 50) /*enough characters to have coords even on this last line*/
      {
         ireturn = 1; /*line has content*/
         LatEOF = 0; /*pretend not EOF, let next pass catch it*/
         Lhitend = 1; /*but store info that EOF was hit*/
         texts[itext] = '\0';
      }
      else
      {
         ireturn = 0; /*not a valid line of any sort*/
         LatEOF = 1;  /* EOF */
         Lhitend = 1;
      }
   }/*bypass for EOF on previous line*/
   else
   {
      LatEOF = 1;  /* EOF from previous call now registered*/
   }
   return(ireturn);
}
/*___Getoneinputline()_______________________________________________________*/

/*****floatfromstr()*********************************************************/
float    floatfromstr(char ins[256])
{
        int            Lstart,n,m,OK;
        char         s[256];
        float        freturn;

   /*970703 now only handles C-strings*/

        freturn = 0.0;
        Lstart = 1;
        n = 0;
        m = 0;
        OK = 1;
        if(ins[0] != '\0')
        {/* not zero length C string */
            while(OK)
            {
                /* start at position 0 of c str */
                if(ins[n] == ' ' && Lstart)  ;
                else if(ins[n] == '0' || ins[n] == '1' || ins[n] == '2' ||
                        ins[n] == '3' || ins[n] == '4' || ins[n] == '5' ||
                        ins[n] == '6' || ins[n] == '7' || ins[n] == '8' ||
                        ins[n] == '9' || ins[n] == '.' || ins[n] == '-'   )
                {
                    s[m] = ins[n];
                    m++;
                    Lstart = 0;
                }
                else
                {
                    s[m] = '\0';
                    OK = 0;
                }
                n++;
            }

/*c*/       sscanf(s,"%f",&freturn);
        }
        return(freturn);
}
/*___floatfromstr()_________________________________________________________*/



