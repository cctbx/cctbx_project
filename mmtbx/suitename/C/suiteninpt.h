/*                        suiteninpt.h                            */

/*only suiteninpt.c defines SUITENINPT */

#ifndef SUITENINPT
extern char *NAListA,*NAListG,*NAListC,*NAListU,*NAListT;
extern float angle[10];
extern char  ptID[10][32];
extern char  basechr[2];
#else
 char *NAListA = ":ADE:  A:A  : Ar:ATP:ADP:AMP:T6A:1MA:RIA:  I:I  :";
 char *NAListG = ":GUA:  G:G  : Gr:GTP:GDP:GMP:GSP:1MG:2MG:M2G:OMG:\
 YG: 7MG:YG :";
 char *NAListC = ":CYT:  C:C  : Cr:CTP:CDP:CMP:5MC:OMC:";
 char *NAListU = ":URA:URI:  U: Ur:U  :UTP:UDP:UMP:5MU:H2U:PSU:4SU:";
 char *NAListT = ":THY:  T:T  : Tr:TTP:TDP:TMP:";
 float angle[10];
 char  ptID[10][32];
 char  basechr[2];
#endif

int   getresidue(void);
int   getsuite(void);   /*070211 alternative to getresidue()*/
int   loadsuite(void);
int   movenewtooldresidue(void);
void  printresidue(char*);
int   clearnewresidue(void);
int   loadnewresidue(void);
int   interpretdanglerecord(void);   /*070525*/
int   interpretresiduerecord(void);
int   interpretsuiterecord(void);   /*070211*/
int   Getoneinputline(void);
float floatfromstr(char[256]); /*from MAGEUTIL.c*/

