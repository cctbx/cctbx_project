/*               suitenout.h                  */

#ifdef  SUITENOUT
#undef  SUITENOUT
#define SUITENOUT
#else
#define SUITENOUT extern
#endif

#define DELTAM   1  /*Ltriage reason flags 070628*/
#define EPSILONM 2
#define ZETAM    3
#define ALPHA    4
#define BETA     5
#define GAMMA    6
#define DELTA    7

/*diagnostic flags*/
SUITENOUT int Ltriage; /*reset for each suite, set to reason, e.g. DELTAM*/
SUITENOUT int Liswannabe;
SUITENOUT int Lcomment; /*for 7D distance forced outlier or closest assignment*/
SUITENOUT char commentstr[16]; /*and Lcomment  070628*/
SUITENOUT int Loutlier; /*accummulates*/
SUITENOUT int LTdeltam,LTdelta,LTepsilon,LTzeta,LTalpha,LTbeta,LTgamma;
              /* LT...flags accummulate for pointmaster definition 070628*/
SUITENOUT int L33out,L32out,L23out,L22out,Ltriageout;           /*groups*/
SUITENOUT int binout[MAXBINS]; /*incl 1--12 named bins */       /*subgroups*/
SUITENOUT int clusterout[MAXBINS][MAXCLST];                     /*lists*/
SUITENOUT double suitenesssum[MAXBINS][MAXCLST],     suitenesssumall;
SUITENOUT int    suitenesscnt[MAXBINS][MAXCLST][12]; 
            /*3rd index for suiteness intervals: 10 at 10ths + 2 extras   */
            /*   extra at zero counts valid suites with suiteness == 0    */
            /*        e.g. 4D distance < 1 but 7D distance > 1            */
            /*   extra at 11   counts valid suites but triaged or outlier */
SUITENOUT int    binnedsuitecountall, reportcountall;
SUITENOUT int    triagecountall;  /*070328*/
SUITENOUT int    Lwannabeout;  /*wannabe in output flag 070429*/

SUITENOUT char temps[1024];

static int nout=1; /*Lsringout residue/suite counter*/

/*prototypes*/
void writeoutput(void);
void binstuffout(int, int);
int  transferout(char*);
void kinemageheader(char*); /*070328 char* */
void kinemagestuffer(char* kinemagestuff[]);  /*070414, 070421*/
void  writesuite(int, int, char*, float, float, char*, char*);
void  clearbinout(void);
void  clearclusterout(void);
void  usageout(void);
void  suitenessaverage(int);

