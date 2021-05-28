/*                       suitenutil.h                     */

#ifdef  SUITENUTIL
#undef  SUITENUTIL
#define SUITENUTIL
#else
#define SUITENUTIL extern
#endif

/* prototypes: */

void  vector7ab(float*, float*, float*);
float hyperellipsoiddist(int,int,int,float*);
void  resetcoordw(float*,int,int,int);
float dotproduct(float*, float*,int);
int   confirmsuite(void);
int   CompArgStr(char*, char*, int);

