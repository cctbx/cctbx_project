/*                        suiteninit.h                            */
/* suitename.h MAXBINS 14  1--12 named bins, bin 0 triaged, bin 13 incomplete*/
/* suitename.h MAXCLST 16  practical, <observed, limit of clusters in a bin*/
         /*clst indexed from 1: bin 33p had 10, index to 11, 070428, +1 070919*/
#define MAXSAT  12  /*maximun number of all satellite clusters, 9 as of 070429*/

#ifdef  SUITENINIT
#undef  DECLARATIONS
#define DEFINITIONS  /*storage and initializations done for suiteninit.c*/
#undef  SUITENINIT
#define SUITENINIT   /*no prefix, so do definition */
#else
#undef  DEFINITIONS
#define DECLARATIONS /*just declare existance for everybody*/
#undef  SUITENINIT
#define SUITENINIT extern  /*extern prefix so just declaration*/
#endif

#ifdef DEFINITIONS  
const struct satinfo
   {  
      char name[3];
      /*satellite repacement widths (when non-zero) */
      float satw[9]; /*satellite width when point between sat and dom clst*/
      float domw[9]; /*dominant clst width for the inbetween case*/
       /*count from 1 for 7 angle widths, use 9 for consistency*/
      char doma[3]; /*dominant name -- bookkeeping 070506*/
   }satinfo[MAXSAT] =
    /*satellite widths -- inbetween -- dom widths*/
    /* sat  0 1 s2 s3 4 s5 6 7 8 0 1 d2 d3 4 d5 6 7 8   dom*/
      {"1m",0,0, 0, 0,0,32,0,0,0,0,0, 0, 0,0,64,0,0,0 ,"1a",
       "1L",0,0,18, 0,0,18,0,0,0,0,0,70, 0,0,70,0,0,0 ,"1a",
       "&a",0,0,20,20,0, 0,0,0,0,0,0,60,60,0, 0,0,0,0 ,"1a",
       "1f",0,0, 0, 0,0,47,0,0,0,0,0, 0, 0,0,65,0,0,0 ,"1c",
       "1[",0,0, 0, 0,0,34,0,0,0,0,0, 0, 0,0,56,0,0,0 ,"1b",
       "4a",0,0,40,40,0, 0,0,0,0,0,0,50,50,0, 0,0,0,0 ,"0a",
       "#a",0,0,26,26,0, 0,0,0,0,0,0,36,36,0, 0,0,0,0 ,"0a", /*070506*/
       "0i",0,0, 0, 0,0,60,0,0,0,0,0, 0, 0,0,60,0,0,0 ,"6n",
       "6j",0,0, 0, 0,0,60,0,0,0,0,0, 0, 0,0,60,0,0,0 ,"6n"};
#endif /*DEFINITIONS*/

#ifdef DECLARATIONS
extern struct satinfo
   {  
      char name[3];
      /*satellite repacement widths (when non-zero) */
      float satw[9]; /*satellite width when point between sat and dom clst*/
      float domw[9]; /*dominant clst width for the inbetween case*/
       /*count from 1 for 7 angle widths, use 9 for consistency*/
      char doma[3]; /*dominant name -- bookkeeping 070506*/
   } satinfo[MAXSAT];
#endif /*DECLARATIONS*/

                  /*bin and clst declared for all files*/
SUITENINIT char   clusteraveragesversion[7];
SUITENINIT struct clusterdefinition 
   {
      char  clustername[3]; /*2char designation of consensus names*/
      int   LOK; /*logical 0 or 1 for actual cluster */
      char  status[8]; /*certain OR wannabe, triaged,outlier,nothing,incompl*/
      char  clustercolor[12]; /*std kinemage color names, 12 char sufficient*/
      char  domsatness[4];  /*dom.inant, sat.ellite, ord.inary, out,tri,inc*/ 
      const struct satelliteinfo* satptr;/*NULL except for satellite clusters 070428*/
      float ang[9];        /*7 count from 1 with chi-1 as 0 and chi as 8*/
   } clst[MAXCLST];

SUITENINIT struct bindefinition 
   {
      char   binname[5]; /*4char designation of bin, trig,33 p, ...  */
      struct clusterdefinition clst[MAXCLST];
   }  bin[MAXBINS];

#ifdef DEFINITIONS  
                     /*bin and clst defined for suiteninit.c*/
char   clusteraveragesversion[7] = {"070506"};

              /* NB: uses while(LOK): need 1,... contiguous clusters */
              /* ?depend on compiler to zero trailing undefined arrays?*/
              /* otherwise need to define trailing dummy LOK=0 cluster*/
              /*BEWARE: of MAXCLST when adding clusters !! (suitename.h)*/
/*angles: chi-1, delta ,epsilon,  zeta ,  alpha,  beta , gamma , delta ,chi */
struct bindefinition bin[MAXBINS] = 
    {
      {/*bin  0 */
         "trig",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "triaged", "white      ", "tri", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  1 */
         "33 p",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "1a", 1, "certain", "yellowtint ", "dom", NULL,
       180, 081.495, 212.250, 288.831, 294.967, 173.990, 053.550, 081.035 ,180},
           {/*  2 */ "1m", 1, "certain", "blue       ", "sat", NULL,
       180, 083.513, 218.120, 291.593, 292.247, 222.300, 058.067, 086.093 ,180},
           {/*  3 */ "1L", 1, "certain", "green      ", "sat", NULL,
       180, 085.664, 245.014, 268.257, 303.879, 138.164, 061.950, 079.457 ,180},
           {/*  4 */ "&a", 1, "certain", "cyan       ", "sat", NULL,
       180, 082.112, 190.682, 264.945, 295.967, 181.839, 051.455, 081.512 ,180},
           {/*  5 */ "7a", 1, "certain", "pink       ", "ord", NULL,
       180, 083.414, 217.400, 222.006, 302.856, 160.719, 049.097, 082.444 ,180},
           {/*  6 */ "3a", 1, "certain", "magenta    ", "ord", NULL,
       180, 085.072, 216.324, 173.276, 289.320, 164.132, 045.876, 084.956 ,180},
           {/*  7 */ "9a", 1, "certain", "hotpink    ", "ord", NULL,
       180, 083.179, 210.347, 121.474, 288.568, 157.268, 049.347, 081.047 ,180},
           {/*  8 */ "1g", 1, "certain", "sea        ", "ord", NULL,
       180, 080.888, 218.636, 290.735, 167.447, 159.565, 051.326, 085.213 ,180},
           {/*  9 */ "7d", 1, "certain", "purple     ", "ord", NULL,
       180, 083.856, 238.750, 256.875, 069.562, 170.200, 052.800, 085.287 ,180},
           {/* 10 */ "3d", 1, "certain", "peach      ", "ord", NULL,
       180, 085.295, 244.085, 203.815, 065.880, 181.130, 054.680, 086.035 ,180},
           {/* 11 */ "5d", 1, "certain", "yellow     ", "ord", NULL,
       180, 079.671, 202.471, 063.064, 068.164, 143.450, 049.664, 082.757 ,180},
           {/* 12 */ "3g", 1, "wannabe", "gray       ", "ord", NULL,
       180, 084.000, 195.000, 146.000, 170.000, 170.000, 052.000, 084.000 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  2 */
         "33 t",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "1e", 1, "certain", "red        ", "ord", NULL,
       180, 080.514, 200.545, 280.510, 249.314, 082.662, 167.890, 085.507 ,180},
           {/*  2 */ "1c", 1, "certain", "gold       ", "dom", NULL,
       180, 080.223, 196.591, 291.299, 153.060, 194.379, 179.061, 083.648 ,180},
           {/*  3 */ "1f", 1, "certain", "lime       ", "sat", NULL,
       180, 081.395, 203.030, 294.445, 172.195, 138.540, 175.565, 084.470 ,180},
           {/*  4 */ "5j", 1, "certain", "sky        ", "ord", NULL,
       180, 087.417, 223.558, 080.175, 066.667, 109.150, 176.475, 083.833 ,180},
           {/*  5 */ "5n", 1, "wannabe", "gray       ", "ord", NULL,
       180, 086.055, 246.502, 100.392, 073.595, 213.752, 183.395, 085.483 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  3 */
         "33 m",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  4 */
         "32 p",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "1b", 1, "certain", "cyan       ", "dom", NULL,
       180, 084.215, 215.014, 288.672, 300.420, 177.476, 058.307, 144.841 ,180},
           {/*  2 */ "1[", 1, "certain", "pink       ", "sat", NULL,
       180, 082.731, 220.463, 288.665, 296.983, 221.654, 054.213, 143.771 ,180},
           {/*  3 */ "3b", 1, "certain", "lilac      ", "ord", NULL,
       180, 084.700, 226.400, 168.336, 292.771, 177.629, 048.629, 147.950 ,180},
           {/*  4 */ "1z", 1, "certain", "peach      ", "ord", NULL,
       180, 083.358, 206.042, 277.567, 195.700, 161.600, 050.750, 145.258 ,180},
           {/*  5 */ "5z", 1, "certain", "purple     ", "ord", NULL,
       180, 082.614, 206.440, 052.524, 163.669, 148.421, 050.176, 147.590 ,180},
           {/*  6 */ "7p", 1, "certain", "sea        ", "ord", NULL,
       180, 084.285, 236.600, 220.400, 068.300, 200.122, 053.693, 145.730 ,180},
           {/*  7 */ "5p", 1, "wannabe", "gray       ", "ord", NULL,
       180, 084.457, 213.286, 069.086, 075.500, 156.671, 057.486, 147.686 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  5 */
         "32 t",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "1t", 1, "certain", "red        ", "ord", NULL,
       180, 081.200, 199.243, 288.986, 180.286, 194.743, 178.200, 147.386 ,180},
           {/*  2 */ "5q", 1, "certain", "yellow     ", "ord", NULL,
       180, 082.133, 204.933, 069.483, 063.417, 115.233, 176.283, 145.733 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  6 */
         "32 m",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "1o", 1, "certain", "sky        ", "ord", NULL,
       180, 083.977, 216.508, 287.192, 297.254, 225.154, 293.738, 150.677 ,180},
           {/*  2 */ "7r", 1, "certain", "lilactint  ", "ord", NULL,
       180, 084.606, 232.856, 248.125, 063.269, 181.975, 295.744, 149.744 ,180},
           {/*  3 */ "5r", 1, "wannabe", "gray       ", "ord", NULL,
       180, 083.000, 196.900, 065.350, 060.150, 138.425, 292.550, 154.275 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  7 */
         "23 p",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "2a", 1, "certain", "cyan       ", "ord", NULL,
       180, 145.399, 260.339, 288.756, 288.444, 192.733, 053.097, 084.067 ,180},
           {/*  2 */ "4a", 1, "certain", "yellow     ", "sat", NULL,
       180, 146.275, 259.783, 169.958, 298.450, 169.583, 050.908, 083.967 ,180},
           {/*  3 */ "0a", 1, "certain", "green      ", "dom", NULL,
       180, 149.286, 223.159, 139.421, 284.559, 158.107, 047.900, 084.424 ,180},
           {/*  4 */ "#a", 1, "certain", "hotpink    ", "sat", NULL,
       180, 148.006, 191.944, 146.231, 289.288, 150.781, 042.419, 084.956 ,180},
           {/*  5 */ "4g", 1, "certain", "greentint  ", "ord", NULL,
       180, 148.028, 256.922, 165.194, 204.961, 165.194, 049.383, 082.983 ,180},
           {/*  6 */ "6g", 1, "certain", "gold       ", "ord", NULL,
       180, 145.337, 262.869, 079.588, 203.863, 189.688, 058.000, 084.900 ,180},
           {/*  7 */ "8d", 1, "certain", "red        ", "ord", NULL,
       180, 148.992, 270.596, 240.892, 062.225, 176.271, 053.600, 087.262 ,180},
           {/*  8 */ "4d", 1, "certain", "sky        ", "ord", NULL,
       180, 149.822, 249.956, 187.678, 080.433, 198.133, 061.000, 089.378 ,180},
           {/*  9 */ "6d", 1, "certain", "orange     ", "ord", NULL,
       180, 146.922, 241.222, 088.894, 059.344, 160.683, 052.333, 083.417 ,180},
           {/* 10 */ "2g", 1, "wannabe", "gray       ", "ord", NULL,
       180, 141.900, 258.383, 286.517, 178.267, 165.217, 048.350, 084.783 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  8 */
         "23 t",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "2h", 1, "certain", "sea        ", "ord", NULL,
       180, 147.782, 260.712, 290.424, 296.200, 177.282, 175.594, 086.565 ,180},
           {/*  2 */ "4n", 1, "certain", "peach      ", "ord", NULL,
       180, 143.722, 227.256, 203.789, 073.856, 216.733, 194.444, 080.911 ,180},
           {/*  3 */ "0i", 1, "certain", "lilactint  ", "sat", NULL,
       180, 148.717, 274.683, 100.283, 080.600, 248.133, 181.817, 082.600 ,180},
           {/*  4 */ "6n", 1, "certain", "lilac      ", "dom", NULL,
       180, 150.311, 268.383, 084.972, 063.811, 191.483, 176.644, 085.600 ,180},
           {/*  5 */ "6j", 1, "certain", "purple     ", "sat", NULL,
       180, 141.633, 244.100, 066.056, 071.667, 122.167, 182.200, 083.622 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin  9 */
         "23 m",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "0k", 1, "wannabe", "gray       ", "ord", NULL,
            180, 149.07, 249.78, 111.52, 278.37, 207.78, 287.82,  86.65 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin 10 */
         "22 p",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "2[", 1, "certain", "sea        ", "ord", NULL,
       180, 146.383, 259.402, 291.275, 291.982, 210.048, 054.412, 147.760 ,180},
           {/*  2 */ "4b", 1, "certain", "gold       ", "ord", NULL,
       180, 145.256, 244.622, 162.822, 294.159, 171.630, 045.900, 145.804 ,180},
           {/*  3 */ "0b", 1, "certain", "red        ", "ord", NULL,
       180, 147.593, 248.421, 112.086, 274.943, 164.764, 056.843, 146.264 ,180},
           {/*  4 */ "4p", 1, "certain", "purple     ", "ord", NULL,
       180, 150.077, 260.246, 213.785, 071.900, 207.638, 056.715, 148.131 ,180},
           {/*  5 */ "6p", 1, "certain", "sky        ", "ord", NULL,
       180, 146.415, 257.831, 089.597, 067.923, 173.051, 055.513, 147.623 ,180},
           {/*  6 */ "2z", 1, "wannabe", "gray       ", "ord", NULL,
       180, 142.900, 236.550, 268.800, 180.783, 185.133, 054.467, 143.350 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin 11 */
         "22 t",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "4s", 1, "certain", "lime       ", "ord", NULL,
       180, 149.863, 247.562, 170.488, 277.938, 084.425, 176.413, 148.087 ,180},
           {/*  2 */ "2u", 1, "wannabe", "gray       ", "ord", NULL,
       180, 143.940, 258.200, 298.240, 279.640, 183.680, 183.080, 145.120 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin 12 */
         "22 m",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "!!", 0, "outlier", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  1 */ "2o", 1, "certain", "hotpink    ", "ord", NULL,
       180, 147.342, 256.475, 295.508, 287.408, 194.525, 293.725, 150.458 ,180},
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
      {/*bin 13 incomplete angles*/
         "inc ",
         {/*clst definitions:name, LOK, status, color, domsat, satinfo, angles*/
           {/*  0 */ "__", 0, "incompl", "white      ", "inc", NULL,
            0,0,0,0,0,0,0,0,0 },
           {/*  N */ "!!", 0, "nothing", "white      ", "out", NULL,
            0,0,0,0,0,0,0,0,0 } /*trailing cluster with LOK==0*/
         }
      },
    };
#endif /*DEFINITIONS*/

                      /*axis weights (effective widths) */
#ifdef DECLARATIONS
              /*general cluster radial widths: */
extern char clusterhalfwidthsversion[7];
extern float  deltamw, epsilonw, zetaw, alphaw, betaw, gammaw, deltaw;
         /*narrower widths for 9 satellite clusters 070328*/
extern float  epsilonsatw, zetasatw, alphasatw, betasatw; /*070328*/
#endif /*DECLARATIONS*/

#ifdef DEFINITIONS
   char clusterhalfwidthsversion[] = "070328";
   float    deltamw  = 28;
   float    epsilonw = 60, epsilonsatw = 50; /*satw 070328*/
   float    zetaw    = 55, zetasatw    = 50; /*satw 070328*/
   float    alphaw   = 50, alphasatw   = 45; /*satw 070328*/
   float    betaw    = 70, betasatw    = 60; /*satw 070328*/
   float    gammaw   = 35;
   float    deltaw   = 28;
#endif /*DEFINITIONS*/

SUITENINIT int    Lgeneralsatw; /*070328 modify general case satellite widths*/


                  /* triage filters (single axis limits) */
#ifdef DECLARATIONS
extern char axeslimitsversion[7];
extern float epsilonmin, epsilonmax;
extern float delta3min , delta3max ;
extern float delta2min , delta2max ;
extern float gammapmin , gammapmax ;
extern float gammatmin , gammatmax ;
extern float gammammin , gammammax ;
extern float alphamin  , alphamax  ;
extern float betamin   , betamax   ;
extern float zetamin   , zetamax   ;
#endif /*DECLARATIONS*/

#ifdef DEFINITIONS
   char  axeslimitsversion[] = "070326"; 
   float epsilonmin = 155, epsilonmax = 310; /*070130*/
   float delta3min  =  60, delta3max  = 105; /* changed by S.J. on 06/06/2011 */
   float delta2min  = 125, delta2max  = 165;
   float gammapmin  =  20, gammapmax  =  95; /*max 070326*/
   float gammatmin  = 140, gammatmax  = 215; /*max 070326*/
   float gammammin  = 260, gammammax  = 335; /*max 070326*/
   float alphamin   =  25, alphamax   = 335;
   float betamin    =  50, betamax    = 290;
   float zetamin    =  25, zetamax    = 335;
#endif /*DEFINITIONS*/

         /*hyperellipsoid power default, override with  -power #.# */
#ifdef DECLARATIONS
extern double power; /*hyperellipsoid power*/
#endif /*DECLARATIONS*/

#ifdef DEFINITIONS
   double power = 3.0; /*3.0 070414; 2.5 070326; 2.0 earlier */
#endif /*DEFINITIONS*/

SUITENINIT int   dominate[13], satellites[13][4];
SUITENINIT float clusterav[13][MAXCLST][8];

                      /*flags*/
SUITENINIT int   Lsuitesin,Lresiduesin,NptIDfields,Nanglefields,Lnewfile,altIDfield; /* added altIDfield and altID variable - S.J. 01/07/2014*/
SUITENINIT int   Lstringout,Lhelpout,Lchangesout,Lkinemageout,LNptIDfields,LaltIDfield,Luseincorrect,LaltID; /* added last four variables as flags to check for correct usage of altIDfield - S.J. 01/07/2014*/
SUITENINIT int   Lreportout,Ltestout; 
SUITENINIT int   Lchart,Lsourout,Letatheta,Ldangle; /*070521,070524,070525*/
SUITENINIT char  NameStr[256],altID[2]; /* last variable added by S.J. 01/07/2014*/
SUITENINIT int   Lsequence, Loverlap, Loneline; /*070409*/
SUITENINIT int   Lwannabe; /*070429*/

                      /*prototypes*/
int   initializations(void); /*070325*/
void  assignweights(int, int, float*, float*);
int   parsecommandline(int*, char**);

