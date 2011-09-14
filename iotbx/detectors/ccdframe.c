/* proteus.c conversion 033.26 */
/* based on difoff.c */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

int ccdata[1024][1024], aduhist[4][10];
unsigned char ccbyte[4096]; unsigned short ccshort[4096];
unsigned int ccint[4096];

main ()
{ int i,j,k, err,m,n, hist[100], sump, summ, kount;
  int swab,ii,jj,kk,mm,nn, nhoriz, nvert, hoff,max,mxflg;
  int nunder=0, n2over=0, n4over=0, lendata=1, lenunder=1,form=0,
    nrows=1024, ncols=1024, zeroff=0, baseoff=0;
        short s,sptp[3000], sptm[3000], delspt[3000];
        unsigned short bs;
        /*      int bx = 1492, by = 1666 */
        int  pd, md, pdo, mdo, zz;
        double dsum, esum, rat;  float x,ave;
      union { unsigned char b[2];  short s; unsigned short us; } u;
      union { unsigned char b[4];  int i; unsigned int ui; } u4;
FILE *fi, *fc, *fo; char filename[136],line[136];

  printf(" Proteus CCD  frame file:"); scanf("%s",filename);
  fi=fopen(filename,"rb");
        if (fi==NULL) { printf(" open error (+) \n"); exit(1); };

  printf(" xdip outfile:"); scanf("%s",filename);
  fo=fopen(filename,"wb");
        if (fo==NULL) { printf(" open error (out) \n"); exit(1); };

  printf(" headeroutfile [] :"); scanf("%s",filename); fc = NULL;
  if (filename[0]!=' ') {  fc=fopen(filename,"w");
        if (fc==NULL) { printf(" open error (out) \n"); exit(1); }; }

  pd = md = pdo = mdo = zz = 0;   dsum = esum = 0.0;
  for (j=0; j<100; j++) hist[j] = 0;

  /* look for format info, etc. in header lines... */
  nhoriz = nvert = 1024;  swab = 1; line[80] = 0;

  /* 15 blocks of 80 char lines == 96 lines */
  for (j=0;j<96;j++) { n=fread(line,1,80,fi); printf("%s\n",line);
  if (strncmp(line,"FORMAT :",8)==0)
    { n = sscanf(line+8," %d",&form);
         printf(">>>>> FORMAT %d \n",form); }
  else if (strncmp(line,"NOVERFL:",8)==0)
    { n = sscanf(line+8," %d %d %d",&nunder,&n2over,&n4over);
      printf(">>>>> %d overflows  %d %d %d\n",n,nunder,n2over,n4over);
      if (form<100) { n4over = nunder; nunder = n2over = 0; } }
  else if (strncmp(line,"NPIXELB:",8)==0)
    { n = sscanf(line+8," %d %d",&lendata,&lenunder);
      printf(">>>>>pixel bytes data %d under %d\n",lendata,lenunder); }
  else if (strncmp(line,"NROWS  :",8)==0)
    { n = sscanf(line+8," %d",&nrows);
      printf(">>>>rows  %d \n",nrows); nvert = nrows; }
 else if (strncmp(line,"NCOLS  :",8)==0)
    { n = sscanf(line+8," %d",&ncols);
      printf(">>>>>columns  %d\n",ncols); nhoriz = ncols; }
 else if (strncmp(line,"NEXP   :",8)==0)
    { n = sscanf(line+8," %d %d %d",&k,&k,&zeroff);
      printf(">>>>>zero offset  %d \n",zeroff); }

  if (fc) { for (k=79; k>10; k--) if (line[k]!=' ') break;
               line[k+1] = 0; fprintf(fc,"%s\n",line); }

  } /* for (j<96) */

# if (0)
  /* zeroes from 7680 to 9712 */
 n = fread( sptp, sizeof(s),(size_t) 1016, fi);
    if (j<10) { /* determine byte swap */ kk = mm = 0; m = 2*nhoriz;
        for (k=0; k<nhoriz; k++) { u.s = sptp[k];
                         if (u.b[0]>8) kk++; if (u.b[1]>8) mm++; }
            printf(" %4d large[0] %4d large[1]\n",kk,mm);
            if (mm<kk) swab = 1; else swab = 0;
            u.s = sptp[500];
            if (swab) { bs = u.b[0]; u.b[0] = u.b[1]; u.b[1] = bs; }
            printf(" u.s %d u.b[0] %d u.b[1] %d \n",u.s,u.b[0],u.b[1]);
            }

# endif
    /* 03x.29 have probably misinterpreted zeroff... should have subtracted */
    /* this isn''t working:     zeroff = -zeroff;  */
    /* may differ between format 86 and format 100 files */
    /* try defining ave = esum - zeroff;  */

if (lendata==1) { for (k=0; k<nrows; k++)
  {   n = fread(ccbyte , lendata, (size_t) ncols, fi);
      if (n != ncols) printf(" short data %d  row %d\n",n,k);
      for (j=0; j<ncols; j++) ccdata[k][j] = ccbyte[j]+zeroff;
  } /* for (k) */  }  /* if (lendata==1) */
 else if (lendata==2) { for (k=0; k<nrows; k++)
  {   n = fread(ccshort , lendata, (size_t) ncols, fi);
      if (n != ncols) printf(" short data %d  row %d\n",n,k);
      for (j=0; j<ncols; j++) ccdata[k][j] = ccshort[j]+zeroff;
  } /* for (k) */  }  /* if (lendata==1) */
 else  printf("input not implemented for length %d\n",lendata);

 n = m = 0;  ii = 255; if (lendata==2) ii = 65535; mxflg = zeroff + ii;
 for (j=0;j<ncols;j++) for (k=0;k<nrows;k++)
   { i = ccdata[k][j]; if (i==zeroff) m++; if (i==mxflg) n++; }

 printf(" %d byte data read under %d over %d\n",lendata,m,n);

 if (form>99) {  /* binary overflow data for FORMAT 100+ */
   nn = nunder; kk = jj = 0; mm = 512;
while (nn>0) { if (nn<512) mm = nn; nn -= 512; ii=0;
      n = fread(ccbyte, lenunder, (size_t) mm, fi);
      if (n != mm) printf(" short under %d at %d %d\n",n,kk,jj);
      while (ii<mm) { if (ccdata[kk][jj]==zeroff)
                           { ccdata[kk][jj] = ccbyte[ii]; ii++;}
                      jj++; if (jj>=ncols) { jj=0; if (kk<nrows-1) kk++; }
     }  /* while (ii<mm) */
} /* while (nn>0) get underflows */
  mm = nunder%16; if (mm>0)  n = fread(ccbyte,lenunder, (size_t) 16-mm, fi);

 printf(" underflows read,  %d\n",mm);

   nn = n2over; kk = jj = 0; mm = 512;  m = zeroff+255;
while (nn>0) { if (nn<512) mm = nn; nn -= 512; ii=0;
      n = fread(ccshort, 2, (size_t) mm, fi);
      if (n != mm) printf(" short over %d  at %d  %d\n",n,kk,jj);
     if (swab) { for (i=0;i<mm; i++)
                    { u.us = ccshort[i];  bs = u.b[0]; u.b[0] = u.b[1];
                      u.b[1] = bs; ccshort[i] = u.us; } }
      while (ii<mm) { if (ccdata[kk][jj]==m)
                           { ccdata[kk][jj] = ccshort[ii] + zeroff; ii++;}
                      jj++; if (jj>=ncols) { jj=0; if (kk<nrows-1) kk++; }
     }  /* while (ii<mm) */
} /* while (nn>0) get 2 byte overflows */
  mm = n2over%8; if (mm>0)  n = fread(ccshort, 2, (size_t) 8-mm, fi);

 printf(" n2over read, extra %d\n", mm);

   nn = n4over; kk = jj = 0; mm = 512;  m = zeroff+65535;
while (nn>0) { if (nn<512) mm = nn; nn -= 512; ii=0;
      n = fread(ccint, 4, (size_t) mm, fi);
      if (n != mm) printf(" int over %d  at %d  %d\n",n,kk,jj);
     if (swab) { for (i=0;i<mm; i++)
                    { u4.ui = ccint[i];
                      printf( "%d %x  %d ",u4.ui, u4.ui, u4.ui+zeroff);
                      bs = u4.b[0]; u4.b[0] = u4.b[3]; u4.b[3] = bs;
                      bs = u4.b[1]; u4.b[1] = u4.b[2]; u4.b[2] = bs;
                      printf( " swapped  %d %x \n",u4.ui,u4.ui);
                      ccint[i] = u4.ui; } }
      while (ii<mm) { if (ccdata[kk][jj]==m)
                           {
                             ccdata[kk][jj] = ccint[ii] + zeroff;
                             printf("n4over at %d %d\n",kk,jj);
                             ii++;}

                      jj++; if (jj>=ncols) { jj=0; if (kk<nrows-1) kk++; }
     }  /* while (ii<mm) */
} /* while (nn>0) get 4 byte overflows */
  mm = n4over%4; if (mm>0)  n = fread(ccshort, 4, (size_t) 4-mm, fi);

 printf(" n4over read\n");
 }  /* if (form>99) */
 else { /* 16 byte ascii overflow records */  line[16] = 0;
   for (j=0; j<n4over; j++) { n = fread(line, 1, (size_t) 16, fi);
       if (n==16) { sscanf(line+9,"%7d",&ii);  line[9] = 0;
                     sscanf(line,"%9d",&k);
        jj = ii%nrows; kk = ii/nrows;
        n = ccdata[kk][jj]; if (n==mxflg) ccdata[kk][jj] = k;
        if (j<20) printf("%s value %10d pixel %10d  [ %d ][ %d ] old %d\n",
                       line,k,ii,jj,nhoriz-kk-1,n); } }
                         /*  old indices      line,k,ii,kk,jj,n); } } */
 } /* if (form<100) */

 hoff = (3000-nhoriz)/2; max = 0;

 dsum = 0.0;
 for (j=0; j<nvert; j++) { summ = 0;
    for (k=0; k<nhoriz; k++) summ += ccdata[j][k];
    dsum += summ; }
 esum = dsum/(nvert*nhoriz);
 ave = esum - zeroff;
 printf("average value %g  esum %g zeroff %d\n",ave,esum,zeroff);

 for (mm=0;mm<4;mm++) for (nn=0;nn<10;nn++) aduhist[mm][nn] = 0;
 for (j=0; j<nvert; j++) for (k=0; k<nhoriz; k++)
  {  x = (ccdata[j][k]-zeroff)/ave;      /* was ccdata/esum */
     if (j<3||j>nvert-4||k<3||k>nhoriz-4) mm = 3; else
     if (j<9||j>nvert-9||k<9||k>nhoriz-9) mm = 2; else mm = 1;
     if (x<0.25) kk = 0; else if (x<0.5) kk = 1; else if (x<0.75) kk = 2;
     else if (x<1.0) kk = 3; else if (x<1.5) kk = 4; else if (x<2.0) kk = 5;
     else if (x<3.0) kk = 6; else if (x<4.0) kk = 7; else kk = 8;
     aduhist[mm][kk]++;  aduhist[0][kk]++; }

 for (mm=0;mm<4;mm++) { for (nn=0;nn<9;nn++) printf(" %6d",aduhist[mm][nn]);
      printf("\n"); }

 x = 4.0;
 while (x>0.0) {
   m = x*ave + zeroff; kount = 0;   /* was x*esum */
 for (j=0; j<nvert; j++) for (k=0; k<nhoriz; k++)
   { kk = j;  jj = nhoriz - k - 1;  /* SMART index permutation */
     if (x>1.0) { if(ccdata[jj][kk]>m) { kount++;
        printf("large pixel %7d  x %4d y %4d\n",ccdata[jj][kk],j,k);}}
     else { if(ccdata[jj][kk]<m) { kount++;
     if (kount<10 || ((j>16&&j<1008)&&(k>16&&k<1008)))
       printf("small pixel %7d  x %4d y %4d\n",ccdata[jj][kk],j,k);}}
   } /* for (j) for (k) scan frame */
  printf("%f kount  %d new factor: ",x,kount);
  scanf("%g",&x); } /* while (x>0.0) */

  for (k=0; k<3000; k++) delspt[k] = zeroff;
if (fo) { for (k=0; k<hoff; k++)
             n = fwrite( delspt, sizeof(s),(size_t) 3000, fo); }

for (j = 0; j < nvert; j++)
{    for (k=0; k<nhoriz; k++) { m = ccdata[j][k]; if (m>max) max = m;
              if (m>32767) m = -m/32 - 1; delspt[k+hoff] = m; }
           if (fo) n = fwrite( delspt, sizeof(s),(size_t) 3000, fo);
        }

  for (k=0; k<3000; k++) delspt[k] = zeroff;
if (fo) { for (k=0; k<hoff; k++)
           n = fwrite( delspt, sizeof(s),(size_t) 3000, fo);
/* pseudo MAC Science trailer */
           n = fwrite( delspt, sizeof(s),(size_t) 512, fo);  }
printf(" max value %d\n",max);

exit(0);

for (j = 0; j < 3000; j++)
        {  n = fread( sptp, sizeof(s),(size_t) 3000, fi);
           if (n != 3000) { err = ferror(fi);
             printf(" (+) read err %d  n %d \n",err,n); /*exit(1);*/ }
           n = fread( sptm, sizeof(s),(size_t) 3000, fc);
           if (n != 3000) { err = ferror(fc);
             printf(" (-) read err %d  n %d \n",err,n); /*exit(1);*/ }
           for (i=0; i < 3000; i++)
           {  n = sptp[i]; if (n<0) n = 32*(65536+n);
              m = sptm[i]; if (m<0) m = 32*(65536+m);
              sump += n; summ += m;  }
           dsum += sump;  esum += summ;  rat = dsum; rat = rat/esum;

           for (i=0; i < 3000; i++)
           {  n = sptp[i]; if (n<0) n = 32*(65536+n);
              m = sptm[i]; if (m<0) m = 32*(65536+m);
              k = n - rat*m;
                if (k>0) { pd++; if (k>99) pdo++; }
/* use absolute value for initial statistics...  */
                if (k<0) { md++; if (k<-99) mdo++;  /* k = -k; */ }
        /*      if (k<99) hist[k]++; else hist[99]++;  */

 k += 100;  if (k<1) k = 0;

                if (k>32767) k = k/32 - 65536;
                if (n==0&&m==0) zz++;
                delspt[i] = k;
           }   /* for i */
           n = fwrite( delspt, sizeof(s),(size_t) 3000, fo);
        }  /* for j */
   printf(" pd %d md %d pdo %d mdo %d zz %d \n",pd,md,pdo,mdo,zz);
  printf(" +sum %g -sum %g ratio (+/-) %g \n",dsum,esum,dsum/esum);
   for (j=0; j<100; j += 10)
   { for (k=0; k<10; k++) printf("%6d ",hist[j+k]); printf("\n"); }

/* transfer one header-trailer */
        n = fread(sptp, 1, 1024, fi);
        n = fwrite(delspt, 1, 1024, fo);
   err = fclose(fi);
   err = fclose(fc);
   err = fclose(fo);
}
