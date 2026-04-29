#
# Ported "AS IS" by P. Afonine, 10-JUN-2025
# Subject to revision for coding complience and formal testing (no tests!)
#
#

#
#   calculating FSC-limit value using 1/d**3 scale histograms     A.Urzhumtsev, 2025-06-10
#
from __future__ import absolute_import, division, print_function
import math

# input parameters

# 5 arrays for the complex Fourier coefficients :
# ss2
# Fx1, Fy1 - real and imaginary part of the Fourier coefficients, set1
# Fx2, Fy2 - real and imaginary part of the Fourier coefficients, set1
# FSC_cut   = 0.143 - FSC cut-off
# nbin_min  = 100   - miminal allower number of reflections per bin
# precision = 0.01  - relative (with respect to the value) precision of the anwser

# output values

# res_val     - estimated resolution corresponding to the FSC cut-off
# res_del     - estimated precision, in A
# res_min     - bottom bound of the estimation interval
# res_max     - upper bound of the estimation interval
# Ncycle      - numer of cycles done
# hist_middle - histogram, when the searched bin is at a middle of the histogram
#               indicated by reason = 'hist'
# hist_final  - histogram at the end of the search; indicated by reason = 'last'
# reason      - termination flag ('none' if failed to find the value)

# each histogram line contains 9 numbers (hist_save write such histograms to file):
# Nbin                     - number of corefficnets per bin)
# ss3min ss3max            - bin resolution limits, in 1/d**3)
# resmin                   - bin resoluion in A, corresponding to ss3max
# FSC                      - mean FSC value)
# Fav/Fmax(1), Fav/Fmax(2) - ratio of the mean to the maximum amplitude in bin,
#                            sets 1 and 2
# Fav/Fav0(1), Fav/Fav0(2) - ratio of the mean amplitudes: in this and in the 1st bin
#                            sets 1 and 2

def run(f1, f2):
  from scitbx.array_family import flex
  Fx1, Fy1, Fx2, Fy2 = [],[], [],[]
  ss2 = list( 1./flex.pow2(f1.d_spacings().data()) )
  for d1,d2 in zip(f1.data(), f2.data()):
    Fx1.append( d1.real )
    Fy1.append( d1.imag )
    Fx2.append( d2.real )
    Fy2.append( d2.imag )
  r = FSC_limit(ss2=ss2, Fx1=Fx1, Fy1=Fy1, Fx2=Fx2, Fy2=Fy2,
                    FSC_cut=0.143, nbin_min=100, precision=0.01)
  res_val, res_del, res_min, res_max, Ncycle, hist_middle, hist_final, reason = r
  return res_val, res_del

def FSC_limit(ss2, Fx1, Fy1, Fx2, Fy2, FSC_cut, nbin_min, precision, verbose=False) :
   if verbose:
     print()
     print('resolution estimation corresponding to the FSC cut-off =',f'{FSC_cut:6.3f}')
     print('relative presition required',26*' ',f'= {precision:6.3f}')
     print('minimal allowed number of coefficients per bin         = ',f'{nbin_min:5}')
     print()

   Ncoef = len(ss2)

   ss3  = [0.0 for icoef in range(Ncoef)]
   Fsq1 = [Fx1[icoef]*Fx1[icoef] + Fy1[icoef]*Fy1[icoef] for icoef in range(Ncoef)]
   Fsq2 = [Fx2[icoef]*Fx2[icoef] + Fy2[icoef]*Fy2[icoef] for icoef in range(Ncoef)]
   Fr12 = [Fx1[icoef]*Fx2[icoef] + Fy1[icoef]*Fy2[icoef] for icoef in range(Ncoef)]

   dhigh       = 0.0
   Ncycle      = 0
   reason      = 'none'
   hist_middle = []
   hist_final  = []

   Ncoef_hist = Ncoef
   ss2min    = ss2[0]
   ss2max    = 0.0
   for icoef in range(Ncoef) :
       ss2coef = ss2[icoef]
       if ss2min > ss2coef : ss2min = ss2coef
       if ss2max < ss2coef : ss2max = ss2coef
       ss3coef = ss2coef * math.sqrt(ss2coef)
       ss3[icoef] = ss3coef
   ss3max = ss2max * math.sqrt(ss2max)

   res_min0 = 1.0 / math.sqrt(ss2max)
   res_max0 = 1.0 / math.sqrt(ss2min)
   if verbose:
     print('total number of coefficients    ',f'{Ncoef:15}')
     print('in the resolution range (A)',f'{res_min0:10.4f}{res_max0:10.4f}')
     print()
     print('In the table below :')
     print('N        - consecutive iteration number')
     print('FSC_res  - estimated resolution for the FSC cut-off')
     print('res_prec - precision of the estimate')
     print('res_prev - higher-resolution estimate bound')
     print('res_next - lower-resolution estimate bound')
     print('Ntotal   - total number of coefficients in the histogram')
     print('Nused    - number of coefficients till the bin with FSC < cut-off')
     print('Ncoef(1) - number of coefficients in the first bin')
     print('Ncoef(F) - number of cofficients in the last bin')
     print('Ncoef(N) - number of cofficients in the bin with FSC < cut-off')
     print('Nbins    - total number of bin')
     print('Fbin     - bin number with FSC < cut-off')
     print('bintype  - type of the interval with FSC < cut-off')
     print()
     print('  N FSC_res res_prec  res_prev  res_next     Ntotal    Nused',
           'Ncoef(1) Ncoef(F) Ncoef(N) Nbins Fbin type')

   while True :

      Ncycle += 1

#     calculate histogram

      hist, shrink = Calc_Hist(ss3, Fsq1, Fsq2, Fr12, ss3max, Ncoef_hist, nbin_min, precision)

#     process histogram

      ss3max, ss3min, ss3_val, dss3, Ncoef_used, ihist, Nint, reason = \
                               Proc_Hist(hist, ss3max, FSC_cut, nbin_min, precision, reason)

      Nhist  = len(hist)
      Ncoef0 = hist[0][0]
      NcoefN = hist[Nhist-1][0]
      res_min = ss3max ** float(-1.0/3.0)
      if ss3min > 0.0 :
         res_max = ss3min ** float(-1.0/3.0)
      else :
         res_max = res_max0
      res_val = ss3_val ** float(-1.0/3.0)
      res_del = (ss3max-dss3) ** float(-1.0/3.0) - res_min
      if verbose:
        print(f'{Ncycle:3} {res_val:7.3f}{res_del:9.3f} {res_min:9.3f} {res_max:9.3f}',
              f'{Ncoef_hist:10}{Ncoef_used:9}{Ncoef0:8} {Nint:8} {NcoefN:8}{Nhist:7}',
              f'{ihist+1:4}',reason)

      Ncoef_hist = Ncoef_used

      if reason == 'hist'                     : hist_middle = hist
      if reason == 'last' or reason == 'none' : hist_final  = hist

      if   reason == 'none' or reason == 'last' :
         break

   if verbose:
     print()
     if reason == 'none' :
        print('resolution estimate is beyond the dataset resolution {res_min0:8.3f}')
     else :
        print(f'resolution estimate {res_val:8.3f} in the bounds {res_min:8.3f} {res_max:8.3f}')
        if (shrink < 1.0) :
           print('warning : too few data to get a required accuracy')

   return res_val, res_del, res_min, res_max, Ncycle, hist_middle, hist_final, reason

#-------------------------------------------

def Proc_Hist(hist, ss3max, FSC_cut, nbin_min, precision, reason) :

   Nhist = len(hist)
   dss3 = ss3max / float(Nhist)

   Nint    = 0
   ss3min  = 0.0
   ss3_val = ss3max

   Ncoef_used = 0

#  check bins

   for ihist in range(0,Nhist) :

       Ncoef_used       += hist[ihist][0]
       hist_value = hist[ihist][4]
       if hist_value < FSC_cut :
          ss3min =  ihist    * dss3
          ss3max = (ihist+1) * dss3

#         calculate the mode of the found interval and inetrpolated cut-off

          if ihist == 0 :
             reason  = 'first'
             ss3_val = dss3
          else :
             if ihist * 4 < Nhist :
                reason = 'begin'
             else :
                if reason == 'none' or reason == 'first' or reason == 'begin' :
                   reason = 'hist'
                else :
                   reason = 'end'

             if ihist == Nhist-1   : reason = 'last'
             hist_prev = hist[ihist-1][4]
             ss3_val = ss3min + dss3 * (hist_prev - FSC_cut) / (hist_prev - hist_value)

          Nint = hist[ihist][0]
          break

   return ss3max, ss3min, ss3_val, dss3, Ncoef_used, ihist, Nint, reason

#-------------------------------------------

def Calc_Hist(ss3, Fsq1, Fsq2, Fr12, ss3max, Ncoef_hist, nbin_min, precision) :

   Ncoef = len(ss3)
   Nhist = int (0.33333333 / precision)
   Nin_bin = float(Ncoef_hist) / float(Nhist)
   shrink = Nin_bin / float(nbin_min)
   if shrink < 1.0 : Nhist = int(Nhist * shrink)


   hist = [[0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for ihist in range(Nhist)]

   dss3 = ss3max / float(Nhist)

#  collect statistics

   for icoef in range(Ncoef) :
       ss3icoef = ss3[icoef]
       ihist   = int(ss3icoef / dss3)
       if ihist < Nhist :
          F1icoef  = Fsq1[icoef]
          F2icoef  = Fsq2[icoef]
          hist[ihist][0] += 1
          hist[ihist][1] += F1icoef
          hist[ihist][2] += F2icoef
          hist[ihist][4] += Fr12[icoef]
          if hist[ihist][7] < F1icoef : hist[ihist][7] = F1icoef
          if hist[ihist][8] < F2icoef : hist[ihist][8] = F2icoef

#  process statistics

   Ncoef0    = hist[0][0]
   Fsq1Mean0 = float ( hist[0][1] / Ncoef0 )
   Fsq2Mean0 = float ( hist[0][2] / Ncoef0 )

   for ihist in range(Nhist) :
       Ncoefi = float (hist[ihist][0])
       Fsq1i  = hist[ihist][1]
       Fsq2i  = hist[ihist][2]
       hist[ihist][1] = ihist     * dss3
       ss3val         = (ihist+1) * dss3
       hist[ihist][2] = ss3val
       hist[ihist][3] = ss3val ** (-1.0 / 3.0)

       if Ncoefi > 0 :
          Fsq1Meani = Fsq1i / Ncoefi
          Fsq2Meani = Fsq2i / Ncoefi
          hist[ihist][4] = hist[ihist][4] / math.sqrt(Fsq1i * Fsq2i)
          hist[ihist][5] = math.sqrt(Fsq1Meani / hist[ihist][7])
          hist[ihist][6] = math.sqrt(Fsq2Meani / hist[ihist][8])
          hist[ihist][7] = math.log10(Fsq1Meani / Fsq1Mean0) / 2.0
          hist[ihist][8] = math.log10(Fsq2Meani / Fsq2Mean0) / 2.0

   return hist, shrink

#-------------------------------------------

def hist_save(FileHist,hist_middle, hist_final) :

   FileOut  = open(FileHist, 'w')

   Nhist = len(hist_middle)
   print('Middle-position histogram',file=FileOut)
   print(' Ncoef   ss3min    ss3max   resmin     FSC    Fav/Fmax-1 Fav/Fmax-2 Fav/Fav0-1 Fav/Fav0-2',
         file=FileOut)

   for ihist in range(Nhist) :
       (Ncoef, ss3min, ss3max, resmin, FSC_val, Fmax1, Fmax2, Faver1, Faver2) = hist_middle[ihist]
       print(f'{Ncoef:6}{ss3min:10.7f}{ss3max:10.7f}{resmin:8.3f}{FSC_val:11.7f}',
             f'{Fmax1:10.7f}{Fmax2:11.7f}{Faver1:11.7f}{Faver2:11.7f}',file=FileOut)

   Nhist = len(hist_final)
   print('',file=FileOut)
   print('Final histogram',file=FileOut)
   print(' Ncoef   ss3min    ss3max   resmin     FSC    Fav/Fmax-1 Fav/Fmax-2 Fav/Fav0-1 Fav/Fav0-2',
         file=FileOut)

   for ihist in range(Nhist) :
       (Ncoef, ss3min, ss3max, resmin, FSC_val, Fmax1, Fmax2, Faver1, Faver2) = hist_final[ihist]
       print(f'{Ncoef:6}{ss3min:10.7f}{ss3max:10.7f}{resmin:8.3f}{FSC_val:11.7f}',
             f'{Fmax1:10.7f}{Fmax2:11.7f}{Faver1:11.7f}{Faver2:11.7f}',file=FileOut)

   return
