from __future__ import absolute_import, division, print_function
from six.moves import range
import math
from libtbx.table_utils import Spreadsheet
from spotfinder.core_toolbox import bin_populations
from spotfinder.array_family import flex

'''Detailed resolution-bin analysis for augmented spotfinder.'''

class spotreporter(Spreadsheet):
  def __init__(self,spots,phil_params,targetResolution=None,wavelength=None,
               max_total_rows=None,fractionCalculator=None,use_binning_of=None):

    if len(spots) == 0: return # No spots, no binwise analysis
    if use_binning_of is None or not hasattr(use_binning_of,"binning"):
      self.original_binning=True
      assert targetResolution is not None and wavelength is not None
      self.wavelength = wavelength

      #Total rows by formula
      Total_rows = int(math.pow(
         targetResolution/spots[len(spots)-1],3.))+1
      if max_total_rows is not None:
        Total_rows = min(Total_rows, max_total_rows)

      Spreadsheet.__init__(self,rows=Total_rows)
      self.addColumn('Limit')
      self.addColumn('Fract')
      self.addColumn('Missing') # fraction res shell remaining after missing cone is subtracted
      self.addColumn('Population',0)

      if fractionCalculator is not None and \
        phil_params.distl.bins.corner:
        last_d_star = 1./flex.min(fractionCalculator.corner_resolutions())
      else:
        last_d_star = 1./spots[-1]
      volume_fraction = (1./Total_rows)*math.pow(last_d_star,3)

      def BinRes(shellnumber):
        return 1./math.pow((shellnumber+1)*volume_fraction,1.0/3.0)

      #this little loop consumes 1.2 seconds of CPU time for a 1400-row spreadsheet
      for xrow in range(0,Total_rows):
        self.Limit[xrow] = BinRes(xrow)
        self.Fract[xrow] = fractionCalculator(self.Limit[xrow])

      self.populate_missing_cone(Total_rows)

    if use_binning_of is not None:
      if not hasattr(use_binning_of,"binning"):
        use_binning_of.binning = self #inject this spreadsheet into caller's space
      else:
        self.original_binning=False
        self.wavelength = wavelength

        Spreadsheet.__init__(self,rows=use_binning_of.binning.S_table_rows)
        self.addColumn('Limit')
        self.addColumn('Fract')
        self.addColumn('Missing')
        self.addColumn('Population',0)

        for xrow in range(0,self.S_table_rows):
          self.Limit[xrow] = use_binning_of.binning.Limit[xrow]
          self.Fract[xrow] = use_binning_of.binning.Fract[xrow]
          self.Missing[xrow] = use_binning_of.binning.Missing[xrow]
        Total_rows = self.S_table_rows

    bp = bin_populations(spots,self.Limit[0])
    for c in range(min(Total_rows,len(bp))): #reconcile inconsistent bin counts
        self.Population[c]=bp[c]

    self.Limit.format = "%.2f"
    self.Fract.format = "%.2f"
    self.Missing.format = "%.2f"
    self.Population.format = "%7d"


  def total_signal(self,fstats,indices):
    #for index in indices:
    #  print index,fstats.master[index].intensity(),flex.sum(fstats.master[index].wts)
    self.addColumn('Integrated',0) #total integrated signal within the shell
    self.Integrated.format = "%.0f"
    self.addColumn('MeanI',0) #mean integrated signal within the shell
    self.MeanI.format = "%.0f"
    self.addColumn('MeanBkg',0) #mean integrated background within the shell
    self.MeanBkg.format = "%.1f"
    self.addColumn('MeanSz',0) #mean pixel count for spots within the shell
    self.MeanSz.format = "%.0f"
    self.addColumn('MnEccen',0) #mean eccentricity spots within the shell
    self.MnEccen.format = "%.3f"
    self.addColumn('MnSkew',0) #mean eccentricity spots within the shell
    self.MnSkew.format = "%.3f"
    index=0
    for row in range(self.S_table_rows):
      row_values = flex.double()
      bkg_values = flex.double()
      sz_values = flex.double()
      eccen_values = flex.double()
      skew_values = flex.double()
      for idx in range(self.Population[row]):
        spot_signal = flex.sum(fstats.master[indices[index]].wts)
        bkg_signal = flex.sum(fstats.master[indices[index]].bkg)
        spot_sz = len(fstats.master[indices[index]].wts)
        spot_eccen = fstats.master[indices[index]].model_eccentricity()
        spot_skew = fstats.master[indices[index]].skewness()
        self.Integrated[row]+=spot_signal
        row_values.append(spot_signal)
        bkg_values.append(bkg_signal)
        sz_values.append(spot_sz)
        eccen_values.append(spot_eccen)
        skew_values.append(spot_skew)
        index+=1
      if len(row_values)>0:
        self.MeanI[row]=flex.mean(row_values)
        self.MeanBkg[row]=flex.mean(bkg_values)
        self.MeanSz[row]=flex.mean(sz_values)
        self.MnEccen[row]=flex.mean(eccen_values)
        self.MnSkew[row]=flex.mean(skew_values)
    self.persist_fstats = fstats
    self.persist_indices = indices

  def populate_missing_cone(self,total_rows):
    from spotfinder.math_support.sphere_formulae import sphere_volume
    from spotfinder.math_support.sphere_formulae import sphere_volume_minus_missing_cone
    last_cumulative_sphere_volume = 0.0
    last_cumulative_corrected_volume = 0.0
    for xrow in range(0,total_rows):
      this_sphere_volume = sphere_volume(1./self.Limit[xrow])
      this_corrected_volume = sphere_volume_minus_missing_cone(1./self.Limit[xrow],self.wavelength)
      shell_volume = this_sphere_volume - last_cumulative_sphere_volume
      shell_corrected_volume = this_corrected_volume - last_cumulative_corrected_volume
      self.Missing[xrow] = shell_corrected_volume/shell_volume
      last_cumulative_sphere_volume = this_sphere_volume
      last_cumulative_corrected_volume = this_corrected_volume

  def background(self,spotfinder):
    self.addColumn('PxlBkgrd') # mean pixel background at center of windows in the shell
    self.PxlBkgrd.format = "%.1f"
    self.addColumn('MnWndwSz') # mean size (in pixels) of the scanbox windows in the shell
    self.MnWndwSz.format = "%.0f"
    for irow in range(self.S_table_rows):
      self.PxlBkgrd[irow] = flex.double()
      self.MnWndwSz[irow] = flex.double()
    sortedidx = flex.sort_permutation(spotfinder.background_resolutions(),reverse=True)
    reverse_resolutions = spotfinder.background_resolutions().select(sortedidx)
    reverse_backgrounds = spotfinder.background_means().select(sortedidx)
    reverse_wndw_sz = spotfinder.background_wndw_sz().select(sortedidx)
    row = 0
    for x in range(len(reverse_resolutions)):
      while reverse_resolutions[x] < self.Limit[row]:
        row += 1
        if row >= self.S_table_rows: break
      if row >= self.S_table_rows: break
      self.PxlBkgrd[row].append(reverse_backgrounds[x])
      self.MnWndwSz[row].append(reverse_wndw_sz[x])
    for irow in range(self.S_table_rows):
      #print irow, len(self.PxlBkgrd[irow])
      #print list(self.PxlBkgrd[irow])
      #print list(self.PxlSigma[irow])
      if len(self.PxlBkgrd[irow])==0: self.PxlBkgrd[irow]=None; continue #No analysis without mean background
      self.PxlBkgrd[irow] = flex.mean(self.PxlBkgrd[irow])
      self.MnWndwSz[irow] = flex.mean(self.MnWndwSz[irow])

  def sigma_analysis(self):
    #assumes that self.total_signal() and self.background() have already been called.
    self.addColumn('MeanIsigI',0) #mean I over sig(I) for spots within the shell
    self.MeanIsigI.format = "%.3f"
    index=0
    for row in range(self.S_table_rows):
      if self.PxlBkgrd[row] == None: continue
      IsigI_values = flex.double()
      for idx in range(self.Population[row]):
        # Use International Tables Vol F (2001) equation 11.2.5.9 (p. 214)
        I_s = flex.sum(self.persist_fstats.master[self.persist_indices[index]].wts)
        I_bg = flex.sum(self.persist_fstats.master[self.persist_indices[index]].bkg)
        m_sz = len(self.persist_fstats.master[self.persist_indices[index]].wts)
        n_sz = self.MnWndwSz[row]
        gain = 1.00

        spot_variance = gain * (I_s + I_bg + (m_sz*m_sz) * (1./n_sz) * self.PxlBkgrd[row])

        IsigI_values.append( I_s / math.sqrt(spot_variance))

        index+=1
      if len(IsigI_values)>0: self.MeanIsigI[row]=flex.mean(IsigI_values)


  def show(self,message,default=['Limit','Population','Missing','Fract','PxlBkgrd','MnWndwSz','Integrated','MeanI',
           'MeanBkg','MeanSz','MeanIsigI','MnEccen','MnSkew']):
    self.summaries=['Population','Integrated']
    self.wtmean=['MeanI','MeanBkg','MeanSz','MeanIsigI','MnEccen','MnSkew']
    self.weights=self.Population
    legend = """Limit: outer edge of (equal-volume) resolution shell (Angstroms)
Population: number of spots in the resolution shell
Missing: fraction of reciprocal space shell recorded in complete rotation, accounting for missing cone
Fract: fraction of shell recorded, accounting for truncated detector corner
PxlBkgrd: mean pixel background at center of scanbox windows in the shell
MnWndwSz: mean scanbox window size in pixels after discarding spot pixels
Integrated: total integrated signal within the shell in analog-digital units above background
MeanI: mean integrated signal for spots within the shell
MeanBkg: mean integrated background for spots within the shell
MeanSz: mean pixel count for spots within the shell
MeanIsigI: mean I over sig(I) for spots within the shell
MnEccen: mean eccentricity for spots within the shell; 0=circle, 1=parabola
MnSkew: mean skewness for spots within the shell; skew:= (maximum-center_of_gravity)/semimajor_axis
"""; not_implemented="[PxlGain: estimate of the average gain for scanbox windows in the shell]"
    if self.original_binning: print(legend)
    print(message+":")
    to_print = [max([self.Population[i] for i in range(j,self.S_table_rows)])>0
                and self.PxlBkgrd[j] != None
           for j in range(self.S_table_rows)]
    self.printTable(default,printed_rows=to_print)

  def vetter(self,cutoff,last):
    redcutoff = 10
    """The problem:  we are given a list of shell populations, pop.
       We are told we are only interested from index 0 thru last.
       We have to find a new slice from 0 to newlast such that there is
       no contiguous stretch of populations < cutoff, which is longer than 10.
    """
    pop = []
    for x in range(0,last+1):
      if self.Population[x]>cutoff: pop.append(1)
      else: pop.append(0)

    requalify=[0,]
    for x in range(1,last+1):
      if pop[x]==0: requalify.append(requalify[x-1]+1)
      else:
         requalify.append(0)
      if requalify[x] == redcutoff:
        break

    pop = pop[0:x+1]
    for x in range(len(pop)-1,-1,-1):
      if pop[x]==1:
        #print "In Vetter with input",last,"output",x
        return x
