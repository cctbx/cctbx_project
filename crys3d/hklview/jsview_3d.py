
from __future__ import absolute_import, division, print_function
from libtbx.math_utils import roundoff
from cctbx.miller import display2 as display
from cctbx.array_family import flex
from scitbx import graphics_utils
from cctbx import miller
from libtbx.utils import Sorry, to_str
from websocket_server import WebsocketServer
import threading, math, sys
from time import sleep
import os.path, time
import libtbx
import webbrowser, tempfile
from six.moves import range




class ArrayInfo:
  def __init__(self, millarr, mprint=sys.stdout.write):
    from iotbx.gui_tools.reflections import get_array_description
    data = millarr.data()
    if (isinstance(data, flex.int)):
      data = [e for e in data if e!= display.inanval]
    if millarr.is_complex_array():
      data = flex.abs(millarr.data())
    data = [e for e in data if not math.isnan(e)]
    self.maxdata =max( data )
    self.mindata =min( data )
    self.maxsigmas = self.minsigmas = None
    if millarr.sigmas() is not None:
      data = millarr.sigmas()
      data = [e for e in data if not math.isnan(e)]
      self.maxsigmas =max( data )
      self.minsigmas =min( data )
    self.minmaxdata = (roundoff(self.mindata), roundoff(self.maxdata))
    self.minmaxsigs = (roundoff(self.minsigmas), roundoff(self.maxsigmas))
    self.labels = self.desc = ""
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    if millarr.info():
      self.labels = millarr.info().label_string()
      self.desc = get_array_description(millarr)
    self.span = ("?" , "?")
    dmin = 0.0
    dmax = 0.0
    try:
      self.span = ( millarr.index_span().min(), millarr.index_span().max())
      dmin = millarr.d_max_min()[1]
      dmax = millarr.d_max_min()[0]
    except Exception as e:
      mprint(to_str(e))
    issymunique = millarr.is_unique_set_under_symmetry()
    self.infotpl = (self.labels, self.desc, millarr.indices().size(), self.span,
     self.minmaxdata, self.minmaxsigs, (roundoff(dmin), roundoff(dmax)), issymunique )
    self.infostr = "%s (%s), %s HKLs: %s, MinMax: %s, MinMaxSigs: %s, d_minmax: %s, SymUnique: %d" %self.infotpl


class hklview_3d:
  def __init__ (self, *args, **kwds) :
    self.settings = kwds.get("settings")
    self.miller_array = None
    self.d_min = None
    self.scene = None
    self.NGLscriptstr = ""
    self.cameratype = "orthographic"
    self.url = ""
    self.binarray = "Resolution"
    self.icolourcol = None
    self.iradiicol = None
    self.iarray = None
    self.isnewfile = False
    self.binvals = []
    self.workingbinvals = []
    self.proc_arrays = []
    self.otherscenes = []
    self.othermaxdata = []
    self.othermindata = []
    self.othermaxsigmas = []
    self.otherminsigmas = []
    self.sceneisdirty = False
    self.matchingarrayinfo = []
    self.match_valarrays = []
    self.binstrs = []
    self.mapcoef_fom_dict = {}
    self.verbose = True
    if kwds.has_key('verbose'):
      self.verbose = kwds['verbose']
    self.mprint = sys.stdout.write
    if 'mprint' in kwds:
      self.mprint = kwds['mprint']
    self.nbin = 0
    self.websockclient = None
    self.lastmsg = ""
    self.StartWebsocket()
    tempdir = tempfile.gettempdir()
    self.hklfname = os.path.join(tempdir, "hkl.htm" )
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)
    if 'htmlfname' in kwds:
      self.hklfname = kwds['htmlfname']
    self.hklfname = os.path.abspath( self.hklfname )
    self.jscriptfname = os.path.join(tempdir, "hkljstr.js")
    if os.path.isfile(self.jscriptfname):
      os.remove(self.jscriptfname)
    if 'jscriptfname' in kwds:
      self.jscriptfname = kwds['jscriptfname']
    self.mprint('Output will be written to \"%s\"\n' \
      'including reference to NGL JavaScript \"%s\"' %(self.hklfname, self.jscriptfname))
    self.hklhtml = r"""
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">

<html>
<head>
<meta charset="utf-8" />
</head>

<body>
<script src="%s" type="text/javascript"></script>

<script src="%s" type="text/javascript"></script>
    """
    self.htmldiv = """
<div id="viewport" style="width:100%; height:100%;"></div>

</body></html>

    """
    self.colourgradientvalues = []
    self.UseOSBrowser = True
    if 'UseOSBrowser' in kwds:
      self.UseOSBrowser = kwds['UseOSBrowser']
    self.viewmtrxelms = None
    self.pendingmessage = None



  def __exit__(self, exc_type, exc_value, traceback):
    # not called unless instantiated with a "with hklview_3d ... " statement
    self.server.shutdown()
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)


  def set_miller_array(self, col, merge=None, details="", proc_arrays=[]):
    self.iarray = col
    if col:
      self.miller_array = proc_arrays[col]
    self.proc_arrays = proc_arrays
    self.merge = merge
    if (self.miller_array is None):
      return
    self.d_min = self.miller_array.d_min()
    array_info = self.miller_array.info()
    self.binvals = [ 1.0/self.miller_array.d_max_min()[0], 1.0/self.miller_array.d_max_min()[1]  ]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % self.miller_array.unit_cell().parameters()
    self.mprint( "Data: %s %s, %d reflections in space group: %s, unit Cell: %s" \
      % (array_info.label_string(), details, self.miller_array.indices().size(), \
          self.miller_array.space_group_info(), uc) )
    #self.construct_reciprocal_space(merge=merge)


  def ExtendMillerArraysUnionHKLs(self):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.match_valarrays = []
    # loop over all miller arrays to find the subsets of hkls common between currently selected
    # miler array and the other arrays. hkls found in the currently selected miller array but
    # missing in the subsets are populated populated with NaN values
    # create miller indices being a superset of available indices in all arrays
    self.mprint("Gathering superset of miller indices")
    superset_array = self.proc_arrays[0]
    for i,validarray in enumerate(self.proc_arrays):
      if i==0:
        continue
      # first match indices in currently selected miller array with indices in the other miller arrays
      matchindices = miller.match_indices(superset_array.indices(), validarray.indices() )
      #print validarray.info().label_string()
      valarray = validarray.select( matchindices.pairs().column(1) )
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      if valarray.anomalous_flag() and not superset_array.anomalous_flag():
        # valarray gets its anomalous_flag from validarray. But it cannot have more HKLs than self.miller_array
        # so set its anomalous_flag to False if self.miller_array is not anomalous data
        valarray._anomalous_flag = False
      if not valarray.anomalous_flag() and superset_array.anomalous_flag():
        # temporarily expand other arrays to anomalous if self.miller_array is anomalous
        valarray = valarray.generate_bijvoet_mates()
      missing = superset_array.lone_set( valarray )
      # insert NAN values for reflections in self.miller_array not found in validarray
      valarray = display.ExtendMillerArray(valarray, missing.size(), missing.indices())
      match_valindices = miller.match_indices(superset_array.indices(), valarray.indices() )
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      match_valarray = valarray.select( match_valindices.pairs().column(1) )
      match_valarray.sort(by_value="packed_indices")
      #match_valarray.set_info(validarray.info() )
      superset_array = match_valarray
      #print "supersetsize:", superset_array.size()
    # now extend each miller array to contain any missing indices from the superset miller array
    for i,validarray in enumerate(self.proc_arrays):
      # first match indices in currently selected miller array with indices in the other miller arrays
      matchindices = miller.match_indices(superset_array.indices(), validarray.indices() )
      #print validarray.info().label_string()
      valarray = validarray.select( matchindices.pairs().column(1) )
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      if valarray.anomalous_flag() and not superset_array.anomalous_flag():
        # valarray gets its anomalous_flag from validarray. But it cannot have more HKLs than self.miller_array
        # so set its anomalous_flag to False if self.miller_array is not anomalous data
        valarray._anomalous_flag = False
      if not valarray.anomalous_flag() and superset_array.anomalous_flag():
        # temporarily expand other arrays to anomalous if self.miller_array is anomalous
        valarray = valarray.generate_bijvoet_mates()
      missing = superset_array.lone_set( valarray )
      # insert NAN values for reflections in self.miller_array not found in validarray
      valarray = display.ExtendMillerArray(valarray, missing.size(), missing.indices())
      match_valindices = miller.match_indices(superset_array.indices(), valarray.indices() )
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      match_valarray = valarray.select( match_valindices.pairs().column(1) )
      match_valarray.sort(by_value="packed_indices")
      match_valarray.set_info(validarray.info() )
      self.match_valarrays.append( match_valarray )


  def construct_reciprocal_space(self, merge=None):
    self.mprint("Constructing HKL scenes")
    #if len(self.match_valarrays)==0:
    #  self.ExtendMillerArraysUnionHKLs()
    self.miller_array = self.match_valarrays[self.iarray]
    if not self.sceneisdirty:
      return
    self.otherscenes = []
    self.othermaxdata = []
    self.othermindata = []
    self.othermaxsigmas = []
    self.otherminsigmas = []
    self.matchingarrayinfo = []
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    foms_array = None
    if self.miller_array.is_complex_array():
      fomcolm = self.mapcoef_fom_dict.get(self.miller_array.info().label_string())
      if fomcolm:
        foms_array = self.proc_arrays[fomcolm].deep_copy()
    #self.scene = display.scene(miller_array=self.miller_array, merge=merge,
    # settings=self.settings, foms_array=foms_array)
    # compute scenes for each of the miller arrays
    for i,match_valarray in enumerate(self.match_valarrays):
      foms = None
      if match_valarray.is_complex_array():
        fomcolm = self.mapcoef_fom_dict.get(match_valarray.info().label_string())
        #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
        if fomcolm:
          foms = self.match_valarrays[fomcolm]
      #commonarray = match_valarray
      commonarray = match_valarray.common_set(self.proc_arrays[self.iarray] )
      commonarray.set_info(match_valarray.info() )
      if i==self.iradiicol or i==self.icolourcol or i==self.iarray:
        bfullprocess = True
      else:
        bfullprocess = False
      otherscene = display.scene(miller_array=commonarray, merge=merge,
        settings=self.settings, foms_array=foms, fullprocessarray=True )
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      if not otherscene.SceneCreated:
        self.mprint("The " + commonarray.info().label_string() + " array was not processed")
        continue
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      # cast any NAN values to 1 of the colours and radii to 0.2 before writing javascript
      b = flex.bool([bool(math.isnan(e[0]) + math.isnan(e[1]) + math.isnan(e[2])) for e in otherscene.colors])
      otherscene.colors = otherscene.colors.set_selected(b, (1.0, 1.0, 1.0))
      b = flex.bool([bool(math.isnan(e)) for e in otherscene.radii])
      otherscene.radii = otherscene.radii.set_selected(b, 0.2)
      self.otherscenes.append( otherscene)
      ainf = ArrayInfo(otherscene.work_array, self.mprint)
      self.othermaxdata.append( ainf.maxdata )
      self.othermindata.append( ainf.mindata )
      self.othermaxsigmas.append(ainf.maxsigmas)
      self.otherminsigmas.append(ainf.minsigmas)
      infostr = ainf.infostr
      self.mprint("%d, %s" %(i, infostr) )
      self.matchingarrayinfo.append(infostr)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.sceneisdirty = False


  def update_settings(self, diffphil) :
    if hasattr(diffphil, "filename"):
      self.ExtendMillerArraysUnionHKLs()
      self.sceneisdirty = True

    if hasattr(diffphil, "spacegroupchoice") or \
      hasattr(diffphil, "mergedata") or \
      hasattr(diffphil, "column") or \
      hasattr(diffphil, "fomcolumn") or \
      hasattr(diffphil, "viewer") and ( hasattr(diffphil.viewer, "show_anomalous_pairs") \
      or hasattr(diffphil.viewer, "show_data_over_sigma") \
      or hasattr(diffphil.viewer, "show_missing") \
      or hasattr(diffphil.viewer, "show_only_missing") \
      or hasattr(diffphil.viewer, "show_systematic_absences") \
      or hasattr(diffphil.viewer, "slice_axis") \
      or hasattr(diffphil.viewer, "slice_mode") \
      or hasattr(diffphil.viewer, "scale") \
      or hasattr(diffphil.viewer, "nth_power_scale_radii") \
      or hasattr(diffphil.viewer, "expand_anomalous") \
      or hasattr(diffphil.viewer, "expand_to_p1")
      ):
        self.sceneisdirty = True
    self.construct_reciprocal_space(merge=self.merge)
    if self.miller_array is None or self.iarray < 0:
      return
    self.scene = self.otherscenes[self.iarray]
    self.DrawNGLJavaScript()
    msg = "Rendered %d reflections\n" % self.scene.points.size()
    return msg


  def UpdateBinValues(self, binvals = [] ):
    if binvals:
      self.binvals = binvals
    else:
      self.binvals = [ 1.0/self.miller_array.d_max_min()[0], 1.0/self.miller_array.d_max_min()[1] ]


  def DrawNGLJavaScript(self):
    if self.miller_array is None :
      self.mprint( "A miller array must be selected for drawing" )
      return
    self.mprint("Composing NGL JavaScript...")
    h_axis = self.scene.axes[0]
    k_axis = self.scene.axes[1]
    l_axis = self.scene.axes[2]
    nrefls = self.scene.points.size()

    Hstararrowstart = roundoff( [-h_axis[0]*100, -h_axis[1]*100, -h_axis[2]*100] )
    Hstararrowend = roundoff( [h_axis[0]*100, h_axis[1]*100, h_axis[2]*100] )
    Hstararrowtxt  = roundoff( [h_axis[0]*102, h_axis[1]*102, h_axis[2]*102] )
    Kstararrowstart = roundoff( [-k_axis[0]*100, -k_axis[1]*100, -k_axis[2]*100] )
    Kstararrowend = roundoff( [k_axis[0]*100, k_axis[1]*100, k_axis[2]*100] )
    Kstararrowtxt  = roundoff( [k_axis[0]*102, k_axis[1]*102, k_axis[2]*102] )
    Lstararrowstart = roundoff( [-l_axis[0]*100, -l_axis[1]*100, -l_axis[2]*100] )
    Lstararrowend = roundoff( [l_axis[0]*100, l_axis[1]*100, l_axis[2]*100] )
    Lstararrowtxt  = roundoff( [l_axis[0]*102, l_axis[1]*102, l_axis[2]*102] )
    # make arrow font size roughly proportional to radius of highest resolution shell
    fontsize = str(1.0 + roundoff(math.pow( max(self.miller_array.index_span().max()), 1.0/3.0)))

    arrowstr = """
    // xyz arrows
    shape.addSphere( [0,0,0] , [ 1, 1, 1 ], 0.3, 'Origo');
    //blue-x
    shape.addArrow( %s, %s , [ 0, 0, 1 ], 0.1);
    //green-y
    shape.addArrow( %s, %s , [ 0, 1, 0 ], 0.1);
    //red-z
    shape.addArrow( %s, %s , [ 1, 0, 0 ], 0.1);

    shape.addText( %s, [ 0, 0, 1 ], %s, 'h');
    shape.addText( %s, [ 0, 1, 0 ], %s, 'k');
    shape.addText( %s, [ 1, 0, 0 ], %s, 'l');
    """ %(str(Hstararrowstart), str(Hstararrowend), str(Kstararrowstart), str(Kstararrowend),
          str(Lstararrowstart), str(Lstararrowend), Hstararrowtxt, fontsize,
          Kstararrowtxt, fontsize, Lstararrowtxt, fontsize)

    # Make colour gradient array used for drawing a bar of colours next to associated values on the rendered html
    mincolourscalar = self.othermindata[self.icolourcol]
    maxcolourscalar = self.othermaxdata[self.icolourcol]
    if self.settings.sigma_color:
      mincolourscalar = self.otherminsigmas[self.icolourcol]
      maxcolourscalar = self.othermaxsigmas[self.icolourcol]
    span = maxcolourscalar - mincolourscalar
    ln = 60
    incr = span/ln
    colourgradarrays = []
    val = mincolourscalar
    colourscalararray =flex.double()
    colourscalararray.append( val )
    for j,sc in enumerate(range(ln)):
      val += incr
      colourscalararray.append( val )
    if self.otherscenes[self.icolourcol].miller_array.is_complex_array():
      # When displaying phases from map coefficients together with fom values
      # compute colour map chart as a function of fom and phase values (x,y axis)
      incr = 360.0/ln
      val = 0.0
      colourscalararray = flex.double()
      colourscalararray.append( val )
      for j in enumerate(range(ln)):
        val += incr
        colourscalararray.append( val )

      fomarrays = []
      if self.otherscenes[self.icolourcol].isUsingFOMs():
        fomln = 50
        fom = 1.0
        fomdecr = 1.0/(fomln-1.0)
      # make fomln fom arrays of size len(colourscalararray) when calling colour_by_phi_FOM
        for j in range(fomln):
          fomarrays.append( flex.double(len(colourscalararray), fom) )
          fom -= fomdecr
        for j in range(fomln):
          colourgradarrays.append( graphics_utils.colour_by_phi_FOM( colourscalararray*(math.pi/180.0), fomarrays[j] ) * 255.0)
      else:
        fomln =1
        fomarrays = [1.0]
        colourgradarrays.append( graphics_utils.colour_by_phi_FOM( colourscalararray*(math.pi/180.0) ) * 255.0)
    else:
      fomln = 1
      fomarrays = [1.0]
      colourgradarrays.append(graphics_utils.color_by_property(
        properties= flex.double(colourscalararray),
        selection=flex.bool( len(colourscalararray), True),
        color_all=False,
        gradient_type= self.settings.color_scheme) * 255.0)

    colors = self.otherscenes[self.icolourcol].colors
    radii = self.otherscenes[self.iradiicol].radii
    points = self.scene.points
    hkls = self.scene.indices
    dres = self.scene.dres
    colstr = self.scene.miller_array.info().label_string()
    data = self.scene.data
    colourlabel = self.otherscenes[self.icolourcol].colourlabel
    fomlabel = self.otherscenes[self.icolourcol].fomlabel
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    assert (colors.size() == radii.size() == nrefls)
    colours = []
    positions = []
    radii2 = []
    spbufttips = []

    self.workingbinvals = []
    if not self.binarray=="Resolution":
      ibinarray= int(self.binarray)
      self.workingbinvals = [ self.othermindata[ibinarray] - 0.1 , self.othermaxdata[ibinarray] + 0.1 ]
      self.workingbinvals.extend( self.binvals )
      self.workingbinvals.sort()
      if self.workingbinvals[0] < 0.0:
         self.workingbinvals.append(0.0)
         self.workingbinvals.sort()
      bindata = self.otherscenes[ibinarray].data
      if self.otherscenes[ibinarray].work_array.is_complex_array():
        bindata = self.otherscenes[ibinarray].ampl
    else:
      self.workingbinvals = self.binvals
      colstr = "dres"
      bindata = 1.0/dres
    self.nbin = len(self.workingbinvals)

    for ibin in range(self.nbin):
      colours.append([]) # colours and positions are 3 x size of data()
      positions.append([])
      radii2.append([])
      spbufttips.append([])

    def data2bin(d):
      for ibin, binval in enumerate(self.workingbinvals):
        if (ibin+1) == self.nbin:
          return ibin
        if d > binval and d <= self.workingbinvals[ibin+1]:
          return ibin
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      raise Sorry("Should never get here")

    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    for i, hklstars in enumerate(points):
      # bin currently displayed data according to the values of another miller array
      ibin = data2bin( bindata[i] )
      spbufttip = 'H,K,L: %s, %s, %s' %(hkls[i][0], hkls[i][1], hkls[i][2])
      spbufttip += '\ndres: %s ' %str(roundoff(dres[i])  )
      spbufttip += '\' + AA + \''
      for j,otherscene in enumerate(self.otherscenes):
        ocolstr = self.proc_arrays[j].info().label_string()
        odata = otherscene.data
        od =""
        if self.proc_arrays[j].is_complex_array():
          if not math.isnan(otherscene.foms[i]):
            od = str(roundoff(otherscene.ampl[i])) + ", " + str(roundoff(otherscene.phases[i])  ) + \
              "\' + DGR + \'" +  ", " + str(roundoff(otherscene.foms[i])  )
          else:
            od = str(roundoff(otherscene.ampl[i])) + ", " + str(roundoff(otherscene.phases[i])  ) + \
              "\' + DGR + \'"
        elif self.proc_arrays[j].sigmas() is not None:
          od = str(roundoff(odata[i]) ) + ", " + str(roundoff(otherscene.sigmas[i]))
        else:
          od = str(roundoff(odata[i]) )
        if not (math.isnan( abs(odata[i]) ) or odata[i] == display.inanval):
          spbufttip += "\n%s: %s" %(ocolstr, od)
      positions[ibin].extend( roundoff(list(hklstars)) )
      colours[ibin].extend( roundoff(list( colors[i] ), 2) )
      radii2[ibin].append( roundoff(radii[i], 2) )
      spbufttips[ibin].append(spbufttip)

    spherebufferstr = """
  ttips = new Array(%d)
  positions = new Array(%d)
  colours = new Array(%d)
  radii = new Array(%d)
  spherebufs = new Array(%d)
    """ %(self.nbin, self.nbin, self.nbin, self.nbin, self.nbin)

    negativeradiistr = ""
    cntbin = 0
    self.binstrs = []
    for ibin in range(self.nbin):
      mstr =""
      nreflsinbin = len(radii2[ibin])
      if (ibin+1) < self.nbin and nreflsinbin > 0:
        bin1= self.workingbinvals[ibin]
        bin2= self.workingbinvals[ibin+1]
        if colstr=="dres":
          bin1= 1.0/self.workingbinvals[ibin]
          bin2= 1.0/self.workingbinvals[ibin+1]
        mstr= "bin:%d, %d reflections with %s in ]%2.2f; %2.2f]" %(cntbin, nreflsinbin, \
                colstr, bin1, bin2)
        self.binstrs.append(mstr)
        self.mprint(mstr, verbose=True)
        spherebufferstr += """
  // %s
  ttips[%d] = %s
  positions[%d] = new Float32Array( %s )
  colours[%d] = new Float32Array( %s )
  radii[%d] = new Float32Array( %s )
  spherebufs[%d] = new NGL.SphereBuffer({
    position: positions[%d],
    color: colours[%d],
    radius: radii[%d],
    picking: ttips[%d],
  })
  //}, { disableImpostor: true }) // if true allows wireframe spheres but does not allow resizing spheres

  shape.addBuffer(spherebufs[%d])
      """ %(mstr, cntbin, str(spbufttips[ibin]).replace('\"', '\''), \
         cntbin, str(positions[ibin]), cntbin, str(colours[ibin]), \
         cntbin, str(radii2[ibin]), cntbin, cntbin, cntbin, cntbin, cntbin, cntbin )

        if self.workingbinvals[ibin] < 0.0:
          negativeradiistr += "spherebufs[%d].setParameters({metalness: 1})\n" %cntbin
        cntbin += 1

    spherebufferstr += """
// create tooltip element and add to the viewer canvas
  tooltip = document.createElement("div");
  Object.assign(tooltip.style, {
    display: "none",
    position: "absolute",
    zIndex: 10,
    pointerEvents: "none",
    backgroundColor: "rgba(255, 255, 255, 0.75)",
    color: "black",
    padding: "0.1em",
    fontFamily: "sans-serif"
  });

  stage.viewer.container.appendChild(tooltip);
  // listen to `hovered` signal to move tooltip around and change its text
  stage.signals.hovered.add(function (pickingProxy) {
    if (pickingProxy && (Object.prototype.toString.call(pickingProxy.picker) === '[object Array]'  )){
      var sphere = pickingProxy.sphere;
      var cp = pickingProxy.canvasPosition;
      tooltip.innerText = pickingProxy.picker[pickingProxy.pid];
      tooltip.style.bottom = cp.y + 7 + "px";
      tooltip.style.left = cp.x + 8 + "px";
      tooltip.style.fontSize = "smaller";
      tooltip.style.display = "block";
    }else{
      tooltip.style.display = "none";
    }
  });

  stage.signals.clicked.add(function (pickingProxy) {
  if (pickingProxy && (Object.prototype.toString.call(pickingProxy.picker) === '[object Array]'  )){
    var innerText = pickingProxy.picker[pickingProxy.pid];
    mysocket.send( innerText);
  }
  });

    """
    colourgradstrs = "colourgradvalarray = new Array(%s)\n" %fomln
    # if displaying phases from map coefficients together with fom values then
    for g,colourgradarray in enumerate(colourgradarrays):
      self.colourgradientvalues = []
      for j,e in enumerate(colourgradarray):
        self.colourgradientvalues.append( [colourscalararray[j], e] )
      self.colourgradientvalues = roundoff( self.colourgradientvalues )

      fom = fomarrays[g]
      colourgradstr = []
      for j,val in enumerate(self.colourgradientvalues):
        vstr = ""
        alpha = 1.0
        gradval = "rgba(%s, %s, %s, %s)" %(val[1][0], val[1][1], val[1][2], alpha)
        if j%10 == 0 or j==len(self.colourgradientvalues)-1 :
          vstr = str( roundoff(val[0], 2) )
        colourgradstr.append([vstr , gradval])

      colourgradstrs += "  colourgradvalarray[%s] = %s\n" %(g, str(colourgradstr) )

    #negativeradiistr = ""
    #for ibin in range(self.nbin):
    #  if self.workingbinvals[ibin] < 0.0:
    #    negativeradiistr += "spherebufs[%d].setParameters({metalness: 1})\n" %ibin



    self.NGLscriptstr = """
// Microsoft Edge users follow instructions on
// https://stackoverflow.com/questions/31772564/websocket-to-localhost-not-working-on-microsoft-edge
// to enable websocket connection

var pagename = location.pathname.substring(1);
var mysocket = new WebSocket('ws://127.0.0.1:7894/');

mysocket.onopen = function (e) {
  mysocket.send('%s now connected via websocket to ' + pagename + '\\n');
};

mysocket.onclose = function (e) {
  mysocket.send('%s now disconnecting from websocket ' + pagename + '\\n');
};

// Log errors to debugger of your browser
mysocket.onerror = function (error) {
  console.log('WebSocket Error ' + error);
};


window.addEventListener( 'resize', function( event ){
    stage.handleResize();
}, false );




var stage;
var shape;
var shapeComp;
var repr;
var AA = String.fromCharCode(197); // short for angstrom
var DGR = String.fromCharCode(176); // short for degree symbol

function createElement (name, properties, style) {
// utility function used in for loop over colourgradvalarray
  var el = document.createElement(name)
  Object.assign(el, properties)
  Object.assign(el.style, style)
  Object.assign(el.style, {
      display: "block",
      position: "absolute",
      color: "black",
      fontFamily: "sans-serif",
      fontSize: "smaller",
    }
  )
  return el
}


function addElement (el) {
// utility function used in for loop over colourgradvalarray
  Object.assign(el.style, {
    position: "absolute",
    zIndex: 10
  })
  stage.viewer.container.appendChild(el)
}


function addDivBox(txt, t, l, w, h, bgcolour='rgba(255.0, 255.0, 255.0, 0.0)') {
  divbox = createElement("div",
  {
    innerText: txt
  },
  {
    backgroundColor: bgcolour,
    color:  'rgba(0.0, 0.0, 0.0, 1.0)',
    top: t.toString() + "px",
    left: l.toString() + "px",
    width: w.toString() + "px",
    height: h.toString() + "px",
  }
  );
  addElement(divbox)
}


var hklscene = function () {
  shape = new NGL.Shape('shape');
  stage = new NGL.Stage('viewport', { backgroundColor: "grey", tooltip:false,
                                      fogNear: 100, fogFar: 100 });
  stage.setParameters( { cameraType: "%s" } );

  %s

  %s

  shapeComp = stage.addComponentFromObject(shape);
  repr = shapeComp.addRepresentation('buffer');
  shapeComp.autoView();
  repr.update()

  // if some radii are negative draw them with wireframe
  %s

  //colourgradvalarrays
  %s

  var ih = 3;
  var topr = 35
  var topr2 = 10
  var lp = 10
  var wp = 40
  var lp2 = lp + wp
  var gl = 3
  var wp2 = gl
  var fomlabelheight = 25
  if (colourgradvalarray.length === 1) {
    wp2 = 15
    fomlabelheight = 0
  }

  var wp3 = wp + colourgradvalarray.length * wp2 + 2

  totalheight = ih*colourgradvalarray[0].length + 35 + fomlabelheight
  // make a white box on top of which boxes with transparent background are placed
  // containing the colour values at regular intervals as well as label legend of
  // the displayed miller array
  addDivBox("", topr2, lp, wp3, totalheight, 'rgba(255.0, 255.0, 255.0, 1.0)');

  // print label of the miller array used for colouring
  addDivBox("%s", topr2, lp, wp, 20);

  if (colourgradvalarray.length > 1) {
    // print FOM label, 1, 0.5 and 0.0 values below colour chart
    fomtop = topr2 + totalheight - 18
    fomlp = lp + wp
    fomwp = wp3
    fomtop2 = fomtop - 13
    // print the 1 number
    addDivBox("1", fomtop2, fomlp, fomwp, 20)
    // print the 0.5 number
    leftp = fomlp + 0.48 * gl * colourgradvalarray.length
    addDivBox("0.5", fomtop2, leftp, fomwp, 20)
    // print the FOM label
    addDivBox("%s", fomtop, leftp, fomwp, 20);
    // print the 0 number
    leftp = fomlp + 0.96 * gl * colourgradvalarray.length
    addDivBox("0", fomtop2, leftp, fomwp, 20)
  }

  for (j = 0; j < colourgradvalarray[0].length; j++) {
    rgbcol = colourgradvalarray[0][j][1];
    val = colourgradvalarray[0][j][0]
    topv = j*ih + topr
    toptxt = topv - 5
    // print value of miller array if present in colourgradvalarray[0][j][0]
    addDivBox(val, toptxt, lp, wp, ih);
  }

  // draw the colour gradient
  for (g = 0; g < colourgradvalarray.length; g++) {
    leftp = g*gl + lp + wp

    // if FOM values are supplied draw colour gradients with decreasing
    // saturation values as stored in the colourgradvalarray[g] arrays
    for (j = 0; j < colourgradvalarray[g].length; j++) {
      rgbcol = colourgradvalarray[g][j][1];
      val = colourgradvalarray[g][j][0]
      topv = j*ih + topr
      addDivBox("", topv, leftp, wp2, ih, rgbcol);
    }
  }

}

document.addEventListener('DOMContentLoaded', function() { hklscene() }, false );


mysocket.onmessage = function (e) {
  //alert('Server: ' + e.data);
  var c
  var alpha
  var si
  mysocket.send('got ' + e.data ); // tell server what it sent us
  try {
    val = e.data.split(",")

    if (val[0] === "alpha") {
      ibin = parseInt(val[1])
      alpha = parseFloat(val[2])
      spherebufs[ibin].setParameters({opacity: alpha})
      stage.viewer.requestRender()
    }

    if (val[0] === "colour") {
      ibin = parseInt(val[1])
      si =  parseInt(val[2])
      colours[ibin][3*si] = parseFloat(val[3])
      colours[ibin][3*si+1] = parseFloat(val[4])
      colours[ibin][3*si+2] = parseFloat(val[5])
      spherebufs[ibin].setAttributes({ color: colours[ibin] })
      stage.viewer.requestRender()
    }

    if (val[0] === "Redraw") {
      stage.viewer.requestRender()
    }

    if (val[0] === "ReOrient") {
      mysocket.send( 'Reorienting ' + pagename );
      sm = new Float32Array(16);
      for (j=0; j<16; j++)
        sm[j] = parseFloat(val[j + 2]) // first 2 are "ReOrient", "NGL\\n"

      var m = new NGL.Matrix4();
      m.fromArray(sm);
      stage.viewerControls.orient(m);
      stage.viewer.requestRender();
    }

    if (val[0] === "Reload") {
    // refresh browser with the javascript file
      cvorient = stage.viewerControls.getOrientation().elements
      msg = String(cvorient)
      mysocket.send('Current vieworientation:\\n, ' + msg );

      mysocket.send( 'Refreshing ' + pagename );
      window.location.reload(true);
    }

    if (val[0] === "Testing") {
      // test something new
      mysocket.send( 'Testing something new ' + pagename );
      var newradii = radii[0].map(function(element) {
        return element*1.5;
      });

      spherebufs[0].setAttributes({
          radius: newradii
        })
      stage.viewer.requestRender()
    }

  }
  catch(err) {
    mysocket.send('error: ' + err );
  }
};




    """ % (self.__module__, self.__module__, self.cameratype, arrowstr, spherebufferstr, \
            negativeradiistr, colourgradstrs, colourlabel, fomlabel)
    if self.jscriptfname:
      with open( self.jscriptfname, "w") as f:
        f.write( self.NGLscriptstr )
    self.ReloadNGL()


  def OnConnectWebsocketClient(self, client, server):
    #if not self.websockclient:
    self.websockclient = client
    self.mprint( "Browser connected:" + str( self.websockclient ), self.verbose )
    #else:
    #  self.mprint( "Unexpected browser connection was rejected" )



  def OnWebsocketClientMessage(self, client, server, message):
    if message != "":
      self.mprint( message, verbose=False)
    self.lastmsg = message
    if "Current vieworientation:" in message:
      # The NGL.Matrix4 with the orientation is a list of floats.
      self.viewmtrxelms = message[ message.find("\n") : ]
      sleep(0.2)
      self.mprint( "Reorienting client after refresh:" + str( self.websockclient ) )
      if not self.isnewfile:
        self.pendingmessage = u"ReOrient, NGL" + self.viewmtrxelms
      self.isnewfile = False


  def WebBrowserMsgQueue(self):
    while True:
        sleep(1)
        if hasattr(self, "pendingmessage") and self.pendingmessage:
          self.SendWebSockMsg(self.pendingmessage)
          self.pendingmessage = None


  def StartWebsocket(self):
    self.server = WebsocketServer(7894, host='127.0.0.1')
    if not self.server:
      raise Sorry("Could not connect socket to web browser")
    self.server.set_fn_new_client(self.OnConnectWebsocketClient)
    self.server.set_fn_message_received(self.OnWebsocketClientMessage)
    self.wst = threading.Thread(target=self.server.run_forever)
    self.wst.daemon = True
    self.wst.start()
    self.msgqueuethrd = threading.Thread(target = self.WebBrowserMsgQueue )
    self.msgqueuethrd.daemon = True
    self.msgqueuethrd.start()


  def SendWebSockMsg(self, msg):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    #print "self.server.clients: ", self.server.clients
    #print "self.websockclient: ",
    if self.websockclient:
      while "Refreshing" in self.lastmsg:
        sleep(0.5)
        self.lastmsg = ""
      self.server.send_message(self.websockclient, msg )
    else:
      self.OpenBrowser()


  def SetOpacity(self, bin, alpha):
    if bin > self.nbin:
      self.mprint( "There are only %d bins of data present" %self.nbin )
      return
    msg = u"alpha, %d, %f" %(bin, alpha)
    self.SendWebSockMsg(msg)


  def RedrawNGL(self):
    self.SendWebSockMsg( u"Redraw, NGL\n" )


  def ReloadNGL(self): # expensive as javascript may be several Mbytes large
    self.mprint("Rendering NGL JavaScript...")
    self.SendWebSockMsg( u"Reload, NGL\n" )


  def OpenBrowser(self):
    NGLlibpath = libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","ngl.js") )
    htmlstr = self.hklhtml %(NGLlibpath, os.path.abspath( self.jscriptfname))
    htmlstr += self.htmldiv
    with open(self.hklfname, "w") as f:
      f.write( htmlstr )
    self.url = "file://" + os.path.abspath( self.hklfname )
    self.mprint( "Writing %s and connecting to its websocket client" %self.hklfname, self.verbose )
    if self.UseOSBrowser:
      webbrowser.open(self.url, new=1)
    self.isnewfile = False


  def TestNewFunction(self):
    self.SendWebSockMsg( u"Testing, NGL\n" )



"""
# python2 code

from websocket_server import WebsocketServer
import threading, math
from time import sleep

nc = {}
def new_client(client, server):
  nc = client
  print "got a new client:", nc

def on_message(client, server, message):
    print message

websocket.enableTrace(True)
server = WebsocketServer(7894, host='127.0.0.1')
server.set_fn_new_client(new_client)
server.set_fn_message_received(on_message)

wst = threading.Thread(target=server.run_forever)
wst.daemon = True
wst.start()

def LoopSendMessages():
  x = 0.0
  i=0
  while server.clients:
    nc = server.clients[0]
    x += 0.2
    alpha =  (math.cos(x) +1.0 )/2.0
    msg = u"alpha, 2, %f" %alpha
    server.send_message(server.clients[0], msg )
    r = (math.cos(x) +1.0 )/2.0
    g = (math.cos(x+1) +1.0 )/2.0
    b = (math.cos(x+2) +1.0 )/2.0
    msg = u"colour, 1, %d, %f, %f, %f" %(i,r,g,b)
    server.send_message(server.clients[0], msg )
    sleep(0.2)



"""


"""

# python3 code


import asyncio
import math
import websockets

async def time(websocket, path):
  x = 0
  for i in range(1000):
    x += 0.2
    alpha =  (math.cos(x) +1.0 )/2.0
    msg = u"alpha, 2, %f" %alpha
    await websocket.send( msg )
    r = (math.cos(x) +1.0 )/2.0
    g = (math.cos(x+1) +1.0 )/2.0
    b = (math.cos(x+2) +1.0 )/2.0
    msg = u"colour, 1, %d, %f, %f, %f" %(i,r,g,b)
    await websocket.send( msg )
    message = await websocket.recv()
    print( message)
    await asyncio.sleep(0.2)

start_server = websockets.serve(time, '127.0.0.1', 7894)
asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()



"""
