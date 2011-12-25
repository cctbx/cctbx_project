from spotfinder.array_family import flex
from spotfinder.applications.wrappers import DistlOrganizer

class Empty: pass

def run_signal_strength(params):
  E = Empty()
  E.argv=['Empty']
  E.argv.append(params.distl.image)
  return run_signal_strength_core(params,E)

def run_signal_strength_core(params,E):
  verbose = params.distl.verbose
  if params.distl.res.inner!=None:
    params.distl_lowres_limit = params.distl.res.inner
  if params.distl.res.outer!=None:
    params.force_method2_resolution_limit = params.distl.res.outer
    params.distl_highres_limit = params.distl.res.outer

  params.distl_force_binning = False
  params.distl_permit_binning = False
  params.wedgelimit = len(E.argv)
  params.spotfinder_header_tests = False
  Org = DistlOrganizer(verbose = True, argument_module=E,
                       phil_params=params)
  Org.printSpots()

  #Image analysis requested by NE-CAT (Point of contact: Craig Ogata)
  for key in Org.S.images.keys():
    # List of spots between specified high- and low-resolution limits
    if Org.S.images[key].has_key('lo_pass_resolution_spots'):
      spots = Org.S.images[key]['lo_pass_resolution_spots']
    elif Org.S.images[key].has_key('inlier_spots'):
      spots = Org.S.images[key]['inlier_spots']
    else:
      spots = []

    saturation = Org.Files.imageindex(key).saturation

    #Total number of spots in this range
    print
    print "Number of focus spots on image #%d within the input resolution range: %d"%(
      key,len(spots))

    signals=flex.double()
    saturations=flex.double()

    #Each spot
    for i,spot in enumerate(spots):
     signals.append(flex.sum(spot.wts))
     saturations.append(flex.max(spot.wts)/saturation)
     if verbose:
      #peak height given in ADC units above local background
      #integrated signal strength given in pixel-ADC units above local background
      print "%2d: Area in pixels=%d Peak=%.1f, Total integrated signal=%.1f (in pixel-ADC units above background)"%(
        i, spot.area(), flex.max(spot.wts), flex.sum(spot.wts))

      #peak signal-to-noise expressed in standard deviations above local background
      print "    Peak signal-to-noise=%.1f"%(spot.intensity())

      #peak height expressed in ADC units, without background subtraction
      image = Org.Files.imageindex(key)
      print "    Peak position x=%4d y=%4d (pixels); pixel value=%5d"%(
        spot.max_pxl_x(), spot.max_pxl_y(),
        image.linearintdata[(spot.max_pxl_x(),spot.max_pxl_y())])

      #Gory detail, looping through each pixel on each spot
      for j,pixel in enumerate(spot.bodypixels):
        print "       body pixel x=%4d y=%4d; pixel value=%5d; ADC height above background=%.1f"%(
          pixel.x,pixel.y,image.linearintdata[(pixel.x,pixel.y)],spot.wts[j])
    if signals.size()>0:
      print "Total integrated signal, pixel-ADC units above local background (just the good Bragg candidates) %d"%(
            flex.sum(flex.double([flex.sum(spot.wts) for spot in Org.S.images[key]['inlier_spots']]))
      )
      print "Signals range from %.1f to %.1f with mean integrated signal %.1f"%(
      flex.min(signals), flex.max(signals), flex.mean(signals) )
      print "Saturations range from %.1f%% to %.1f%% with mean saturation %.1f%%"%(
      100.*flex.min(saturations), 100.*flex.max(saturations), 100.*flex.mean(saturations) )

  if params.distl.pdf_output != None:
    #later, put this in a separate module so reportlab is not imported unless requested
    from labelit.publications.sublattice.sublattice_pdf import SublatticePDF,graphic
    from labelit.publications.sublattice.sublattice_pdf import PointTransform
    class genPDF(SublatticePDF):
      def make_image_plots_detail(self):
         params.pdf_output.window_fraction=1.0
         params.pdf_output.window_offset_x=0.0
         params.pdf_output.window_offset_y=0.0
         params.pdf_output.markup_inliers=True
         couple=(params.pdf_output.window_offset_x,
                 params.pdf_output.window_offset_y)
         #instead of self.R.setTransform, which requires pickled spotfinder:
         self.R.T = PointTransform()
         self.R.S = self.R.spotfinder
         self.R.T.setImage(spotfinder=self.R.S,subwindow_origin=couple,commands=params)
         self.R.title(self.image_name)
         #try:
         pil_image = graphic(filein = self.image_name,
                             couple = couple,
                             commands = params)
         self.R.image(pil_image)
         #except:
         #  print "failure, file %s"%self.filename
         if params.pdf_output.markup_inliers:
           self.R.show_ellipse(
           image_number=self.R.spotfinder.images.keys()[0],
           tags = ['goodspots','spots_non-ice','hi_pass_resolution_spots',
                    'spots_unimodal'],
           detail=True)
         self.R.c.showPage()
         return self
    pdf = genPDF(params.distl.pdf_output)
    pdf.filename = params.distl.pdf_output
    pdf.image_name = params.distl.image
    pdf.set_spotfinder(Org.S)
    pdf.make_image_plots_detail()
  return Org
