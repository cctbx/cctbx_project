import math
from scitbx.array_family import flex
from scitbx.matrix import col

def get_model_ref_limits(self,raw_image,spotfinder,imageindex,inputai,
      spot_prediction_limiting_resolution):
    return 0. # for now, disable the special model refinement limits as they do not help

    body_pixel_reference = flex.double()
    delx=[]
    dely=[]
    tt=flex.double()
    d_ang = flex.double()
    for spot in spotfinder.images[imageindex]["goodspots"]:

      spot_center = col((spot.ctr_mass_x(),spot.ctr_mass_y(),0.0))*raw_image.pixel_size
      offs = self.detector_origin+spot_center
      radius_mm = math.hypot(offs[0],offs[1])
      spot_two_theta_rad = math.atan(radius_mm/inputai.getBase().distance)
      d_ang.append(  inputai.wavelength / (2.*math.sin (spot_two_theta_rad/2.)) )
      tt.append(spot_two_theta_rad)
      for pxl in spot.bodypixels:
        body_pixel_reference.append(pxl.y + 0.5)
        body_pixel_reference.append(pxl.x + 0.5)

    obs_limited_d = max ( d_ang.select(flex.sort_permutation(d_ang))[0],
                          spot_prediction_limiting_resolution )

    # Get optimal mosaicity to match the observations using uniform sampling method
    reference_mosaicity = inputai.getMosaicity()

    corr_best = None
    mos_best = None
    pred_two_theta_rad_best = None
    pred_full_set_best = None
    mos_min,mos_max = (0.001,1.5)
    mos_low = mos_min + 0.
    mos_high = mos_max + 0.
    steps = math.pow(mos_high/mos_low,1./20.)
    from bpcx_regression.use_case_1_2.compare_bpcx_xds import cc,cc_slope
    for i in xrange(20):
      mos_test = mos_low * math.pow(steps,i)



      #for the moment, assume image center = 0 (always true for xfel stills)
      inputai.setMosaicity(mos_test)
      pred = inputai.predict_all( 0.0,obs_limited_d )
      pred_two_theta_rad = inputai.getOrientation().unit_cell().two_theta(inputai.hklpredict(),
        wavelength=inputai.wavelength)

      obs_tt_sort = tt.select(flex.sort_permutation(tt))
      pred_tt_sort = pred_two_theta_rad.select(flex.sort_permutation(pred_two_theta_rad))
      print "len obs",len(obs_tt_sort)
      print "len pred",len(pred_tt_sort)
      if not len(pred_tt_sort)>len(obs_tt_sort):continue
      conversion = len(pred_tt_sort)/len(obs_tt_sort)
      select_pred_tt_sort = flex.double()
      for x in xrange(len(obs_tt_sort)):
        select_pred_tt_sort.append(pred_tt_sort[int(x * conversion)])
      print "len select pred",len(select_pred_tt_sort)
      corr_test = cc(obs_tt_sort, select_pred_tt_sort)
      print "test mosaicity",mos_test,"Correlation",corr_test
      if mos_best==None:
        mos_best=mos_test;corr_best=corr_test;pred_two_theta_rad_best=select_pred_tt_sort
        pred_full_set_best=pred_tt_sort
      if corr_test > corr_best:
        mos_best=mos_test;corr_best=corr_test;pred_two_theta_rad_best=select_pred_tt_sort
        pred_full_set_best=pred_tt_sort
    print "best mosaicity",mos_best

    # got best model mosaicity, now go back again and match up predictions to observations
    # get get prediction set with same length and two theta values as the observation set
    modified_pred_x=flex.double()
    modified_pred_y=flex.double()
    ipred = 0
    for iobs in xrange(len(obs_tt_sort)):
      try:
        while pred_full_set_best[ipred]<obs_tt_sort[iobs]: ipred+=1
      except IndexError:
        break
      modified_pred_x.append(ipred)
      modified_pred_y.append(obs_tt_sort[iobs])

    # now do some linear least squares fitting of the modified pred so it fits low-resolution obs
    #  x = low resolution obs_x
    #  y = low resolution modified_pred_x

    #N = int(len(obs_tt_sort)*0.75) # take the first half of obs
    #X = flex.double(xrange(N)).as_numpy_array()
    #Y = flex.double(modified_pred_x[:N]).as_numpy_array()
    #from labelit.diffraction.limits import linear_fit
    #K,B,corr = linear_fit(X,Y)
    #print K,B,corr

    ccc=[];cccx=[];
    for x in xrange(int(len(obs_tt_sort)/3),len(obs_tt_sort)):
      cccx.append(x)
      ccc.append(cc(xrange(x),modified_pred_x[0:x]))
      #print x,cc(xrange(x),modified_pred_x[0:x])
      corr,slope = cc_slope(xrange(x),modified_pred_x[0:x])
      scaled_modified_pred_x = modified_pred_x/slope

      second_modified_pred_x=flex.double()
      second_modified_pred_y=flex.double()
      ipred = 0
      for ixx in xrange(len(obs_tt_sort)):
        try:
          while scaled_modified_pred_x[ipred]<ixx: ipred+=1
        except IndexError:
          break
        second_modified_pred_x.append(ixx)
        second_modified_pred_y.append(modified_pred_y[ipred])

      diff = second_modified_pred_y - obs_tt_sort[0:len(second_modified_pred_y)]
      stats = flex.mean_and_variance( second_modified_pred_y - obs_tt_sort[0:len(second_modified_pred_y)] )

      rmsd = stats.unweighted_sample_standard_deviation()
      mean = stats.mean()
      print x,mean,rmsd,-diff[x]/rmsd
      if -diff[x]/rmsd > 0.5:  # ad hoc cutoff seems like reasonable criterion for resolution cutoff
        print "twotheta rad cutoff",obs_tt_sort[x]
        print "angstrom cutoff",inputai.wavelength / (2.*math.sin (  obs_tt_sort[x] /2.))
        break

    if False:
      from matplotlib import pyplot as plt
      #plt.plot(delx,dely,"r.")
      #plt.show()
      plt.plot(xrange(len(tt)),obs_tt_sort,"r.")
      plt.plot(flex.double(xrange(len(pred_two_theta_rad_best))),pred_two_theta_rad_best,"g.")
      #plt.plot(flex.double(xrange(len(pred_full_set_best))),pred_full_set_best,"b.")
      plt.plot(cccx, ccc,"b.")
      plt.plot(scaled_modified_pred_x, modified_pred_y,"b.")
      #plt.plot(modified_pred_x/B,modified_pred_y,"y.")
      plt.plot(second_modified_pred_x, second_modified_pred_y,"r+")
      plt.show()

    model_refinement_limiting_resolution = inputai.wavelength / (2.*math.sin (  obs_tt_sort[x] /2.))
    return model_refinement_limiting_resolution
