import math
from scitbx.array_family import flex

# This analysis is meant for CXI CSPAD only

def spot_shape_verbose(rawdata,beam_center_pix,indexed_pairs,spotfinder_observations,
      distance_mm, mm_per_pixel, hkllist, unit_cell, wavelength_ang):
    """requires:
    rawdata -- a 2d flex.int() with raw data
    beam_center_pix -- a scitbx.matrix.col() with beam xy in pixels
    indexed_pairs -- custom data structure
    spotfinder_observations -- spotfinder results
    """

    #--------------------------------------------- work on spot shape -------------
    from scitbx.matrix import col
    plotradial = flex.double()
    plotazim = flex.double()
    plotresol = flex.double()
    bodyx = flex.int()
    bodyy = flex.int()

    domain_sizes = flex.double()
    rot_mosaicities_deg = flex.double()
    implied_mosaicities_deg = flex.double()
    bandpasses = flex.double()
    unit_cell_deltas = flex.double()
    for ipidx,item in enumerate(indexed_pairs):
      #Not sure if xbeam & ybeam need to be swapped; values are very similar!!!
      radial, azimuthal = spotfinder_observations[item["spot"]].get_radial_and_azimuthal_size(
        beam_center_pix[0], beam_center_pix[1])

      model_center = col((spotfinder_observations[item["spot"]].ctr_mass_x(),
                          spotfinder_observations[item["spot"]].ctr_mass_y()))
      spot_vec_px = model_center - beam_center_pix
      spot_rvec = spot_vec_px.normalize()
      spot_avec = spot_rvec.rotate_2d(angle=90.,deg=True)

      from rstbx.new_horizons.pixel_spread import fwhm_2d_response
      RADIAL = fwhm_2d_response(rawdata,projection_vector = spot_rvec,spotfinder_spot = spotfinder_observations[item["spot"]])
      AZIMUT = fwhm_2d_response(rawdata,projection_vector = spot_avec,spotfinder_spot = spotfinder_observations[item["spot"]])


      overloaded = False
      for point in spotfinder_observations[item["spot"]].bodypixels:
        pixel_value = rawdata[(point.x,point.y)]
        if pixel_value > 5000.: #for XFEL only
          overloaded=True
      if not overloaded:
        # print out some properties:
        # first, simple model where all photon energy is assigned to the middle of the pixel
        print "%5d Johan radial %6.2f px, Johan azimuthal %6.2f px"%(ipidx, radial, azimuthal),
        Miller = hkllist[item["pred"]]
        resolution = unit_cell.d(Miller)
        print "resolution %6.2f Miller %s"%(resolution,str(Miller))

        fwhm_azi = AZIMUT.fwhm_pix()
        #small angle approximation
        rotational_mosaicity_deg = (fwhm_azi/spot_vec_px.length())*180./math.pi

        crystal_to_spot_ray_pix = math.hypot(distance_mm/mm_per_pixel,spot_vec_px.length())
        scherrer_fwhm_rad = fwhm_azi/crystal_to_spot_ray_pix
        scherrer_fwhm_deg = (scherrer_fwhm_rad)*180./math.pi
        if scherrer_fwhm_deg > 0.:
          scherrer_domain_size = wavelength_ang / (fwhm_azi/crystal_to_spot_ray_pix)
          domain_sizes.append(scherrer_domain_size)
        else:
          scherrer_domain_size = 0.

        print "      FWHM  radial %6.2f px, FWHM  azimuthal %6.2f px"%(RADIAL.fwhm_pix(), fwhm_azi)

        print "      Scherrer diffracted ray divergence FWHM %7.3f deg. Domain size %8.0f Angstrom"%(
          scherrer_fwhm_deg,scherrer_domain_size)

        print "          isotropic rotational FWHM mosaicity %7.3f deg"%rotational_mosaicity_deg

        spot_vec_mm = spot_vec_px.length() * mm_per_pixel
        radial_hwhm_mm = (RADIAL.fwhm_pix() * mm_per_pixel)/2.
        two_theta_low = math.atan( (spot_vec_mm+radial_hwhm_mm)/distance_mm )
        two_theta_high = math.atan( (spot_vec_mm-radial_hwhm_mm)/distance_mm )
        mosaicity_fwhm_rad = 0.5 * ( two_theta_low - two_theta_high )

        Elow_Ehigh_ratio = math.sin( two_theta_high/2. ) / math.sin ( two_theta_low/2. )
        bandpass_fwhm = 1. - Elow_Ehigh_ratio

        print "     radial divergence implied FWHM mosaicity %7.3f deg bandpass FWHM = %8.5f"%(
          mosaicity_fwhm_rad * 180/math.pi, bandpass_fwhm)


        Delta_a_over_a = scherrer_fwhm_rad * resolution / wavelength_ang

        print "    diffracted ray divergence FWHM (delta a)/a %8.5f"%( Delta_a_over_a )

        if resolution < 4.0:
          rot_mosaicities_deg.append(rotational_mosaicity_deg)
          implied_mosaicities_deg.append(mosaicity_fwhm_rad * 180/math.pi)
          bandpasses.append(bandpass_fwhm)
          unit_cell_deltas.append(Delta_a_over_a)
        print

    if len(domain_sizes) > 10:
      print "  IMAGE      Scherrer domain size, lower bound = %8.0F Angstrom"%flex.min(domain_sizes)
    if len(rot_mosaicities_deg) > 10:
      print "  IMAGE rotational FWHM mosaicity, upper bound = %7.3f deg. mean = %7.3f deg."%(
        flex.max(rot_mosaicities_deg), flex.mean(rot_mosaicities_deg))
      print "  IMAGE    implied FWHM mosaicity, upper bound = %7.3f deg. mean = %7.3f deg."%(
        flex.max(implied_mosaicities_deg), flex.mean(implied_mosaicities_deg))
      print "  IMAGE             FWHM bandpass, upper bound =  %8.5f"%(
        flex.max(bandpasses))
      print "  IMAGE   azimuthal FWHM deltaa/a, upper bound =  %8.5f from radial: %8.5f"%(
        flex.max(unit_cell_deltas), flex.max(implied_mosaicities_deg)*math.pi/180.)


    #-------------------- done with spot shape
