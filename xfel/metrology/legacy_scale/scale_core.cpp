#include <xfel/metrology/legacy_scale/parameters.h>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>

xfel_legacy::parameter::parameter_array::parameter_array(
  int const& ndata,int const& kdim,farray d):
  parameters(d),
  gradients(ndata*kdim,scitbx::af::init_functor_null<double>()),
  curvatures(ndata*kdim,scitbx::af::init_functor_null<double>()),
  ndata(ndata),
  local_flag(false)
{
  SCITBX_ASSERT(d.size()==ndata*kdim);
}

xfel_legacy::parameter::parameter_array
xfel_legacy::parameter::organizer_base::register_array(
  std::string const& tag, int const& ndata, int const& kdim, farray d){
  _P[ tag ] = parameter_array(ndata,kdim,d);
  return _P[ tag ];
}

xfel_legacy::parameter::parameter_array
xfel_legacy::parameter::organizer_base::register_local_array(
  std::string const& tag, int const& ndata, int const& kdim, farray d){
  _P[ tag ] = parameter_array(ndata,kdim,d);
  _P[ tag ].local_flag=true;
  return _P[ tag ];
}

xfel_legacy::farray
xfel_legacy::parameter::organizer_base::as_x_array()const{
  farray result;
  typedef std::map<std::string,parameter_array>::const_iterator citer;
  for (citer M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    const double* bare_ptr = &(M->second.parameters.const_ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      result.push_back(bare_ptr[count]);
    }
  }
  return result;
}

xfel_legacy::farray
xfel_legacy::parameter::organizer_base::get_gradient_array()const{
  farray result;
  typedef std::map<std::string,parameter_array>::const_iterator citer;
  for (citer M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    const double* bare_ptr = &(M->second.gradients.const_ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      //std::cout<<M->first<<" "<<bare_ptr[count]<<std::endl;
      result.push_back(bare_ptr[count]);
    }
  }
  return result;
}

void
xfel_legacy::parameter::organizer_base::set_gradient_array(
  std::string const& aname, farray d)
{
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    if (aname == M->first) {
      double* bare_ptr = &(M->second.gradients.ref()[0]);

      SCITBX_ASSERT(d.size()==M->second.size());
      for (int ix = 0; ix < d.size(); ++ix){
        *bare_ptr++ = d[ix];
      }
    }
  }
}

xfel_legacy::farray
xfel_legacy::parameter::organizer_base::get_curvature_array()const{
  farray result;
  typedef std::map<std::string,parameter_array>::const_iterator citer;
  for (citer M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    const double* bare_ptr = &(M->second.curvatures.const_ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      result.push_back(bare_ptr[count]);
    }
  }
  // rest of this is for debugging only and can be removed.
  bool is_negative=false;
  for (int i = 0; i<result.size(); ++i){
    if (result[i] < 0.){is_negative=true;break;}
  }
  if (is_negative){
    for (citer M = _P.begin(); M != _P.end(); ++M){
      if (M->second.local_flag) {continue;}
      std::cout<<M->first<<std::endl;
      const double* bare_ptr = &(M->second.curvatures.const_ref()[0]);
      for (int count = 0; count < M->second.size(); ++count){
        std::cout<<"item "<<count<<" value "<<(bare_ptr[count])<<std::endl;
      }
    }
  }
  return result;
}

void
xfel_legacy::parameter::organizer_base::from_x_array(farray const& newx){
  farray result;
  const double* bare_x = &(newx.const_ref()[0]);
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    double* bare_ptr = &(M->second.parameters.ref()[0]);
    for (int count = 0; count < M->second.size(); ++count){
      *bare_ptr++ = *bare_x++;
    }
  }
}

void
xfel_legacy::parameter::organizer_base::initialize_gradients_curvatures(){
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    M->second.gradients = farray(M->second.parameters.size());
    M->second.curvatures = farray(M->second.parameters.size());
  }
}

void
xfel_legacy::parameter::organizer_base::rezero_gradients_curvatures(){
  typedef std::map<std::string,parameter_array>::iterator iter;
  for (iter M = _P.begin(); M != _P.end(); ++M){
    if (M->second.local_flag) {continue;}
    double* grad_ptr = &(M->second.gradients.ref()[0]);
    double* curv_ptr = &(M->second.curvatures.ref()[0]);
    for (int i = 0; i < M->second.gradients.size(); ++i){
      *grad_ptr++ = 0.;
      *curv_ptr++ = 0.;
    }
  }
}

void
xfel_legacy::algorithm::mark5_iteration::set_refined_origins_to_c(vec3array d){
  frame_origins = d;
}

double
xfel_legacy::algorithm::mark5_iteration::compute_target(
    farray tox, farray toy, farray spotcx, farray spotcy,
    farray spotfx, farray spotfy,
    iarray master_tiles,iarray frames,
    vec3array partial_r_partial_distance){

    model_calcx=farray(spotcx.size(),scitbx::af::init_functor_null<double>());
    model_calcy=farray(spotcx.size(),scitbx::af::init_functor_null<double>());

    int ntiles = _P["tile_trans"].ndata;
    SCITBX_ASSERT(tox.size()==ntiles);
    SCITBX_ASSERT(toy.size()==ntiles);

    functional = 0.;
    sine = std::vector<double>();
    cosine = std::vector<double>();
    farray rotations = _P["tile_rot"].parameters;
    for (int tidx=0; tidx < ntiles; ++tidx){
      cosine.push_back(std::cos(rotations[tidx]*scitbx::constants::pi_180));
      sine.push_back(std::sin(rotations[tidx]*scitbx::constants::pi_180));
    }

    const double* translations    = &(_P["tile_trans"].parameters.const_ref()[0]);
    double* translation_gradients = &(_P["tile_trans"].gradients.ref()[0]);
    double* rotation_gradients    = &(_P["tile_rot"].gradients.ref()[0]);
    double* frame_trans_gradients = &(_P["frame_trans"].gradients.ref()[0]);
    double* frame_dist_gradients  = &(_P["frame_distance"].gradients.ref()[0]);
    double* frame_rotz_gradients  = &(_P["frame_rotz"].gradients.ref()[0]);
    double* half_mos_gradients    = &(_P["half_mosaicity_rad"].gradients.ref()[0]);
    double* translation_curvatures= &(_P["tile_trans"].curvatures.ref()[0]);
    double* rotation_curvatures   = &(_P["tile_rot"].curvatures.ref()[0]);
    double* frame_trans_curvatures= &(_P["frame_trans"].curvatures.ref()[0]);
    double* frame_dist_curvatures = &(_P["frame_distance"].curvatures.ref()[0]);
    double* frame_rotz_curvatures = &(_P["frame_rotz"].curvatures.ref()[0]);
    double* half_mos_curvatures   = &(_P["half_mosaicity_rad"].curvatures.ref()[0]);

    int nframes = _P["frame_trans"].ndata;
    for (int ridx=0; ridx < spotcx.size(); ++ridx){
      int itile = master_tiles[ridx];
      int frame_param_no = frames[ridx];
      /*printf("ridx %5d frame %5d  partmos %8.5f %8.5f %8.5f \n",ridx, frame_param_no,
        vecc._V["half_mosaicity_rad"].gradients[ridx][0],
        vecc._V["half_mosaicity_rad"].gradients[ridx][1],
        vecc._V["half_mosaicity_rad"].gradients[ridx][2]

        );
      */
      calc_minus_To_x = spotcx[ridx] - tox[itile];
      calc_minus_To_y = spotcy[ridx] - toy[itile];

      rotated_o_x = calc_minus_To_x * cosine[itile]
                  - calc_minus_To_y * sine[itile];
      rotated_o_y = calc_minus_To_x * sine[itile]
                  + calc_minus_To_y * cosine[itile];

      model_calcx[ridx] = rotated_o_x + (tox[itile] + translations[2*itile]);
      model_calcy[ridx] = rotated_o_y + (toy[itile] + translations[2*itile+1]);

      partial_partial_theta_x = -calc_minus_To_x * sine[itile]
                               - calc_minus_To_y * cosine[itile];
      partial_partial_theta_y =  calc_minus_To_x * cosine[itile]
                               - calc_minus_To_y * sine[itile];

      partial_sq_theta_x = -calc_minus_To_x * cosine[itile]
                          + calc_minus_To_y * sine[itile];
      partial_sq_theta_y = -calc_minus_To_x * sine[itile]
                          - calc_minus_To_y * cosine[itile];

      double delx = model_calcx[ridx] - spotfx[ridx];
      double dely = model_calcy[ridx] - spotfy[ridx];
      double delrsq(delx*delx + dely*dely);
      functional += delrsq; // sum of square differences

      vec3 part_r_part_d = partial_r_partial_distance[ridx];
      double rotated_part_r_c_x = part_r_part_d[0] * cosine[itile]
                                - part_r_part_d[1] * sine[itile];
      double rotated_part_r_c_y = part_r_part_d[0] * sine[itile]
                                + part_r_part_d[1] * cosine[itile];

      translation_gradients[2*itile]  += 2. *  delx;
      translation_gradients[2*itile+1]+= 2. *  dely;

      if (frame_param_no < nframes){
        frame_trans_gradients[2*frame_param_no]  += 2. *  delx;
        frame_trans_gradients[2*frame_param_no+1]+= 2. *  dely;
        frame_dist_gradients[frame_param_no] += 2. * (
          delx * rotated_part_r_c_y + dely * rotated_part_r_c_x
        );//  SOMETHING IS ROTTEN IN DENMARK
      }

      rotation_gradients[itile] += scitbx::constants::pi_180 * 2.* (
        delx * partial_partial_theta_x +
        dely * partial_partial_theta_y
      );

      translation_curvatures[2*itile] += 2.;
      translation_curvatures[2*itile+1] += 2.;
      if (frame_param_no < nframes){
        frame_trans_curvatures[2*frame_param_no] += 2.;
        frame_trans_curvatures[1+2*frame_param_no] += 2.;
        frame_dist_curvatures[frame_param_no] +=2 * (
          rotated_part_r_c_x * rotated_part_r_c_x
        + rotated_part_r_c_y * rotated_part_r_c_y
        );
      }

      rotation_curvatures[itile] += scitbx::constants::pi_180 *
                                    scitbx::constants::pi_180 * 2. * (
        ( partial_partial_theta_x*partial_partial_theta_x +
          partial_partial_theta_y*partial_partial_theta_y ) +
        ( delx*partial_sq_theta_x + dely*partial_sq_theta_y )
      );

      // XXX make this map lookup more efficient
      if (_P["frame_rotz"].size()>0 && frame_param_no < nframes){
        //rotz angle is in radians; no need for pi_180 prefactor
        double origin_x = frame_origins[frame_param_no][0];
        double origin_y = frame_origins[frame_param_no][1];
        //
        double d_rcth_d_th_x = -(model_calcy[ridx] - origin_y);
        double d_rcth_d_th_y =   model_calcx[ridx] - origin_x;
        frame_rotz_gradients[frame_param_no] += 2. * (
          delx * d_rcth_d_th_x + dely * d_rcth_d_th_y );

        frame_rotz_curvatures[frame_param_no] += 2. * (
          ( delx * (-d_rcth_d_th_y) + dely * d_rcth_d_th_x) +
          ( d_rcth_d_th_x * d_rcth_d_th_x + d_rcth_d_th_y * d_rcth_d_th_y)
        );

      }

      if ( !_P["half_mosaicity_rad"].local_flag && _P["half_mosaicity_rad"].size()>0 &&
          frame_param_no < nframes){
        half_mos_gradients[frame_param_no] += 2. * (
          delx * vecc._V["half_mosaicity_rad"].gradients[ridx][1] +
          dely * vecc._V["half_mosaicity_rad"].gradients[ridx][0] );
          //  SOMETHING IS ROTTEN IN DENMARK (swap x & y for correct derivative)

        //first order approximation to curvature suggested by David Waterman
        half_mos_curvatures[frame_param_no] += 2. * (
          vecc._V["half_mosaicity_rad"].gradients[ridx].length_sq() // v.v dot product

          //These two lines are probably wrong.  Switch 0 & 1? vecc._V source incorrect?
          + (delx * vecc._V["half_mosaicity_rad"].curvatures[ridx][1])
          + (dely * vecc._V["half_mosaicity_rad"].curvatures[ridx][0])
        );

      }
    }
    return functional;
}

double
xfel_legacy::algorithm::mark5_iteration::compute_functional_only(
    farray tox, farray toy, farray spotcx, farray spotcy,
    farray spotfx, farray spotfy,
    iarray master_tiles,iarray frames,
    vec3array partial_r_partial_distance){

    model_calcx=farray(spotcx.size(),scitbx::af::init_functor_null<double>());
    model_calcy=farray(spotcx.size(),scitbx::af::init_functor_null<double>());

    int ntiles = _P["tile_trans"].ndata;
    SCITBX_ASSERT(tox.size()==ntiles);
    SCITBX_ASSERT(toy.size()==ntiles);

    functional = 0.;
    sine = std::vector<double>();
    cosine = std::vector<double>();
    farray rotations = _P["tile_rot"].parameters;
    for (int tidx=0; tidx < ntiles; ++tidx){
      cosine.push_back(std::cos(rotations[tidx]*scitbx::constants::pi_180));
      sine.push_back(std::sin(rotations[tidx]*scitbx::constants::pi_180));
    }

    const double* translations    = &(_P["tile_trans"].parameters.const_ref()[0]);

    for (int ridx=0; ridx < spotcx.size(); ++ridx){
      int itile = master_tiles[ridx];

      calc_minus_To_x = spotcx[ridx] - tox[itile];
      calc_minus_To_y = spotcy[ridx] - toy[itile];

      rotated_o_x = calc_minus_To_x * cosine[itile]
                  - calc_minus_To_y * sine[itile];
      rotated_o_y = calc_minus_To_x * sine[itile]
                  + calc_minus_To_y * cosine[itile];

      model_calcx[ridx] = rotated_o_x + (tox[itile] + translations[2*itile]);
      model_calcy[ridx] = rotated_o_y + (toy[itile] + translations[2*itile+1]);

      double delx = model_calcx[ridx] - spotfx[ridx];
      double dely = model_calcy[ridx] - spotfy[ridx];
      double delrsq(delx*delx + dely*dely);
      functional += delrsq; // sum of square differences
    }
    return functional;
}

xfel_legacy::vec3array
xfel_legacy::algorithm::mark5_iteration::uncorrected_detector_to_laboratory_frame(
    farray tox, farray toy,
    farray spotfx, farray spotfy,
    iarray master_tiles) const {

    /* purpose:  take observed positions on one frame as input.
       using currect detector tile metrology, convert these positions to
       the laboratory frame, i.e., correct for tile translations and
       rotations.
    */
    vec3array returnarray(spotfx.size(),scitbx::af::init_functor_null<vec3 >());

    std::map<std::string,parameter::parameter_array>::const_iterator findarray;
    findarray = _P.find("tile_trans");
    int ntiles = (findarray->second).ndata;
    SCITBX_ASSERT(tox.size()==ntiles);
    SCITBX_ASSERT(toy.size()==ntiles);
    const double* translations    = &((findarray->second).parameters.const_ref()[0]);

    std::vector<double> sine;
    std::vector<double> cosine;
    findarray = _P.find("tile_rot");
    farray rotations = (findarray->second).parameters;
    for (int tidx=0; tidx < ntiles; ++tidx){
      cosine.push_back(std::cos(rotations[tidx]*scitbx::constants::pi_180));
      sine.push_back(std::sin(rotations[tidx]*scitbx::constants::pi_180));
    }


    for (int ridx=0; ridx < spotfx.size(); ++ridx){
      int itile = master_tiles[ridx];

      // Explpicitly code a swap in X & Y here. Seems to be a result of legacy.

      double rotdx = spotfy[ridx] - translations[2*itile]   - tox[itile];
      double rotdy = spotfx[ridx] - translations[2*itile+1] - toy[itile];

      double unrotx = rotdx * cosine[itile]
                    + rotdy * sine[itile];
      double unroty = rotdx * -sine[itile]
                    + rotdy * cosine[itile];

      returnarray[ridx] = vec3(unrotx + tox[itile], unroty + toy[itile], 0.);

    }

    return returnarray;
}
