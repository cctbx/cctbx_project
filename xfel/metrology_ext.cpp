#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/constants.h>
#include <scitbx/math/mean_and_variance.h>
#include <scitbx/vec3.h>
#include <cctbx/miller.h>

#include <vector>
#include <map>

using namespace boost::python;

namespace xfel {

struct mark2_iteration {
  typedef scitbx::af::shared<double> farray;
  typedef scitbx::af::shared<int> iarray;
  farray values;
  farray tox;
  farray toy;
  farray spotcx;
  farray spotcy;
  farray spotfx;
  farray spotfy;
  iarray master_tiles;
  double functional;
  farray gradients_,curvatures_;

  farray model_calcx, model_calcy;
  double calc_minus_To_x, calc_minus_To_y;
  double rotated_o_x, rotated_o_y;
  double partial_partial_theta_x, partial_partial_theta_y;
  double partial_sq_theta_x, partial_sq_theta_y;

  farray sine,cosine;

  mark2_iteration(){}
  mark2_iteration(farray values, farray tox, farray toy, farray spotcx, farray spotcy,
                  farray spotfx, farray spotfy,
                  iarray master_tiles):
    values(values),tox(tox),toy(toy),spotcx(spotcx),spotcy(spotcy),
    master_tiles(master_tiles),
    model_calcx(spotcx.size(),scitbx::af::init_functor_null<double>()),
    model_calcy(spotcx.size(),scitbx::af::init_functor_null<double>()),
    functional(0.),
    gradients_(3*64),
    curvatures_(3*64)
  {
    SCITBX_ASSERT(tox.size()==64);
    SCITBX_ASSERT(toy.size()==64);
    SCITBX_ASSERT(values.size()==3*64);
    for (int tidx=0; tidx < 64; ++tidx){
      cosine.push_back(std::cos(values[128+tidx]*scitbx::constants::pi_180));
      sine.push_back(std::sin(values[128+tidx]*scitbx::constants::pi_180));
    }

    for (int ridx=0; ridx < spotcx.size(); ++ridx){
      int itile = master_tiles[ridx];
      calc_minus_To_x = spotcx[ridx] - tox[itile];
      calc_minus_To_y = spotcy[ridx] - toy[itile];

      rotated_o_x = calc_minus_To_x * cosine[itile] - calc_minus_To_y * sine[itile];
      rotated_o_y = calc_minus_To_x * sine[itile] +   calc_minus_To_y * cosine[itile];

      model_calcx[ridx] = rotated_o_x + (tox[itile] + values[2*itile]);
      model_calcy[ridx] = rotated_o_y + (toy[itile] + values[2*itile+1]);

      partial_partial_theta_x = -calc_minus_To_x * sine[itile] - calc_minus_To_y * cosine[itile];
      partial_partial_theta_y =  calc_minus_To_x * cosine[itile] - calc_minus_To_y * sine[itile];

      partial_sq_theta_x = -calc_minus_To_x * cosine[itile] + calc_minus_To_y * sine[itile];
      partial_sq_theta_y = -calc_minus_To_x * sine[itile] - calc_minus_To_y * cosine[itile];

      double delx = model_calcx[ridx] - spotfx[ridx];
      double dely = model_calcy[ridx] - spotfy[ridx];
      double delrsq(delx*delx + dely*dely);
      functional += delrsq; // sum of square differences

      gradients_[2*itile]  += 2. *  delx;
      gradients_[2*itile+1]+= 2. *  dely;

      gradients_[128+itile] += scitbx::constants::pi_180 * 2.* (
        delx * partial_partial_theta_x +
        dely * partial_partial_theta_y
      );

      curvatures_[2*itile] += 2.;
      curvatures_[2*itile+1] += 2.;

      curvatures_[128+itile] += scitbx::constants::pi_180 * scitbx::constants::pi_180 * 2. * (
        ( partial_partial_theta_x*partial_partial_theta_x +
          partial_partial_theta_y*partial_partial_theta_y ) +
        ( delx*partial_sq_theta_x + dely*partial_sq_theta_y )
      );
    }
  }
  mark2_iteration(farray values, farray tox, farray toy, farray spotcx, farray spotcy,
                  farray spotfx, farray spotfy,
                  iarray master_tiles,iarray frames,int const& nframes,
                  bool const& inheritable):
    values(values),tox(tox),toy(toy),spotcx(spotcx),spotcy(spotcy),
    master_tiles(master_tiles),
    model_calcx(spotcx.size(),scitbx::af::init_functor_null<double>()),
    model_calcy(spotcx.size(),scitbx::af::init_functor_null<double>()),
    functional(0.){}

  mark2_iteration(farray values, farray tox, farray toy, farray spotcx, farray spotcy,
                  farray spotfx, farray spotfy,
                  iarray master_tiles,iarray frames,int const& nframes):
    values(values),tox(tox),toy(toy),spotcx(spotcx),spotcy(spotcy),
    master_tiles(master_tiles),
    model_calcx(spotcx.size(),scitbx::af::init_functor_null<double>()),
    model_calcy(spotcx.size(),scitbx::af::init_functor_null<double>()),
    functional(0.),
    gradients_(3*64+2*nframes),
    curvatures_(3*64+2*nframes)
  {
    SCITBX_ASSERT(tox.size()==64);
    SCITBX_ASSERT(toy.size()==64);
    SCITBX_ASSERT(values.size()==3*64+2*nframes);
    for (int tidx=0; tidx < 64; ++tidx){
      cosine.push_back(std::cos(values[128+tidx]*scitbx::constants::pi_180));
      sine.push_back(std::sin(values[128+tidx]*scitbx::constants::pi_180));
    }

    for (int ridx=0; ridx < spotcx.size(); ++ridx){
      int itile = master_tiles[ridx];
      int frame_param_no = frames[ridx];
      calc_minus_To_x = spotcx[ridx] - tox[itile];
      calc_minus_To_y = spotcy[ridx] - toy[itile];

      rotated_o_x = calc_minus_To_x * cosine[itile] - calc_minus_To_y * sine[itile];
      rotated_o_y = calc_minus_To_x * sine[itile] +   calc_minus_To_y * cosine[itile];

      model_calcx[ridx] = rotated_o_x + (tox[itile] + values[2*itile]);
      model_calcy[ridx] = rotated_o_y + (toy[itile] + values[2*itile+1]);

      partial_partial_theta_x = -calc_minus_To_x * sine[itile] - calc_minus_To_y * cosine[itile];
      partial_partial_theta_y =  calc_minus_To_x * cosine[itile] - calc_minus_To_y * sine[itile];

      partial_sq_theta_x = -calc_minus_To_x * cosine[itile] + calc_minus_To_y * sine[itile];
      partial_sq_theta_y = -calc_minus_To_x * sine[itile] - calc_minus_To_y * cosine[itile];

      double delx = model_calcx[ridx] - spotfx[ridx];
      double dely = model_calcy[ridx] - spotfy[ridx];
      double delrsq(delx*delx + dely*dely);
      functional += delrsq; // sum of square differences

      gradients_[2*itile]  += 2. *  delx;
      gradients_[2*itile+1]+= 2. *  dely;
      if (frame_param_no < nframes){
        gradients_[192+2*frame_param_no]  += 2. *  delx;
        gradients_[193+2*frame_param_no]  += 2. *  dely;
      }

      gradients_[128+itile] += scitbx::constants::pi_180 * 2.* (
        delx * partial_partial_theta_x +
        dely * partial_partial_theta_y
      );

      curvatures_[2*itile] += 2.;
      curvatures_[2*itile+1] += 2.;
      if (frame_param_no < nframes){
        curvatures_[192+2*frame_param_no] += 2.;
        curvatures_[193+2*frame_param_no] += 2.;
      }

      curvatures_[128+itile] += scitbx::constants::pi_180 * scitbx::constants::pi_180 * 2. * (
        ( partial_partial_theta_x*partial_partial_theta_x +
          partial_partial_theta_y*partial_partial_theta_y ) +
        ( delx*partial_sq_theta_x + dely*partial_sq_theta_y )
      );
    }
  }

  double f(){ return functional; }
  farray gradients(){ return gradients_; }
  farray curvatures(){ return curvatures_; }
};

struct mark3_collect_data{
  //adapt all-frame data to individual-frame parameter refinement
  typedef scitbx::af::shared<double> farray;
  typedef scitbx::af::shared<int> iarray;
  typedef scitbx::af::shared<cctbx::miller::index<> > marray;
  typedef scitbx::af::shared<bool> barray;
  marray HKL;
  std::map<int,int> frame_first_index, frame_match_count;
  farray result_model_cx,result_model_cy;
  barray result_flags;
  scitbx::af::shared<scitbx::vec3<double> > result_part_distance;

  mark3_collect_data(){}
  mark3_collect_data(iarray frame_id, marray indices):
    HKL(indices),
    result_model_cx(frame_id.size(),scitbx::af::init_functor_null<double>()),
    result_model_cy(frame_id.size(),scitbx::af::init_functor_null<double>()),
    result_flags(frame_id.size(),scitbx::af::init_functor_null<bool>()),
    result_part_distance(frame_id.size(),scitbx::af::init_functor_null<scitbx::vec3<double> >())
  {
    SCITBX_ASSERT(frame_id.size()==indices.size());

    for (int idx=0; idx < frame_id.size(); ++idx){
      int iframe = frame_id[idx];
      if (frame_first_index.find(iframe)==frame_first_index.end()){
        frame_first_index[iframe]=idx;
        frame_match_count[iframe]=1;
      } else {
        SCITBX_ASSERT(
          frame_first_index[iframe]+frame_match_count[iframe] == idx);// each frame all contiguous
        frame_match_count[iframe]+=1;
      }
    }
  }

  int
  get_first_index(int const& frame_id) const {
    return frame_first_index.find(frame_id)->second;
  }

  marray
  frame_indices(int const& frame_id)const{
    marray result;
    int last = frame_first_index.find(frame_id)->second + frame_match_count.find(frame_id)->second;
    for (int idx=frame_first_index.find(frame_id)->second; idx < last; ++idx){
      result.push_back(HKL[idx]);
    }
    return result;
  }

  barray
  selection(int const& frame_id)const{
    barray result(result_model_cx.size());
    int last = frame_first_index.find(frame_id)->second + frame_match_count.find(frame_id)->second;
    for (int idx=frame_first_index.find(frame_id)->second; idx < last; ++idx){
      result[idx] = true;
    }
    return result;
  }

  void collect(scitbx::af::shared<scitbx::vec3<double> > hi_E_limit,
               scitbx::af::shared<scitbx::vec3<double> > lo_E_limit,
               barray observed_flag,
               int const& frame_id){
    int first = frame_first_index.find(frame_id)->second;
    for (int im = 0; im < hi_E_limit.size(); ++im){
      result_model_cx[first + im] = (hi_E_limit[im][1] + lo_E_limit[im][1])/2.;
      result_model_cy[first + im] = (hi_E_limit[im][0] + lo_E_limit[im][0])/2.;
      result_flags[first + im] = observed_flag[im];
      SCITBX_ASSERT (observed_flag[im]); // no current support for masked-out spots
    }
  }
  void collect_mean_position(scitbx::af::shared<scitbx::vec3<double> > mean_position,
               barray observed_flag,
               int const& frame_id){
    int first = frame_first_index.find(frame_id)->second;
    for (int im = 0; im < mean_position.size(); ++im){
      result_model_cx[first + im] = mean_position[im][1];
      result_model_cy[first + im] = mean_position[im][0];
      result_flags[first + im] = observed_flag[im];
      SCITBX_ASSERT (observed_flag[im]); // no current support for masked-out spots
    }
  }
  void collect_distance(scitbx::af::shared<scitbx::vec3<double> > part_distance,
               int const& frame_id){
    int first = frame_first_index.find(frame_id)->second;
    for (int im = 0; im < part_distance.size(); ++im){
      result_part_distance[first + im] = part_distance[im];
    }
  }
};

struct mark5_iteration: public mark2_iteration {
  typedef scitbx::af::shared<scitbx::vec3<double> > vec3array;
  typedef scitbx::vec3<double>                      vec3;

  mark5_iteration(farray values, farray tox, farray toy, farray spotcx, farray spotcy,
                  farray spotfx, farray spotfy,
                  iarray master_tiles,iarray frames,int const& nframes,
                  vec3array partial_r_partial_distance):
    mark2_iteration(values, tox, toy, spotcx, spotcy, spotfx, spotfy,
                    master_tiles, frames, nframes, true)
  {
    gradients_ = farray(3*64+3*nframes);
    curvatures_ = farray(3*64+3*nframes);
    SCITBX_ASSERT(tox.size()==64);
    SCITBX_ASSERT(toy.size()==64);
    SCITBX_ASSERT(values.size()==3*64+3*nframes);
    for (int tidx=0; tidx < 64; ++tidx){
      cosine.push_back(std::cos(values[128+tidx]*scitbx::constants::pi_180));
      sine.push_back(std::sin(values[128+tidx]*scitbx::constants::pi_180));
    }

    for (int ridx=0; ridx < spotcx.size(); ++ridx){
      int itile = master_tiles[ridx];
      int frame_param_no = frames[ridx];

      calc_minus_To_x = spotcx[ridx] - tox[itile];
      calc_minus_To_y = spotcy[ridx] - toy[itile];

      rotated_o_x = calc_minus_To_x * cosine[itile] - calc_minus_To_y * sine[itile];
      rotated_o_y = calc_minus_To_x * sine[itile] +   calc_minus_To_y * cosine[itile];

      model_calcx[ridx] = rotated_o_x + (tox[itile] + values[2*itile]);
      model_calcy[ridx] = rotated_o_y + (toy[itile] + values[2*itile+1]);

      partial_partial_theta_x = -calc_minus_To_x * sine[itile] - calc_minus_To_y * cosine[itile];
      partial_partial_theta_y =  calc_minus_To_x * cosine[itile] - calc_minus_To_y * sine[itile];

      partial_sq_theta_x = -calc_minus_To_x * cosine[itile] + calc_minus_To_y * sine[itile];
      partial_sq_theta_y = -calc_minus_To_x * sine[itile] - calc_minus_To_y * cosine[itile];

      double delx = model_calcx[ridx] - spotfx[ridx];
      double dely = model_calcy[ridx] - spotfy[ridx];
      double delrsq(delx*delx + dely*dely);
      functional += delrsq; // sum of square differences

      vec3 part_r_part_d = partial_r_partial_distance[ridx];
      double rotated_part_r_c_x = part_r_part_d[0] * cosine[itile] - part_r_part_d[1] * sine[itile];
      double rotated_part_r_c_y = part_r_part_d[0] * sine[itile] +   part_r_part_d[1] * cosine[itile];

      //if (frame_param_no==10){
      //        printf("%5d %8.4f %8.4f", itile, part_r_part_d[0], part_r_part_d[1]);
      //        printf("%5d %8.4f %8.4f\n", itile, rotated_part_r_c_x, rotated_part_r_c_y);
      //}

      gradients_[2*itile]  += 2. *  delx;
      gradients_[2*itile+1]+= 2. *  dely;
      if (frame_param_no < nframes){
        gradients_[192+2*frame_param_no]  += 2. *  delx;
        gradients_[193+2*frame_param_no]  += 2. *  dely;
        gradients_[192 + 2*nframes + frame_param_no] += 2. * (
          delx * rotated_part_r_c_y + dely * rotated_part_r_c_x
        );//  SOMETHING IS ROTTEN IN DENMARK
      }

      gradients_[128+itile] += scitbx::constants::pi_180 * 2.* (
        delx * partial_partial_theta_x +
        dely * partial_partial_theta_y
      );

      curvatures_[2*itile] += 2.;
      curvatures_[2*itile+1] += 2.;
      if (frame_param_no < nframes){
        curvatures_[192+2*frame_param_no] += 2.;
        curvatures_[193+2*frame_param_no] += 2.;
        curvatures_[192 + 2*nframes + frame_param_no] +=2 * (
          rotated_part_r_c_x * rotated_part_r_c_x + rotated_part_r_c_y * rotated_part_r_c_y
        );
      }

      curvatures_[128+itile] += scitbx::constants::pi_180 * scitbx::constants::pi_180 * 2. * (
        ( partial_partial_theta_x*partial_partial_theta_x +
          partial_partial_theta_y*partial_partial_theta_y ) +
        ( delx*partial_sq_theta_x + dely*partial_sq_theta_y )
      );
    }

  }
};

namespace boost_python { namespace {

  void
  metrology_init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;

    class_<mark2_iteration>("mark2_iteration",no_init)
      .def(init< >())
      .def(init<mark2_iteration::farray, mark2_iteration::farray, mark2_iteration::farray,
                mark2_iteration::farray, mark2_iteration::farray,
                mark2_iteration::farray, mark2_iteration::farray,
                mark2_iteration::iarray >(
        (arg_("values"),arg("tox"),arg_("toy"),arg_("spotcx"),arg_("spotcy"),
         arg_("spotfx"),arg_("spotfy"),
         arg_("master_tiles"))))
      .def(init<mark2_iteration::farray, mark2_iteration::farray, mark2_iteration::farray,
                mark2_iteration::farray, mark2_iteration::farray,
                mark2_iteration::farray, mark2_iteration::farray,
                mark2_iteration::iarray, mark2_iteration::iarray, int const& >(
        (arg_("values"),arg("tox"),arg_("toy"),arg_("spotcx"),arg_("spotcy"),
         arg_("spotfx"),arg_("spotfy"),
         arg_("master_tiles"), arg_("frames"), arg_("nframes"))))
      .def("f",&mark2_iteration::f)
      .def("gradients",&mark2_iteration::gradients)
      .def("curvatures",&mark2_iteration::curvatures)
      .add_property("model_calcx", make_getter(&mark2_iteration::model_calcx, rbv()))
      .add_property("model_calcy", make_getter(&mark2_iteration::model_calcy, rbv()))
    ;
    class_<mark3_collect_data>("mark3_collect_data",no_init)
      .def(init<mark3_collect_data::iarray,mark3_collect_data::marray>())
      .def("get_first_index",&mark3_collect_data::get_first_index)
      .def("frame_indices",&mark3_collect_data::frame_indices)
      .def("selection",&mark3_collect_data::selection)
      .def("collect",&mark3_collect_data::collect)
      .def("collect_mean_position",&mark3_collect_data::collect_mean_position)
      .add_property("cx", make_getter(&mark3_collect_data::result_model_cx, rbv()))
      .add_property("cy", make_getter(&mark3_collect_data::result_model_cy, rbv()))
      .add_property("flags", make_getter(&mark3_collect_data::result_flags, rbv()))
      .def("collect_distance",&mark3_collect_data::collect_distance)
      .add_property("part_distance", make_getter(&mark3_collect_data::result_part_distance, rbv()))
    ;
    class_<mark5_iteration, bases<mark2_iteration> >("mark5_iteration",no_init)
      .def(init<mark5_iteration::farray, mark5_iteration::farray, mark5_iteration::farray,
                mark5_iteration::farray, mark5_iteration::farray,
                mark5_iteration::farray, mark5_iteration::farray,
                mark5_iteration::iarray, mark5_iteration::iarray, int const&,
                mark5_iteration::vec3array >(
        (arg_("values"),arg("tox"),arg_("toy"),arg_("spotcx"),arg_("spotcy"),
         arg_("spotfx"),arg_("spotfy"),
         arg_("master_tiles"), arg_("frames"), arg_("nframes"), arg_("part_distance"))))
    ;
}
}}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xfel_metrology_ext)
{
  xfel::boost_python::metrology_init_module();

}
