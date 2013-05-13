#ifndef XFEL_VECTOR_COLLECTION_H
#define XFEL_VECTOR_COLLECTION_H

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <cctbx/miller.h>

namespace xfel {
typedef scitbx::af::shared<double> farray;
typedef scitbx::af::shared<int>    iarray;
typedef scitbx::af::shared<scitbx::vec3<double> > vec3array;
typedef scitbx::vec3<double>                      vec3;
typedef scitbx::af::shared<cctbx::miller::index<> > marray;

namespace parameter {

struct vector_array {
  vec3array gradients;
  vec3array curvatures;
};

struct vector_collection{
  //Take vector gradients from all reflections on one frame;
  // add them to a list of gradients of all reflections on all frames.
  typedef scitbx::af::shared<cctbx::miller::index<> > marray;
  marray HKL;
  std::map<int,int> frame_first_index, frame_match_count;
  std::map<std::string,vector_array> _V;

  inline
  vector_collection(){}

  inline
  vector_collection(iarray frame_id, marray indices):
    HKL(indices)
  {
    SCITBX_ASSERT(frame_id.size()==indices.size());

    for (int idx=0; idx < frame_id.size(); ++idx){
      int iframe = frame_id[idx];
      if (frame_first_index.find(iframe)==frame_first_index.end()){
        frame_first_index[iframe]=idx;
        frame_match_count[iframe]=1;
      } else {
        SCITBX_ASSERT(
          frame_first_index[iframe]+frame_match_count[iframe] == idx);
          // each frame all contiguous
        frame_match_count[iframe]+=1;
      }
    }
  }

  inline
  vector_array
  register_tag(std::string const& tag){
    if (_V.find(tag) == _V.end()){
      _V[tag]=vector_array();
      _V[tag].gradients =
         vec3array(HKL.size(),scitbx::af::init_functor_null<vec3>());
      _V[tag].curvatures=
         vec3array(HKL.size(),scitbx::af::init_functor_null<vec3>());
      return _V[tag];
    }
    return vector_array();
  }

  inline
  void collect_vector_information( std::string const& tag,
    vec3array vector_input,
    int const& frame_id){
    SCITBX_ASSERT(_V.find(tag) != _V.end());

    int first = frame_first_index.find(frame_id)->second;
    for (int im = 0; im < vector_input.size(); ++im){
      // make this more efficient with pointers later
      _V[tag].gradients[first + im] = vector_input[im];
    }
  }
};
}
} //namespace xfel
#endif// XFEL_VECTOR_COLLECTION_H
