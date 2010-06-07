#ifndef SINGLEMASK_H
#define SINGLEMASK_H

#include <string>
#include <cmath>

#include <scitbx/error.h>//SCITBX_CHECK_POINT
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/constants.h>

#include <spotfinder/core_toolbox/libdistl.h>

namespace af = scitbx::af;
namespace spotfinder {
namespace distltbx {

inline int iround(double const& x)
    {if (x < 0) return static_cast<int>(x-0.5);
      return static_cast<int>(x+.5);}

typedef Distl::spot spot_t;
typedef af::shared<spot_t> spot_list_t;
typedef spot_t::point_list_t::const_iterator pt_iterator;
typedef std::vector<unsigned short> coord_list;
typedef coord_list::const_iterator c_iterator;

struct SingleMask {
  std::vector<short> xwidth;
  std::vector<short> ywidth;

  int x,y;

  inline SingleMask (spot_list_t master, af::shared<int> selection,
    int minimum_spot_count){
    //SCITBX_CHECK_POINT;
    //SCITBX_EXAMINE(master.size());
    //SCITBX_EXAMINE(selection.size());

    for (int si = 0; si< selection.size(); ++si) {
      //SCITBX_EXAMINE(selection[si]);
      std::vector<unsigned short> body_x_values;
      std::vector<unsigned short> body_y_values;
      spot_t* T = &master[selection[si]];
      pt_iterator B = T->bodypixels.begin();
      for (;B!=T->bodypixels.end(); ++B) {
        body_x_values.push_back(B->x);
        body_y_values.push_back(B->y);
        //SCITBX_EXAMINE(B->x);
        //SCITBX_EXAMINE(B->y);
      }
unsigned short xmin = *(std::min_element<c_iterator>(body_x_values.begin(),
                                                    body_x_values.end()));
unsigned short xmax = *(std::max_element<c_iterator>(body_x_values.begin(),
                                                    body_x_values.end()));
unsigned short ymin = *(std::min_element<c_iterator>(body_y_values.begin(),
                                                    body_y_values.end()));
unsigned short ymax = *(std::max_element<c_iterator>(body_y_values.begin(),
                                                    body_y_values.end()));
      xwidth.push_back(xmax - xmin + 3); //full width + 2 border pixels
      ywidth.push_back(ymax - ymin + 3); //full width + 2 border pixels
    }

    typedef std::vector<short>::iterator v_it;
    //for (v_it A= xwidth.begin(); A!=xwidth.end(); ++A){
    //  SCITBX_EXAMINE(*A);
    //}
    std::sort(xwidth.begin(),xwidth.end());
    std::sort(ywidth.begin(),ywidth.end());
    //for (v_it BB= xwidth.begin(); BB!=xwidth.end(); ++BB){
    //  SCITBX_EXAMINE(*BB);
    //}
    SCITBX_ASSERT(xwidth.size()> minimum_spot_count);
    int SampleSize = std::max(minimum_spot_count, int(xwidth.size()/10));
    double xmean = 0;
    double ymean = 0;
    v_it A= xwidth.end()-SampleSize;
    v_it B= ywidth.end()-SampleSize;
    for (; A!=xwidth.end(); ++A,++B){
      xmean += *A; ymean += *B;
    }
    //SCITBX_EXAMINE(xmean/SampleSize);
    //SCITBX_EXAMINE(ymean/SampleSize);
    x = iround(xmean/SampleSize);
    y = iround(ymean/SampleSize);
  }
};


} //namespace

} //namespace

#endif //SINGLEMASK_H
