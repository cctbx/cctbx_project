#include <iostream>
#include <spotfinder/core_toolbox/libdistl.h>
#include <scitbx/vec2.h>

using namespace Distl;

void
spot_base::find_weighted_center(image_rawdata_t const& pixelvalue,
                                flag_array_t const& pixelismaxima,
                                float_array_t const& pixellocalmean){
  std::cout<<
  "spot_base::find_weighted_center DEPRECATION WARNING--contact authors"<<
  std::endl;
}

void
spot::find_weighted_center(image_rawdata_t const& pixelvalue,
                           flag_array_t const& pixelismaxima,
                           float_array_t const& pixellocalmean){

  scitbx::af::shared<scitbx::vec2<double> > pts;
  for (point_list_t::const_iterator q=bodypixels.begin();
       q!=bodypixels.end(); ++q) {
    pts.push_back(scitbx::vec2<double>(q->x, q->y));
    wts.push_back(static_cast<double>(pixelvalue[q->x][q->y])
                 - pixellocalmean[q->x][q->y]);
    bkg.push_back(pixellocalmean[q->x][q->y]);
  }
  model_ellipse(pts.const_ref(),wts.const_ref());
  //show_summary(pixelvalue,pixellocalmean);
}

void
spot::show_summary(image_rawdata_t const& pixelvalue,
                   double_array_t const& pixellocalmean) {
  std::cout<<"This spot: peakx "<<peak.x<<" peaky "<<peak.y<<std::endl;
  for (spot::point_list_t::const_iterator q=bodypixels.begin();
                     q!=bodypixels.end(); q++) {
    std::cout<<"  body x "<<q->x<<" y "<<q->y<<" height "<<
    (double)pixelvalue[q->x][q->y] - pixellocalmean[q->x][q->y]<<std::endl;
  }
  show_axes();
}
