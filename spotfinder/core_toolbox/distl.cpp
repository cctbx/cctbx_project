#include <boost/tokenizer.hpp>
#include <spotfinder/core_toolbox/distl.h>
#include <scitbx/vec3.h>


namespace af = scitbx::af;
namespace di = spotfinder::distltbx;
//this function was developed but never used; deprecate soon
Distl::spot::spot(af::flex_int::const_iterator intptr) :
  peakintensity(0.0),
  nmaxima(0),
  p_gotaxes(false)
{
  peak = Distl::point(*intptr,*(intptr+1));
  intptr+=2;
  int bodylen = *intptr; intptr++;
  bodypixels.reserve(bodylen);
  for (int i=0; i<bodylen; ++i, intptr+=2){
    bodypixels.push_back(Distl::point(*intptr,*(intptr+1)));}
}

void
di::w_spot::p_getaxes() {
//DEPRECATE SOON---DEPRECATE SOON
  if (p_gotaxes) return;

  // iteration to find average position
  af::tiny<double, 2> average_pixel_position = af::tiny<double, 2>(0.,0.);
  for (int i = 0; i < bodypixels.size() ; ++i) {
    average_pixel_position[0]+=bodypixels[i].x;
    average_pixel_position[1]+=bodypixels[i].y;
  }
  average_pixel_position[0]/=bodypixels.size();
  average_pixel_position[1]/=bodypixels.size();

  //iteration to find major axis
  std::vector<double> sq_distance;
  for (int i = 0; i < bodypixels.size() ; ++i) {
    sq_distance.push_back(
     (bodypixels[i].x - average_pixel_position[0]) *
     (bodypixels[i].x - average_pixel_position[0]) +
     (bodypixels[i].y - average_pixel_position[1]) *
     (bodypixels[i].y - average_pixel_position[1]));
  }
  std::vector<double>::const_iterator ptrmajor =
    std::max_element(sq_distance.begin(),sq_distance.end());
  double c_majoraxis = 2.0*std::sqrt(*ptrmajor);
  size_t posmajor = ptrmajor - sq_distance.begin();
  //printf("posmajor %7d\n",posmajor);

  //double check to verify the major axis vector:
  double dx = bodypixels[posmajor].x - average_pixel_position[0];
  double dy = bodypixels[posmajor].y - average_pixel_position[1];
  //double d_majoraxis = 2.0*std::sqrt( dx*dx + dy*dy );
  //printf("C_major %7.1f   D_major %7.1f\n",c_majoraxis,d_majoraxis);

  //iteration to find minor axis
  std::vector<double> abs_cross_product;
  for (int i = 0; i < bodypixels.size() ; ++i) {
    abs_cross_product.push_back( std::fabs(
       dx * (bodypixels[i].y - average_pixel_position[1]) -
       (bodypixels[i].x - average_pixel_position[0]) * dy
    ) );
  }
  std::vector<double>::const_iterator ptrminor =
    std::max_element(abs_cross_product.begin(),abs_cross_product.end());
  double c_minoraxis = 4.0*(*ptrminor)/c_majoraxis;
  //printf("C_major %7.1f   C_minor %7.1f\n",c_majoraxis,c_minoraxis);

  if (c_majoraxis >= c_minoraxis) {
    p_majoraxis = c_majoraxis;
    p_minoraxis = c_minoraxis;
  } else {
    p_majoraxis = c_minoraxis;
    p_minoraxis = c_majoraxis;
  }

  p_gotaxes=true;
}

// major and minor axes given in units of pixels
double
di::w_spot::get_majoraxis() {
  p_getaxes();
  return p_majoraxis;
}

double
di::w_spot::get_minoraxis() {
  p_getaxes();
  return p_minoraxis;
}

void di::w_spot::setstate(Distl::point const& pk){
  nmaxima=0;
  p_gotaxes=false;
  peak=pk;
}


/*
Commentary on DISTL, relevant to the use of DISTL for virus images.

1) the function pxlclassify considers background boxes in sequence across
   the detector face, but then when it comes to the final edge, it
   calculates the final box boundaries differently.  This makes for clumsy
   code, in particular the pxlclassify_scanbox() is called in 4 places.
   Duplicated code, harder to maintain.  Better to change the function so
   it establishes a single loop with an adjusted-size box so it only
   needs to call scanbox once.

2) consecutively scanned boxes overlap in 2 pixels.  This was probably
   intentional to avoid edge-related discontinuties; but it also has an
   unwanted consequence.  The pixelintensity[x][y] value is both used and
   modified by pxlclassify_scanbox().  Therefore, for the two pixels that
   overlap with the previous box, we are using
   pixelintensity[x][y] values modified in this cycle, while all other pixels
   have values from the last cycle.  It would be cleaner to implement boxes
   without overlaps.

3) on the same topic, it would help avoid discontinuties to determine a
   background plane using International Tables, volume F, p. 213.

4) the as-published gamma-I values don't really make sense:
   a) in cycle 1, all initial pixelintensity[x][y] values are supposed to
      be zero (but see 2 above) so it doesn't mean anything to have a
      bgupperint for this cycle.
   b) in cycle 2, we remove some box pixels from the background calculation.
      Then in cycle 3, the published procedure increases the bgupperint
      cutoff and recalculates the mask, having the effect of adding some
      pixels back to the background.  But this seems wrong; once a pixel
      is removed from background (because it may have signal in it) it
      probably shouldn't be added back.  So the bgupperint cutoff should
      remain constant across cycles.

5) regarding #3, for the HK97 image (with many spots) we are spending an
   enormous amount of time increasing the background-box size and retesting
   it to attain 2/3 of spots in background.  It is possible that by
   applying the background plane correction, we will attain the 2/3 threshhold
   much sooner.  No--just a modest improvement when the plane correction is
   applied to the pixelvalue[x][y].  Second try: the plane correction must
   also be applied when calculating the boxstd(); otherwise the estimate
   of standard deviation will be generally too high.  No:  no improvement in
   the occurrence of background-box increases (even though the sd is lower).
   There are dramatically more spots chosen around the water ring.

   Actions:
   -add a DetectorImageBase.debug_write() command to the iotbx.
   -add a w_Distl.mod_data() function to get a custom-written modified dataarray

6) Conclusions:
   - the get_underload() function is at fault; it is choosing a lower cutoff
   of 4035, which is much too high for reasonable interpretation of the
   virus diffraction.  Most true background is masked out.
   Change the code so it does a sanity check, and never masks out more than
   10% of the pixels as underloads (and thus ~90% of pixels are potential
   background pixels).

7) tackle some performance issues.  the function spotlist_to_flexint() takes
   as much CPU as the libdistl calls.  Solution: alleviate by coding in C++.

8) In tnear, there are numerous boost.python function calls for each spot,
   allowing for spot filtering.  This becomes prohibitive for large unit
   cells.  Solution: recode spot filtering tests in C++ code.
   Also: implement the SpotManager to give a more general way of subsetting
   the spot list.

   Tuesday goals:
     tweak parameters for this image--bg cutoff to 1.5, d1 to 2.5
       & see where I am.
     still don't understand why the spots bleed into each other when bg
       cutoff is at 1.5.  I thought the criteria for growing spots had to
       do with d1, not bg2.
*/

di::w_Distl::w_Distl(std::string optionstring, bool report_overloads){
/*
"       -s2     Smallest acceptible spot area. Spots with area smaller than this value are ignored.\n\n"
"               Recommended values: [3,5]\n\n"
"       -s3     Spot base area.\n\n"
"               Summary of spot shape, strength, etc. are based on spots no smaller than this size.\n\n"
"       -s7     Spot area upper bound factor. If area exceeds\n"
"                   median + (95th prctile - 5th prctile) * factor,\n"
"               the spot is eliminated.\n\n"
"               Recommended value: [2, 5]\n\n"
"       -s8     Spot peak intensity upper bound factor. If peak intensity exceeds \n"
"                   median + (95th prctile - 5th prctile) * factor,\n"
"               the spot is eliminated.\n\n"
"               Recommended value: [2, 10]\n\n"
"   Ice-Ring Detection Parameters:\n\n"
"       -i3     Intensity percentile as a measure of ice-ring strength.\n\n"
"                               Recommended value: [0.1, 0.3]\n\n"
        -d1     Diffraction lower intensity; a lower bound for finding maxima, default 3.5
    More parameters added Feb. 2006 for virus work
        -bx0,1,2 Scanboxsize integer values for cycles 1,2,and 3.
        -bg0,1,2 Bgupperint cutoff values for cycles 1,2,and 3.
*/
    finder.spotbasesize = 10; // distl initializes it as 16, but the LABELIT default is 10.
    finder.bgupperint[0] = 1.5; //See note 4a above; value shouldn't be used
    finder.bgupperint[1] = 1.5; //See note 4b above
    finder.bgupperint[2] = 1.5; //See note 4b above
    if (report_overloads) {finder.report_overloads = true;}
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    boost::char_separator<char> sep(" ");
    tokenizer tok(optionstring,sep);
    for (tokenizer::iterator tok_iter = tok.begin();
         tok_iter!=tok.end(); ++tok_iter){
      if        (*tok_iter == "-s2") {
        finder.spotarealowcut = atoi((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-s3") {
        finder.spotbasesize = atoi((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-s7") {
        finder.spotareamaxfactor = atof((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-s8") {
        finder.spotpeakintmaxfactor = atof((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-i3") {
        finder.icering_strengthprctile = atof((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-d1") {
        finder.difflowerint = atof((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-bx0") {
        finder.scanboxsize[0] = atoi((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-bx1") {
        finder.scanboxsize[1] = atoi((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-bx2") {
        finder.scanboxsize[2] = atoi((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-bg0") {
        finder.bgupperint[0] = atof((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-bg1") {
        finder.bgupperint[1] = atof((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-bg2") {
        finder.bgupperint[2] = atof((*(++tok_iter)).c_str());
      } else if (*tok_iter == "-ro") {
        set_resolution_outer( atof((*(++tok_iter)).c_str()) );
      }
    }
}

void
di::w_Distl::set_resolution_outer(const double& newvalue)
{
  // only meaningful if the resolution value is a positive number.
  SCITBX_ASSERT(newvalue>0.0);
  finder.resolution_outer = newvalue;
  //SCITBX_EXAMINE(finder.resolution_outer);
}

void
di::w_Distl::setspotimg(const double& pixel_size,
                        const double& distance,
                        const double& wavelength,
                        const double& beamx,
                        const double& beamy,
                        af::flex_int const& intdata,
                        const int& peripheral_margin,
                        const double& saturation )
{
  finder.overloadvalue = saturation;
  //Take default angle settings because we don't use them anyway.
  finder.set_imageheader(pixel_size, distance, wavelength, 0., 0.,
                         beamx, beamy);

  int ncols = intdata.accessor().all()[0];
  int nrows = intdata.accessor().all()[1];

  finder.imgmargin = peripheral_margin;
  finder.set_imagedata(intdata.begin(),ncols,nrows);
}

void
di::w_Distl::set_tiling(const string& vendortype)
{
  if (vendortype=="Pilatus-6M") {
    finder.tiling = Distl::ptr_tiling(new Distl::scanbox_tiling_pilatus6M(
      finder.firstx, finder.lastx, finder.firsty, finder.lasty));
  } else if (vendortype=="Pilatus-2M") {
    finder.tiling = Distl::ptr_tiling(new Distl::scanbox_tiling_pilatus2M(
      finder.firstx, finder.lastx, finder.firsty, finder.lasty));
  } else {
    finder.tiling = Distl::ptr_tiling(new Distl::scanbox_tiling(
      finder.firstx, finder.lastx, finder.firsty, finder.lasty));
  }
}

void
di::w_Distl::set_tiling(af::flex_int const& explicit_tiling,int const& peripheral_margin)
{
  finder.tiling = Distl::ptr_tiling(new Distl::scanbox_tiling_explicit(
      explicit_tiling, peripheral_margin));
}

af::flex_double
di::w_Distl::Z_data()
{
  int nrows=finder.pixelvalue.ny;
  af::flex_double z(af::flex_grid<>(finder.pixelvalue.nx,nrows));

  double* begin = z.begin();

  for (int x=0; x<finder.pixelvalue.nx; x++) {
    for (int y=0; y<nrows; y++){
        //capture DISTL's Z-function;
        *begin++ = finder.pixelintensity[x][y];
    }
  }
  return z;
}

af::flex_int
di::w_Distl::mod_data()
{
  int nrows=finder.pixelvalue.ny;
  af::flex_int z(af::flex_grid<>(finder.pixelvalue.nx,nrows));

  int* begin = z.begin();

  for (int x=0; x<finder.pixelvalue.nx; x++) {
    for (int y=0; y<nrows; y++){
        //example mod_data function; pixel value conditional on some property
        if (finder.pixelintensity[x][y]<1.5 &&
        finder.pixelvalue[x][y]>finder.underloadvalue &&
        finder.pixelvalue[x][y]<finder.overloadvalue
        ){
            *begin = 10;
        }else{
            *begin = finder.pixelvalue[x][y];
        }
        begin++;
    }
  }
  return z;
}

void
di::w_Distl::finish_analysis(){

  // Now done as python calls to individual methods-->finder.process();

  //get the spots out of std::list and in to af::shared

  spots.reserve(finder.spots.size());
  list<Distl::spot>::const_iterator position = finder.spots.begin();
  list<Distl::spot>::const_iterator end = finder.spots.end();
  for (; position!=end;++position ) {

    //It turns out that minimum spot area is an extremely important
    //  filter without which all sorts of junk are reported.
    if ((*position).area() > finder.spotbasesize) {
      spots.push_back(*position);
    }
  }

  icerings.reserve(finder.icerings.size());
  vector<Distl::icering>::const_iterator iposition = finder.icerings.begin();
  vector<Distl::icering>::const_iterator iend = finder.icerings.end();

  for (; iposition!=iend;++iposition ) {
    icerings.push_back(*iposition);
  }

}

bool
di::w_Distl::isIsolated(const w_spot& spot, const double& mmradius) const {
  //Future:  1) find cases that use this algorithm
  //         2) change the algorithm so it relies on elliptical modelling
  //            rather than on the borderpixels data structure
  double pixelradius = mmradius / finder.pixel_size;
  typedef scitbx::vec2<double>                 vpoint;
  af::shared<Distl::point>::const_iterator sptr;
  af::shared<Distl::point>::const_iterator send;

  //vpoint focusspot(spot.ctr_mass_x(),spot.ctr_mass_y());
  vpoint focusspot(spot.max_pxl_x(),spot.max_pxl_y());

  spot_list_t::const_iterator p = spots.begin();
  spot_list_t::const_iterator e = spots.end();
  for (; p!=e; ++p) {
      //vpoint target((*p).ctr_mass_x(),(*p).ctr_mass_y());
      vpoint target((*p).max_pxl_x(),(*p).max_pxl_y());
      vpoint diff = target - focusspot;
      double dist = std::sqrt(diff*diff);
      if ( dist > pixelradius ) {continue;}

      // Consider border pixels of the focus spot
      vpoint bisector = 0.45*diff;
      double bisectorsq = bisector*bisector;
      sptr = spot.borderpixels.begin();
      send = spot.borderpixels.end();
      for (;sptr!=send; ++sptr) {
        vpoint borderpt((*sptr).x,(*sptr).y);
        vpoint bordervec = borderpt - focusspot;
        if ( bordervec*bisector > bisectorsq ) {return false;}
      }

      // Consider border pixels of the target spot
      bisector = 0.55*diff;
      bisectorsq = bisector*bisector;
      sptr = (*p).borderpixels.begin();
      send = (*p).borderpixels.end();
      for (;sptr!=send; ++sptr) {
        vpoint borderpt((*sptr).x,(*sptr).y);
        vpoint bordervec = borderpt - focusspot;
        if ( bordervec*bisector < bisectorsq ) {return false;}
      }
  }
  return true;
}
