/************************************************************************
                        Copyright 2003
                              by
                 The Board of Trustees of the
               Leland Stanford Junior University
                      All rights reserved.
                       Disclaimer Notice
     The items furnished herewith were developed under the sponsorship
 of the U.S. Government.  Neither the U.S., nor the U.S. D.O.E., nor the
 Leland Stanford Junior University, nor their employees, makes any war-
 ranty, express or implied, or assumes any liability or responsibility
 for accuracy, completeness or usefulness of any information, apparatus,
 product or process disclosed, or represents that its use will not in-
 fringe privately-owned rights.  Mention of any product, its manufactur-
 er, or suppliers shall not, nor is it intended to, imply approval, dis-
 approval, or fitness for any particular use.  The U.S. and the Univer-
 sity at all times retain the right to use and disseminate the furnished
 items for any purpose whatsoever.                       Notice 91 02 01
   Work supported by the U.S. Department of Energy under contract
   DE-AC03-76SF00515; and the National Institutes of Health, National
   Center for Research Resources, grant 2P41RR01209.
                       Permission Notice
 Permission is hereby granted, free of charge, to any person obtaining a
 copy of this software and associated documentation files (the "Software"),
 to deal in the Software without restriction, including without limitation
 the rights to use, copy, modify, merge, publish, distribute, sublicense,
 and/or sell copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following conditions:
 The above copyright notice and this permission notice shall be included
 in all copies or substantial portions of the Software.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTA-
 BILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
 EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 THE USE OR OTHER DEALINGS IN THE SOFTWARE.
************************************************************************/

/*
 * This library is used to analyze images of crystal diffraction patterns.
 * Major functions include:
 *  1) Calculate an "intensity" score for each pixel,
 *     based on the characteristics if local background.
 *  2) Mark out diffraction spots, indicating number of peaks within the spot.
 *  3) Estimate a resolution boundary within (larger than) which
 *     a major fraction of the diffraction spots reside.
 *  4) Summarize the diffraction spots on a variety of measures.
 *  5) Accurately locate and mark out ice-rings.
 *
 * The library consists of two files: libdiffspot.h, libdiffspot.cc
 *
 * Developed by
 *  Zepu Zhang,    zpzhang@stanford.edu
 *  Ashley Deacon, adeacon@slac.stanford.edu
 *  and others.
 *
 * July 2001 - May 2003.
 */


#ifndef LIBDISTL_H
#define LIBDISTL_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>
#include <numeric>
// functional
// ctime
#define DISTL_CHECKPT {std::cout<<"checkpoint at line "<<__LINE__<<std::endl;}
#define DISTL_EXAMINE(A) {std::cout<<"variable "<<#A<<": "<<A<<std::endl;}

#include <spotfinder/core_toolbox/spot_types.h>
#include <spotfinder/core_toolbox/spot_shapes.h>
#include <iotbx/detectors/scanbox.h>

using namespace std;

namespace Distl {

typedef std::vector< std::vector< bool > > flag_array_t;
typedef std::vector< std::vector< double > > double_array_t;
typedef std::vector< std::vector< float > > float_array_t;

template <class T>
class constmat {

private:
        vector<const T*> data;

public:
        size_t nx;
        size_t ny;

        typedef T  data_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;

        typedef pointer iterator;
        typedef const_pointer const_iterator;

        iterator begin() { return data[0]; }
        iterator end() { return data[0] + nx*ny; }
        const_iterator begin() const { return data[0]; }
        const_iterator end() const {return data[0] + nx*ny; }

        //pointer operator[](size_t n) {return data[n]; }
        const_pointer operator[](size_t n) const { return data[n]; }

        size_t size() const { return nx*ny; }

        constmat() {
                nx = 0;
                ny = 0;
                data = vector<const T*>();
        }

        constmat(const T* mat, size_t n1, size_t n2) {
                nx = n1;
                ny = n2;
                data = vector<const T*>(nx);
                for (int x = 0; x < nx; x++)
                        data[x] = mat + x*ny;
        };

        void clear() {
                data.clear();
        }
};

typedef constmat<int> image_rawdata_t;

struct point {
  int x;
  int y;
  int value;

  point(): x(0), y(0) {;}
  point(const int xx, const int yy): x(xx), y(yy) {;}
  point(const int xx, const int yy, const int vv): x(xx), y(yy), value(vv) {;}
};


struct spot_base {
  typedef Distl::list_types<point>::list_t point_list_t;
  point_list_t bodypixels;    // area is bodypixels.size()
  point_list_t borderpixels;  // perimeter is borderpixels.size()
  // spot doesn't include the pixels in 'borderpixels'

  point_list_t maximas;

  double centerx;
  double centery;



  point peak;


  double peakresol;

  // added by Qingping
  double m_PixelValueSum; // total intensity over all pixel values

  //Methods
  int area() const {return bodypixels.size();}//number of pixels in spot
  double total_intensity() const {return m_PixelValueSum;}
  virtual
  void find_weighted_center(image_rawdata_t const&, flag_array_t const&,
                            float_array_t const&);
  virtual
  ~spot_base(){}
};

struct spot: public spot_base, public spot_shapes {

private:
  double p_majoraxis;
  double p_minoraxis;
  bool p_gotaxes;
  void p_getaxes();
public:
  double peakintensity;
  double peakheight;
  scitbx::af::shared<double> wts; // pixel-by-pixel ADC units above local background
  scitbx::af::shared<double> bkg; // pixel-by-pixel local background
  int nmaxima;         // number of local maximas
  inline spot():p_gotaxes(false){}
  spot(scitbx::af::flex_int::const_iterator);/* Special method to reconstruct
      the spot from a flex int.  Meant to be used for unpickling.
      The const_iterator is intentionally meant to be incremented by
      this constructor, so it can be used by the calling code to unpickle
      the next spot.
  */
  inline int max_pxl_x() const {return peak.x;}
  inline int max_pxl_y() const {return peak.y;}
  inline double intensity() const {return peakintensity;}
  inline int perimeter() const {return borderpixels.size();}
  inline double get_shape() const {return 4 * area() *
         scitbx::constants::pi / (perimeter() * perimeter());}
  double shape() const {return get_shape();}
  //approximate the spot by an ellipse
  double get_majoraxis();
  double get_minoraxis();
  //

  scitbx::vec2<double> get_radial_and_azimuthal_size(double, double);

  double resolution;//added Feb. 2006 for C++ manipulation of large-cell problems
  void setstate(Distl::point const&);//only used in unpickling
  virtual
  void find_weighted_center(image_rawdata_t const&, flag_array_t const&,
                            float_array_t const&);
  void show_summary(image_rawdata_t const&,double_array_t const&);
  inline double skewness() const {
    //A general indication of skewness, not a rigorous calculation
    //Formula gives max-pixel to center-of-mass distance as fraction of semi-major axis
    //Important:integer truncation prevents overestimation due to discrete image sampling
    scitbx::vec2<double>
      diffvec(max_pxl_x() - ctr_mass_x(), max_pxl_y() - ctr_mass_y());
    int difflength_trunc(static_cast<int>(diffvec.length()));
    return difflength_trunc / a();
  }
  virtual
  ~spot(){}
};


struct icering {
        // squared radius of inner and outter bounds
        double lowerr2;
        double upperr2;

        // resolution of inner and outter bounds
        double lowerresol;
        double upperresol;

        double strength;  // between 0 and 1.
        int npx;          // number of pixels on the ring.
};

struct background_plane_stats;

enum detector_shape { UNKNOWN, SQUARE, CIRCLE, RECTANGULAR_PIXEL };

inline detector_shape get_image_geometry( image_rawdata_t const& pixels){

  if (pixels.nx < 100 || pixels.ny < 100) { return UNKNOWN; }

  // simple heuristic to guess the image geometry.  If a small square region in
  // the upper-left of the image has a constant pixel value, then the active
  // area is likely an inscribed circle.  Otherwise, image is assumed square.

  image_rawdata_t::data_type reference = pixels [50][50];
  for (int x = 50; x < 100; ++x) {
    for (int y = 50; y < 100; ++y) {
      if (pixels [x][y] != reference) { return SQUARE; }
    }
  }
  return CIRCLE;
}

class diffimage {
public:
        int get_underload() const;

        void pxlclassify();
        void pxlclassify_scanbox(const int, const int, const int, const int,
                                                                const double);

        void search_icerings();

        void search_maximas();

        void search_spots();
        void search_overloadpatches();
private:
        void search_border_spot(const int, const int, spot&,
                                                        vector< vector<bool> >&);
        void search_border_overload(const int, const int, spot&,
                                                        vector< vector<bool> >&);
        template <typename DataType, typename CutoffType>
        void search_border_generic
        (const int, const int, spot&, vector< vector<bool> >&,
         DataType const&, CutoffType const&);
        bool pixelisonice(const int x, const int y) const;

        void screen_spots();
public:
        void imgresolution();
private:
        // Utilities
        double resol_to_r2(const double) const;
        double r2_to_resol(const double) const;
        double xy2resol(const double, const double) const;

        //! Convert x,y pixel coordinates to resolution
        /*! Gives exact same answer as xy2resol() to 12 decimal places.
            Prove the identity???
            Future: libdistl should implement two-theta/detector tilt.
         */
        double xy2resol_exact_normal(const double, const double) const;
        background_plane_stats* BP;
public:
        // Image facts

        double pixel_size;       // (mm)
        double distance;         // (mm)
        double wavelength;       // (angstroms)
        double resolb;           // (pixel_size/distance)^2
        double osc_start;
        double osc_range;
        double beam_center_x;    // (mm)
        double beam_center_y;    // (mm)
        int beam_x; // x index of beam center on image
        int beam_y; // y ...

        constmat<int> pixelvalue;
        // Each element vector of pixelvalue contains the values of
        // one column of pixels, from top to bottom.
        // Element vectors of pixelvalue are columns from left
        // to right.

        Distl::ptr_tiling tiling; //definition of image active area & scanbox tiling.

        // Processing parameters and options

        // Define processing region on the image.
        int firstx;
        int lastx;
        int firsty;
        int lasty;
        // Image geometry; important to exclude periphery from circular image
        detector_shape image_geometry;

        int underloadvalue;
        int overloadvalue;    // >= OVERLOAD: overloaded pixel
        bool report_overloads; // ignore the cutoff when spot searching
        int imgmargin;
        // Specifies margin width to be ignored from processing.

        const int npxclassifyscan; // = 3
        int scanboxsize[3];
        // border length of square boxes scanned in determining
        // local background and pixel intensities
        double bgupperint[3];
        // upper intensity threshold for bg pixel.
        // intensity < intensity_bguppercutoff &&
        // value > UNDERLOAD
        // => bg pixel

        // "intensity_bguppercutoff"
        // is set at the beginning of each scan, and used in the scan,
        // so the code works the same for each scan, except for different
        // parameters; but the code doesn't need to know how many
        // scans have been done.

        double difflowerint;
        // lower intensity threshold for diff pixel.
        // intensity > intensity_difflowercutoff &&
        // value < OVERLOAD
        // => diffract

        int iceringwidth;
        // Width of rings used while checking for ice-rings.
        double iceresolmax;
        double iceresolmin;
        // Maximum  and minimum resolution of the region checked
        // for ice-rings.

        const int nicecutoff; // = 2
        // number of elements in icering_cutoff and prctile.
        double icering_cutoff[2];
        // cutoff value for the ice-ring intensity.
        double icering_prctile[2];
        // intensity percentiles corresponding to cutoff values
        double icering_strengthprctile;
        // intensity percentile checked to assess continuity and strength of ice-ring

        int spotarealowcut;
        // Spot area lower bound.
        double spotareamaxfactor;
        // If spot area exceeds  median + (95th prctile - 5th prctile) * maxfactor,
        // it is suspicious and may cause trouble. Eliminate the spot/.
        double spotpeakintmaxfactor;
        // If spot peak intensity exceeds  median + (95th prctile - 5th prctile) * maxfactor,
        // it is out of the place and may cause trouble. Eliminate the spot/.
        //

        int imgresolringsize;
        // number of rings in resolution determination
        double imgresolringpow;
        // power of resol used in defining the rings in resolution determination.
        //double imgresolringsmoothspan;

        int spotbasesize;
        // Spots of size no smaller than this are used for summarizing
        // shape, intensity, etc.


        // Processing results and summaries

        vector< vector<float> > pixelintensity;
        // (PixelValue - LocalBackground) / LocalStandardError
        // Same location arrangement as in pixelvalue.

        vector< vector<float> > pixellocalmean;

        list<point> maximas;
        list<spot> spots;
        list<spot> overloadpatches;
        vector<icering> icerings;

        vector<double> imgresol_ringresol;
        vector<double> imgresol_ringresolpow;

        double imgresol;            // image resolution limit

        // Public functions

        diffimage();
        ~diffimage();

        void set_imageheader(const double, const double, const double,
                                        const double, const double, const double, const double);
        void cleardata();
        void set_imagedata(const int* const, const int, const int);

        int process();
        double resolution_outer; //initialize with negative value; >0 means impose a limit
        scitbx::af::shared<double> scanbox_background_resolutions;
        scitbx::af::shared<double> scanbox_background_means;
        scitbx::af::shared<double> scanbox_background_wndw_sz;
};



template<class T>
int percentiles(const vector<T>& data, const vector<double>& x, vector<T>& prctile);

}


#endif
