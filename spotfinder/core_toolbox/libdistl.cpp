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
 * The library consists of two files: libdistl.h, libdistl.cc.
 *
 * Help information is obtained by running the program
 * with no command-line parameters.
 *
 * Developed by
 *  Zepu Zhang,    zpzhang@stanford.edu
 *  Ashley Deacon, adeacon@slac.stanford.edu
 *  and others.
 *
 * July 2001 - May 2003.
 */

#include <spotfinder/core_toolbox/libdistl.h>
//#include <boost/timer.hpp>
//#include <omptbx/omp_or_stubs.h>
#include <stack> // to implement std::stack fix for stack overflow, see search_border_spot

using namespace Distl;

diffimage::diffimage(): report_overloads(false), npxclassifyscan(3), nicecutoff(2)
{
        // Only processing parameters not bound to any specific image
        // are specified here.
        // Image properties are initialized when image data is initialized.

    overloadvalue = 65535;

    //imgmargin = 20;
    imgmargin = 20;

    scanboxsize[0] = 101;
    scanboxsize[1] = 51;
    scanboxsize[2] = 51;

    bgupperint[0] = 1.5;
    bgupperint[1] = 2.0;
    bgupperint[2] = 2.5;

    difflowerint = 3.5;

    iceringwidth = 4;
    iceresolmin = 1.0;
    iceresolmax = 8.0;

    icering_prctile[0] = 0.45;
    icering_prctile[1] = 0.80;
        icering_strengthprctile = 0.20;

    icering_cutoff[0] = 0.0;
    icering_cutoff[1] = 1.5;

    spotarealowcut = 4;
        spotareamaxfactor = 5.0;
        spotpeakintmaxfactor = 10.0;

        spotbasesize = 16;

    imgresolringsize = 20;
    imgresolringpow = -3.0;
    resolution_outer = -1.;
}

diffimage::~diffimage() { cleardata(); }

inline int iround(double const& x) {
  if (x < 0) {return static_cast<int>(x-0.5);}
  return static_cast<int>(x+.5);
}

void diffimage::set_imageheader(const double pxsize, const double dist, const double wavelen,
                const double oscstart, const double oscrange,
                const double beamctrx, const double beamctry)
{


        cleardata();

        // In setimageheader and setimagedata
        // eventually all image property values will be reset,
        // except for processing parameters.
        // This is to avoid problems while multiple files
        // are processed in a row.

        pixel_size = pxsize;
        distance = dist;
        wavelength = wavelen;





        iceresolmin = std::max(wavelength/2.0,iceresolmin); //respect sampling theorem
        resolb = pxsize * pxsize / dist / dist;
        osc_start = oscstart;
        osc_range = oscrange;
        beam_center_x = beamctrx;
        beam_center_y = beamctry;
        // The origin point wrt which the beam center is located is not clear.
        // Here it is assumed that beam_center_x and beam_center_y use
        // the same coord system as X and Y do, i.e.,
        // beam_center_x: top->bottom
        // beam_center_y: left->right

        //beam_x/beam_y are assigned as integers to support
        // array indexing in subsequent steps.  However, it is dangerous
        // to use static_cast<int> because this produces a truncation
        // instead of a proper rounding to the nearest integer.
        // Use of the new function iround() corrects this.

        beam_x = iround(beam_center_x / pixel_size);
        beam_y = iround(beam_center_y / pixel_size);
}

void diffimage::cleardata(){
        for (int i = 0; i < pixelintensity.size(); i++)
                pixelintensity[i].clear();
        pixelintensity.clear();
        for (int i = 0; i < pixellocalmean.size(); i++)
                pixellocalmean[i].clear();
        pixellocalmean.clear();

        imgresol_ringresol.clear();
        imgresol_ringresolpow.clear();

        maximas.clear();
        spots.clear();
        overloadpatches.clear();
        icerings.clear();

        imgresol = 0;
}

void diffimage::set_imagedata(const int* const data, const int ncol, const int nrow)
{
        // DATA are stored by column, from left to right.
        // Type is int.

        pixelvalue = constmat<int>(data, ncol, nrow);

        firstx = imgmargin;
        lastx = pixelvalue.nx-1-imgmargin;
        firsty = imgmargin;
        lasty = pixelvalue.ny-1-imgmargin;
}


int diffimage::get_underload() const
{
//printf("spotarealowcut(-s2) %3d, spotbasesize(-s3) %3d bgupperint[2] %4.2f difflowerint %4.2f spotareamaxfactor %4.2f\n",spotarealowcut, spotbasesize, bgupperint[2], difflowerint, spotareamaxfactor);
  scanbox_tiling_pilatus6M* possible_pilatus =
    dynamic_cast<scanbox_tiling_pilatus6M*>(&(*tiling));

  if (possible_pilatus != NULL){
    //std::cout<<"THIS IS A PILATUS"<<std::endl;
    return -1;
  } else {
    //std::cout<<"THIS IS NOT A PILATUS"<<std::endl;
  }
        // *****************************************************************
        // Determine the upperbound value for underloaded pixels, i.e.,
        // pixels with value <= UNDERLOAD are considered underloaded, like
        // those on the border or blocked by the beam stick.
        //
        // This function checks the whole image, not restricted by
        // 'firstx', 'lastx', 'firsty', 'lasty'.
        // *****************************************************************

        int ncols = pixelvalue.nx;
        int nrows = pixelvalue.ny;
        if (ncols < 100 || nrows < 100) {
          // No underload or offset determination for small test cases.
          return 0;
        }
        int crosswid = 40;
        int cornerfrac = 5;
        int cornerwidth = pixelvalue.ny / cornerfrac;
        int np = 4 * cornerwidth * cornerwidth +
                 crosswid * (nrows + ncols) - crosswid*crosswid;

        vector<int> px(np);
        vector<int>::iterator px_itr = px.begin();

        // Four corners.
        for (int anchor = 0; //linux performance: first iteration is very fast
                 anchor <= pixelvalue.nx - cornerwidth; //2nd consumes more time
                 anchor += pixelvalue.nx - cornerwidth) {
          for (int x = anchor; x < anchor + cornerwidth; x++) {
                  px_itr = std::copy(pixelvalue[x],
                                     pixelvalue[x] + cornerwidth,
                                     px_itr);
                  px_itr = std::copy(pixelvalue[x] + pixelvalue.ny - cornerwidth,
                                     pixelvalue[x]+ pixelvalue.ny,
                                     px_itr);
          }
        }

        // Central cross.
        for (int icol=0; icol<ncols; icol++) {
            px_itr =  std::copy(pixelvalue[icol] + nrows/2 - crosswid/2,
                                pixelvalue[icol] + nrows/2 + crosswid/2,
                                px_itr);
        }
        for (int icol=ncols/2-crosswid/2; icol<ncols/2+crosswid/2; icol++) {
            px_itr =  std::copy(pixelvalue[icol],
                                pixelvalue[icol] + nrows/2 - crosswid/2,
                                px_itr);
            px_itr =  std::copy(pixelvalue[icol] + nrows/2 + crosswid/2,
                                pixelvalue[icol] + nrows,
                                px_itr);
        }


        // Check 'nq' percentiles of array 'px', identify the value that is
        // several consecutive percentiles, which indicates that value happens
        // a lot.
        int nq = 30;
        vector<int> q(nq);

        for (int i=0; i<nq; i++) {
                nth_element(px.begin()+i*np/nq, px.begin()+np*(i+1)/(nq+1), px.end());
                q[i] = px[np*(i+1)/(nq+1)];
        }

        int ul = q[nq/4];

        if (q[nq-1]==q[nq-2] && q[nq-2]==q[nq-3]) {
                nth_element(px.begin(), px.begin()+np*(nq+2)/(nq+3), px.end());
                ul = px[np*(nq+2)/(nq+3)];
        } else {
                for (int i=nq-2; i>1; i--)
                        if (q[i]==q[i-1] && q[i-1]==q[i-2]) {
                                ul = q[i+1]; break;
                        }
        }

/* NKS.  Add sanity check to make sure that no more than 10% of pixels
   are masked out. The most straightforward way is to copy the data to a new
   container and use the nth_element algorithm. This is the only way to
   insure that true background pixels are not masked out as underloads.
   Two concerns: 1) nth_element might be time-consuming, so only use every
   4th pixel. 2) want the algorithm to work equally well on circular &
   square plates, so filter out everything but the inscribed circle.
*/
        double cutof_frac = 0.1; /*arbitrary choice: no more than
                                 10% of pixels should be masked as underloads*/
        double imagearea = pixelvalue.size();
        int speedup = 2;
        std::vector<int> sdata; sdata.reserve((int)imagearea/(speedup*speedup));
        if (image_geometry==CIRCLE) {
          for (int x=0; x<pixelvalue.nx; x+=speedup) {
            double half_chord = std::sqrt((imagearea/4-(ncols/2-x)*(ncols/2-x)));
            for (int y=(nrows/2)-(int)half_chord+1;
                     y<(nrows/2)+(int)half_chord-1; y+=speedup){
              sdata.push_back(pixelvalue[x][y]);
            }
          }
        } else {
          for (int x=0; x<pixelvalue.nx; x+=speedup) {
            for (int y=0; y<nrows; y+=speedup){
              sdata.push_back(pixelvalue[x][y]);
            }
          }
        }
        int n_cutoff = (int)((cutof_frac) * sdata.size());
        nth_element(sdata.begin(),sdata.begin()+n_cutoff,sdata.end());
        int n_cutoff_value = *(sdata.begin()+n_cutoff);

        return std::min(n_cutoff_value,ul);
}

int diffimage::process()
{
//boost::timer T = boost::timer();
        // 50-pixel corner sample to detect circular image geometry
        image_geometry = get_image_geometry(pixelvalue);
        underloadvalue = get_underload();
//std::cout<<"underloads..."<<T.elapsed()<<std::endl;T.restart();
        pxlclassify();
//std::cout<<"classify..."<<T.elapsed()<<std::endl;T.restart();
        search_icerings();
//std::cout<<"icerings..."<<T.elapsed()<<std::endl;T.restart();
        search_maximas();
//std::cout<<"maximas..."<<T.elapsed()<<std::endl;T.restart();
        search_spots();
//std::cout<<"spots..."<<T.elapsed()<<std::endl;T.restart();
        search_overloadpatches();
//std::cout<<"overloadpatches..."<<T.elapsed()<<std::endl;
        imgresolution();

        return 0;
}


void diffimage::pxlclassify()
{
        // *********************************************
        // Classify pixels into underloaded, background,
        // diffracted, overloaded, etc.
        //
        // Scan the image box by box.
        // **********************************************************
        pixelintensity = float_array_t(
                           pixelvalue.nx, vector<float>(pixelvalue.ny, 0.0));
        pixellocalmean = float_array_t(
                           pixelvalue.nx, vector<float>(pixelvalue.ny, 0.0));
        for (int i=0; i<npxclassifyscan; i++) {
                scanbox_background_resolutions.clear();
                scanbox_background_means.clear();
                scanbox_background_wndw_sz.clear();
                if (scanboxsize[i] == 0)
                        continue;

                int boxsize = scanboxsize[i];

                int x_box = boxsize;
                int y_box = boxsize;

                //#pragma omp parallel for
                /*
                for (int xstart = firstx;xstart <= (lastx - x_box); xstart += x_box){
                  int xend = xstart + x_box - 1;
                  for (int ystart = firsty;ystart <= (lasty - y_box); ystart += y_box) {
                    int yend = ystart + y_box - 1;
                    pxlclassify_scanbox(xstart, xend, ystart, yend, bgupperint[i]);
                  }
                }
                */

                tiling->reset();
                for (interval_ptr xs = tiling->x_tiles(x_box); xs != tiling->x_end(); ++xs){
                  //if (i==0) {printf("x = %4d to %4d inclusive\n",xs->first,xs->last);}
                  for (interval_ptr ys = tiling->y_tiles(y_box); ys != tiling->y_end(); ++ys) {
                    pxlclassify_scanbox(xs->first, xs->last, ys->first, ys->last, bgupperint[i]);
                  }
                }
        }
}

struct backplane {
 public:
  int boxnbg;
  double boxmean, boxvar, boxstd;
  double Sum_x,Sum_x2;
  backplane(){
  //    cout<<"base constructor"<<endl;
clear(); }
  virtual void accumulate (const int&, const int&, const int& px){
    //cout<<"base accumulate"<<endl;
    Sum_x += px;
    Sum_x2 += (double)(px)*px;
    boxnbg++;
  }
  virtual void clear(){
    boxnbg = 0;
    boxmean = boxvar = boxstd = 0.;
    Sum_x = Sum_x2 = 0.;
  }
  virtual void reinitialize(const int&, const int&){
    //cout<<"base reinitialize"<<endl;
    clear();
  }
  virtual void finish(){
    //cout<<"base finish"<<endl;
    boxmean = Sum_x/boxnbg;
    boxvar = Sum_x2/boxnbg - boxmean*boxmean;
    boxstd = std::sqrt(boxvar);
  }
  virtual inline double localmean (const int&, const int&) {
      //cout<<"base mean "<<boxmean<<endl;
    return boxmean; }

  // Trivial calculation of detector gain; re-check derivation before commenting in
  //inline double gain () const { //Suggested by Mosflm manual: digitized value = GAIN * Equiv # of photons
  //  double sample_average = Sum_x/boxnbg;
  //  double sample_variance = (1./boxnbg)*(Sum_x2 - sample_average*sample_average);
  //  return std::sqrt(sample_variance) / sample_average; // is mean_squared / mean
  //}

  virtual ~backplane(){}
};

struct corrected_backplane: public backplane {
 private:
  int Sum_p2,Sum_pq,Sum_p,Sum_q2,Sum_q;
  double Sum_xp,Sum_xq;
  int xstart,ystart;
  double a,b,c;
  std::vector<int> rho_cache;
  std::vector<int> p_cache;
  std::vector<int> q_cache;
  double rmsd;
  double p,q; //temporary values
 public:
  corrected_backplane(const int& xst, const int& yst):
    xstart(xst),ystart(yst) {
    //cout<<"corrected constructor"<<endl;
    clear();
  }
  void reinitialize(const int& xst, const int& yst){
    //cout<<"corrected reinitialize"<<endl;
    xstart = xst; ystart = yst; clear();
  }
  inline void clear(){
    backplane::clear();
    Sum_p2=0;Sum_pq=0;Sum_p=0;Sum_q2=0;Sum_q=0;Sum_xp=0;Sum_xq=0;
    rho_cache.clear();
    p_cache.clear();
    q_cache.clear();
    rmsd=0.;
  }
  inline void accumulate(const int& x, const int& y, const int& px){
    //cout<<"corrected accumulate"<<endl;
    backplane::accumulate(x,y,px);
    int p = x-xstart;
    int q = y-ystart;
    Sum_p2+=p*p;
    Sum_pq+=p*q;
    Sum_p+=p;
    Sum_q2+=q*q;
    Sum_q+=q;
    Sum_xp+=px*p;
    Sum_xq+=px*q;
    rho_cache.push_back(px);p_cache.push_back(p);q_cache.push_back(q);
  }
  inline void finish(){
    //cout<<"corrected finish"<<endl;
    scitbx::mat3<double> rossmann(Sum_p2,Sum_pq,Sum_p,
                          Sum_pq,Sum_q2,Sum_q,
                          Sum_p,Sum_q,boxnbg);
    scitbx::vec3<double> obs(Sum_xp,Sum_xq,Sum_x);
    scitbx::mat3<double> rinv = rossmann.inverse();
    //scitbx::vec3<double> abc = rossmann.inverse()*obs;
    //a=abc[0]; b= abc[1]; c=abc[2];
    a = rinv[0]*Sum_xp + rinv[1]*Sum_xq +rinv[2]*Sum_x;
    b = rinv[3]*Sum_xp + rinv[4]*Sum_xq +rinv[5]*Sum_x;
    c = rinv[6]*Sum_xp + rinv[7]*Sum_xq +rinv[8]*Sum_x;
    for (int v=0; v<boxnbg; ++v){
      double bgobs_bgplane = rho_cache[v] - a*p_cache[v] - b*q_cache[v] -c;
      rmsd +=  bgobs_bgplane*bgobs_bgplane;
    }
    rmsd = std::sqrt(rmsd/boxnbg);
    boxstd=rmsd;
  }
  inline double localmean(const int&x, const int&y){
    return a*(x-xstart)+b*(y-ystart)+c;
  }
};

void diffimage::pxlclassify_scanbox(const int xstart, const int xend,
                                                                        const int ystart, const int yend,
                                    const double intensity_bguppercutoff)
{
  // ******************************************************************************
  // Cutoff based on an imposed resolution limit.  This is just an ad-hoc fix.
  //  Known issues:  1) This does not take two-theta offsets or detector tilt into account.
  //
  // ******************************************************************************
  if (resolution_outer > 0.0) {
    bool any_point_within_resol_limit = false;
    for (int x = xstart; x <= xend; x+=xend-xstart){
      for (int y = ystart; y <= yend; y+=yend-ystart){
        if ( xy2resol(x,y) > resolution_outer ) { any_point_within_resol_limit = true; }
      }
    }
    if (!any_point_within_resol_limit) { return; }
  }
  // ******************************************************************************
  // if any corner of the box is outside of the image active area, then do not
  // scan the box
  // ******************************************************************************
  if (image_geometry==CIRCLE) {
    int iradius = pixelvalue.nx / 2;
    int iradius_sq = iradius*iradius;
    for (int x = xstart; x <= xend; x+=xend-xstart){
      for (int y = ystart; y <= yend; y+=yend-ystart){
        if ( (x-iradius)*(x-iradius) + (y-iradius)*(y-iradius) > iradius_sq )
          { return; }
      }
    }
  }
  // ******************************************************************************
  // Calculate mean and std of all background pixels in the box.
  // Assign the calculated mean and std to each pixel in the box as its background.
  // Classify each pixel in the box based on this background.
  //
  // Expand the box if number of background pixels falls short of the required
  // amount.
  // ******************************************************************************

  const int speedup = 2;  // Sample from the window. No need to be exhaustive.
  const bool use_plane_correction = true;
  backplane* bp;
  if (use_plane_correction) { bp = new corrected_backplane(xstart,ystart); }
  else { bp = new backplane(); }

  const int nbgmin = (xend - xstart + 1) * (yend - ystart + 1) * 2/3
      /(speedup*speedup);
  // Making 'nbgmin' reasonably small compared with the box size
  // is to avoid frequent need to expand the box due to shortage
  // of background pixels.

  for (int x=xstart; x<=xend; x += speedup){
    for (int y=ystart; y<=yend; y += speedup){
      int px = pixelvalue[x][y];
      if (pixelintensity[x][y]<intensity_bguppercutoff &&
          px > underloadvalue &&
          (report_overloads || px < overloadvalue)) {
        bp->accumulate(x,y,px);
      }
    }
  }


  int xbox = xend - xstart + 1;
  int ybox = yend - ystart + 1;
  int halfincrease = min(xbox, ybox) / 3;
  int x0 = xstart;
  int y0 = ystart;

  while (bp->boxnbg<nbgmin) {
    // Usually this loop is skipped because
    // usually boxnbg > nbgmin.
    xbox += 2*halfincrease;
    if (xbox > 1000) {
    delete bp;
    throw scitbx::error("pxlclassify scanbox function: cannot distinguish signal/background");
    } //do not permit scanbox to increase without limit
    ybox += 2*halfincrease;
    x0 = min(max(0,x0-halfincrease), static_cast<int>(pixelvalue.nx - xbox));
    y0 = min(max(0,y0-halfincrease), static_cast<int>(pixelvalue.ny - ybox));
    // Whole image, rather than confined by firstx, firsty, lastx, lasty,
    // is used for scan.

    bp->clear();
    for (int x=x0; x<x0+xbox; x += speedup){
      for (int y=y0; y<y0+ybox; y += speedup){
        int px = pixelvalue[x][y];
        if (pixelintensity[x][y]<intensity_bguppercutoff &&
            px > underloadvalue &&
            (report_overloads || px < overloadvalue)) {
          bp->accumulate(x,y,px);
        }
      }
    }

  }

  bp->finish();
  //BP->evaluate_plane(xstart,xend,ystart,yend);

  // Calculate pixel intensity.
  //
  // Intensity is calculated and assigned to the original box,
  // instead of the expansion, if any, to the box.

  //avoid a divide-by-zero if the signal is flat
  double multfactor = bp->boxstd == 0. ? 0. : 1.0 / bp->boxstd;

  for (int x=xstart; x<=xend; x++) {
    for (int y=ystart; y<=yend; y++) {
      double bxmean = bp->localmean(x,y);
      pixelintensity[x][y] = (pixelvalue[x][y] - bxmean) * multfactor;
      pixellocalmean[x][y] = bxmean;
          //pixellocalstd[x][y] = boxstd;
    }
  }
  scanbox_background_resolutions.push_back(xy2resol( int((xstart+xend)/2.),int((ystart+yend)/2.) ));
  scanbox_background_means.push_back( bp->localmean(int((xstart+xend)/2.),int((ystart+yend)/2.)));
  scanbox_background_wndw_sz.push_back( bp->boxnbg );
  delete bp;
}



void diffimage::search_icerings()
{
        // *************************
        // Detect ice-rings.
        // *************************

        double icermin, icermax;

        int nrings;
        // ring index increases as distance from ring to image center increases.
        // allocate one extra ring to avoid potential problems caused by weird (x,y) values.

        vector< vector<float> > ringpxints;
        // pixel intensities on each ring.

        bool checksquare = false;
        // If only a square area symmetrically surrounding the beam center is to be checked,
        // set checksquare to TRUE.
        // If non-symmetrical, possibly rectangular, are surrounding the beam center
        // is acceptible, set checksquare to FALSE.

        if (checksquare) {
                int deltaxy = min(min(beam_x-firstx, lastx-beam_x), min(beam_y-firsty, lasty-beam_y));
                icermax = sqrt(2.0*deltaxy*deltaxy);
                if (r2_to_resol(icermax*icermax) < iceresolmin)
                        icermax = sqrt(resol_to_r2(iceresolmin));

                icermin = sqrt(resol_to_r2(iceresolmax));

                nrings = static_cast<int>(icermax / iceringwidth + 1);
                ringpxints.resize(nrings);

                for (int ringidx=0; ringidx<nrings; ringidx++)
                        ringpxints[ringidx].reserve((ringidx+1)*iceringwidth*8*iceringwidth);
                        // Reserve enough space for pixels on each ring.


                for (int dx=static_cast<int>(0.7*icermax); dx>0; dx--) {
                        for (int dy=static_cast<int>(0.7*icermax); dy>0; dy--) {
                                // 0.7 = 1/sqrt(2)
                                double r = sqrt(static_cast<double>(dx*dx + dy*dy));
                                int ringidx = static_cast<int>(r / iceringwidth);
                                // Calculate ring index for the upper-left 1/4 of the image,
                                // take pixels on the other 3/4 image by symmetry.

                                ringpxints[ringidx].push_back(pixelintensity[beam_x-dx][beam_y-dy]);
                                ringpxints[ringidx].push_back(pixelintensity[beam_x-dx][beam_y+dy]);
                                ringpxints[ringidx].push_back(pixelintensity[beam_x+dx][beam_y-dy]);
                                ringpxints[ringidx].push_back(pixelintensity[beam_x+dx][beam_y+dy]);
                        }
                }
        } else {
                double dx1 = beam_x - firstx;
                double dx2 = lastx - beam_x;
                double dy1 = beam_y - firsty;
                double dy2 = lasty - beam_y;
                dx1 = dx1 * dx1;
                dx2 = dx2 * dx2;
                dy1 = dy1 * dy1;
                dy2 = dy2 * dy2;
                icermax = max(max(dx1+dy1, dx1+dy2), max(dx2+dy1, dx2+dy2));
                double temp = resol_to_r2(iceresolmin);
                if (temp > icermax)
                        icermax = sqrt(icermax);
                else
                        icermax = sqrt(temp);

                icermin = sqrt(resol_to_r2(iceresolmax));

                nrings = static_cast<int>(icermax / iceringwidth + 1);
                ringpxints.resize(nrings);

                for (int ringidx=0; ringidx<nrings; ringidx++)
                        ringpxints[ringidx].reserve((ringidx+1)*iceringwidth*8*iceringwidth);
                        // Reserve enough space for pixels on each ring.

                //#pragma omp parallel for
                for (int x=firstx; x<=lastx; x++) {
                        double dx = x - beam_x;
                        for (int y=firsty; y<=lasty; y++) {
                                double dy = y - beam_y;
                                double r = sqrt(dx*dx + dy*dy);
                                int ringidx = static_cast<int>(r / iceringwidth);
                                if (ringidx < ringpxints.size()) {
                                        ringpxints[ringidx].push_back(pixelintensity[x][y]);
                                }
                        }
                }

        }
        int iringmin = static_cast<int>(icermin / iceringwidth);
        // Inner-most ring to consider for ice-ring checking.

        int iringmax = nrings - 2;
        // Outer-most ring to consider for ice-ring checking.

        vector<bool> ringice(nrings, false);
        // Keep track of whether or not a ring is an ice-ring

        icerings.reserve(5);

        ///////////////////////////////////////////
        // Record the ice-rings found.
        // Contiguous ice-rings are combined to be counted as one.
        ///////////////////////////////////////////
        int lasticering = -3;
        for (int ringidx=iringmin; ringidx<=iringmax; ringidx++) {
                bool passed = true;
                double prctints[2];

                for (int i=0; i<2; i++) {
                        int prc = static_cast<int>(ringpxints[ringidx].size()
                                                         * icering_prctile[i]);
                        if (prc==0) { passed = false; break;} // not enough points in ring to analyze
                        nth_element(ringpxints[ringidx].begin(),
                                                ringpxints[ringidx].begin()+prc,
                                                ringpxints[ringidx].end());
                        prctints[i] = ringpxints[ringidx][prc];
                        if(prctints[i] < icering_cutoff[i]) {
                                // Specified intensity percentile too small -- not ice.
                                passed = false;
                                break;
                        }
                }

                if (passed) {
                        double r1 = iceringwidth*ringidx;
                        double r2 = r1 + iceringwidth;
                        int prc = static_cast<int>(ringpxints[ringidx].size() * icering_strengthprctile);
                        nth_element(ringpxints[ringidx].begin(),
                                                ringpxints[ringidx].begin()+prc,
                                                ringpxints[ringidx].end());
                        double icestrength = ringpxints[ringidx][prc];
                        if (ringidx > lasticering+1) {
                                // a new ice-ring, not continuous with the previous one.
                                icerings.push_back(icering());
                                icerings.back().lowerr2 = r1*r1;
                                icerings.back().upperr2 = r2*r2;
                                icerings.back().upperresol = r2_to_resol(r1 * r1);
                                icerings.back().lowerresol = r2_to_resol(r2 * r2);
                                icerings.back().strength = icestrength;
                        } else {
                                icerings.back().upperr2 = r2*r2;
                                icerings.back().lowerresol = r2_to_resol(r2 * r2);
                                icerings.back().strength = max(icerings.back().strength, icestrength);
                        }
                        lasticering = ringidx;
                }
        }

}



void diffimage::search_maximas()
{
        // ********************************
        // Search for local maximas.
        // ********************************


        // A local maxima is a pixel which is a diffraction pixel,
        // is not overloaded, and whose value is not smaller than
        // any of its 8 neighbors.

        int neighboramount = min(4, spotarealowcut-1);
        int diffnum = neighboramount - 1;

        for (int x=firstx+1; x<lastx; x++) {
                for (int y=firsty+1; y<lasty; y++) {

                        if (pixelvalue[x][y] >= overloadvalue) {
                                maximas.push_front(point(x,y));

                        } else {
                                double pxint = pixelintensity[x][y];
                                int pxv = pixelvalue[x][y];

                                if (pxint > difflowerint &&
                                        pxv >  underloadvalue &&
                                        pxv >= pixelvalue[x-1][y-1] &&
                                        pxv >= pixelvalue[x][y-1] &&
                                        pxv >= pixelvalue[x+1][y-1] &&
                                        pxv >= pixelvalue[x-1][y] &&
                                        pxv >= pixelvalue[x+1][y] &&
                                        pxv >= pixelvalue[x-1][y+1] &&
                                        pxv >= pixelvalue[x][y+1] &&
                                        pxv >= pixelvalue[x+1][y+1]) {
                                        int goodneighbors =
                                                (pixelintensity[x-1][y-1] > difflowerint) +
                                                (pixelintensity[x][y-1] > difflowerint) +
                                                (pixelintensity[x+1][y-1] > difflowerint) +
                                                (pixelintensity[x-1][y] > difflowerint) +
                                                (pixelintensity[x+1][y] > difflowerint) +
                                                (pixelintensity[x-1][y+1] > difflowerint) +
                                                (pixelintensity[x][y+1] > difflowerint) +
                                                (pixelintensity[x+1][y+1] > difflowerint);
                                        if (goodneighbors > diffnum) {
                                                maximas.push_front(point(x,y));
                                        }
                                }
                        }
                }
        }
}


void diffimage::search_spots()
{
        // ************************************************
        // Search for spots.
        // Record area and calculate intensity of each spot.
        //
        // Every maxima belongs to a unique spot; but
        // one spot could have multiple maximas, depending
        // on whether neighborhoods of two maximas are
        // connnected, in which case both maximas are
        // enclosed in one single spot.
        // ************************************************

        const double PI = 3.14159265;

        int nx = pixelvalue.nx;
        int ny = pixelvalue.ny;
        vector< vector<bool> > pixelvisited(nx, vector<bool>(ny, false));

        for (list<point>::iterator p=maximas.begin(); p!=maximas.end(); p++) {
                if (!pixelvisited[p->x][p->y]) {
                        if (report_overloads || pixelvalue[p->x][p->y] < overloadvalue) {
                                spots.push_back(spot());
                                search_border_spot(p->x, p->y, spots.back(), pixelvisited);

                                // Discard this spot if any of its pixel lies on an ice-ring.
                                for (spot::point_list_t::const_iterator q=spots.back().borderpixels.begin();
                                                q!=spots.back().borderpixels.end(); q++) {
                                        if (pixelisonice(q->x, q->y)) {
                                                spots.pop_back();
                                                break;
                                        }
                                }
                        }
                }
        }

        // Calculate several basic properties.

        vector< vector<bool> > pixelismaxima(nx, vector<bool>(ny, false));
        for (list<point>::const_iterator p=maximas.begin(); p!=maximas.end(); p++)
                pixelismaxima[p->x][p->y] = true;

        for (list<spot>::iterator p=spots.begin(); p!=spots.end(); p++) {
                p->nmaxima = 0;
                for (spot::point_list_t::const_iterator q=p->bodypixels.begin();
                     q!=p->bodypixels.end(); q++) {
                        if (pixelismaxima[q->x][q->y]) {
                          p->maximas.push_back(*q); //this line is solely to produce imagemagick output 12/06
                          p->nmaxima ++;
                          if (p->nmaxima == 1)
                            p->peak = *q;
                          else
                            if (pixelvalue[q->x][q->y] > pixelvalue[p->peak.x][p->peak.y])
                              p->peak = *q;
                        }
                }

                p->peakintensity = pixelintensity[p->peak.x][p->peak.y];
        }

        // Screen on spots based on size and peak intensity.
        screen_spots();

        // Calculate more properties.
        for (list<spot>::iterator p=spots.begin(); p!=spots.end(); p++) {
                // Spot size is guaranteed to exceed 1.

                p->find_weighted_center(pixelvalue,pixelismaxima,pixellocalmean);

                p->peakresol = xy2resol(p->peak.x, p->peak.y);

                p->peakheight = pixelvalue[p->peak.x][p->peak.y] - pixellocalmean[p->peak.x][p->peak.y];

                double xoffset = p->peak.x - beam_x;
                double yoffset = p->peak.y - beam_y;

                if (yoffset==0) {
                        if (xoffset < 0)
                                p->angle = PI * 0.5;
                        else if (xoffset > 0)
                                p->angle = PI * 1.5;
                        else
                                p->angle = 0;
                } else {
                        p->angle = atan(static_cast<double>(-xoffset)/yoffset);
                        if (yoffset<0)
                                p->angle = p->angle + PI;
                        if (p->angle < 0)
                                p->angle = 2*PI + p->angle;
                }
        }
}


void diffimage::search_border_spot(const int x, const int y, spot& curspot,
                                   vector< vector<bool> >& pixelvisited)
{
  // **********************************************************
  // Searches for contiguous diffraction pixels to be included
  // in one spot.
  // ************************************************************
  search_border_generic<float_array_t, double> (
    x,y,curspot,pixelvisited,pixelintensity,difflowerint);
  return;
}

void diffimage::search_border_overload(const int x, const int y, spot& curspot,
                                   vector< vector<bool> >& pixelvisited)
{
  search_border_generic<constmat<int>, int> (
    x,y,curspot,pixelvisited,pixelvalue,overloadvalue);
  return;
}

template <typename DataType, typename CutoffType>
void diffimage::search_border_generic
  (const int x,const int y,spot& curspot,vector< vector<bool> >& pixelvisited,
   DataType const& pixeldata, CutoffType const& cutoff)
{ std::stack<point> Q;
  Q.push(point(x,y));
  while (!Q.empty()) {
    point pt = Q.top();
    Q.pop();
    if (pixelvisited[pt.x][pt.y]) {continue;}
    if (pt.x < firstx || pt.x > lastx || pt.y < firsty || pt.y > lasty) {continue;}
    pixelvisited[pt.x][pt.y] = true;
    if (pixeldata[pt.x][pt.y]<=cutoff) {
      // Hit a nighboring pixel, which does not belong to
      // this spot. Record border pixels and return.
      curspot.borderpixels.push_back(pt);
      // From this step it is clear that a border pixel
      // is one that is bordering a spot, but does not
      // belong to the spot, i.e, is not counted in the
      // spot area.
      // Perimeter is number of such border pixels.
    } else {
      // An interior pixel on the spot.
      // Label current pixel and search on.
      curspot.bodypixels.push_back(pt);
      Q.push(point(pt.x-1,pt.y-1));
      Q.push(point(pt.x  ,pt.y-1));
      Q.push(point(pt.x+1,pt.y-1));
      Q.push(point(pt.x-1,pt.y  ));
      Q.push(point(pt.x+1,pt.y  ));
      Q.push(point(pt.x-1,pt.y+1));
      Q.push(point(pt.x  ,pt.y+1));
      Q.push(point(pt.x+1,pt.y+1));
    }
  }
  return;
}


void diffimage::search_overloadpatches()
{
        int nx = pixelvalue.nx;
        int ny = pixelvalue.ny;
        vector< vector<bool> > pixelvisited(nx, vector<bool>(ny, false));

        for (list<point>::iterator p=maximas.begin(); p!=maximas.end(); p++) {
                if (!pixelvisited[p->x][p->y]) {
                        if (pixelvalue[p->x][p->y] >= overloadvalue) {
                                // If maxima is overloaded, search for a overloaded patch.
                                overloadpatches.push_back(spot());
                                search_border_overload(p->x, p->y, overloadpatches.back(), pixelvisited);

                                // Remove this maxima from maxima list.
                                list<point>::iterator pp = p;
                                p--;
                                maximas.erase(pp);
                        }
                }
        }
}


bool diffimage::pixelisonice(const int x, const int y) const
{
        // ******************************************************************
        // Check whether or not a pixel in on an ice-ring.
        // This function is called only when ice-ring does exist on the image.
        // ******************************************************************
        double resol = xy2resol(x, y);

        for (vector<icering>::const_iterator q=icerings.begin(); q!=icerings.end(); q++) {
                if (resol >= q->lowerresol) {
                        if (resol <= q->upperresol)
                                return (true);
                        else
                                return (false);
                }
        }

        return (false);
}


void diffimage::screen_spots()
{
        // Eliminate bad spots according to chosen criteria.

        for (list<spot>::iterator p=spots.begin(); p!=spots.end(); p++) {
                if (p->bodypixels.size()<spotarealowcut) {
                        list<spot>::iterator q = p;
                        p--;
                        spots.erase(q);
                        continue;
                }
        }

        vector<int> spotarea;
        vector<double> peakint;
        spotarea.reserve(spots.size());
        peakint.reserve(spots.size());

        for (list<spot>::iterator p=spots.begin(); p!=spots.end(); p++) {
                spotarea.push_back(p->bodypixels.size());
                peakint.push_back(p->peakintensity);
        }

        vector<double> prctX(3);
        prctX[0] = 0.05;
        prctX[1] = 0.50;
        prctX[2] = 0.95;

        vector<int> areaprcts(3);
        vector<double> peakintprcts(3);

        percentiles<int>(spotarea,prctX,areaprcts);
        percentiles<double>(peakint,prctX,peakintprcts);

        int areaupper = static_cast<int>(areaprcts[1] + spotareamaxfactor * (areaprcts[2] - areaprcts[0]));
        double peakintupper = peakintprcts[1] +
                spotpeakintmaxfactor * (peakintprcts[2] - peakintprcts[0]);


        for (list<spot>::iterator p=spots.begin(); p!=spots.end(); p++) {
                if (p->peakintensity>peakintupper || p->bodypixels.size()>areaupper) {
                        list<spot>::iterator q = p;
                        p--;
                        spots.erase(q);
                        continue;
                }
        }

}

void diffimage::imgresolution()
{
        const double PI = 3.14159265;

        if (spots.size()<=10) {
                imgresol = 10.0;
                return;
        } else {
                imgresol = 2.0;
        }

        int nimgresolrings = min(40,max(10,static_cast<int>(spots.size() / imgresolringsize)));


        /////////////////////////////////////////
        // Calulate corner modifiers
        // to be used to adjust amount of spots.
        /////////////////////////////////////////


        // Resolution and corner factor of each spot.
        vector<double> spotresol;
        vector<int> spotcornerfactor;
        spotresol.reserve(spots.size());
        spotcornerfactor.reserve(spots.size());

        int ndirections = 15;
        // Number of spots in each direction.
        vector<int> nanglespot(ndirections,0);

        // Largest radius with complete circle on the image
        int commonr = min(min(beam_x-firstx, lastx-beam_x),
                                        min(beam_y-firsty, lasty-beam_y));
        double unitangle = 2*PI / ndirections;
        double rupper = beam_x - firstx;
        double rlower = lastx - beam_x;
        double rleft = beam_y - firsty;
        double rright = lasty - beam_y;


        for (list<spot>::const_iterator p = spots.begin();
                        p != spots.end(); p++) {
                spotresol.push_back(p->peakresol);

                double dx = p->peak.x - beam_x;
                double dy = p->peak.y - beam_y;
                double r = sqrt(dx*dx + dy*dy);
                if (r < commonr) {
                        spotcornerfactor.push_back(1);

                        // spots in corners are ignored in directional assessment.
                        int idx = static_cast<int>(p->angle/unitangle);
                        idx = max(min(ndirections-1, idx), 0);
                        nanglespot[idx] ++;
                } else {
                        double goodang = 2.0 * PI;

                        if (r > rupper)
                        goodang -= 2.0 * acos(rupper/r);
                        if (r > rright)
                        goodang -= 2.0 * acos(rright/r);
                        if (r > rlower)
                        goodang -= 2.0 * acos(rlower/r);
                        if (r > rleft)
                        goodang -= 2.0 * acos(rleft/r);

                        if (goodang < 1e-2)
                                goodang = 1e-2;

                        spotcornerfactor.push_back(static_cast<int>(2.0*PI/goodang));
                }
        }


        // Duplicate spots whose cornerfactor is greater than 1.
        int nspots = accumulate(spotcornerfactor.begin(),spotcornerfactor.end(),0);
        if (nspots > spots.size()) {
                // Copy corner spots to make total number of spots the expected amount.
                spotresol.reserve(nspots);
                for (int sptidx=spots.size()-1; sptidx>=0; sptidx--) {
                        if (spotcornerfactor[sptidx]>1)
                                spotresol.insert(spotresol.begin()+sptidx,
                                                spotcornerfactor[sptidx]-1,spotresol[sptidx]);
                }
        }
        if (spotresol.size() < nspots) {
                cout << "Warning: spot amount inconsistent!\n";
                nspots = spotresol.size();
        }



        // resolution at the outer border of each ring.
        imgresol_ringresol = vector<double>(nimgresolrings);
        imgresol_ringresolpow = vector<double>(nimgresolrings);
        int idx0 = 0;
        for (int ringidx=0; ringidx<nimgresolrings; ringidx++) {
                int idx1 = static_cast<int>(ringidx * nspots / nimgresolrings);
                idx1 = min(nspots-1, max(idx1, 0));

                nth_element(spotresol.begin()+idx0, spotresol.begin()+idx1, spotresol.end());

                // Spots should have been ordered by reverse order of resolution.
                // So do it the following way.
                imgresol_ringresol[nimgresolrings-1-ringidx] = spotresol[idx1];
                imgresol_ringresolpow[nimgresolrings-1-ringidx]
                        = pow(imgresol_ringresol[nimgresolrings-1-ringidx], imgresolringpow);

                idx0 = idx1 + 1;
        }



        // Check directional distribution of spots.
        // They are expected to be located uniformly in all directions.
        // Ignore those spots in corners.
        // Comment this whole section out because it leads to a division-by-zero
        // error in the case where the direct beam is off the detector face.
        // The variable imgangularanomaly is not used by LABELIT anyway
        /*
        double nanglespotexpected
                = static_cast<double>(accumulate(nanglespot.begin(),nanglespot.end(),0)) / ndirections;
        imgangularanomaly = 0;
        for (vector<int>::const_iterator p=nanglespot.begin(); p!=nanglespot.end(); p++) {
                double bias;
                if (*p > nanglespotexpected)
                        bias = 1 - nanglespotexpected/(*p);
                else
                        bias = 1 - (*p)/nanglespotexpected;
                if (bias > imgangularanomaly)
                        imgangularanomaly = bias;
        }
        */

        /*
        vector<double> ringresolpowsmooth(nimgresolrings);

        vector<double> ringindices(nimgresolrings);
        for (int ringidx=0; ringidx<nimgresolrings; ringidx++)
                ringindices[ringidx] = static_cast<double>(ringidx);
        ksmooth(ringindices, ringresolpow, ringindices, ringresolpowsmooth, smoothspan);
        */



        // Slope of line connecting 1st and last rings.
        vector<double> slope(nimgresolrings);
        int maxslopeidx = 0;
        for (int ringidx=1; ringidx<nimgresolrings; ringidx++) {
                slope[ringidx] = (imgresol_ringresolpow[ringidx] - imgresol_ringresolpow[0]) / ringidx;
                if (slope[ringidx] > slope[maxslopeidx])
                        maxslopeidx = ringidx;
        }
        //int maxslopeidx = distance(slope.begin(), max_element(slope.begin(),slope.end()));

        double a = slope[maxslopeidx];
        double b = imgresol_ringresolpow[0];
        vector<double> gap(maxslopeidx+1);
        int bendidx = 0;
        for (int ringidx=0; ringidx<=maxslopeidx; ringidx++) {
                gap[ringidx] = a*ringidx+b - imgresol_ringresolpow[ringidx];
                if (gap[ringidx] > gap[bendidx])
                        bendidx = ringidx;
        }
        //int bendidx = distance(gap.begin(), max_element(gap.begin(),gap.end()));

        double temp1 = 0;
        double temp2 = 0;
        for (int ringidx=0; ringidx<=maxslopeidx; ringidx++) {
                temp1 += gap[ringidx];
                temp2 += gap[ringidx] * gap[ringidx];
        }
        temp1 /= (maxslopeidx + 1.0);
        temp2 /= (maxslopeidx + 1.0);
        double gapstd = sqrt(temp2 - temp1*temp1);


        double gapmax = gap[bendidx];
        for (int ringidx=bendidx; ringidx<=maxslopeidx; ringidx++) {
                if (gap[ringidx] >= gapmax - 0.5*gapstd)
                        bendidx = ringidx;
                else
                        break;
        }

        imgresol = imgresol_ringresol[bendidx];
}



double diffimage::r2_to_resol(const double r2) const
{
  // ********************************************
  // Determine resolution from squared radius.
  // This function is the reverse of resol_to_r2.
  // ********************************************

  // Calculation of resolution:
  //
  // d = lambda / (2*sin(theta))
  // sin(theta) = sqrt( 1/2 - 1/(2 * sqrt(1+(h/L)^2)) )
  // lambda: wavelength
  // h: distance from spot to beam center
  //    h^2 = ( (x-beam_x)^2 + (y-beam_y)^2 ) * pixel_size^2
  // L: distance from detector to crystal
  // ====>
  //    d = lambda / sqrt( 2 - 2/sqrt(1 + h^2/L^2) )
  //      = wavelength / sqrt( 2 - 2/sqrt(1 + pixel_size^2/L^2 *
  //                      ( (x-beam_x)^2 + (y-beam_y)^2 ) ) )
  //    b = pixel_size^2 / L^2
  //    c = (x-beam_x)^2 + (y-beam_y)^2
  //    d = sqrt(1 + c * b)
  //    two_sin_theta = sqrt(2 - 2/d)
  //    resolution = wavelength / two_sin_theta

  // resol = wavelength / sqrt(2 - 2/sqrt(1 + b*r^2) )
  // wavelength on the order of 1, b on the order of 1e-7.

  // for b = 1.6e-7, wavelength = 1,
  // 1/d^3 is approx proportional to r^2.8

  // d = wavelength / sqrt(2 - 2/(1 + b*r*r))
  // 1 / sqrt(1 + b*r*r) approx. 1 + b*r*r*(-1/2)
  // ==> d approx. wavelength/sqrt(b*r*r) = A / r
  // So resolution is approx. proportional to 1/r,
  // where r is distance from pixel to beam center.



  //static double resolb = pow(pixel_size / distance, 2);

  double d = sqrt(1.0 + resolb * r2);
  double two_sin_theta = sqrt(2 - 2/d);

  if (two_sin_theta<10e-8)
    return(10e8);

  return wavelength/two_sin_theta;

}


double diffimage::xy2resol(const double x, const double y) const
{
        double r2 = (x - beam_x)*(x - beam_x) + (y - beam_y)*(y - beam_y);
        return r2_to_resol(r2);
}

double diffimage::xy2resol_exact_normal(const double x, const double y) const
{
        double lin_radius = pixel_size *
                            std::sqrt((x - beam_x)*(x - beam_x) + (y - beam_y)*(y - beam_y));
        double two_sin_theta = 2. * std::sin(0.5*std::atan(lin_radius/distance));
        if (two_sin_theta<10e-8) {return(10e8);}
        return wavelength/two_sin_theta;
}


double diffimage::resol_to_r2(const double resol) const
{
  // ********************************************
  // Determine squared radius from resolution.
  // This function is the reverse of r2_to_resol.
  // ********************************************

  // static double b = pow(pixel_size / distance, 2);

  double d, two_sin_theta;

  two_sin_theta = wavelength / resol;
  // When this function is called properly,
  // 'resol' should be reasonably larger than 0.

  d = 2.0 / (2 - two_sin_theta * two_sin_theta);
  /* Exact formulae:
  double arg = std::tan(2.*std::asin(two_sin_theta/2.0));
  double should_be = arg*arg/resolb;
  */
  return (d*d - 1)/resolb;
}


template<class T>
int Distl::percentiles(const vector<T>& data, const vector<double>& x, vector<T>& prctile)
{
// X is a sorted vector in [0,1] specifying the percentages to be calculated.

        //if (!is_sorted(x.begin(),x.end()))
        //      return 1;

        vector<T> datacopy = data;
        int n = datacopy.size();

        int idx0 = 0;

        for (int i=0; i<x.size(); i++) {
                prctile[i] = 0; //initialize expecting default
                if (datacopy.size()==0) {continue;} //protect against empty array
                int idx1 = static_cast<int>((n-1) * x[i]);
                if (idx1 > n-1) idx1 = n - 1;
                if (idx1 < 0)  idx1 = 0;

                nth_element(datacopy.begin()+idx0, datacopy.begin()+idx1, datacopy.end());
                prctile[i] = datacopy[idx1];

                idx0 = idx1 + 1;
        }

        return 0;
}
