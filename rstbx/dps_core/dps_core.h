#ifndef RSTBX_DPS_CORE_H
#define RSTBX_DPS_CORE_H

#include <boost/shared_ptr.hpp>

#include <scitbx/array_family/flex_types.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>
#include <cctbx/crystal_orientation.h>

#include <rstbx/dps_core/directional_fft.h>
#include <rstbx/dps_core/spotclass.h>

namespace af = scitbx::af;
namespace rstbx {

typedef af::shared<scitbx::vec3<double> >    pointlistmm; //memory management
typedef boost::shared_ptr<Directional_FFT>   fftptr;
typedef cctbx::crystal_orientation           Orientation;
typedef af::shared<SpotClass>                dps_statuslist;
typedef af::shared<cctbx::miller::index<> >  hkllistmm;   //memory management

//! Autoindexing of observed reciprocal space vectors using DPS algorithm
/*! Reference:
    I. Steller, R. Bolotovsky, M.G. Rossmann
    An Algorithm for Automatic Indexing of Oscillation Images using
    Fourier Analysis.  J. Appl. Cryst. (1997) 30, 1036-1040.
 */
class dps_core {
 protected:
  double granularity; //! number of bins between lattice planes (n in Rossmann)
  double amax;        //! max cell in Angstroms
  pointlistmm xyzdata;//! original spots, reciprocal orthogonal coords xyz,
                      //!   expressed in inverse Angstroms
  af::shared<Direction> hemisphere_solutions;
  Orientation orientation;
  bool outliers_marked;// Boolean for indicating if outliers have been marked
 public:
  dps_statuslist status;  // scratch pad for spot status
  dps_statuslist Estatus; // scratch pad for spot status
  pointlistmm obsdata;// original spots, reciprocal skew coords expressed
                      // in the basis of the reciprocal unit cell, with the
                      // simplifying assumption that image is a still
                      // photograph collected at center of oscillation range;
  pointlistmm hkldata;// the whole number portion of obsdata
  void reset_spots(); // classify all as GOOD

  dps_core();
  void setMaxcell(const double &);
  void setXyzData(pointlistmm inputdata){ xyzdata = inputdata; }
  int getXyzSize() const {return xyzdata.size();}
  fftptr fft_factory(const Direction&) const;
  void setSolutions(af::shared<Direction>);
  af::shared<Direction> getSolutions() const;
  inline int n_candidates() const { return hemisphere_solutions.size(); }
  inline Direction candidate(const int &i) const {
                                    return hemisphere_solutions[i];}
  void setOrientation(const Orientation&);
  void set_orientation_direct_matrix(const matrix&);
  void set_orientation_reciprocal_matrix(const matrix&);
  inline Orientation getOrientation() const {return orientation;}
  virtual void classify_spots(); // find good, bad & ugly
  double rmsdev() const; // rmsd (real minus nearest whole),
                         // based on skew coordinates of GOOD spots;
  pointlistmm observed() // original spots, reciprocal skew coords expressed
                  const; // in the basis of the reciprocal unit cell, with the
                         // simplifying assumption that image is a still
                         // photograph collected at center of oscillation range;
                         // also zero mosaicity
  hkllistmm hklobserved() const; // the whole number portion of observed()
  hkllistmm hklobserved(const pointlistmm&) const; // choose an observed()
};

} //namespace
#endif //RSTBX_DPS_CORE_H
