#ifndef RSTBX_SIMPLE_INT_HPP
#define RSTBX_SIMPLE_INT_HPP
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#include <map>
#include <annlib_adaptbx/ann_adaptor.h>

#define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

namespace rstbx { namespace integration {

  template <typename NumType=int>
  class fast_less_than
  {
    public:
      //! This fast comparison function is implemented as operator().
      bool operator()(scitbx::vec2<NumType> const& h1,
                      scitbx::vec2<NumType> const& h2) const
      {
        for(std::size_t i=0;i<2;i++) {
          if (h1[i] < h2[i]) return true;
          if (h1[i] > h2[i]) return false;
        }
        return false;
      }
  };

  struct simple_integration {
    double pixel_size;
    scitbx::vec2<int> detector_size;
    typedef std::map<scitbx::vec2<int>, bool,fast_less_than<> > mask_t;
    typedef scitbx::af::shared<mask_t>         masks_t;
    masks_t ISmasks, BSmasks;
    scitbx::af::shared<scitbx::vec2<int> > detector_xy_draft;
    int FRAME;
    const int MAXOVER; //number of nearest neighbors
    double nbr_cutoff_sq;

    simple_integration(): MAXOVER(6){}

    void set_pixel_size(double const& pxsz) {pixel_size=pxsz;}

    void set_detector_size(int const& x, int const& y)
     {detector_size = scitbx::vec2<int>(x,y);}

    void set_frame(int const& f) {FRAME=f;}

    void set_nbr_cutoff_sq(double const& b) {nbr_cutoff_sq=b;}

    void append_ISmask(scitbx::af::shared<int> mask){
      mask_t newmask;
      for (int i=0; i < mask.size(); i+=2){
        newmask[scitbx::vec2<int>(mask[i],mask[i+1])]=true;
      }
      ISmasks.push_back(newmask);
    }

    scitbx::af::shared<int> get_bsmask(int const& i){
      scitbx::af::shared<int> return_value;
      for (mask_t::const_iterator k=BSmasks[i].begin();
                               k != BSmasks[i].end(); ++k){
        return_value.push_back(k->first[0]);
        return_value.push_back(k->first[1]);
      }
      return return_value;
    }

    scitbx::af::shared<scitbx::vec2<double> >
    safe_background(
      scitbx::af::shared<scitbx::vec3<double> > predicted,
      scitbx::af::shared<scitbx::vec2<double> > corrections,
      annlib_adaptbx::AnnAdaptor const& OS_adapt,
      scitbx::af::shared<int > flex_sorted
      ){
      int guard_width_sq =10;
      for (int i=0; i<predicted.size(); ++i){
        scitbx::vec3<double> pred = predicted[i];
        double predX = pred[0]/pixel_size;
        double predY = pred[1]/pixel_size;
        scitbx::vec2<double>correction = corrections[i];
        mask_t const& I_S_mask = ISmasks[i];
        // now consider the background
        mask_t B_S_mask;
        int i_bs = 0;
        scitbx::vec2<int> spot_position(
          round(predX + correction[0]),
          round(predY + correction[1]) );
        detector_xy_draft.push_back(spot_position);

        //insert a test to make sure spot is within FRAME
        if (spot_position[0] > FRAME && spot_position[1] > FRAME &&
            spot_position[0] < detector_size[0] - FRAME &&
            spot_position[1] < detector_size[1] - FRAME){

          mask_t spot_keys = I_S_mask;
          int base_spot_size = spot_keys.size();

          //Look for potential overlaps
          for (int n=0; n<MAXOVER; ++n){
            double distance = OS_adapt.distances[i*MAXOVER+n];
            if (distance < nbr_cutoff_sq){
              mask_t const& other_mask = ISmasks[ OS_adapt.nn[i*MAXOVER+n] ];
              for (mask_t::const_iterator k=other_mask.begin();
                k != other_mask.end(); ++k){
                spot_keys[k->first]=true;
              }
            }
          }

          for (int isort=0; isort<flex_sorted.size(); isort+=2){
            scitbx::vec2<int> increment(flex_sorted[isort],flex_sorted[isort+1]);
            scitbx::vec2<int> candidate_bkgd = spot_position + increment;
            if (spot_keys.find(candidate_bkgd)==spot_keys.end()){
             //eliminate if in guard region
             bool guard = false;
             for (mask_t::const_iterator key=spot_keys.begin();
               key != spot_keys.end(); ++key){
               int dx = candidate_bkgd[0]-(key->first)[0];
               int dy = candidate_bkgd[1]-(key->first)[1];
               if (dx*dx + dy*dy < guard_width_sq){
                 guard = true;
                 break;
               }
             }
             if (guard){ continue; }
             i_bs += 1;
             B_S_mask[candidate_bkgd] = true;
            }
            if (i_bs == base_spot_size){break;}
          }
        }
        BSmasks.push_back(B_S_mask);
      }
      scitbx::af::shared<scitbx::vec2<double> > return_detector_xy_draft;
      for (int i=0; i<detector_xy_draft.size(); ++i){
        return_detector_xy_draft.push_back(
          scitbx::vec2<double>(
            detector_xy_draft[i][0],detector_xy_draft[i][1]));
      }
      return return_detector_xy_draft;
    }

  };

}} // namespace rstbx::simage

#endif// RSTBX_SIMPLE_INT_HPP
