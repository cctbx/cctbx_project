#ifndef IOTBX_DET_CONTEXT_CONV_H
#define IOTBX_DET_CONTEXT_CONV_H

namespace iotbx {
namespace detectors {
namespace context {

namespace af = scitbx::af;

class spot_xy_convention{
  double W1,W2,pxlsz;
  int index;
public:
  spot_xy_convention(){ /*default*/  }
  spot_xy_convention(const double& width1,const double& width2,const double& px,const int& i){
    W1 = width1;
    W2 = width2;
    pxlsz = px; //pixel size in mm
    index = i;
    if (index==1 || index==3 || index==5 || index==7){
      SCITBX_ASSERT(W1==W2); // index swapping not permitted with rectangular image
    }
  }
  //based on code in python detectors/convention.py
  template <typename vec_t>
  af::tiny<vec_t,2 >  call(const scitbx::vec3<double>* xy) const {
    switch (index) {
      case 0: return af::tiny<vec_t,2 >(vec_t((*xy)[0]/pxlsz),vec_t((*xy)[1]/pxlsz));
      case 1: return af::tiny<vec_t,2 >(vec_t((*xy)[1]/pxlsz),vec_t((*xy)[0]/pxlsz));
      case 2: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[0]/pxlsz),vec_t((*xy)[1]/pxlsz));
      case 3: return af::tiny<vec_t,2 >(vec_t((*xy)[1]/pxlsz),vec_t(W1 - (*xy)[0]/pxlsz));
      case 4: return af::tiny<vec_t,2 >(vec_t((*xy)[0]/pxlsz),vec_t(W2 - (*xy)[1]/pxlsz));
      case 5: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[1]/pxlsz),vec_t((*xy)[0]/pxlsz));
      case 6: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[0]/pxlsz),vec_t(W2 - (*xy)[1]/pxlsz));
      case 7: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[1]/pxlsz),vec_t(W1 - (*xy)[0]/pxlsz));
      default: break;
    }
    throw;
  }
  //based on code in python detectors/convention.py
  //this set of overloads gives a one-to-one mapping of input to transformed pixels
  template <typename vec_t>
  af::tiny<vec_t,2 >  call(const scitbx::vec2<int>* xy) const {
    switch (index) {
      case 0: return af::tiny<vec_t,2 >(vec_t((*xy)[0]),vec_t((*xy)[1]));
      case 1: return af::tiny<vec_t,2 >(vec_t((*xy)[1]),vec_t((*xy)[0]));
      case 2: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[0] - 1),vec_t((*xy)[1]));
      case 3: return af::tiny<vec_t,2 >(vec_t((*xy)[1]),vec_t(W1 - (*xy)[0] - 1));
      case 4: return af::tiny<vec_t,2 >(vec_t((*xy)[0]),vec_t(W2 - (*xy)[1] - 1));
      case 5: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[1] - 1),vec_t((*xy)[0]));
      case 6: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[0] - 1),vec_t(W2 - (*xy)[1]- 1));
      case 7: return af::tiny<vec_t,2 >(vec_t(W1 - (*xy)[1] - 1),vec_t(W1 - (*xy)[0]- 1));
      default: break;
    }
    throw;
  }
};

} //namespace context
} //namespace detectors
} //namespace iotbx

#endif //IOTBX_DET_CONTEXT_CONV_H
