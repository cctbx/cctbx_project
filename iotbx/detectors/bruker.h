#ifndef IOTBX_DETECTORS_BRUKER_H
#define IOTBX_DETECTORS_BRUKER_H
#include <scitbx/array_family/boost_python/flex_fwd.h>

/* brukccd.c from ccd2ipf.c, Stan Swanson, TAMU, 25 April 2006 */
/* proteus.c conversion 033.26 */
/* based on difoff.c */
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

#include <vector>
#include <scitbx/array_family/flex_types.h>

namespace iotbx { namespace detectors {

struct bruker{
  void v_alloc();
  bruker(std::string);
  int v_read(const char*);

  std::vector<std::vector<int> > ccdata;
  std::vector<unsigned char> v_ccbyte;
  std::vector<unsigned short> v_ccshort;
  std::vector<unsigned int> v_ccint;

  unsigned char* ccbyte;
  unsigned short* ccshort;
  unsigned int* ccint;

  int maxpixel,saturate;
  double pixsizemm,distance,delta,wavelen,centerx,centery,oscrange,twoth;
  double osc_start;
  scitbx::af::flex_int linearintdata();
};

}} // namespace iotbx::detectors
#endif // IOTBX_DETECTORS_BRUKER_H
