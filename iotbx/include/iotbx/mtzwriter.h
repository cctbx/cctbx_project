#ifndef MTZWRITER_H
#define MTZWRITER_H

#include <iostream>
#include <exception>
#include <string>

#include <sys/stat.h> // a cure for errors in irix compile
#include "cmtzlib.h"
#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#include <cctbx/hendrickson_lattman.h>
#include <cctbx/uctbx.h>
#include <iotbx/cppmtz.h>
#include <cctbx/sgtbx/space_group.h>

namespace af = scitbx::af;

namespace iotbx{ namespace mtz {

class MtzWriter {
private:
  CMtz::MTZ* mtz;
  CMtz::MTZXTAL* onextal;
  CMtz::MTZSET*  oneset;
  void safe_ccp4_lwrefl(const float*, CMtz::MTZCOL **,
           const int, const int);
public:
  MtzWriter();
  ~MtzWriter();

  void setTitle(const std::string&);
  void setSpaceGroup(const cctbx::sgtbx::space_group&, const std::string&);
  void oneCrystal(const std::string&,const std::string&,
                  const cctbx::uctbx::unit_cell&);
  void oneDataset(const std::string&,const double&);
  void addColumn(
    const std::string& name,
    char type_code,
    af::const_ref<cctbx::miller::index<> > const& miller_indices,
    af::const_ref<double> const& data);
  void write(const std::string&);
};
}} //namespaces

#endif /* MTZWRITER_H*/
