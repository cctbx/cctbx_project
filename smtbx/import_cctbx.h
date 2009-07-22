#ifndef SMTBX_IMPORT_CCTBX_H
#define SMTBX_IMPORT_CCTBX_H

namespace cctbx {
  namespace adptbx {}
  namespace uctbx {}
  namespace sgtbx {}
  namespace xray {}
  namespace miller {}
}

/// Small Molecule Toolbox namespace
namespace smtbx {
  namespace adptbx = cctbx::adptbx;
  namespace uctbx = cctbx::uctbx;
  namespace sgtbx = cctbx::sgtbx;
  namespace xray = cctbx::xray;
  namespace miller = cctbx::miller;
  using cctbx::cartesian;
  using cctbx::fractional;
}

#endif
