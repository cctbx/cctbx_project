#ifndef FEM_COMMON_HPP
#define FEM_COMMON_HPP

#include <fem/close_chain.hpp>
#include <fem/file_positioning_chain.hpp>
#include <fem/inquire_chain.hpp>
#include <fem/open_chain.hpp>
#include <fem/utils/misc.hpp>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

namespace fem {

  struct common
  {
    fem::io io;
  };

  typedef boost::any cmn_sve;

  using utils::one_time_flag;

} // namespace fem

#define FEM_CMN_SVE(FUNC) \
  if (cmn.FUNC##_sve.empty()) { \
    cmn.FUNC##_sve = boost::shared_ptr<FUNC##_save>( \
      new FUNC##_save); \
  } \
  FUNC##_save& sve = *boost::any_cast< \
    boost::shared_ptr<FUNC##_save> >(cmn.FUNC##_sve).get()

#define FEM_CMN_SVE_DYNAMIC_PARAMETERS(FUNC) \
  if (cmn.FUNC##_sve.empty()) { \
    cmn.FUNC##_sve = boost::shared_ptr<FUNC##_save>( \
      new FUNC##_save(cmn.dynamic_parameters)); \
  } \
  FUNC##_save& sve = *boost::any_cast< \
    boost::shared_ptr<FUNC##_save> >(cmn.FUNC##_sve).get()

#endif // GUARD
