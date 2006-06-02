#ifndef CCTBX_CRYSTAL_WORKAROUNDS_BPL_H
#define CCTBX_CRYSTAL_WORKAROUNDS_BPL_H

#if defined(__sgi) && !defined(__GNUC__)

namespace cctbx { namespace crystal {

  namespace work_around_edg_typename_handling {
#if defined(CCTBX_CRYSTAL_DIRECT_SPACE_ASU_H)
    typedef direct_space_asu::asu_mappings<>
      asu_mappings_default;
    typedef boost::shared_ptr<direct_space_asu::asu_mappings<> >
      boost_shared_ptr_asu_mappings_default;
#endif
#if defined(CCTBX_CRYSTAL_PAIR_TABLES_H)
    typedef pair_asu_table<>
      pair_asu_table_default;
    typedef boost::shared_ptr<pair_asu_table<> >
      boost_shared_ptr_pair_asu_table_default;
#endif
  }

}} // namespace cctbx::crystal

#endif

#endif // CCTBX_CRYSTAL_WORKAROUNDS_BPL_H
