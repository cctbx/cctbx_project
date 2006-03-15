#ifndef CCTBX_XRAY_REFINEMENT_FLAGS_H
#define CCTBX_XRAY_REFINEMENT_FLAGS_H

namespace cctbx { namespace xray {

  struct refinement_flags
  {
    static const unsigned site_bit =      0x00000001U;
    static const unsigned u_iso_bit =     0x00000002U;
    static const unsigned u_aniso_bit =   0x00000004U;
    static const unsigned occupancy_bit = 0x00000008U;
    static const unsigned fp_bit =        0x00000010U;
    static const unsigned fdp_bit =       0x00000020U;
    static const unsigned tan_u_iso_bit = 0x00000040U;

    unsigned bits;
    int param;

    refinement_flags(
      bool site=false,
      bool u_iso=false,
      bool u_aniso=false,
      bool occupancy=false,
      bool fp=false,
      bool fdp=false,
      bool tan_u_iso=false,
      int param_=0)
    :
      bits(0),
      param(param_)
    {
      set_site(site);
      set_u_iso(u_iso);
      set_u_aniso(u_aniso);
      set_occupancy(occupancy);
      set_fp(fp);
      set_fdp(fdp);
      set_tan_u_iso(tan_u_iso);
    }

#define CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(attr) \
    bool \
    attr() const { return bits & attr##_bit; } \
\
    void \
    set_##attr(bool state) \
    { \
      if (state) bits |= attr##_bit; \
      else       bits &= ~attr##_bit; \
    }

    CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(site)
    CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(u_iso)
    CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(u_aniso)
    CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(occupancy)
    CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(fp)
    CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(fdp)
    CCTBX_XRAY_REFINEMENT_FLAGS_GET_SET(tan_u_iso)
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_REFINEMENT_FLAGS_H
