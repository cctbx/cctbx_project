#ifndef CCTBX_XRAY_GRADIENT_FLAGS_H
#define CCTBX_XRAY_GRADIENT_FLAGS_H

namespace cctbx { namespace xray {

  struct gradient_flags
  {
    gradient_flags(
      bool site_,
      bool u_iso_,
      bool u_aniso_,
      bool occupancy_,
      bool fp_,
      bool fdp_)
    :
      site(site_),
      u_iso(u_iso_),
      u_aniso(u_aniso_),
      occupancy(occupancy_),
      fp(fp_),
      fdp(fdp_)
    {}

    bool site;
    bool u_iso;
    bool u_aniso;
    bool occupancy;
    bool fp;
    bool fdp;

    bool
    all_false() const
    {
      return !(site || u_iso || u_aniso || occupancy || fp || fdp);
    }

    gradient_flags
    adjust(bool anisotropic_flag) const
    {
      gradient_flags result = *this;
      if (anisotropic_flag) result.u_iso = false;
      if (!anisotropic_flag) result.u_aniso = false;
      return result;
    }
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_GRADIENT_FLAGS_H
