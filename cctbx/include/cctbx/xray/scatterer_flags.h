#ifndef CCTBX_XRAY_SCATTERER_FLAGS_H
#define CCTBX_XRAY_SCATTERER_FLAGS_H

#include <scitbx/array_family/ref.h>

namespace cctbx { namespace xray {

  struct scatterer_flags
  {
    static const unsigned use_bit =                      0x00000001U;
    static const unsigned use_u_iso_bit =                0x00000002U;
    static const unsigned use_u_aniso_bit =              0x00000004U;
    static const unsigned grad_site_bit =                0x00000008U;
    static const unsigned grad_u_iso_bit =               0x00000010U;
    static const unsigned grad_u_aniso_bit =             0x00000020U;
    static const unsigned grad_occupancy_bit =           0x00000040U;
    static const unsigned grad_fp_bit =                  0x00000080U;
    static const unsigned grad_fdp_bit =                 0x00000100U;
    static const unsigned curv_site_site_bit =           0x00000200U;
    static const unsigned curv_site_u_iso_bit =          0x00000400U;
    static const unsigned curv_site_u_aniso_bit =        0x00000800U;
    static const unsigned curv_site_occupancy_bit =      0x00001000U;
    static const unsigned curv_site_fp_bit =             0x00002000U;
    static const unsigned curv_site_fdp_bit =            0x00004000U;
    static const unsigned curv_u_iso_u_iso_bit =         0x00008000U;
    static const unsigned curv_u_iso_u_aniso_bit =       0x00010000U;
    static const unsigned curv_u_iso_occupancy_bit =     0x00020000U;
    static const unsigned curv_u_iso_fp_bit =            0x00040000U;
    static const unsigned curv_u_iso_fdp_bit =           0x00080000U;
    static const unsigned curv_u_aniso_u_aniso_bit =     0x00100000U;
    static const unsigned curv_u_aniso_occupancy_bit =   0x00200000U;
    static const unsigned curv_u_aniso_fp_bit =          0x00400000U;
    static const unsigned curv_u_aniso_fdp_bit =         0x00800000U;
    static const unsigned curv_occupancy_occupancy_bit = 0x01000000U;
    static const unsigned curv_occupancy_fp_bit =        0x02000000U;
    static const unsigned curv_occupancy_fdp_bit =       0x04000000U;
    static const unsigned curv_fp_fp_bit =               0x08000000U;
    static const unsigned curv_fp_fdp_bit =              0x10000000U;
    static const unsigned curv_fdp_fdp_bit =             0x20000000U;
    static const unsigned tan_u_iso_bit =                0x40000000U;

    unsigned bits;
    int param;

    scatterer_flags(
      unsigned bits_=use_bit,
      int param_=0)
    :
      bits(bits_),
      param(param_)
    {}

    scatterer_flags(
      bool use,
      bool use_u_iso=false,
      bool use_u_aniso=false,
      bool grad_site=false,
      bool grad_u_iso=false,
      bool grad_u_aniso=false,
      bool grad_occupancy=false,
      bool grad_fp=false,
      bool grad_fdp=false,
      bool curv_site_site=false,
      bool curv_site_u_iso=false,
      bool curv_site_u_aniso=false,
      bool curv_site_occupancy=false,
      bool curv_site_fp=false,
      bool curv_site_fdp=false,
      bool curv_u_iso_u_iso=false,
      bool curv_u_iso_u_aniso=false,
      bool curv_u_iso_occupancy=false,
      bool curv_u_iso_fp=false,
      bool curv_u_iso_fdp=false,
      bool curv_u_aniso_u_aniso=false,
      bool curv_u_aniso_occupancy=false,
      bool curv_u_aniso_fp=false,
      bool curv_u_aniso_fdp=false,
      bool curv_occupancy_occupancy=false,
      bool curv_occupancy_fp=false,
      bool curv_occupancy_fdp=false,
      bool curv_fp_fp=false,
      bool curv_fp_fdp=false,
      bool curv_fdp_fdp=false,
      bool tan_u_iso=false,
      int  param_=0)
    :
      bits(0),
      param(param_)
    {
      set_use(use);
      set_use_u_iso(use_u_iso);
      set_use_u_aniso(use_u_aniso);
      set_grad_site(grad_site);
      set_grad_u_iso(grad_u_iso);
      set_grad_u_aniso(grad_u_aniso);
      set_grad_occupancy(grad_occupancy);
      set_grad_fp(grad_fp);
      set_grad_fdp(grad_fdp);
      set_curv_site_site(curv_site_site);
      set_curv_site_u_iso(curv_site_u_iso);
      set_curv_site_u_aniso(curv_site_u_aniso);
      set_curv_site_occupancy(curv_site_occupancy);
      set_curv_site_fp(curv_site_fp);
      set_curv_site_fdp(curv_site_fdp);
      set_curv_u_iso_u_iso(curv_u_iso_u_iso);
      set_curv_u_iso_u_aniso(curv_u_iso_u_aniso);
      set_curv_u_iso_occupancy(curv_u_iso_occupancy);
      set_curv_u_iso_fp(curv_u_iso_fp);
      set_curv_u_iso_fdp(curv_u_iso_fdp);
      set_curv_u_aniso_u_aniso(curv_u_aniso_u_aniso);
      set_curv_u_aniso_occupancy(curv_u_aniso_occupancy);
      set_curv_u_aniso_fp(curv_u_aniso_fp);
      set_curv_u_aniso_fdp(curv_u_aniso_fdp);
      set_curv_occupancy_occupancy(curv_occupancy_occupancy);
      set_curv_occupancy_fp(curv_occupancy_fp);
      set_curv_occupancy_fdp(curv_occupancy_fdp);
      set_curv_fp_fp(curv_fp_fp);
      set_curv_fp_fdp(curv_fp_fdp);
      set_curv_fdp_fdp(curv_fdp_fdp);
      set_tan_u_iso(tan_u_iso);
    }

#define CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(attr) \
    bool \
    attr() const { return bits & attr##_bit; } \
\
    scatterer_flags& \
    set_##attr(bool state) \
    { \
      if (state) bits |= attr##_bit; \
      else       bits &= ~attr##_bit; \
      return *this; \
    }

    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(use)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(use_u_iso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(use_u_aniso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(grad_site)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(grad_u_iso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(grad_u_aniso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(grad_occupancy)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(grad_fp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(grad_fdp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_site_site)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_site_u_iso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_site_u_aniso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_site_occupancy)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_site_fp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_site_fdp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_iso_u_iso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_iso_u_aniso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_iso_occupancy)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_iso_fp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_iso_fdp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_aniso_u_aniso)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_aniso_occupancy)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_aniso_fp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_u_aniso_fdp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_occupancy_occupancy)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_occupancy_fp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_occupancy_fdp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_fp_fp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_fp_fdp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(curv_fdp_fdp)
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(tan_u_iso)

    void set_use_u(bool iso, bool aniso)
    {
      set_use_u_iso(iso);
      set_use_u_aniso(aniso);
    };

    void set_use_u(bool state)
    {
      set_use_u_iso(state);
      set_use_u_aniso(!state);
    };

  };

  class scatterer_grad_flags_counts {
    public:
      unsigned site;
      unsigned u_iso;
      unsigned u_aniso;
      unsigned occupancy;
      unsigned fp;
      unsigned fdp;
      unsigned tan_u_iso;
      unsigned use_u_iso;
      unsigned use_u_aniso;

      scatterer_grad_flags_counts() {}

      template <typename ScattererType>
      scatterer_grad_flags_counts(
                        scitbx::af::const_ref<ScattererType> const& scatterers)
      :
        site(0),
        u_iso(0),
        u_aniso(0),
        occupancy(0),
        fp(0),
        fdp(0),
        tan_u_iso(0),
        use_u_iso(0),
        use_u_aniso(0)
      {
        for(std::size_t i=0;i<scatterers.size();i++) {
            ScattererType const& sc = scatterers[i];
            if(sc.flags.use()) {
              if (sc.flags.grad_site()) site = site + 3;
              if (sc.flags.grad_u_iso() && sc.flags.use_u_iso()) u_iso++;
              if (sc.flags.grad_u_aniso() &&
                  sc.flags.use_u_aniso()) u_aniso = u_aniso+6;
              if (sc.flags.grad_occupancy()) occupancy++;
              if (sc.flags.grad_fp()) fp++;
              if (sc.flags.grad_fdp()) fdp++;
              if (sc.flags.tan_u_iso()) tan_u_iso++;
              if(sc.flags.use_u_iso()) use_u_iso++;
              if(sc.flags.use_u_aniso()) use_u_aniso++;
            }
        }
      }

      unsigned
      n_parameters() const { return site+u_iso+u_aniso+occupancy+fp+fdp; }
  };

  template <typename ScattererType>
  void
  set_scatterer_grad_flags(
                    scitbx::af::ref<ScattererType> const& scatterers,
                    bool site      = false,
                    bool u_iso     = false,
                    bool u_aniso   = false,
                    bool occupancy = false,
                    bool fp        = false,
                    bool fdp       = false,
                    bool tan_u_iso = false,
                    int  param     = 0)
  {
    for(std::size_t i=0;i<scatterers.size();i++) {
        ScattererType& sc = scatterers[i];
        if(sc.flags.use()) {
          sc.flags.set_grad_site(site);
          if (sc.flags.use_u_iso())
            sc.flags.set_grad_u_iso(u_iso);
          else
            sc.flags.set_grad_u_iso(false);
          if (sc.flags.use_u_aniso())
            sc.flags.set_grad_u_aniso(u_aniso);
          else
            sc.flags.set_grad_u_aniso(false);
          sc.flags.set_grad_occupancy(occupancy);
          sc.flags.set_grad_fp(fp);
          sc.flags.set_grad_fdp(fdp);
          if (sc.flags.use_u_iso()) sc.flags.set_tan_u_iso(tan_u_iso);
          sc.flags.param = param;
        }
    }
  }

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERER_FLAGS_H
