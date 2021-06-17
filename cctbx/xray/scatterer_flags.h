#ifndef CCTBX_XRAY_SCATTERER_FLAGS_H
#define CCTBX_XRAY_SCATTERER_FLAGS_H

#include <scitbx/array_family/ref.h>
#include <scitbx/sym_mat3.h>
#include <cctbx/import_scitbx_af.h>
#include <cctbx/error.h>

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
    static const unsigned use_fp_fdp_bit =               0x80000000U;

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
      bool use_fp_fdp=false,
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
      set_use_fp_fdp(use_fp_fdp);
    }

    /// Whether for each corresponding pair of bits a and b from bits
    // and other.bits, one has a => b (logical implication)
    bool implies (scatterer_flags const& other) {
      return (bits & ~other.bits) == 0;
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
    CCTBX_XRAY_SCATTERER_FLAGS_GET_SET(use_fp_fdp)

    //! An exception is thrown if u_iso and u_aniso are both true.
    bool use_u_iso_only() const
    {
      bool result = use_u_iso();
      if (result) {
        if (use_u_aniso()) {
          throw std::runtime_error(
            "scatterer.flags.u_iso_only(): u_iso and u_aniso both true.");
        }
      }
      else {
        if (!use_u_aniso()) {
          throw std::runtime_error(
            "scatterer.flags.u_iso_only(): u_iso and u_aniso both false.");
        }
      }
      return result;
    }

    //! An exception is thrown if u_iso and u_aniso are both true.
    bool use_u_aniso_only() const
    {
      bool result = use_u_aniso();
      if (result) {
        if (use_u_iso()) {
          throw std::runtime_error(
            "scatterer.flags.u_aniso_only(): u_iso and u_aniso both true.");
        }
      }
      else {
        if (!use_u_iso()) {
          throw std::runtime_error(
            "scatterer.flags.u_aniso_only(): u_iso and u_aniso both false.");
        }
      }
      return result;
    }

    void set_use_u_iso_only()
    {
      set_use_u_iso(true);
      set_use_u_aniso(false);
    };

    void set_use_u_aniso_only()
    {
      set_use_u_iso(false);
      set_use_u_aniso(true);
    };

    void set_use_u(bool iso, bool aniso)
    {
      set_use_u_iso(iso);
      set_use_u_aniso(aniso);
    };

    void
    set_grads(bool state)
    {
      set_grad_site(state);
      set_grad_u_iso(state);
      set_grad_u_aniso(state);
      set_grad_occupancy(state);
      set_grad_fp(state);
      set_grad_fdp(state);
    }
  };

  class grad_flags_counts_core
  {
    public:
      unsigned site;
      unsigned u_iso;
      unsigned u_aniso;
      unsigned u_anharmonic;
      unsigned occupancy;
      unsigned fp;
      unsigned fdp;
      unsigned tan_u_iso;
      unsigned use_u_iso;
      unsigned use_u_aniso;
      unsigned use_fp_fdp;

      grad_flags_counts_core()
        : site(0),
          u_iso(0),
          u_aniso(0),
          u_anharmonic(0),
          occupancy(0),
          fp(0),
          fdp(0),
          tan_u_iso(0),
          use_u_iso(0),
          use_u_aniso(0),
          use_fp_fdp(0)
      {}

      void process(scatterer_flags const &f)
      {
        if(f.use()) {
          if (f.grad_site()) site = site + 3;
          if (f.grad_u_iso() && f.use_u_iso()) u_iso++;
          if (f.grad_u_aniso() && f.use_u_aniso()) u_aniso = u_aniso+6;
          if (f.grad_occupancy()) occupancy++;
          if (f.grad_fp()) fp++;
          if (f.grad_fdp()) fdp++;
          if (f.tan_u_iso()) tan_u_iso++;
          if (f.use_u_iso()) use_u_iso++;
          if (f.use_u_aniso()) use_u_aniso++;
          if (f.use_fp_fdp()) use_fp_fdp++;
        }
      }

      unsigned
      n_parameters() const {
        return site + u_iso + u_aniso + u_anharmonic + occupancy + fp + fdp;
      }
  };


  class scatterer_grad_flags_counts : public grad_flags_counts_core
  {
    public:
      scatterer_grad_flags_counts() {}

      template <typename ScattererType>
      scatterer_grad_flags_counts(
                        af::const_ref<ScattererType> const& scatterers)
        : grad_flags_counts_core()
      {
        for(std::size_t i=0;i<scatterers.size();i++) {
          process(scatterers[i].flags);
          if (scatterers[i].flags.grad_u_aniso() &&
            scatterers[i].flags.use_u_aniso() &&
            scatterers[i].anharmonic_adp)
          {
            u_anharmonic += 25;
          }
        }
      }
  };


  class grad_flags_counts : public grad_flags_counts_core
  {
    public:
      grad_flags_counts(af::const_ref<scatterer_flags> const &flags) {
        for (std::size_t i=0; i < flags.size(); ++i) process(flags[i]);
      }
  };

  /// Set the flags of each given scatterers to the given values
  /** The logic is as follow:
        -# if a scatterer is not used (the 'use' flag is false), no flag is set;
        -# the grad_u_iso flag (as well as the tan_u_iso flag)
             -# is set to the specified value only if the use_u_iso flag is set;
             -# otherwise it is set to false;
        -# the same logic applies for grad_u_aniso and use_u_aniso;
        -# the use_fp_fdp flag
             -# is set iff either grad_fp or grad_fdp is set;
             -# otherwise it is left unchanged.

   */
  template <typename ScattererType>
  void
  set_scatterer_grad_flags(
                    af::ref<ScattererType> const& scatterers,
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
           if(sc.flags.use_u_iso()) {
              sc.flags.set_grad_u_iso(u_iso);
              CCTBX_ASSERT(sc.u_iso != -1.0);
           }
           else {
              sc.flags.set_grad_u_iso(false);
           }
           if(sc.flags.use_u_aniso()) {
              sc.flags.set_grad_u_aniso(u_aniso);
              CCTBX_ASSERT(
                     sc.u_star != scitbx::sym_mat3<double>(-1,-1,-1,-1,-1,-1));
           }
           else {
              sc.flags.set_grad_u_aniso(false);
           }
           sc.flags.set_grad_occupancy(occupancy);
           sc.flags.set_use_fp_fdp(fp || fdp);
           sc.flags.set_grad_fp(fp);
           sc.flags.set_grad_fdp(fdp);
           if(sc.flags.use_u_iso()) sc.flags.set_tan_u_iso(tan_u_iso);
           sc.flags.param = param;
        }
    }
  }

  template <typename ScattererType>
  void
  flags_set_grads(
    af::ref<ScattererType> const& self,
    bool state)
  {
    for(std::size_t i=0;i<self.size();i++) {
      self[i].flags.set_grads(state);
    }
  }

#define CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD(attr) \
  template <typename ScattererType> \
  void \
  flags_set_grad_##attr( \
    af::ref<ScattererType> const& self, \
    af::const_ref<std::size_t> const& iselection) \
  { \
    for(std::size_t i=0;i<iselection.size();i++) { \
      std::size_t i_seq = iselection[i]; \
      CCTBX_ASSERT(i_seq < self.size()); \
      self[i_seq].flags.set_grad_##attr(true); \
    } \
  }

#define CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD_U(attr) \
  template <typename ScattererType> \
  void \
  flags_set_grad_##attr( \
    af::ref<ScattererType> const& self, \
    af::const_ref<std::size_t> const& iselection) \
  { \
    for(std::size_t i=0;i<iselection.size();i++) { \
      std::size_t i_seq = iselection[i]; \
      CCTBX_ASSERT(i_seq < self.size()); \
      scatterer_flags& f = self[i_seq].flags; \
      CCTBX_ASSERT(f.use_##attr()); \
      f.set_grad_##attr(true); \
    } \
  }

  CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD(site)
  CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD_U(u_iso)
  CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD_U(u_aniso)
  CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD(occupancy)
  CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD(fp)
  CCTBX_XRAY_SCATTERERS_FLAGS_SET_GRAD(fdp)

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_SCATTERER_FLAGS_H
