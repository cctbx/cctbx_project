#pragma once

#include <cctbx/xray/thickness.h>
#include <smtbx/ED/utils.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/ED/dyn_calculator.h>
#include <smtbx/ED/n_beam.h>
#include <smtbx/refinement/least_squares_fc.h>
#include <cctbx/xray/fc_correction.h>
#include <boost/scoped_ptr.hpp>

namespace smtbx { namespace ED {
  using namespace smtbx::refinement::least_squares;

  template <typename FloatType>
  struct beam_group_profiler : public N_beam_shared_data_base<FloatType> {
    ED_UTIL_TYPEDEFS;
    N_NEAM_SHARED_DATA_TYPEDES;
    typedef N_beam_shared_data_base<FloatType> base_t;

    beam_group_profiler(const BeamGroup<FloatType>& beam_group,
      f_calc_function_base_t& f_calc_function,
      const sgtbx::space_group& space_group,
      const uctbx::unit_cell& unit_cell,
      const RefinementParams<FloatType>& params,
      const cctbx::xray::thickness<FloatType>& thickness)
      : base_t(f_calc_function, space_group, unit_cell, thickness, params),
      beam_group(beam_group),
      mat_type(params.getMatrixType())
    {
      if (mat_type >= 100) {
        mat_type -= 100;
      }

      af::shared<BeamGroup<FloatType> > beam_groups;
      beam_groups.push_back(beam_group);
      base_t::init(beam_groups);

      base_t::build();
    }

    struct process_beam_group_profile {
      process_beam_group_profile(const beam_group_profiler& parent,
        const BeamGroup<FloatType>& beam_group)
        : parent(parent),
        beam_group(beam_group)
      {}

      // N-Beam specific
      af::shared<FloatType> process(const miller::index<> &h, const af::shared<FloatType> &angles) {
        af::shared<FloatType> rv(angles.size());
        try {
          bool floating = parent.params.isNBeamFloating();
          dyn_calculator_n_beam<FloatType> n_beam_dc(
            parent, beam_group, parent.mat_type);
          if (!floating) {
            FloatType da = beam_group.get_diffraction_angle(h, parent.K);
            n_beam_dc.init(h, da, parent.Fcs_kin, parent.mi_lookup);
          }
          else {
            n_beam_dc.init(h, angles);
          }
          for (size_t ai = 0; ai < angles.size(); ai++) {
            mat3_t R = beam_group.get_R(angles[ai]);
            if (floating) {
              n_beam_dc.init(h, R, parent.Fcs_kin, parent.mi_lookup);
            }
            rv[ai] = std::norm(n_beam_dc.calc_amp(R));
          }
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
        return rv;
      }

      af::shared<FloatType> process(const mat3_t& RMf) {
        af::shared<FloatType> Is(beam_group.beams.size());
        try {
          /* a LOT of overhead here!!! may need to change the logic to speed up
            as the matrix rebuilt for each angle
            */
          if (parent.params.getBeamN() > 2) {
            bool floating = parent.params.isNBeamFloating();
            for (size_t bi = 0; bi < beam_group.beams.size(); bi++) {
              const miller::index<>& h = beam_group.beams[bi].h;
              dyn_calculator_n_beam<FloatType> n_beam_dc(
                parent, beam_group, parent.mat_type);
              FloatType da = beam_group.get_diffraction_angle(h, parent.K);
              mat3_t da_R = beam_group.get_R(da);
              if (floating) {
                n_beam_dc.init(h, RMf, parent.Fcs_kin, parent.mi_lookup);
              }
              else {
                n_beam_dc.init(h, da_R, parent.Fcs_kin, parent.mi_lookup);
              }
              Is[bi] = std::norm(n_beam_dc.calc_amp(RMf));
            }
          }
          else {
            for (size_t bi = 0; bi < beam_group.beams.size(); bi++) {
              const miller::index<>& h = beam_group.beams[bi].h;
              int ii = parent.mi_lookup.find_hkl(h);
              complex_t Fc = ii != -1 ? parent.Fcs_kin[ii] : 0;
              Is[bi] = std::norm(utils<FloatType>::calc_amp_2beam(
                h, Fc,
                parent.thickness.value,
                parent.K, RMf, beam_group.get_N()));
            }
          }
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
        return Is;
      }

      FloatType process_incident(const mat3_t& RMf) {
        try {
          FloatType I_sum = 0;
          if (parent.params.getBeamN() > 2) {
            bool floating = parent.params.isNBeamFloating();
            for (size_t bi = 0; bi < beam_group.beams.size(); bi++) {
              const miller::index<> &h = beam_group.beams[bi].h;
              dyn_calculator_n_beam<FloatType> n_beam_dc(
                parent, beam_group, parent.mat_type);
              if (floating) {
                n_beam_dc.init(h, RMf, parent.Fcs_kin, parent.mi_lookup);
              }
              else {
                FloatType da = beam_group.get_diffraction_angle(h, parent.K);
                n_beam_dc.init(h, da, parent.Fcs_kin, parent.mi_lookup);
              }
              I_sum += std::norm(n_beam_dc.calc_amp(RMf, 0));
            }
          }
          else {
            for (size_t bi = 0; bi < beam_group.beams.size(); bi++) {
              const miller::index<>& h = beam_group.beams[bi].h;
              int ii = parent.mi_lookup.find_hkl(h);
              complex_t Fc = ii != -1 ? parent.Fcs_kin[ii] : 0;
              I_sum += std::norm(utils<FloatType>::calc_amp_2beam(
                h, Fc,
                parent.thickness.value,
                parent.K, RMf, beam_group.get_N(), 0));
            }
          }
          return I_sum / beam_group.beams.size();
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
        // suppress warning
        return FloatType();
      }
      boost::shared_ptr< dyn_calculator_n_beam<FloatType> > n_beam_dc;
      const beam_group_profiler& parent;
      const BeamGroup<FloatType>& beam_group;
      FloatType thickness;
      boost::scoped_ptr<smtbx::error> exception_;
    };

    af::shared<PeakProfilePoint<FloatType> > build_group_profile(
      const af::shared<FloatType>& angles, bool inc_incident)
    {
      af::shared<PeakProfilePoint<FloatType> > rv(
        af::reserve((angles.size() + (inc_incident ? 1 : 0))
          * beam_group.beams.size()));
      for (size_t ai = 0; ai < angles.size(); ai++) {
        FloatType ang = angles[ai];
        mat3_t R = beam_group.get_R(ang);
        process_beam_group_profile fp(*this, beam_group);
        af::shared<FloatType> Is = fp.process(R);
        if (fp.exception_) {
          throw* fp.exception_.get();
        }
        if (inc_incident) {
          process_beam_group_profile ibp(*this, beam_group);
          FloatType I = ibp.process_incident(R);
          rv.push_back(PeakProfilePoint<FloatType>(I, 0.0, ang, this->Kl));
        }
        for (size_t bi = 0; bi < beam_group.beams.size(); bi++) {
          const miller::index<>& h = beam_group.beams[bi].h;
          cart_t g = R * cart_t(h[0], h[1], h[2]);
          FloatType Sg = utils<FloatType>::calc_Sg(g, this->K);
          rv.push_back(PeakProfilePoint<FloatType>(Is[bi], Sg, ang, (this->K + g).length()));
        }
      }
      return rv;
    }

    af::shared<PeakProfilePoint<FloatType> > build_reflection_profile(
      const miller::index<> &h,
      const af::shared<FloatType>& angles)
    {
      af::shared<PeakProfilePoint<FloatType> > rv(af::reserve(angles.size()));
      process_beam_group_profile fp(*this, beam_group);
      af::shared<FloatType> Is = fp.process(h, angles);
      if (fp.exception_) {
        throw* fp.exception_.get();
      }
      cart_t h_c(h[0], h[1], h[2]);
      for (size_t ai = 0; ai < angles.size(); ai++) {
        FloatType ang = angles[ai];
        cart_t g = beam_group.get_R(ang) * h_c;
        FloatType Sg = utils<FloatType>::calc_Sg(g, this->K);
        rv.push_back(PeakProfilePoint<FloatType>(Is[ai], Sg, ang, (this->K + g).length()));
      }
      return rv;
    }

    af::shared<PeakProfilePoint<FloatType> > build_incident_profile(
      const af::shared<FloatType>& angles)
    {
      af::shared<PeakProfilePoint<FloatType> > rv(
        af::reserve(angles.size()));
      for (size_t ai = 0; ai < angles.size(); ai++) {
        FloatType ang = angles[ai];
        process_beam_group_profile fp(*this, beam_group);
        FloatType I = fp.process_incident(beam_group.get_R(ang));
        if (fp.exception_) {
          throw* fp.exception_.get();
        }
        rv.push_back(PeakProfilePoint<FloatType>(I, 0.0, ang, this->Kl));
      }
      return rv;
    }

    const BeamGroup<FloatType>& beam_group;
    int mat_type;
    cmat_t A;
  };
}} // end of smtbx::ED
