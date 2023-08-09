#pragma once
#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/ED/utils.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/refinement/constraints/reparametrisation.h>
#include <scitbx/array_family/selections.h>

namespace smtbx {
  namespace ED {
    using namespace refinement::least_squares;
    template <typename FloatType>
    struct a_frame_processor {
      ED_UTIL_TYPEDEFS;
      a_frame_processor(const FrameInfo<FloatType>& frame,
        const cmat_t& Ugs,
        const cart_t& K,
        const cctbx::xray::thickness<FloatType>& thickness,
        // kinematic dFc_dP
        const af::shared<cmat_t>& Ds_kin,
        bool calc_grad)
        : frame(frame),
        Ugs(Ugs),
        K(K),
        thickness(thickness),
        Ds_kin(Ds_kin),
        calc_grad(calc_grad)
      {}

      void operator()() {
        process(frame.RMf, frame.normal);
      }

      virtual void process(const mat3_t& RM, const cart_t& N) = 0;
      virtual void process_1(size_t i, const mat3_t& RM, const cart_t& N) = 0;

      const FrameInfo<FloatType>& frame;
      cmat_t Ugs;
      const cart_t& K;
      const cctbx::xray::thickness<FloatType>& thickness;
      af::shared<cmat_t> Ds_kin;
      bool calc_grad;
      boost::scoped_ptr<smtbx::error> exception_;
      // results
      af::shared<complex_t> CIs;
      mat_t D_dyn;
    };

    template <typename FloatType>
    struct process_frame_2013 : public a_frame_processor<FloatType> {
      ED_UTIL_TYPEDEFS;
      typedef a_frame_processor<FloatType> fp_t;
      process_frame_2013(FrameInfo<FloatType>& frame,
        const cmat_t& Ugs,
        const cart_t& K,
        const cctbx::xray::thickness<FloatType>& thickness,
        // kinematic dFc_dP
        const af::shared<cmat_t>& Ds_kin,
        bool calc_grad)
        : fp_t(frame,
          Ugs, K, thickness,
          Ds_kin,
          calc_grad),
        strong_indices(af::select(frame.indices.const_ref(),
          frame.strong_beams.const_ref()))
      {}

      void process(const mat3_t& RM, const cart_t& N) {
        try {
          size_t beam_n = fp_t::frame.strong_measured_beams.size();
          af::shared<FloatType> ExpDen;
          cmat_t A = fp_t::Ugs.deep_copy();
          utils<FloatType>::build_eigen_matrix_2013(A, strong_indices,
            fp_t::K,
            RM, N, ExpDen);
          if (fp_t::calc_grad) {
            fp_t::CIs = utils<FloatType>::calc_amps_2013_ext(A,
              fp_t::Ds_kin, ExpDen,
              fp_t::thickness.value,
              fp_t::thickness.grad,
              fp_t::D_dyn,
              fp_t::frame.strong_measured_beams.size());
          }
          else {
            fp_t::CIs = utils<FloatType>::calc_amps_2013(A,
              ExpDen, fp_t::thickness.value,
              fp_t::frame.strong_measured_beams.size());
          }
        }
        catch (smtbx::error const& e) {
          fp_t::exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          fp_t::exception_.reset(new smtbx::error(e.what()));
        }
      }
      
      void process_1(size_t idx,const mat3_t& RM, const cart_t& N) {
        try {
          size_t beam_n = fp_t::frame.strong_measured_beams.size();
          af::shared<FloatType> ExpDen;
          cmat_t A = fp_t::Ugs.deep_copy();
          utils<FloatType>::build_eigen_matrix_2013(A, strong_indices,
            fp_t::K,
            RM, N, ExpDen);
          if (fp_t::calc_grad) {
            fp_t::CIs.resize(1);
            fp_t::CIs[0] = utils<FloatType>::calc_amps_2013_ext_1(A,
              fp_t::Ds_kin, ExpDen,
              fp_t::thickness.value,
              fp_t::thickness.grad,
              fp_t::D_dyn,
              idx);
          }
          else {
            fp_t::CIs.resize(1);
            fp_t::CIs[0] = utils<FloatType>::calc_amps_2013_1(A,
              ExpDen, fp_t::thickness.value,
              idx);
          }
        }
        catch (smtbx::error const& e) {
          fp_t::exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          fp_t::exception_.reset(new smtbx::error(e.what()));
        }
      }
      af::shared<miller::index<> > strong_indices;
    };
  }
}
