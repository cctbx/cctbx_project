#pragma once
#include <cctbx/xray/thickness.h>

#include <smtbx/ED/utils.h>
#include <smtbx/ED/frame_processor.h>

namespace smtbx {
  namespace ED {
    template <typename FloatType>
    struct frame_integrator {
      ED_UTIL_TYPEDEFS;
      typedef a_frame_processor<FloatType> frame_processor_t;
      frame_integrator(frame_processor_t* p,
        FloatType angle, FloatType step, int mode=0)
        : processor_(p),
        processor(p),
        beam_n(processor->frame.strong_measured_beams.size()),
        angle(angle),
        step(step),
        mode(mode)
      {
        Is.resize(beam_n);
        if (processor->calc_grad) {
          size_t d_T_off = processor->Ds_kin.size();
          D_dyn.resize(
            af::mat_grid(beam_n, d_T_off + (processor->thickness.grad ? 1 : 0)));
        }
      }

      void compute1(bool align) {
        size_t n_cols = D_dyn.accessor().n_columns();
        if (align) {
          const FrameInfo<FloatType>& frame = processor->frame;
          for (size_t i = 0; i < beam_n; i++) {
            size_t beam_idx = frame.strong_measured_beams[i];
            std::pair<mat3_t, cart_t> r = frame.compute_RMf_N(
              frame.beams[beam_idx].diffraction_angle);
            processor->process(r.first, r.second);
            Is[i] = std::norm(processor->CIs[i]);
            if (processor->calc_grad) {
              std::copy(&processor->D_dyn(i, 0), &processor->D_dyn(i, n_cols), &D_dyn(i, 0));
            }
          }
        }
        else {
          (*processor)();
          if (processor->exception_) {
            exception_.swap(processor->exception_);
            return;
          }
          for (size_t i = 0; i < beam_n; i++) {
            Is[i] = std::norm(processor->CIs[i]);
            if (processor->calc_grad) {
              for (size_t j = 0; j < n_cols; j++) {
                D_dyn(i, j) = processor->D_dyn(i, j);
              }
            }
          }
        }
      }

      void integrate() {
        const cart_t& K = processor->K;
        size_t n_cols = D_dyn.accessor().n_columns();
        int steps = round(angle / step);
        const FrameInfo<FloatType>& frame = processor->frame;
        mat_t D_dyn1;
        for (size_t i = 0; i < frame.strong_measured_beams.size(); i++) {
          size_t beam_idx = frame.strong_measured_beams[i];
          miller::index<> h = frame.indices[beam_idx];
          FloatType da = frame.get_diffraction_angle(h, -K[2]);
          FloatType sg_a = frame.Sg_to_angle(0.01, h, -K[2]);
          //int steps = std::abs(round((da - sg_a) / step));
          FloatType I1 = -1, K_g_l1;
          for (int st = -steps; st <= steps; st++) {
            std::pair<mat3_t, cart_t> r = frame.compute_RMf_N(
              da + st * step);
            processor->process_1(i, r.first, r.second);
            if (processor->exception_) {
              exception_.swap(processor->exception_);
              break;
            }
            // sum up intensities and derivatives
            if (I1 >= 0) {
              cart_t K_g = r.first * cart_t(h[0], h[1], h[2]) + K;
              FloatType K_g_l = K_g.length();
              FloatType I = std::norm(processor->CIs[0]);
              FloatType d_st = std::abs(K_g_l1 - K_g_l) / 2;
              Is[i] += (I1 + I) * d_st;
              I1 = I;
              K_g_l1 = K_g_l;
              if (processor->calc_grad) {
                for (size_t j = 0; j < n_cols; j++) {
                  D_dyn(i, j) += (D_dyn1(0, j) + processor->D_dyn(0, j)) * d_st;
                }
              }
              D_dyn1 = processor->D_dyn.deep_copy();
            }
            else {
              D_dyn1 = processor->D_dyn.deep_copy();
              cart_t K_g = r.first * cart_t(h[0], h[1], h[2]) + K;
              K_g_l1 = K_g.length();
              I1 = std::norm(processor->CIs[0]);
            }
          }
        }
      }

      void operator ()() {
        if (mode == 0) {
          integrate();
        }
        else {
          compute1(mode > 1);
        }
      }

      frame_processor_t* processor;
      boost::shared_ptr<frame_processor_t> processor_;
      size_t beam_n;
      FloatType angle, step;
      int mode;
      // output
      af::shared<FloatType> Is;
      mat_t D_dyn;
      boost::scoped_ptr<smtbx::error> exception_;
    };
  }
}
