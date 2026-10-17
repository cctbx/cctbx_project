#pragma once
#include <cctbx/xray/thickness.h>

#include <smtbx/ED/utils.h>
#include <smtbx/ED/beam_group_processor.h>

namespace smtbx { namespace ED
{
  template <typename FloatType>
  struct beam_group_integrator {
    ED_UTIL_TYPEDEFS;
    typedef beam_group_processor<FloatType> beam_group_processor_t;
    beam_group_integrator(beam_group_processor_t *p,
      const af::shared<FloatType>& angles)
      : processor(p),
      angles(angles),
      beam_n(processor->beam_group.strong_measured_beams.size())
    {
      Is.resize(beam_n);
      if (processor->calc_grad) {
        size_t d_T_off = processor->Ds_kin.size();
        D_dyn.resize(
          af::mat_grid(beam_n, d_T_off + (processor->thickness.grad ? 1 : 0)));
      }
    }

    void integrate_1() {
      const cart_t &K = processor->get_K();
      size_t n_cols = D_dyn.accessor().n_columns();
      const BeamGroup<FloatType>& beam_group = processor->beam_group;
      mat_t D_dyn1;
      for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
        size_t beam_idx = beam_group.strong_measured_beams[i];
        miller::index<> h = beam_group.indices[beam_idx];
        //FloatType da = beam_group.get_diffraction_angle(h, -K[2]);
        //FloatType sg_a = beam_group.Sg_to_angle(0.01, h, -K[2]);
        //int steps = std::abs(round((da - sg_a) / step));
        FloatType I1 = -1, K_g_l1;
        //angles = beam_group.get_angles_Sg(h, -K[2], 0.02, 0.001);
        for (size_t ai = 0;  ai < angles.size(); ai++) {
          mat3_t R = beam_group.get_R(angles[ai]);
          processor->process_1(i+1, R);
          if (processor->exception_) {
            exception_.swap(processor->exception_);
            break;
          }
          // sum up intensities and derivatives
          if (I1 >= 0) {
            cart_t K_g = R * cart_t(h[0], h[1], h[2]) + K;
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
            cart_t K_g = R * cart_t(h[0], h[1], h[2]) + K;
            K_g_l1 = K_g.length();
            I1 = std::norm(processor->CIs[0]);
          }
        }
      }
    }

    void integrate() {
      const cart_t& K = processor->get_K();
      size_t n_cols = D_dyn.accessor().n_columns();
      const BeamGroup<FloatType>& beam_group = processor->beam_group;
      af::shared<FloatType> Is1(beam_n), K_g_ls(beam_n);
      mat_t D_dyn1;
      bool second_step = false;
      for (size_t a = 0; a < angles.size(); a++) {
        mat3_t R = beam_group.get_R(angles[a]);
        processor->process(R);
        if (processor->exception_) {
          exception_.swap(processor->exception_);
          break;
        }
        // sum up intensities and derivatives
        if (second_step) {
          for (size_t ai = 0; ai < beam_n; ai++) {
            miller::index<> h = beam_group.indices[beam_group.strong_measured_beams[ai]];
            cart_t K_g = R * cart_t(h[0], h[1], h[2]) + K;
            FloatType K_g_l = K_g.length();
            FloatType I = std::norm(processor->CIs[ai]);
            FloatType d_st = std::abs(K_g_ls[ai] - K_g_l) / 2;
            Is[ai] += (Is1[ai] + I) * d_st;
            Is1[ai] = I;
            K_g_ls[ai] = K_g_l;

            if (processor->calc_grad) {
              for (size_t j = 0; j < n_cols; j++) {
                D_dyn(ai, j) += (D_dyn(ai, j) + processor->D_dyn(ai, j)) * d_st;
              }
            }
          }
        }
        else {
          D_dyn1 = processor->D_dyn.deep_copy();
          for (size_t ai = 0; ai < beam_n; ai++) {
            miller::index<> h = beam_group.indices[beam_group.strong_measured_beams[ai]];
            cart_t K_g = R * cart_t(h[0], h[1], h[2]) + K;
            K_g_ls[ai] = K_g.length();
            Is1[ai] = std::norm(processor->CIs[ai]);
          }
          second_step = true;
        }
      }
    }

    void operator ()() {
      integrate();
    }
    beam_group_processor_t* processor;
    af::shared<FloatType> angles;
    boost::shared_ptr<beam_group_processor_t> processor_;
    size_t beam_n;
    // output
    af::shared<FloatType> Is;
    mat_t D_dyn;
    boost::scoped_ptr<smtbx::error> exception_;
  };
}}
