#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_N_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_N_H

#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/import_scitbx_af.h>
#include <scitbx/vec3.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/refinement/constraints/reparametrisation.h>
#include <scitbx/array_family/versa_matrix.h>
#include <fast_linalg/lapacke.h>

namespace smtbx {  namespace refinement {
namespace least_squares {
  using namespace smtbx::ED;
  using namespace refinement::constraints;

  template <typename FloatType>
  struct ed_n_shared_data {
    typedef std::complex<FloatType> complex_t;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef boost::shared_ptr< f_calc_function_base_t> f_calc_function_base_ptr_t;
    typedef builder_base<FloatType> data_t;
    typedef af::shared<const BeamInfo<FloatType>*> beam_at;
    typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;
    typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;
    typedef boost::shared_ptr<lookup_t> lookup_ptr_t;
    typedef cctbx::xray::fc_correction<FloatType> fc_correction_t;
    typedef boost::shared_ptr< fc_correction_t> fc_correction_ptr_t;

    ed_n_shared_data(reparametrisation const& reparamn,
      f_calc_function_base_t& f_calc_function,
      fc_correction_t const& fc_cr,
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      scitbx::mat3<FloatType> const& UB,
      af::shared<BeamInfo<FloatType> > const& beams,
      cctbx::xray::thickness<FloatType> const& thickness,
      // wavelength, epsilon, maxSg, F000
      af::shared<FloatType> const& params,
      bool compute_grad,
      bool do_build = true)
      : reparamn(reparamn),
      f_calc_function(f_calc_function),
      space_group(space_group),
      fc_cr(fc_cr),
      Kl(params[0]),
      Fc2Ug(params[1]),
      eps(params[2]),
      UB(UB),
      thickness(thickness),
      compute_grad(compute_grad)
    {
      for (size_t i = 0; i < beams.size(); i++) {
        typename std::map<int, beam_at>::iterator itr = frame_beams.find(beams[i].parent->id);
        beam_at* fb;
        if (itr == frame_beams.end()) {
          fb = &frame_beams.insert(beam_me(beams[i].parent->id, beam_at())).first->second;
        }
        else {
          fb = &itr->second;
        }
        fb->push_back(&beams[i]);
      }
      // build lookups for each frame + collect all indices and they diffs
      af::shared<miller::index<> > all_indices;
      size_t offset = 0;
      // treat equivalents independently inside the frames
      sgtbx::space_group P1("P 1");
      for (typename std::map<int, beam_at>::iterator i = frame_beams.begin();
        i != frame_beams.end(); i++)
      {
        frame_offsets.insert(std::make_pair(i->first, offset));
        offset += i->second.size();
        af::shared<miller::index<> > fis(i->second.size());
        for (size_t hi = 0; hi < i->second.size(); hi++) {
          fis[hi] = i->second[hi]->index;
          all_indices.push_back(fis[hi]);
          for (size_t hj = hi+1; hj < i->second.size(); hj++) {
            all_indices.push_back(fis[hi] - i->second[hj]->index);
          }
        }
        frame_indices.insert(std::make_pair(i->first, fis));
        lookup_ptr_t mi_l(new lookup_t(
          fis.const_ref(),
          P1,
          false));
        frame_lookups.insert(std::make_pair(i->first, mi_l));
      }
      beam_n = offset;
      // a tricky way of getting unique only...
      mi_lookup = lookup_t(
        all_indices.const_ref(),
        space_group,
        anomalous_flag);
      indices = mi_lookup.get_unique();
      mi_lookup = lookup_t(
        indices.const_ref(),
        space_group,
        anomalous_flag);
      if (do_build) {
        build();
      }
    }
    complex_t calc_one_h(miller::index<> const& h) const {
      f_calc_function.compute(h, boost::none, 0, false);
      complex_t fc = f_calc_function.get_f_calc();
      FloatType fc_k = fc_cr.compute(h, f_calc_function.get_observable(), false);
      if (fc_k != 1) {
        fc *= std::sqrt(fc_k);
      }
      return fc;
    }

    /* could be used for both - values and derivatives calculation for the given
    set of beams of a frame. The result is written to Fcs starting at the given
    offset
    */
    struct process_frame {
      process_frame(ed_n_shared_data const& parent,
        const beam_at& beams,
        // source kinematic Fcs, Fc, Fc_+e, Fc_-e
        const af::shared<complex_t> Fcs_k,
        af::shared<complex_t>& Fcs,
        size_t Fc_offset)
        : parent(parent),
        beams(beams),
        thickness(parent.thickness.value),
        Fcs_k(Fcs_k),
        Fcs(Fcs),
        Fc_offset(Fc_offset)
      {}

      void operator()() const {
        const FrameInfo<FloatType>& frame = *beams[0]->parent;
        const FloatType Kl = parent.Kl, Fc2Ug = parent.Fc2Ug;
        const cart_t K = cart_t(0, 0, -Kl);
        // Projection of K onto normal of frame normal, K*frame.normal
        const FloatType Kn = -frame.normal[2] * Kl;
        // A will contain first row and column the initial kinematical structure factor,
        // diagonally the excitation error and off diagonally the relative miller
        // index based structure factor, as defined in 12.11 a&b in order to be evaluated
        // in the eigenvalue problem, defined in 12.10
        using namespace fast_linalg;
        const size_t n_beams = beams.size() + 1; // g0+
        af::versa<complex_t, af::mat_grid> A(af::mat_grid(n_beams, n_beams));
        af::shared<FloatType> M(n_beams);
        M[0] = 1; //for g0
        for (size_t i = 1; i < n_beams; i++) {
          miller::index<> h_i = beams[i - 1]->index;
          int ii = parent.mi_lookup.find_hkl(h_i);
          complex_t Fc_i = 0;
          if (ii == -1) {
            if (!parent.space_group.is_sys_absent(h_i)) {
              SMTBX_ASSERT(ii >= 0)(h_i.as_string());
            }
          }
          else {
            Fc_i = Fcs_k[ii];
          }
          cart_t g_i = frame.RMf * cart_t(h_i[0], h_i[1], h_i[2]);
          FloatType s = (Kl * Kl - (K + g_i).length_sq());
          FloatType i_den = std::sqrt(1. / (1 + g_i * frame.normal / Kn));
          A(i, i) = s * i_den * i_den;

          A(i, 0) = Fc2Ug * Fc_i * i_den;
          A(0, i) = std::conj(A(i, 0));

          M[i] = i_den;
          for (size_t j = i + 1; j < n_beams; j++) {
            miller::index<> h_j = beams[j - 1]->index;
            cart_t g_j = frame.RMf * cart_t(h_j[0], h_j[1], h_j[2]);
            miller::index<> h_i_m_j = h_i - h_j;
            int i_m_j = parent.mi_lookup.find_hkl(h_i_m_j);
            complex_t Fc_i_m_j = 0;
            if (i_m_j == -1) {
              if (!parent.space_group.is_sys_absent(h_i_m_j)) {
                SMTBX_ASSERT(i_m_j >= 0)(h_i_m_j.as_string());
              }
            }
            else {
              Fc_i_m_j = Fcs_k[i_m_j];
            }
            FloatType j_den = std::sqrt(1. / (1 + g_j * frame.normal / Kn));
            A(j, i) = Fc2Ug * Fc_i_m_j * i_den * j_den;
            A(i, j) = std::conj(A(j, i));
          }
        }
/*
        // eigenvalues, this is for generic matrix
        af::shared<complex_t> ev(n_beams);
        // right eigenvectors
        af::versa<complex_t, af::c_grid<2> > eV(af::c_grid<2>(n_beams, n_beams));
        lapack_int info = geev(LAPACK_ROW_MAJOR, 'N', 'V', n_beams,
          A.begin(), n_beams, ev.begin(), 0, n_beams, eV.begin(), n_beams);
        SMTBX_ASSERT(!info)(info);
        af::versa<complex_t, af::mat_grid>
          B(af::mat_grid(n_beams, n_beams), *eV.begin());
        // invert eV now
        {
          af::shared<lapack_int> pivots(n_beams);
          info = getrf(LAPACK_ROW_MAJOR, n_beams, n_beams, eV.begin(),
            n_beams, pivots.begin());
          SMTBX_ASSERT(!info)(info);
          info = getri(LAPACK_ROW_MAJOR, n_beams, eV.begin(),
            n_beams, pivots.begin());
          SMTBX_ASSERT(!info)(info);
        }
*/
        // eigenvalues, for Hermitian - the matrix built above
        af::shared<FloatType> ev(n_beams);
        lapack_int info = heev(LAPACK_ROW_MAJOR, 'V', LAPACK_UPPER, n_beams,
          A.begin(), n_beams, ev.begin());
        SMTBX_ASSERT(!info)(info);
        // heev replaces A with column-wise eigenvectors
        af::versa<complex_t, af::mat_grid> &B = A;
        af::versa<complex_t, af::mat_grid> eV = af::matrix_transpose(A.const_ref());

        const complex_t exp_k(0, scitbx::constants::pi * thickness / Kn);
        //const complex_t exp_k(0, scitbx::constants::pi * thickness);
        // B = B*diag(exp(2*pi*thickness*ev/(2*Kn))
        //af::versa<complex_t, af::mat_grid> P(af::mat_grid(n_beams, n_beams));
        for (size_t i = 0; i < n_beams; i++) {
          complex_t m = std::exp(ev[i] * exp_k);
          for (size_t j = 0; j < n_beams; j++) {
            B(j, i) *= m;
          }
        }
        // P = B*eV^-1
        af::versa<complex_t, af::mat_grid> P = af::matrix_multiply(B.const_ref(), eV.const_ref());
        //P = M*P*M^-1, m - diagonal, M[0] = 1, so start with 1
        for (size_t i = 1; i < n_beams; i++) {
          for (size_t j = 1; j < n_beams; j++) {
            P(i, j) *= M[i]; // M is on the left - apply to rows
            P(j, i) /= M[i]; // M^-1 is on the right - > apply to cols
          }
        }
        af::shared<complex_t> u(n_beams);
        u[0] = complex_t(1, 0);
        /* compute up now and assign to Fcs using at Fc_offset+exited_indices[i]
        * (could save a lot useless calculations if this is used first from the right)
        */
        af::shared<complex_t> up = af::matrix_multiply(P.const_ref(), u.const_ref());
        for (size_t i = 1; i < n_beams; i++) {
          Fcs[Fc_offset + i - 1] = up[i]; // / SfacToVolts; // P(i, 0);
          // is this really needed? Eq 4, 2013
          /*
          miller::index<> h_i = beams[i - 1]->index;
          cart_t g_i = frame.RMf * cart_t(h_i[0], h_i[1], h_i[2]);
          cart_t g_i_p = g_i + frame.normal * (g_i * frame.normal);
          for (size_t j = 1; j < n_beams; j++) {
            if (j == i) {
              continue;
            }
            miller::index<> h_j = beams[j - 1]->index;
            cart_t g_j = frame.RMf * cart_t(h_j[0], h_j[1], h_j[2]);
            cart_t g_j_p = g_j + frame.normal * (g_j * frame.normal);
            if ((g_i_p - g_j_p).length_sq() < 1e-3) {
              Fcs[Fc_offset + i - 1] += up[j];
            }
          }
          */
        }
      }
      ed_n_shared_data const& parent;
      const beam_at& beams;
      FloatType thickness;
      const af::shared<complex_t>& Fcs_k;
      size_t Fc_offset;
      af::shared<complex_t>& Fcs;
    };

    void process_frames_mt(af::shared<complex_t>& Fcs_,
      af::shared<complex_t> const& Fcs_k,
      int thread_count = -1)
    {
      if (thread_count < 0) {
        thread_count = builder_base<FloatType>::get_available_threads();
      }
      boost::thread_group pool;
      typedef boost::shared_ptr<process_frame> frame_processor_t;
      typename std::map<int, beam_at>::iterator f_itr = frame_beams.begin();
      size_t Fc_offset = 0;
      for (size_t fi = 0; fi < frame_beams.size(); fi += thread_count) {
        size_t t_end = std::min(thread_count, (int)(frame_beams.size() - fi));
        if (t_end == 0) {
          break;
        }
        std::vector<frame_processor_t> accumulators;
        for (int thread_idx = 0; thread_idx < t_end; thread_idx++) {
          frame_processor_t pf(
            new process_frame(*this,
              f_itr->second,
              Fcs_k,
              Fcs_,
              Fc_offset)
          );
          accumulators.push_back(pf);
          pool.create_thread(boost::ref(*pf));
          Fc_offset += f_itr->second.size();
          f_itr++;
        }
        pool.join_all();
      }
    }

    void build() {
      // Generate Fcs at current position
      {
        Fcs.resize(beam_n);
        af::shared<complex_t> Fc_cur(indices.size());
        for (size_t ih = 0; ih < indices.size(); ih++) {
          Fc_cur[ih] = calc_one_h(indices[ih]);
        }
        process_frames_mt(Fcs, Fc_cur);
        if (!compute_grad) {
          return;
        }
      }
      af::shared<parameter*> params = reparamn.independent();
      size_t param_n = 0;
      for (size_t i = 0; i < params.size(); i++) {
        param_n += params[i]->components().size();
      }
      design_matrix.resize(af::c_grid<2>(beam_n, param_n));
      // Generate Arrays for storing positive and negative epsilon Fcs for numerical
      // generation of gradients
      af::shared<complex_t> Fcs_p(beam_n), Fcs_m(beam_n), Fc_eps(indices.size());
      FloatType t_eps = 2 * eps;
      for (size_t i = 0, n = 0; i < params.size(); i++) {
        af::ref<double> x = params[i]->components();
        asu_parameter* cp = dynamic_cast<asu_parameter*>(params[i]);
        for (size_t j = 0; j < x.size(); j++, n++) {
          x[j] += eps;
          if (cp != 0) {
            cp->store(reparamn.unit_cell());
          }
          for (size_t i_h = 0; i_h < indices.size(); i_h++) {
            Fc_eps[i_h] = calc_one_h(indices[i_h]);
          }
          //Generate Fcs at x+eps
          process_frames_mt(Fcs_p, Fc_eps);

          x[j] -= t_eps;
          if (cp != 0) {
            cp->store(reparamn.unit_cell());
          }
          for (size_t i_h = 0; i_h < indices.size(); i_h++) {
            Fc_eps[i_h] = calc_one_h(indices[i_h]);
          }
          //Generate Fcs at x-eps
          process_frames_mt(Fcs_m, Fc_eps);
          
          // compute grads
          for (size_t bi = 0; bi < beam_n; bi++) {
            complex_t grad_fc = (Fcs_p[bi] - Fcs_m[bi]) / t_eps;
            design_matrix(bi, n) = 2 * (
              Fcs[bi].real() * grad_fc.real() +
              Fcs[bi].imag() * grad_fc.imag());
          }

          // reset to original value
          x[j] += eps;
          if (cp != 0) {
            cp->store(reparamn.unit_cell());
          }
        }
      }
    }
  

    /* computes the position of the given miller index of the given
    frame in the uniform arrays
    */
    size_t find_hkl(int frame_id, miller::index<> const& h) const {
      typename std::map<int, lookup_ptr_t>::const_iterator i = frame_lookups.find(frame_id);
      SMTBX_ASSERT(i != frame_lookups.end());
      std::map<int, size_t>::const_iterator off = frame_offsets.find(frame_id);
      long hi = i->second->find_hkl(h);
      SMTBX_ASSERT(hi >= 0);
      return hi + off->second;
    }

    reparametrisation const& reparamn;
    f_calc_function_base_t& f_calc_function;
    cctbx::xray::fc_correction<FloatType> const& fc_cr;
    sgtbx::space_group const& space_group;
    af::shared<miller::index<> > indices;
    FloatType Kl, Fc2Ug, eps;
    scitbx::mat3<FloatType> UB;
    size_t beam_n;
    // a map of beams by frame id
    std::map<int, beam_at> frame_beams;
    /* to lookup an index in particular frame, have to keep a copy of the
    indices (not have to be a map!)
    */
    std::map<int, af::shared<miller::index<> > > frame_indices;
    typename std::map<int, lookup_ptr_t> frame_lookups;
    std::map<int, size_t> frame_offsets;

    cctbx::xray::thickness<FloatType> const& thickness;
    bool compute_grad;
    // newly-calculated, aligned by frames
    af::shared<complex_t> Fcs;
    // 
    af::versa<FloatType, af::c_grid<2> > design_matrix;
    lookup_t mi_lookup;
  };

  template <typename FloatType>
  class f_calc_function_ed_n : public f_calc_function_base<FloatType> {
  public:
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef ed_n_shared_data<FloatType> data_t;
    typedef scitbx::vec3<FloatType> cart_t;

    f_calc_function_ed_n(ed_n_shared_data<FloatType> const& data)
      : data(data),
      index(~0)
    {}
    f_calc_function_ed_n(f_calc_function_ed_n const& other)
      : data(other.data),
      index(~0)
    {}

    virtual void compute(
      miller::index<> const& h,
      boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
      twin_fraction<FloatType> const* fraction = 0,
      bool compute_grad = true)
    {
      SMTBX_ASSERT(fraction->tag >= 0);
      index = data.find_hkl(fraction->tag, h);
      Fsq = std::norm(data.Fcs[index]);
    }

    virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
      return boost::shared_ptr<f_calc_function_base_t>(
        new f_calc_function_ed_n(*this));
    }

    virtual FloatType get_observable() const {
      return Fsq;
    }
    virtual std::complex<FloatType> get_f_calc() const {
      return data.Fcs[index];
    }
    virtual af::const_ref<FloatType> get_grad_observable() const {
      typedef af::versa_plain<FloatType> one_dim_type;
      typedef typename one_dim_type::accessor_type one_dim_accessor_type;
      one_dim_accessor_type a(data.design_matrix.accessor().n_columns());
      return af::const_ref<FloatType>(&data.design_matrix(index, 0), a);
    }

    virtual bool raw_gradients() const { return false; }
  private:
    ed_n_shared_data<FloatType> const& data;
    size_t index;
    mutable FloatType Fsq;
  };

}}}

#endif // GUARD
