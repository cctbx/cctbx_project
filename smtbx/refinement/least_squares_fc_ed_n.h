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
    typedef builder_base<FloatType> data_t;
    typedef af::shared<const BeamInfo<FloatType>*> beam_at;
    typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;
    typedef miller::lookup_utils::lookup_tensor<FloatType> lookup_t;

    ed_n_shared_data(data_t const& k_data,
      reparametrisation const &reparamn,
      bool anomalous_flag,
      f_calc_function_base_t & f_calc_function,
      cctbx::xray::fc_correction<FloatType> const& fc_cr,
      sgtbx::space_group const& space_group,
      scitbx::mat3<FloatType> const& UB,
      af::shared<BeamInfo<FloatType> > const& beams,
      cctbx::xray::thickness<FloatType> const& thickness,
      // wavelength, epsilon, maxSg, F000
      af::shared<FloatType> const& params,
      bool compute_grad,
      bool do_build=true)
      : k_data(k_data),
      reparamn(reparamn),
      eps(params[1]),
      f_calc_function(f_calc_function),
      space_group(space_group),
      fc_cr(fc_cr),
      indices(k_data.reflections().indices()),
      wavelength(params[0]),
      F000(params[3]),
      UB(UB),
      thickness(thickness),
      maxSg(params[2]),
      compute_grad(compute_grad)
    {
      mi_lookup = lookup_t(
        indices.const_ref(),
        space_group,
        anomalous_flag);
      for (size_t i = 0; i < beams.size(); i++) {
        std::map<int, beam_at>::iterator itr = frame_beams.find(beams[i].parent->id);
        beam_at *fb;
        if (itr == frame_beams.end()) {
          fb = &frame_beams.insert(beam_me(beams[i].parent->id, beam_at())).first->second;
        }
        else {
          fb = &itr->second;
        }
        fb->push_back(&beams[i]);
      }
      size_t offset = 0;
      for (std::map<int, beam_at>::iterator i = frame_beams.begin();
        i != frame_beams.end(); i++)
      {
        frame_offsets.insert(std::make_pair(i->first, offset));
        offset += i->second.size();
        af::shared<miller::index<> > fis(i->second.size());
        for (size_t hi = 0; hi < i->second.size(); hi++) {
          fis[hi] = i->second[hi]->index;
        }
        frame_indices.insert(std::make_pair(i->first, fis));
        lookup_t mi_l = lookup_t(
          fis.const_ref(),
          space_group,
          anomalous_flag);
        frame_lookups.insert(std::make_pair(i->first, mi_l));
      }
      beam_n = offset;
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
      process_frame(ed_n_shared_data const & parent,
        const beam_at& beams,
        // source kinematic Fcs, Fc, Fc_+e, Fc_-e
        const af::shared<complex_t> Fcs_k,
        af::shared<complex_t>& Fcs,
        size_t Fc_offset)
        : parent(parent),
        beams(beams),
        thickness(thickness),
        Fcs_k(Fcs_k),
        Fcs(Fcs),
        Fc_offset(Fc_offset)
      {}
      void operator()() const {
        // Builds (up)rime based on input vector u according to equation 12.16
        // based on P matrix as defined in 12.17
        const FrameInfo<FloatType>& frame = *beams[0]->parent;
        const FloatType Kl = 1.0 / parent.wavelength,
          Kl_sq = Kl * Kl;
        const cart_t K = cart_t(0, 0, -Kl);
        // Projection of K onto normal of frame
        const FloatType Kn = -frame.normal[2] * Kl;
        // This is the constant factor for eigenvalues used in eq. 12.17
        const FloatType exp_k = scitbx::constants::two_pi * thickness / Kn;
        // find "excited" beams
        af::shared<const BeamInfo<FloatType>*> excited;
        std::vector<size_t> excited_indices;
        // will carry Fcs.real() of the excited beams
        af::shared<FloatType> u;
        // u(0) = U000
        u.push_back(parent.F000);
        for (size_t i = 0; i < beams.size(); i++) {
          miller::index<> h = beams[i]->index;
          // g vector of miller index according to definition on page 276
          cart_t g = frame.RM * parent.UB * cart_t(h[0], h[1], h[2]);
          // Sg is the distance from Bragg condition
          FloatType Sg = (Kl * Kl - scitbx::fn::pow2((K + g).length())) / (2 * Kl);
          bool is_excited = std::abs(Sg) < parent.maxSg;
          int ii = parent.mi_lookup.find_hkl(h);
          if (ii == -1) {
            if (!parent.space_group.is_sys_absent(h)) {
              SMTBX_ASSERT(ii >= 0)(h.as_string());
            }
          }
          if (is_excited) {
            excited.push_back(beams[i]);
            excited_indices.push_back(i);
            u.push_back(ii < 0 ? 0 : Fcs_k[ii].real());
          }
          else {
            // FK: I am not entirely sure that they not maybe simply ignored, that is
            // removed from the refinement completely, if not considered "excited".
            // See also Acta Cryst A 71, 2015, p 240 left bottom
            Fcs[Fc_offset + i] = ii < 0 ? 0 : Fcs_k[ii];
          }
        }

        // A will contain first row and column the initial kienmatical structure factor,
        // diagonally the excitation error and off diagonally the relative miller
        // index based structure factor, as defined in 12.11 a&b in order to be evaluated
        // in the eigenvalue problem, defined in 12.10
        af::versa<FloatType, af::mat_grid> A(
          af::mat_grid(excited.size()+1, excited.size()+1));

        A(0, 0) = 0;
        for (size_t i = 1; i <= excited.size(); i++) {
          miller::index<> h_i = excited[i - 1]->index;
          int ii = parent.mi_lookup.find_hkl(h_i);
          if (ii == -1) {
            if (!parent.space_group.is_sys_absent(h_i)) {
              SMTBX_ASSERT(ii >= 0)(h_i.as_string());
            }
          }
          FloatType Fc_i = ii < 0 ? 0 : Fcs_k[ii].real();
          cart_t g_i = frame.RM * parent.UB * cart_t(h_i[0], h_i[1], h_i[2]);
          FloatType s = (Kl * Kl - scitbx::fn::pow2((K + g_i).length()))/(2*Kl);
          FloatType den = 1./(2 * Kl * g_i * frame.normal);
          A(i, i) = s;
          A(0, i) = Fc_i*den;
          A(i, 0) = A(0, i);
          for (size_t j = i + 1; j <= excited.size(); j++) {
            miller::index<> h_i_m_j = h_i - excited[j - 1]->index;
            cart_t g_i_m_j = frame.RM * parent.UB * cart_t(h_i_m_j[0], h_i_m_j[1], h_i_m_j[2]);
            /*FK: I think this check is what we need to ensure in the F_calcs... 
                  we need to provide all F_calcs in the range of -2*(hkl)_max to 2*(hkl)_max 
                  to be able to fill the matrix. This will be a requirement to k_data 
                  that is given to ed_n_shared_data*/
            int i_m_j = parent.mi_lookup.find_hkl(h_i_m_j);
            if (i_m_j == -1) {
              if (!parent.space_group.is_sys_absent(h_i_m_j)) {
                SMTBX_ASSERT(i_m_j >= 0)(h_i_m_j.as_string());
              }
            }
            FloatType Fc_i_m_j = i_m_j < 0 ? 0 : Fcs_k[i_m_j].real();
            A(i, j) = Fc_i_m_j * den;
            A(j, i) = A(i, j);
          }
        }
        //THIS CRIES FOR PARALLELISATION!
        // This is the eigenvalue system defined in 12.10
        scitbx::matrix::eigensystem::real_symmetric<FloatType> es(A.const_ref());
        // These gamma are the eigenvalues from eq. 12.10
        // FK: At some point I want to look at the es.values() to see whether they are generate
        af::shared<FloatType> gamma = es.values();
        // This will carry the exponential factors of diagonal matrix in equation 12.17
        af::shared<FloatType> exp_L(excited.size());
        for (size_t i = 0; i < excited.size(); i++) {
          exp_L[i] = std::exp(gamma[i] * exp_k);
        }
        // C is a matrix consisting of the eigenvectors of A.
        af::versa<FloatType, af::mat_grid> C = es.vectors();
        // I do not think we need to transpose, if we later use matrix_tanspose_multiply
        // OR is the vectors oriented wrongly?
        //C.ref().transpose_square_in_place();
        
        // Generate P(n+1,n+1) = C * I*gamma * C^(-1),
        // where gamma is a vector with the eigenvalues, I is a diagonal identy matrix
        // and C^(-1), due to the symmetry and "realness" of C, is the transpose of C
        // and n is the number of measured, excited beams.
        af::versa<FloatType, af::c_grid<2> > P(
          af::c_grid<2>(excited.size()+1, excited.size()+1),
          *af::matrix_transpose_multiply_diagonal_multiply_as_packed_u(
            C.const_ref(), exp_L.const_ref()).begin()); 

        //for (size_t i = 0; i < excited.size(); i++) {
        //  for (size_t j = 0; j < excited.size(); j++) {
        //    exp_M1(i, j) *= M[i]; // M is on the left - apply to rows
        //    exp_M1(j, i) /= M[i]; // M^-1 is on the right - > apply to cols
        //  }
        //}
        /* compute up now and assign to Fcs using at Fc_offset+exited_indices[i] */
        af::shared<FloatType> up = af::matrix_multiply(P.const_ref(), u.const_ref());
        for (size_t i = 1; i < up.size(); i++) {
          Fcs[Fc_offset + excited_indices[i - 1]] = up[i];
        }
      }
      ed_n_shared_data const& parent;
      const beam_at& beams;
      FloatType thickness;
      const af::shared<complex_t>& Fcs_k;
      size_t Fc_offset;
      af::shared<complex_t>& Fcs;
    };

    void process_frames_mt(af::shared<complex_t> &Fcs_,
      af::shared<complex_t> const& Fcs_k,
      int thread_count=-1)
    {
      if (thread_count < 0) {
        thread_count = builder_base<FloatType>::get_available_threads();
      }
      boost::thread_group pool;
      typedef boost::shared_ptr<process_frame> frame_processor_t;
      std::map<int, beam_at>::iterator f_itr = frame_beams.begin();
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
      af::shared<parameter*> params = reparamn.independent();
      size_t param_n = 0;
      // Generate Arrays for storing positive and negative epsilon Fcs for numerical
      // generation of gradients
      if (compute_grad) {
        for (size_t i = 0; i < params.size(); i++) {
          param_n += params[i]->components().size();
        }
        Fc_plus_eps.resize(param_n);
        Fc_minus_eps.resize(param_n);
        FloatType t_eps = 2 * eps;
        for (size_t i = 0, n = 0; i < params.size(); i++) {
          af::ref<double> x = params[i]->components();
          asu_parameter* cp = dynamic_cast<asu_parameter*>(params[i]);
          for (size_t j = 0; j < x.size(); j++, n++) {
            Fc_plus_eps[n].resize(indices.size());
            x[j] += eps;
            if (cp != 0) {
              cp->store(reparamn.unit_cell());
            }
            for (size_t i_h = 0; i_h < indices.size(); i_h++) {
              Fc_plus_eps[n][i_h] = calc_one_h(indices[i_h]);
            }
            Fc_minus_eps[n].resize(indices.size());
            x[j] -= t_eps;
            if (cp != 0) {
              cp->store(reparamn.unit_cell());
            }
            for (size_t i_h = 0; i_h < indices.size(); i_h++) {
              Fc_minus_eps[n][i_h] = calc_one_h(indices[i_h]);
            }
            // reset to original value
            x[j] += eps;
            if (cp != 0) {
              cp->store(reparamn.unit_cell());
            }
          }
        }
      }
      // Generate Fcs at current position
      Fcs.resize(beam_n);
      process_frames_mt(Fcs, k_data.f_calc());
      if (compute_grad) {
        af::shared<complex_t> Fcs_n(beam_n), Fcs_p(beam_n);
        design_matrix.resize(af::c_grid<2>(param_n, beam_n));
        for (size_t i = 0; i < param_n; i++) {
          //Generate Fcs at x-eps
          process_frames_mt(Fcs_n, Fc_minus_eps[i]);
          //Generate Fcs at x+eps
          process_frames_mt(Fcs_p, Fc_plus_eps[i]);
          for (size_t j = 0; j < beam_n; j++) {
            design_matrix(j, i) = (Fcs_p[j].real() - Fcs_n[j].real()) / eps;
          }
        }
      }
    }

    /* computes the position of the given miller index of the given
    frame in the uniform arrays
    */
    size_t find_hkl(int frame_id, miller::index<> const& h) const {
      std::map<int, lookup_t>::const_iterator i = frame_lookups.find(frame_id);
      SMTBX_ASSERT(i == frame_lookups.end());
      std::map<int, size_t>::const_iterator off = frame_offsets.find(frame_id);
      long hi = i->second.find_hkl(h);
      SMTBX_ASSERT(hi < 0);
      return hi + off->second;
    }

    data_t const& k_data;
    reparametrisation const& reparamn;
    FloatType eps;
    f_calc_function_base_t& f_calc_function;
    cctbx::xray::fc_correction<FloatType> const& fc_cr;
    sgtbx::space_group const& space_group;
    af::shared<miller::index<> > indices;
    FloatType wavelength, F000;
    scitbx::mat3<FloatType> UB;
    size_t beam_n;
    // a map of beams by frame id
    std::map<int, beam_at> frame_beams;
    /* to lookup an index in particular frame, have to keep a copy of the
    indices (not have to be a map!)
    */
    std::map<int, af::shared<miller::index<> > > frame_indices;
    std::map<int, lookup_t> frame_lookups;
    std::map<int, size_t> frame_offsets;

    cctbx::xray::thickness<FloatType> const& thickness;
    FloatType maxSg;
    bool compute_grad;
    // Fc for each variation in independent param's component
    af::shared<af::shared<complex_t> > Fc_plus_eps,
      Fc_minus_eps;
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
      : data(data)
    {}
    f_calc_function_ed_n(f_calc_function_ed_n const& other)
      : data(other.data)
    {}

    virtual void compute(
      miller::index<> const& h,
      boost::optional<std::complex<FloatType> > const& f_mask = boost::none,
      twin_fraction<FloatType> const* fraction = 0,
      bool compute_grad = true)
    {
      frame_index = fraction->tag;
      index = data.find_hkl(frame_index, h);
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
    int frame_index;
    size_t index;
    mutable FloatType Fsq;
  };

}}}

#endif // GUARD
