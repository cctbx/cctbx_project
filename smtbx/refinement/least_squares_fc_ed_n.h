#ifndef SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_N_H
#define SMTBX_REFINEMENT_LEAST_SQUARES_FC_ED_N_H

#include <cctbx/xray/thickness.h>

#include <smtbx/refinement/least_squares_fc.h>
#include <smtbx/refinement/constraints/reparametrisation.h>
#include <smtbx/ED/ed_data.h>

namespace smtbx {  namespace refinement {
namespace least_squares {
  using namespace smtbx::ED;
  using namespace refinement::constraints;

  template <typename FloatType>
  struct ed_n_shared_data {
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef boost::shared_ptr< f_calc_function_base_t> f_calc_function_base_ptr_t;
    typedef builder_base<FloatType> data_t;
    typedef af::shared<const BeamInfo<FloatType>*> beam_at;
    typedef std::pair<int, af::shared<const BeamInfo<FloatType>*> > beam_me;
    typedef boost::shared_ptr<lookup_t> lookup_ptr_t;
    typedef cctbx::xray::fc_correction<FloatType> fc_correction_t;
    typedef boost::shared_ptr< fc_correction_t> fc_correction_ptr_t;

    ed_n_shared_data(const reparametrisation& reparamn,
      f_calc_function_base_t& f_calc_function,
      const fc_correction_t& fc_cr,
      const sgtbx::space_group& space_group,
      bool anomalous_flag,
      af::shared<BeamGroup<FloatType> > beam_groups,
      const cctbx::xray::thickness<FloatType>& thickness,
      const RefinementParams<FloatType> &params,
      bool compute_grad,
      bool do_build = true)
      : reparamn(reparamn),
      f_calc_function(f_calc_function),
      space_group(space_group),
      fc_cr(fc_cr),
      params(params),
      Kvac(params.getKl_vac()),
      Kl(params.getKl()),
      Fc2Ug(params.getFc2Ug()),
      eps(params.getEpsilon()),
      mat_type(params.getMatrixType()),
      beam_groups(beam_groups),
      thickness(thickness),
      compute_grad(compute_grad),
      thread_n(params.getThreadN())
    {
      if (mat_type >= 100) {
        mat_type -= 100;
        use_n_beam = true;
      }
      else {
        use_n_beam = false;
      }
      // build lookups for each beam_group + collect all indices and they diffs
      af::shared<miller::index<> > all_indices;
      size_t offset = 0;
      // treat equivalents independently inside the beam_groups
      sgtbx::space_group P1("P 1");
      for (size_t i=0; i < beam_groups.size(); i++) {
        BeamGroup<FloatType>& beam_group = beam_groups[i];
        beam_groups_map.insert(std::make_pair(beam_group.id, &beam_group));
        beam_group.offset = offset;
        offset += beam_group.strong_measured_beams.size();
        for (size_t hi = 0; hi < beam_group.indices.size(); hi++) {
          all_indices.push_back(beam_group.indices[hi]);
          for (size_t hj = hi+1; hj < beam_group.indices.size(); hj++) {
            all_indices.push_back(beam_group.indices[hi] - beam_group.indices[hj]);
            all_indices.push_back(beam_group.indices[hj] - beam_group.indices[hi]);
          }
        }
        lookup_ptr_t mi_l(new lookup_t(
          af::select(beam_group.indices.const_ref(),
            beam_group.strong_measured_beams.const_ref()).const_ref(),
          P1,
          true));
        beam_group_lookups.insert(std::make_pair(beam_groups[i].id, mi_l));
      }
      beam_n = offset;
      // a tricky way of getting unique only...
      mi_lookup = lookup_t(
        all_indices.const_ref(),
        P1,
        anomalous_flag);
      indices = mi_lookup.get_unique();
      mi_lookup = lookup_t(
        indices.const_ref(),
        P1,
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
    set of beams of a beam_group. The result is written to Fcs starting at the given
    offset
    */
    struct process_beam_group {
      process_beam_group(ed_n_shared_data const& parent,
        BeamGroup<FloatType> &beam_group,
        // source kinematic Fcs, Fc, Fc_+e, Fc_-e
        const af::shared<complex_t> &Fcs_k,
        af::shared<FloatType>& Is,
        af::shared<complex_t>& CIs, bool use_offset)
        : parent(parent),
        beam_group(beam_group),
        thickness(parent.thickness.value),
        Fcs_k(Fcs_k),
        Is(Is),
        CIs(CIs),
        offset(use_offset ? beam_group.offset : 0)
      {}

      void operator()() {
        try {
          const cart_t K = beam_group.geometry->Kl_as_K(parent.Kl);
          if (parent.use_n_beam) {
            int beam_n = parent.params.getBeamN();
            cmat_t A;
            af::shared<miller::index<> > indices, strong_indices;
            boost::shared_ptr<a_dyn_calculator<FloatType> > dc;
            if (beam_n > 2) {
              strong_indices = af::select(beam_group.indices.const_ref(),
                beam_group.strong_beams.const_ref());
            }
            af::shared<FloatType> angles = beam_group.get_int_angles(
              parent.Kl,
              parent.params.getIntSpan(),
              parent.params.getIntStep(),
              parent.params.getIntPoints(),
              !parent.params.isAngleInt());
            for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
              size_t beam_idx = beam_group.strong_measured_beams[i];
              miller::index<> h = beam_group.indices[beam_idx];
              complex_t Fc;
              FloatType I1 = -1, g1 = -1, I_sum=0;
              
              FloatType da = beam_group.get_diffraction_angle(h, parent.Kl);
              if (beam_n > 2) {
                std::pair<mat3_t, cart_t> FI = beam_group.compute_RMf_N(da);
                indices = utils<FloatType>::build_Ug_matrix_N(A,
                  Fcs_k, parent.mi_lookup,
                  strong_indices, K, h, FI.first, beam_n,
                  parent.params.useNBeamSg(), parent.params.getNBeamWght());
                dc = dyn_calculator_factory<FloatType>(parent.mat_type)
                  .make(indices, K, thickness);
              }
              else {
                int ii = parent.mi_lookup.find_hkl(h);
                Fc = ii != -1 ? Fcs_k[ii] : 0;
              }

              for (size_t ai = 0; ai < angles.size(); ai++) {
                std::pair<mat3_t, cart_t> r = beam_group.compute_RMf_N(angles[ai]);
                cart_t K_g = r.first * cart_t(h[0], h[1], h[2]) + K;
                FloatType I;
                if (beam_n == 2) {
                  I = std::norm(utils<FloatType>::calc_amp_2beam(
                    h, Fc, thickness, K,
                    r.first, r.second));
                }
                else {
                  I = std::norm(dc->reset(A, r.first, r.second).calc_amps_1(0));
                }
                FloatType g = K_g.length();
                if (g1 > 0) {
                  I_sum += (I + I1)*std::abs(g-g1);
                }
                g1 = g;
                I1 = I;
              }
              Is[offset + i] = I_sum;
            }
            return;
          }
          // modified Ug/Felix UNUSED for now
          if (parent.mat_type == 4) {
            cmat_t A;
            af::shared<miller::index<> >
              s = af::select(beam_group.indices.const_ref(), beam_group.strong_beams.const_ref());

            utils<FloatType>::build_Ug_matrix(
              A, Fcs_k,
              parent.mi_lookup,
              s
            );

            utils<FloatType>::modify_Ug_matrix(
              A, s, Fcs_k,
              parent.mi_lookup,
              af::select(beam_group.indices.const_ref(), beam_group.weak_beams.const_ref()),
              af::select(beam_group.excitation_errors.const_ref(), beam_group.weak_beams.const_ref()));

            af::shared<FloatType> ExpDen;
            utils<FloatType>::build_eigen_matrix_modified(
              A, K, beam_group.geometry->get_normal(),
              af::select(beam_group.gs.const_ref(), beam_group.strong_beams.const_ref()),
              af::select(beam_group.excitation_errors.const_ref(), beam_group.strong_beams.const_ref()),
              ExpDen);

            af::shared<complex_t> amps =
              utils<FloatType>::calc_amps_modified(A, ExpDen, thickness,
                beam_group.strong_measured_beams.size());
            
            for (size_t i = 0; i < amps.size(); i++) {
              Is[offset + i] = std::norm(amps[i]);
              if (CIs.size() != 0) {
                CIs[offset + i] = amps[i];
              }
            }
            return;
          }
          cmat_t A;
          af::shared<miller::index<> >
            strong_indices = af::select(beam_group.indices.const_ref(),
              beam_group.strong_beams.const_ref());
          utils<FloatType>::build_Ug_matrix(
            A, Fcs_k,
            parent.mi_lookup,
            strong_indices
          );
          boost::shared_ptr<a_dyn_calculator<FloatType> > dc =
            dyn_calculator_factory<FloatType>(parent.mat_type)
            .make(strong_indices, K, thickness);
          af::shared<FloatType> angles =
            beam_group.get_int_angles(parent.Kl, parent.params.getIntSpan(),
              parent.params.getIntStep(),
              parent.params.getIntPoints(),
              !parent.params.isAngleInt());
          for (size_t i = 0; i < beam_group.strong_measured_beams.size(); i++) {
            size_t beam_idx = beam_group.strong_measured_beams[i];
            miller::index<> h = beam_group.indices[beam_idx];
            int ii = parent.mi_lookup.find_hkl(h);
            complex_t Fc = ii != -1 ? Fcs_k[ii] : 0;
            FloatType I1 = -1, g1 = -1;
            FloatType da = beam_group.get_diffraction_angle(h, parent.Kl);
            for (size_t ai = 0; ai < angles.size(); ai++) {
              std::pair<mat3_t, cart_t> r = beam_group.compute_RMf_N(angles[ai]);
              cart_t K_g = r.first * cart_t(h[0], h[1], h[2]) + K;
              FloatType I = std::norm(dc->reset(A, r.first, r.second).calc_amps_1(i));
              FloatType g = K_g.length();
              if (g1 >= 0) {
                Is[offset + i] += (I + I1) * std::abs(g - g1) / 2;
              }
              g1 = g;
              I1 = I;
            }
          }
        }
        catch (smtbx::error const& e) {
          exception_.reset(new smtbx::error(e));
        }
        catch (std::exception const& e) {
          exception_.reset(new smtbx::error(e.what()));
        }
      }
      ed_n_shared_data const& parent;
      BeamGroup<FloatType>& beam_group;
      FloatType thickness;
      const af::shared<complex_t>& Fcs_k;
      af::shared<FloatType>& Is;
      af::shared<complex_t>& CIs;
      size_t offset;
      boost::scoped_ptr<smtbx::error> exception_;
    };

    af::shared<FloatType> process_beam_group_id(int beam_group_id,
      const af::shared<complex_t> &Fcs_k)
    {
      typename std::map<int, BeamGroup<FloatType>*>::const_iterator fi =
        beam_groups_map.find(beam_group_id);
      SMTBX_ASSERT(fi != beam_groups_map.end());
      BeamGroup<FloatType> &beam_group= *fi->second;
      af::shared<FloatType> Is_(beam_group.beams.size());
      af::shared<complex_t> CIs_;
      process_beam_group(*this, beam_group, Fcs_k, Is_, CIs_,
        false)(); // no offset - write to Is_ at 0;
      return Is_;
    }

    void process_beam_groups_mt(af::shared<FloatType> &Is_,
      af::shared<complex_t>& CIs_,
      af::shared<complex_t> const& Fcs_k)
    {
      if (thread_n < 0) {
        thread_n = builder_base<FloatType>::get_available_threads();
      }
      boost::thread_group pool;
      typedef boost::shared_ptr<process_beam_group> beam_group_processor_t;
      size_t to = 0;
      for (size_t fi = 0; fi < beam_groups.size(); fi += thread_n) {
        size_t t_end = std::min(thread_n, (int)(beam_groups.size() - fi));
        if (t_end == 0) {
          break;
        }
        std::vector<beam_group_processor_t> accumulators;
        for (int thread_idx = 0; thread_idx < t_end; thread_idx++) {
          beam_group_processor_t pf(
            new process_beam_group(*this,
              beam_groups[to],
              Fcs_k,
              Is_,
              CIs_, true)
          );
          accumulators.push_back(pf);
          pool.create_thread(boost::ref(*pf));
          to++;
        }
        pool.join_all();
        for (int thread_idx = 0; thread_idx < t_end; thread_idx++) {
          if (accumulators[thread_idx]->exception_) {
            throw* accumulators[thread_idx]->exception_.get();
          }
        }
      }
    }

    void build() {
      // Generate Fcs at current position
      {
        Is.resize(beam_n);
        Fcs.resize(beam_n);
        if (Fcs_k.size() != indices.size()) {
          Fcs_k.resize(indices.size());
          for (size_t ih = 0; ih < indices.size(); ih++) {
            Fcs_k[ih] = calc_one_h(indices[ih]) * Fc2Ug;
          }
          // expand uniq Fc to beam_group indices
          size_t offset = 0;
          for (size_t i = 0; i < beam_groups.size(); i++) {
            const af::shared<miller::index<> >& fidx = beam_groups[i].indices;
            size_t measured = beam_groups[i].strong_measured_beams.size();
            for (size_t i = 0; i < measured; i++) {
              long idx = mi_lookup.find_hkl(fidx[i]);
              Fcs[offset + i] = Fcs_k[idx];
            }
            offset += measured;
          }
        }
        // replacing dummy with Fc will allow to collect complex amplitudes
        af::shared<complex_t> dummy;
        std::fill(Is.begin(), Is.end(), 0);
        process_beam_groups_mt(Is, dummy, Fcs_k);
        if (!compute_grad) {
          return;
        }
      }
      af::shared<parameter*> params = reparamn.independent();
      af::shared<asu_parameter*> p_owners = reparamn.independent_owners(params);
        size_t param_n = 0;
      for (size_t i = 0; i < params.size(); i++) {
        param_n += params[i]->components().size();
      }
      design_matrix.resize(af::mat_grid(beam_n, param_n));
      // Generate Arrays for storing positive and negative epsilon Fcs for numerical
      // generation of gradients
      af::shared<complex_t> Fc_eps(indices.size()), dummy;
      af::shared<FloatType> Is_p(beam_n), Is_m(beam_n);
      FloatType t_eps = 2 * eps;
      for (size_t i = 0, n = 0; i < params.size(); i++) {
        parameter* param = params[i];
        af::ref<double> x = param->components();
        asu_parameter* cp = p_owners[i];
        for (size_t j = 0; j < x.size(); j++, n++) {
          x[j] += eps;
          if (cp != 0) {
            cp->evaluate(reparamn.unit_cell());
            cp->store(reparamn.unit_cell());
          }
          for (size_t i_h = 0; i_h < indices.size(); i_h++) {
            Fc_eps[i_h] = calc_one_h(indices[i_h]) * Fc2Ug;
          }
          //Generate Fcs at x+eps
          std::fill(Is_p.begin(), Is_p.end(), 0);
          process_beam_groups_mt(Is_p, dummy, Fc_eps);

          x[j] -= t_eps;
          if (cp != 0) {
            cp->evaluate(reparamn.unit_cell());
            cp->store(reparamn.unit_cell());
          }
          for (size_t i_h = 0; i_h < indices.size(); i_h++) {
            Fc_eps[i_h] = calc_one_h(indices[i_h]) * Fc2Ug;
          }
          //Generate Fcs at x-eps
          std::fill(Is_m.begin(), Is_m.end(), 0);
          process_beam_groups_mt(Is_m, dummy, Fc_eps);
          
          // compute grads
          for (size_t bi = 0; bi < beam_n; bi++) {
            FloatType grad_I = (Is_p[bi] - Is_m[bi]) / t_eps;
            design_matrix(bi, n) = grad_I;
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
    beam_group in the uniform arrays
    */
    size_t find_hkl(int beam_group_id, miller::index<> const& h) const {
      typename std::map<int, lookup_ptr_t>::const_iterator i =
        beam_group_lookups.find(beam_group_id);
      SMTBX_ASSERT(i != beam_group_lookups.end());
      typename std::map<int, BeamGroup<FloatType> *>::const_iterator fi =
        beam_groups_map.find(beam_group_id);
      long hi = i->second->find_hkl(h);
      SMTBX_ASSERT(hi >= 0);
      return hi + fi->second->offset;
    }

    const reparametrisation& reparamn;
    f_calc_function_base_t& f_calc_function;
    const cctbx::xray::fc_correction<FloatType>& fc_cr;
    const sgtbx::space_group& space_group;
    af::shared<miller::index<> > indices;
    RefinementParams<FloatType> params;
    FloatType Kvac, Kl, Fc2Ug, eps;
    int mat_type;
    size_t beam_n;
    /* to lookup an index in particular beam_group, have to keep a copy of the
    indices
    */
    typename std::map<int, lookup_ptr_t> beam_group_lookups;
    typename std::map<int, BeamGroup<FloatType>*> beam_groups_map;
    af::shared<BeamGroup<FloatType> > beam_groups;
    const cctbx::xray::thickness<FloatType>& thickness;
    bool compute_grad;
    // newly-calculated, aligned by beam_groups
    af::shared<complex_t> Fcs, Fcs_k;
    af::shared<FloatType> Is;
    // 
    af::versa<FloatType, af::mat_grid> design_matrix;
    lookup_t mi_lookup;
    int thread_n;
    bool use_n_beam;
  };

  template <typename FloatType>
  class f_calc_function_ed_n : public f_calc_function_base<FloatType> {
  public:
    ED_UTIL_TYPEDEFS;
    typedef f_calc_function_base<FloatType> f_calc_function_base_t;
    typedef ed_n_shared_data<FloatType> data_t;

    f_calc_function_ed_n(const ed_n_shared_data<FloatType>& data)
      : data(data),
      index(~0)
    {}
    f_calc_function_ed_n(const f_calc_function_ed_n& other)
      : data(other.data),
      index(~0)
    {}

    virtual void compute(
      const miller::index<>& h,
      const boost::optional<std::complex<FloatType> >& f_mask = boost::none,
      twin_fraction<FloatType> const* fraction = 0,
      bool compute_grad = true)
    {
      SMTBX_ASSERT(fraction->tag >= 0);
      index = data.find_hkl(fraction->tag, h);
    }

    virtual boost::shared_ptr<f_calc_function_base_t> fork() const {
      return boost::shared_ptr<f_calc_function_base_t>(
        new f_calc_function_ed_n(*this));
    }

    virtual FloatType get_observable() const {
      return data.Is[index];
    }
    virtual std::complex<FloatType> get_f_calc() const {
      return data.Fcs[index];
    }
    virtual af::const_ref<complex_t> get_grad_f_calc() const {
      SMTBX_NOT_IMPLEMENTED();
      throw 1;
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
  };

}}}

#endif // GUARD
