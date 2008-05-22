#ifndef CCTBX_XRAY_MULTITHREADED_STRUCTURE_FACTORS_DIRECT_H
#define CCTBX_XRAY_MULTITHREADED_STRUCTURE_FACTORS_DIRECT_H

#include <cctbx/xray/structure_factors_direct.h>
#include <boost/thread/thread.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/type_traits.hpp>
#include <cctbx/import_scitbx_af.h>
#include <vector>

namespace cctbx { namespace xray { namespace structure_factors {

  namespace details {

    template <class CosSinType,
              class ScattererType=scatterer<>,
              class FormFactorType=eltbx::xray_scattering::gaussian>
    struct delayed_direct : public direct<ScattererType, FormFactorType>
    {
      typedef direct<ScattererType, FormFactorType> base_t;
      typedef typename base_t::scattering_type_registry_t
              scattering_type_registry_t;

      CosSinType const *cos_sin;
      uctbx::unit_cell const *unit_cell;
      sgtbx::space_group const *space_group;
      af::const_ref<miller::index<> > miller_indices;
      af::const_ref<ScattererType> const *scatterers;
      scattering_type_registry_t const *scattering_type_registry;

      delayed_direct(
        CosSinType const& cos_sin_,
        uctbx::unit_cell const& unit_cell_,
        sgtbx::space_group const& space_group_,
        af::const_ref<miller::index<> > const& miller_indices_,
        af::const_ref<ScattererType> const& scatterers_,
        scattering_type_registry_t const& scattering_type_registry_)
      : cos_sin(&cos_sin_),
        unit_cell(&unit_cell_),
        space_group(&space_group_),
        miller_indices(miller_indices_),
        scatterers(&scatterers_),
        scattering_type_registry(&scattering_type_registry_)
      {}

      void operator()() {
        compute(*cos_sin, *unit_cell, *space_group, miller_indices,
                *scatterers, *scattering_type_registry);
      }
    };

  } // namespace details


  template <class ScattererType=scatterer<>,
            class FormFactorType=eltbx::xray_scattering::gaussian>
  class multithreaded_direct
  {
    public:
      typedef ScattererType scatterer_type;
      typedef typename ScattererType::float_type float_type;
      typedef xray::generic_scattering_type_registry<FormFactorType>
              scattering_type_registry_t;

      af::shared<std::complex<float_type> > f_calc() const {
        return f_calc_;
      }

      multithreaded_direct() {}

      multithreaded_direct(
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        scattering_type_registry_t const& scattering_type_registry,
        unsigned n_threads)
      {
        math::cos_sin_exact<float_type> cos_sin;
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_type_registry, n_threads);
      }

      template<class CosSinType>
      multithreaded_direct(
        CosSinType const& cos_sin,
        uctbx::unit_cell const& unit_cell,
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<ScattererType> const& scatterers,
        scattering_type_registry_t const& scattering_type_registry,
        unsigned n_threads)
      {
        compute(cos_sin, unit_cell, space_group, miller_indices,
                scatterers, scattering_type_registry, n_threads);
      }

    protected:
      template<class CosSinType>
      void compute(CosSinType const& cos_sin,
                   uctbx::unit_cell const& unit_cell,
                   sgtbx::space_group const& space_group,
                   af::const_ref<miller::index<> > const& miller_indices,
                   af::const_ref<ScattererType> const& scatterers,
                   scattering_type_registry_t const& scattering_type_registry,
                   unsigned n_threads)
      {
        typedef details::delayed_direct<CosSinType,
                                        ScattererType,
                                        FormFactorType> direct_t;
        std::size_t normal_chunk_size = miller_indices.size()/n_threads;
        std::size_t last_chunk_size = miller_indices.size()
                                        - (n_threads - 1)*normal_chunk_size;
        boost::thread_group workers;
        std::vector<direct_t> direct_chunks;
        for(unsigned i_thread=0; i_thread < n_threads; ++i_thread) {
          std::size_t first = i_thread*normal_chunk_size;
          std::size_t chunk_size = i_thread < n_threads - 1 ? normal_chunk_size
                                                            : last_chunk_size;
          af::const_ref<miller::index<> > miller_chunk(&miller_indices[first],
                                                       chunk_size);
          direct_t direct_chunk(cos_sin, unit_cell, space_group, miller_chunk,
                                scatterers, scattering_type_registry);
          direct_chunks.push_back(direct_chunk);
          workers.add_thread(new boost::thread(direct_chunk));
        }
        workers.join_all();
        for(unsigned i_thread=0; i_thread < n_threads; ++i_thread) {
          af::const_ref<std::complex<float_type> >
            f_calc_chunk = direct_chunks[i_thread].f_calc().ref();
          f_calc_.extend(f_calc_chunk.begin(), f_calc_chunk.end());
        }
      }

      af::shared<std::complex<float_type> > f_calc_;
  };


}}} // cctbx::xray::structure_factors

#endif // GUARD
