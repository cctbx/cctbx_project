#ifndef SCITBX_CUBICLES_H
#define SCITBX_CUBICLES_H

#include <scitbx/math/utils.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/misc_functions.h>
#include <map>
#include <cstdio>

namespace scitbx {

  namespace detail {

    template <typename FloatType>
    void
    throw_show_cubicle_dimensions(
      const char* message,
      vec3<FloatType> const& space_span,
      FloatType const& cubicle_edge,
      af::c_grid<3, unsigned> const& n_cubicles,
      unsigned long max_number_of_bytes)
    {
      char buf[512];
      std::sprintf(buf,
        "%s\n"
        "  This may be due to unreasonable parameters:\n"
        "    cubicle_edge=%.6g\n"
        "    space_span=(%.6g,%.6g,%.6g)\n"
        "    n_cubicles=(%u,%u,%u)\n"
        "    max_number_of_bytes=%lu",
        message,
        cubicle_edge,
        space_span[0], space_span[1], space_span[2],
        n_cubicles[0], n_cubicles[1], n_cubicles[2],
        max_number_of_bytes);
      throw std::runtime_error(buf);
    }

    inline
    unsigned long
    cubicles_max_memory_allocation(unsigned long number_of_bytes, bool set)
    {
      static unsigned long n_bytes = 0;
      if (set) n_bytes = number_of_bytes;
      return n_bytes;
    }

  } // namespace detail

  inline
  void
  cubicles_max_memory_allocation_set(unsigned long number_of_bytes)
  {
    detail::cubicles_max_memory_allocation(number_of_bytes, true);
  }

  inline
  unsigned long
  cubicles_max_memory_allocation_get()
  {
    return detail::cubicles_max_memory_allocation(0, false);
  }

  template <typename CubicleContentType, typename FloatType=double>
  struct cubicles
  {
    typedef math::float_int_conversions<FloatType, int> fic;

    //! Default constructor. Some data members are not initialized!
    cubicles() {}

    cubicles(
      vec3<FloatType> const& space_min_,
      vec3<FloatType> const& space_span,
      FloatType const& cubicle_edge_,
      FloatType const& epsilon)
    :
      space_min(space_min_),
      cubicle_edge(cubicle_edge_*(1+epsilon))
    {
      init_finalize(space_span, cubicle_edge_, epsilon);
    }

    cubicles(
      af::const_ref<vec3<FloatType> > const& sites_cart,
      FloatType const& cubicle_edge_,
      FloatType const& epsilon)
    :
      cubicle_edge(cubicle_edge_*(1+epsilon))
    {
      vec3<FloatType> space_span;
      if (sites_cart.size() == 0) {
        space_min.fill(0);
        space_span.fill(cubicle_edge_);
      }
      else {
        space_min = sites_cart[0];
        space_span = sites_cart[0];
        for(std::size_t i=0;i<sites_cart.size();i++) {
          for(unsigned j=0;j<3;j++) {
            math::update_min(space_min[j], sites_cart[i][j]);
            math::update_max(space_span[j], sites_cart[i][j]);
          }
        }
        space_span -= space_min;
      }
      init_finalize(space_span, cubicle_edge_, epsilon);
    }

    void
    init_finalize(
      vec3<FloatType> const& space_span,
      FloatType const& cubicle_edge_,
      FloatType const& epsilon)
    {
      SCITBX_ASSERT(epsilon > 0);
      SCITBX_ASSERT(epsilon < 0.01);
      SCITBX_ASSERT(cubicle_edge > 0);
      af::c_grid<3, unsigned> n_cubicles;
      for(std::size_t i=0;i<3;i++) {
        n_cubicles[i] = static_cast<unsigned>(
          std::max(1, fic::iceil(space_span[i] / cubicle_edge)));
      }
      unsigned long max_alloc = cubicles_max_memory_allocation_get();
      if (math::unsigned_product_leads_to_overflow(
            n_cubicles.begin(), 3)) {
        detail::throw_show_cubicle_dimensions(
          "Excessive number of cubicles:",
          space_span,
          cubicle_edge_,
          n_cubicles,
          max_alloc);
      }
      if (max_alloc != 0) {
        double estimated_memory_allocation = static_cast<double>(
          n_cubicles.size_1d());
        estimated_memory_allocation *= static_cast<double>(
          sizeof(CubicleContentType));
        if (estimated_memory_allocation > static_cast<double>(max_alloc)) {
          detail::throw_show_cubicle_dimensions(
            "Estimated memory allocation for cubicles exceeds"
            " max_number_of_bytes:",
            space_span,
            cubicle_edge_,
            n_cubicles,
            max_alloc);
        }
      }
      try {
        memory.resize(n_cubicles);
      }
      catch(const std::bad_alloc&) {
        detail::throw_show_cubicle_dimensions(
          "Not enough memory for cubicles:",
          space_span,
          cubicle_edge_,
          n_cubicles,
          max_alloc);
      }
      ref = memory.ref();
    }

    template <typename SiteType>
    vec3<unsigned>
    i_cubicle(
      SiteType const& site) const
    {
      vec3<FloatType> delta = site - space_min;
      vec3<unsigned> result;
      for(unsigned i=0;i<3;i++) {
        int j = fic::ifloor(delta[i] / cubicle_edge);
        // compensate for rounding errors
        if      (j < 0) j = 0;
        else if (j >= ref.accessor()[i]) j = ref.accessor()[i]-1;
        result[i] = static_cast<unsigned>(j);
      }
      return result;
    }

    template <typename SiteType>
    vec3<int>
    j_cubicle(
      SiteType const& site) const
    {
      vec3<FloatType> delta = site - space_min;
      vec3<int> result;
      for(unsigned i=0;i<3;i++) {
        result[i] = fic::ifloor(delta[i] / cubicle_edge);
      }
      return result;
    }

    std::map<long, long>
    cubicle_size_counts() const
    {
      std::map<long, long> result;
      for(std::size_t i=0;i<ref.size();i++) {
        result[ref[i].size()]++;
      }
      return result;
    }

    vec3<FloatType> space_min;
    FloatType cubicle_edge;
    af::versa<CubicleContentType, af::c_grid<3, unsigned> > memory;
    af::ref<CubicleContentType, af::c_grid<3, unsigned> > ref;
  };

} // namespace scitbx

#endif // SCITBX_CUBICLES_H
