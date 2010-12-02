#ifndef SMTBX_REFINEMENT_COVARIANCE_H
#define SMTBX_REFINEMENT_COVARIANCE_H

#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/packed_matrix.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/sparse/matrix.h>

#include <cctbx/import_scitbx_af.h>
#include <cctbx/xray/parameter_map.h>
#include <cctbx/xray/scatterer.h>
#include <cctbx/uctbx.h>

namespace cctbx { namespace covariance {

  namespace af=scitbx::af;

  template<typename FloatType>
  af::versa<FloatType, af::packed_u_accessor>
  extract_covariance_matrix_for_sites(
    af::const_ref<std::size_t> const &i_seqs,
    af::const_ref<FloatType, af::packed_u_accessor> const &matrix,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<FloatType> >
      const &parameter_map)
  {
    CCTBX_ASSERT(matrix.size()
      == parameter_map.n_parameters()*(parameter_map.n_parameters()+1)/2);
    af::versa<FloatType, af::packed_u_accessor> result(3*i_seqs.size());
    for (std::size_t i=0; i<i_seqs.size(); i++) {
      for (std::size_t j=i; j<i_seqs.size(); j++) {
        std::size_t i_seq = i_seqs[i];
        std::size_t j_seq = i_seqs[j];
        cctbx::xray::parameter_indices ids_i = parameter_map[i_seq];
        cctbx::xray::parameter_indices ids_j = parameter_map[j_seq];
        CCTBX_ASSERT(ids_i.site> -1);
        CCTBX_ASSERT(ids_j.site> -1);
        for (std::size_t ii=0; ii<3; ii++) {
          for (std::size_t jj=0; jj<3; jj++) {
            if (i==j && ii>jj) continue;
            if (ids_i.site > ids_j.site) {
              result(i*3+ii,j*3+jj) = matrix(ids_j.site+jj, ids_i.site+ii);
            }
            else {
              result(i*3+ii,j*3+jj) = matrix(ids_i.site+ii, ids_j.site+jj);
            }
          }
        }
      }
    }
    return result;
  }

  template<typename FloatType>
  af::versa<FloatType, af::packed_u_accessor>
  extract_covariance_matrix_for_u_aniso(
    std::size_t i_seq,
    af::const_ref<FloatType, af::packed_u_accessor> const &matrix,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<FloatType> >
      const &parameter_map)
  {
    CCTBX_ASSERT(matrix.size()
      == parameter_map.n_parameters()*(parameter_map.n_parameters()+1)/2);
    af::versa<FloatType, af::packed_u_accessor> result(6);
    cctbx::xray::parameter_indices ids = parameter_map[i_seq];
    CCTBX_ASSERT(ids.u_aniso > -1);
    for (std::size_t i=0; i<6; i++) {
      for (std::size_t j=i; j<6; j++) {
        result(i,j) = matrix(ids.u_aniso+i, ids.u_aniso+j);
      }
    }
    return result;
  }

  template<typename FloatType>
  FloatType
  variance_for_u_iso(
    std::size_t i_seq,
    af::const_ref<FloatType, af::packed_u_accessor> const &matrix,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<FloatType> >
      const &parameter_map)
  {
    CCTBX_ASSERT(matrix.size()
      == parameter_map.n_parameters()*(parameter_map.n_parameters()+1)/2);
    cctbx::xray::parameter_indices ids = parameter_map[i_seq];
    CCTBX_ASSERT(ids.u_iso > -1);
    return matrix(ids.u_iso, ids.u_iso);
  }

  template<typename FloatType>
  scitbx::sparse::matrix<FloatType>
  covariance_orthogonalization_matrix(
    cctbx::uctbx::unit_cell const &unit_cell,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<FloatType> >
      const &parameter_map)
  {
    std::size_t n_params = parameter_map.n_parameters();
    std::size_t n_scatterers = parameter_map.n_scatterers();
    scitbx::sparse::matrix<FloatType> O(n_params, n_params);
    scitbx::mat3<FloatType> orth = unit_cell.orthogonalization_matrix();
    for (std::size_t i=0; i<n_scatterers; i++) {
      cctbx::xray::parameter_indices ids = parameter_map[i];
      if (ids.site > -1) {
        for (std::size_t j=0; j<3; j++) {
          for (std::size_t k=j; k<3; k++) {
            O(ids.site+j, ids.site+k) = orth(j, k);
          }
        }
      }
    }
    return O;
  }

  template<typename FloatType>
  af::versa<FloatType, af::packed_u_accessor>
  orthogonalize_covariance_matrix(
    af::const_ref<FloatType, af::packed_u_accessor> const &matrix,
    cctbx::uctbx::unit_cell const &unit_cell,
    cctbx::xray::parameter_map<cctbx::xray::scatterer<FloatType> >
      const &parameter_map)
  {
    CCTBX_ASSERT(matrix.size()
      == parameter_map.n_parameters()*(parameter_map.n_parameters()+1)/2);
    scitbx::sparse::matrix<FloatType> O = covariance_orthogonalization_matrix(
      unit_cell, parameter_map);
    return O.this_times_symmetric_times_this_transpose(matrix);
  }

}} // cctbx::covariance

#endif
