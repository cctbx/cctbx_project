#ifndef SCITBX_MATRIX_TENSORS_H
#define SCITBX_MATRIX_TENSORS_H

#include <vector>
#include <cctbx/sgtbx/rot_mx.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>

/*
Note that in multithreaded environment the librray has to be initialised in the
main thread before executing any threads accessing these objects
For this scitbx::matrix::tensors::initialise<FloatType> function is implemented
here. If you are doing memory leak hunting - also consider using the finalise
method.
*/
namespace scitbx { namespace matrix { namespace tensors {
  namespace utils {
    static size_t factorial(const size_t& a) {
      size_t b = 1;
      for (size_t i = 2; i <= a; i++) {
        b *= i;
      }
      return b;
    }
    static size_t calc_multiplicity(const size_t *idx, size_t sz) {
      size_t mult = 1;
      for (size_t i = 0; i < 3; i++) {
        mult *= factorial(idx[i]);
      }
      return factorial(sz) / mult;
    }
  } // tensor utils

  using cctbx::sgtbx::rot_mx;

  template <class heir_t, typename FloatType>
  class tensor_base {
    heir_t &self() { return *(heir_t*)this; }
    const heir_t &self() const { return *(heir_t*)this; }

  protected:
    /* https://en.wikipedia.org/wiki/Heap%27s_algorithm
    initialises equivalent map indices
    */
    static void init_index(size_t k, std::vector<int> &idx, size_t linear_index) {
      if (k == 1) {
        heir_t::get_linear_index_(idx) = linear_index;
      }
      else {
        init_index(k - 1, idx, linear_index);
        for (size_t i = 0; i < k - 1; i++) {
          if ((k & 1) == 1) { // odd?
            std::swap(idx[0], idx[k - 1]);
          }
          else {
            std::swap(idx[i], idx[k - 1]);
          }
          init_index(k - 1, idx, linear_index);
        }
      }
    }
    /* initialises the rank-D index to linear index map and the indices
    multiplicity
    */
    static void init_map_m() {
      const std::vector<std::vector<int> > &indices = heir_t::get_indices();
      for (size_t i = 0; i < indices.size(); i++) {
        std::vector<int> idx = indices[i];
        init_index(idx.size(), idx, i);
        size_t mps[3] = { 0, 0, 0 };
        for (size_t j = 0; j < idx.size(); j++) {
          mps[idx[j]]++;
        }
        get_multiplicity_()[i] = utils::calc_multiplicity(mps, heir_t::rank());
      }
    }

    static std::vector<size_t> &get_multiplicity_() {
      static std::vector<size_t> multiplicity(heir_t::size());
      return multiplicity;
    }

  public:
    const FloatType & operator [](size_t i) const { return self().data_[i]; }

    FloatType & operator [](size_t i) { return self().data_[i]; }

    static size_t get_linear_idx(const std::vector<int> &idx) {
      return heir_t::get_linear_index_(idx);
    }

    size_t get_multiplicity(size_t i) const {
      return get_multiplicity_()[i];
    }

    static const std::vector<std::vector<int> > &get_indices() {
      return heir_t::get_indices_();
    }
    /* this might be required in a multithreaded environment!
    */
    static void initialise() {
      heir_t::get_map();
    }

    template <typename NumType>
    FloatType sum_up(const scitbx::vec3<NumType> &h) const {
      FloatType r = 0;
      const std::vector<std::vector<int> > &indices = get_indices();
      for (size_t i = 0; i < indices.size(); i++) {
        const std::vector<int> &idx = indices[i];
        FloatType prod_h = 1;
        for (size_t j = 0; j < heir_t::rank(); j++) {
          prod_h *= h[idx[j]];
        }
        r += prod_h * get_multiplicity(i) *
          self().data_[get_linear_idx(idx)];
      }
      return r;
    }

    template <typename NumType>
    af::shared<FloatType> gradient_coefficients(
      const scitbx::vec3<NumType> &h) const
    {
      af::shared<FloatType> r(heir_t::size());
      const std::vector<std::vector<int> > &indices = get_indices();
      for (size_t i = 0; i < indices.size(); i++) {
        const std::vector<int> &idx = indices[i];
        FloatType prod_h = 1;
        for (size_t j = 0; j < heir_t::rank(); j++) {
          prod_h *= h[idx[j]];
        }
        r[i] = prod_h * get_multiplicity(i);
      }
      return r;
    }

    const af::shared<FloatType> &data() const {
      return self().data_;
    }
  };

  template <typename FloatType = double>
  class tensor_rank_2 : public tensor_base<tensor_rank_2<FloatType>, FloatType> {
    typedef tensor_base<tensor_rank_2<FloatType>, FloatType> parent_t;
    friend class tensor_base<tensor_rank_2<FloatType>, FloatType>;
    af::shared<FloatType> data_;
    /*
    should be to be identical with higher order
    0 0, 0 1, 0 2, 1 1, 1 2, 2 2
    but keep in cctbx format for simpler testing:
    0 0, 1 1, 2 2, 0 1, 0 2, 1 2
    */
    static size_t **& get_map_() {
      static size_t ** map = 0;
      return map;
    }

    static size_t **& get_map() {
      size_t **& r = get_map_();
      if (r == 0) {
        r = build_map();
        /* generic procedure
        parent_t::init_map_m();
        */
        // override to match cctbx
        r[0][0] = 0; r[0][1] = 3; r[0][2] = 4;
        r[1][0] = 3; r[1][1] = 1; r[1][2] = 5;
        r[2][0] = 4; r[2][1] = 5; r[2][2] = 2;
        parent_t::get_multiplicity_()[0] = 1;
        parent_t::get_multiplicity_()[1] = 1;
        parent_t::get_multiplicity_()[2] = 1;
        parent_t::get_multiplicity_()[3] = 2;
        parent_t::get_multiplicity_()[4] = 2;
        parent_t::get_multiplicity_()[5] = 2;
      }
      return r;
    }

    static size_t &get_linear_index_(const std::vector<int> &idx) {
      return get_map()[idx[0]][idx[1]];
    }

    static size_t ** build_map() {
      size_t ** map = new size_t *[3];
      for (size_t i = 0; i < 3; i++) {
        map[i] = new size_t[3];
      }
      return map;
    }

    static const std::vector<std::vector<int> > &get_indices_() {
      static std::vector<std::vector<int> > indices;
      if (indices.empty()) {
        indices.resize(size());
        for (size_t i = 0; i < size(); i++) {
          indices[i].resize(2);
        }
        indices[0][0] = 0; indices[0][1] = 0;
        indices[1][0] = 1; indices[1][1] = 1;
        indices[2][0] = 2; indices[2][1] = 2;
        indices[3][0] = 0; indices[3][1] = 1;
        indices[4][0] = 0; indices[4][1] = 2;
        indices[5][0] = 1; indices[5][1] = 2;

        /* overriding to match cctbx
          for (int i = 0, idx = 0; i < 3; i++) {
            for (int j = i; j < 3; j++, idx++) {
              indices[idx].resize(2);
              indices[idx][0] = i;
              indices[idx][1] = j;
            }
          }
          */
      }
      return indices;
    }

   public:
    tensor_rank_2()
      : data_(6)
    {}

    tensor_rank_2(const af::shared<FloatType> &data)
      : data_(data)
    {
      SCITBX_ASSERT(data_.size() == size());
    }

    const FloatType & operator ()(size_t i, size_t j) const {
      return data_[get_map()[i][j]];
    }
    FloatType & operator ()(size_t i, size_t j) {
      return data_[get_map()[i][j]];
    }

    static af::shared <FloatType> get_transform(const std::vector<int> &idx,
      const rot_mx &rm)
    {
      tensor_rank_2 result;
      int i = idx[0], j = idx[1];
      for (int r = 0; r < 3; r++) {
        for (int s = 0; s < 3; s++) {
          result(s, r) += rm(i, r)*rm(j, s);
        }
      }
      return result.data_;
    }

    static size_t rank() { return 2; }
    static size_t size() { return 6; }

    static size_t linearise(size_t i, size_t j) {
      return get_map()[i][j];
    }

    static void cleanup() {
      size_t **map = get_map_();
      if (map != 0) {
        get_map_() = 0;
        for (size_t i = 0; i < 3; i++) {
          delete[] map[i];
        }
        delete map;
      }
    }
  }; // class scitbx::matrx::tensors::tensor_rank_2

  template <typename FloatType = double>
  class tensor_rank_3 : public tensor_base<tensor_rank_3<FloatType>, FloatType>{
    typedef tensor_base<tensor_rank_3<FloatType>, FloatType> parent_t;
    friend class tensor_base<tensor_rank_3<FloatType>, FloatType>;
    af::shared<FloatType> data_;
    /*
    0 0 0, 0 0 1, 0 0 2, 0 1 1, 0 1 2
    0 2 2, 1 1 1, 1 1 2, 1 2 2, 2 2 2
    */
  protected:
    static size_t ***& get_map_() {
      static size_t *** map = 0;
      return map;
    }

    static size_t ***& get_map() {
      size_t ***& r = get_map_();
      if (r == 0) {
        r = build_map();
        parent_t::init_map_m();
      }
      return r;
    }

    static size_t &get_linear_index_(const std::vector<int> &idx) {
      return get_map()[idx[0]][idx[1]][idx[2]];
    }

    static size_t *** build_map() {
      size_t *** map = new size_t **[3];
      for (size_t i = 0; i < 3; i++) {
        map[i] = new size_t*[3];
        for (size_t j = 0; j < 3; j++) {
          map[i][j] = new size_t[3];
        }
      }
      return map;
    }

    static const std::vector<std::vector<int> > &get_indices_() {
      static std::vector<std::vector<int> > indices;
      if (indices.empty()) {
        indices.resize(size());
        for (int i = 0, idx = 0; i < 3; i++) {
          for (int j = i; j < 3; j++) {
            for (int k = j; k < 3; k++, idx++) {
              indices[idx].resize(rank());
              indices[idx][0] = i;
              indices[idx][1] = j;
              indices[idx][2] = k;
            }
          }
        }
      }
      return indices;
    }
  public:
    tensor_rank_3()
      : data_(10)
    {}

    tensor_rank_3(const af::shared<FloatType> &data)
      : data_(data)
    {
      SCITBX_ASSERT(data_.size() == size());
    }

    const FloatType & operator ()(size_t i, size_t j, size_t k) const {
      return data_[get_map()[i][j][k]];
    }
    FloatType & operator ()(size_t i, size_t j, size_t k) {
      return data_[get_map()[i][j][k]];
    }

    static af::shared <FloatType> get_transform(const std::vector<int> &idx,
      const rot_mx &rm)
    {
      tensor_rank_3 result;
      int i = idx[0], j = idx[1], k = idx[2];
      for (int r = 0; r < 3; r++) {
        for (int s = 0; s < 3; s++) {
          for (int t = 0; t < 3; t++) {
            result(r, s, t) += rm(i, r)*rm(j, s)*rm(k, t);
          }
        }
      }
      return result.data_;
    }

    static size_t rank() { return 3; }
    static size_t size() { return 10; }

    static size_t linearise(size_t i, size_t j, size_t k) {
      return get_map()[i][j][k];
    }

    static void cleanup() {
      size_t ***map = get_map_();
      if (map != 0) {
        get_map_() = 0;
        for (size_t i = 0; i < 3; i++) {
          for (size_t j = 0; j < 3; j++) {
            delete[] map[i][j];
          }
          delete[] map[i];
        }
        delete map;
      }
    }
  }; // class scitbx::matrix::tensors::tensor_rank_3

  template <typename FloatType = double>
  class tensor_rank_4 : public tensor_base<tensor_rank_4<FloatType>, FloatType> {
    typedef tensor_base<tensor_rank_4<FloatType>, FloatType> parent_t;
    friend class tensor_base<tensor_rank_4<FloatType>, FloatType>;
    af::shared<FloatType> data_;
    /*
    0 0 0 0, 0 0 0 1, 0 0 0 2, 0 0 1 1, 0 0 1 2
    0 0 2 2, 0 1 1 1, 0 1 1 2, 0 1 2 2, 0 2 2 2
    1 1 1 1, 1 1 1 2, 1 1 2 2, 1 2 2 2, 2 2 2 2
    */
    static size_t ****& get_map_() {
      static size_t **** map = 0;
      return map;
    }

    static size_t ****& get_map() {
      size_t ****& r = get_map_();
      if (r == 0) {
        r = build_map();
        parent_t::init_map_m();
      }
      return r;
    }

    static size_t &get_linear_index_(const std::vector<int> &idx) {
      return get_map()[idx[0]][idx[1]][idx[2]][idx[3]];
    }

    static size_t **** build_map() {
      size_t **** map = new size_t ***[3];
      for (size_t i = 0; i < 3; i++) {
        map[i] = new size_t **[3];
        for (size_t j = 0; j < 3; j++) {
          map[i][j] = new size_t *[3];
          for (size_t k = 0; k < 3; k++) {
            map[i][j][k] = new size_t[3];
          }
        }
      }
      return map;
    }

    static const std::vector<std::vector<int> > &get_indices_() {
      static std::vector<std::vector<int> > indices;
      if (indices.empty()) {
        indices.resize(size());
        for (int i = 0, idx = 0; i < 3; i++) {
          for (int j = i; j < 3; j++) {
            for (int k = j; k < 3; k++) {
              for (int l = k; l < 3; l++, idx++) {
                indices[idx].resize(4);
                indices[idx][0] = i;
                indices[idx][1] = j;
                indices[idx][2] = k;
                indices[idx][3] = l;
              }
            }
          }
        }
      }
      return indices;
    }

  public:
    tensor_rank_4()
      : data_(15)
    {}

    tensor_rank_4(const af::shared<FloatType> &data)
      : data_(data)
    {
      SCITBX_ASSERT(data_.size() == size());
    }

    const FloatType & operator ()(size_t i, size_t j, size_t k, size_t l) const {
      return data_[get_map()[i][j][k][l]];
    }
    FloatType & operator ()(size_t i, size_t j, size_t k, size_t l) {
      return data_[get_map()[i][j][k][l]];
    }

    static af::shared <FloatType> get_transform(const std::vector<int> &idx,
      const rot_mx &rm)
    {
      tensor_rank_4 result;
      int i = idx[0], j = idx[1], k = idx[2], l = idx[3];
      for (int r = 0; r < 3; r++) {
        for (int s = 0; s < 3; s++) {
          for (int t = 0; t < 3; t++) {
            for (int u = 0; u < 3; u++) {
              result(r, s, t, u) += rm(i, r)*rm(j, s)*rm(k, t)*rm(l, u);
            }
          }
        }
      }
      return result.data_;
    }

    static size_t rank() { return 4; }
    static size_t size() { return 15; }

    static size_t linearise(size_t i, size_t j, size_t k, size_t l) {
      return get_map()[i][j][k][l];
    }

    static void cleanup() {
      size_t ****map = get_map_();
      if (map != 0) {
        get_map_() = 0;
        for (size_t i = 0; i < 3; i++) {
          for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
              delete[] map[i][j][k];
            }
            delete[] map[i][j];
          }
          delete[] map[i];
        }
        delete map;
      }
    }
  }; // class scitbx::matrix::tensors::tensor_rank_4

  template <typename FloatType>
  static void initialise() {
    tensor_rank_2<FloatType>::initialise();
    tensor_rank_3<FloatType>::initialise();
    tensor_rank_4<FloatType>::initialise();
  }

  template <typename FloatType>
  static void finalise() {
    tensor_rank_2<FloatType>::cleanup();
    tensor_rank_3<FloatType>::cleanup();
    tensor_rank_4<FloatType>::cleanup();
  }
}}} // namespace scitbx::matrix::tensors

#endif // SCITBX_MATRIX_TENSORS_H
