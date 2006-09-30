#ifndef SCITBX_FFTPACK_DETAIL_REF_H
#define SCITBX_FFTPACK_DETAIL_REF_H

namespace scitbx { namespace fftpack { namespace detail {

  template <class ElementType>
  class ref_2d_tp
  {
    public:
      ref_2d_tp(ElementType* Start,
                std::size_t nx,
                std::size_t)
        : start_(Start),
          nx_(nx) {}
      ElementType&
      operator()(std::size_t ix, std::size_t iy)
      {
        return start_[iy * nx_ + ix];
      }
    private:
      ElementType* start_;
      std::size_t nx_;
  };

  template <class ElementType>
  class ref_3d_tp
  {
    public:
      ref_3d_tp(ElementType* Start,
                std::size_t nx,
                std::size_t ny,
                std::size_t)
        : start_(Start),
          nx_(nx),
          ny_(ny) {}
      ElementType&
      operator()(std::size_t ix, std::size_t iy, std::size_t iz)
      {
        return start_[(iz * ny_ + iy) * nx_ + ix];
      }
    private:
      ElementType* start_;
      std::size_t nx_;
      std::size_t ny_;
  };

}}} // namespace scitbx::fftpack::detail

#endif // SCITBX_FFTPACK_DETAIL_REF_H
