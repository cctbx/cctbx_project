#include <scitbx/lbfgsb/raw.h>

namespace scitbx { namespace lbfgsb {

  template <typename FloatType=double>
  class minimizer
  {
    public:
      minimizer() {}

      minimizer(
        int const& n,
        int const& m,
        af::shared<FloatType> l,
        af::shared<FloatType> u,
        af::shared<int> nbd,
        FloatType const& factr,
        FloatType const& pgtol,
        int const& iprint)
      :
        n_(n),
        m_(m),
        l_(l),
        u_(u),
        nbd_(nbd),
        factr_(factr),
        pgtol_(pgtol),
        wa_(2*m*n+4*n+12*m*m+12*m, FloatType(0)),
        iwa_(3*n, int(0)),
        task_("START"),
        iprint_(iprint),
        lsave_(4, false),
        isave_(44, int(0)),
        dsave_(29, FloatType(0))
      {
        SCITBX_ASSERT(l.size() == n);
        SCITBX_ASSERT(u.size() == n);
        SCITBX_ASSERT(nbd.size() == n);
      }

      std::string
      task() const { return task_; }

      void
      request_stop()
      {
        task_ = "STOP: QUICK";
      }

      void
      request_stop_with_restore()
      {
        task_ = "STOP: CPU";
      }

      std::string
      process(
        af::ref<FloatType> const& x,
        FloatType const& f,
        af::ref<FloatType> const& g)
      {
        f_ = f;
        raw::setulb(
          n_,
          m_,
          raw::ref1<FloatType>(x),
          raw::ref1<FloatType>(l_.ref()),
          raw::ref1<FloatType>(u_.ref()),
          raw::ref1<int>(nbd_.ref()),
          f_,
          raw::ref1<FloatType>(g),
          factr_,
          pgtol_,
          raw::ref1<FloatType>(wa_.ref()),
          raw::ref1<int>(iwa_.ref()),
          task_,
          iprint_,
          csave_,
          raw::ref1<bool>(lsave_.ref()),
          raw::ref1<int>(isave_.ref()),
          raw::ref1<FloatType>(dsave_.ref()));
        return task_;
      }

    protected:
      int n_;
      int m_;
      af::shared<FloatType> l_;
      af::shared<FloatType> u_;
      af::shared<int> nbd_;
      FloatType f_;
      FloatType factr_;
      FloatType pgtol_;
      af::shared<FloatType> wa_;
      af::shared<int> iwa_;
      std::string task_;
      int iprint_;
      std::string csave_;
      af::shared<bool> lsave_;
      af::shared<int> isave_;
      af::shared<FloatType> dsave_;
  };

}} // namespace scitbx::lbfgsb
