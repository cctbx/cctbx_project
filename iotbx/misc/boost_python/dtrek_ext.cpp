#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cctbx/sgtbx/space_group.h>
#include <cctbx/error.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>

namespace {

  namespace af = scitbx::af;

  class anomalous_combined
  {
    public:
      anomalous_combined(
        cctbx::sgtbx::space_group const& space_group,
        af::const_ref<cctbx::miller::index<> > const& miller_indices,
        af::const_ref<double> const& data_plus,
        af::const_ref<double> const& sigmas_plus,
        af::const_ref<double> const& data_minus,
        af::const_ref<double> const& sigmas_minus)
      {
        CCTBX_ASSERT(data_plus.size() == miller_indices.size());
        CCTBX_ASSERT(sigmas_plus.size() == miller_indices.size());
        CCTBX_ASSERT(data_minus.size() == miller_indices.size());
        CCTBX_ASSERT(sigmas_minus.size() == miller_indices.size());
        for(std::size_t i=0;i<miller_indices.size();i++) {
          miller_indices_.push_back(miller_indices[i]);
          if (!space_group.is_centric(miller_indices[i])) {
            miller_indices_.push_back(-miller_indices[i]);
            data_.push_back(data_plus[i]);
            data_.push_back(data_minus[i]);
            sigmas_.push_back(sigmas_plus[i]);
            sigmas_.push_back(sigmas_minus[i]);
          }
          else {
            double sum_w_data = 0;
            double sum_w = 0;
            if (sigmas_plus[i] > 0) {
              double w = 1 / (sigmas_plus[i]*sigmas_plus[i]);
              sum_w_data += w * data_plus[i];
              sum_w += w;
            }
            if (sigmas_minus[i] > 0) {
              double w = 1 / (sigmas_minus[i]*sigmas_minus[i]);
              sum_w_data += w * data_minus[i];
              sum_w += w;
            }
            if (sum_w == 0) {
              data_.push_back(data_plus[i]);
              sigmas_.push_back(sigmas_plus[i]);
            }
            else {
              data_.push_back(sum_w_data / sum_w);
              sigmas_.push_back(1/std::sqrt(sum_w));
            }
          }
        }
      }

      af::shared<cctbx::miller::index<> >
      miller_indices() const { return miller_indices_; }

      af::shared<double>
      data() const { return data_; }

      af::shared<double>
      sigmas() const { return sigmas_; }

    private:
      af::shared<cctbx::miller::index<> > miller_indices_;
      af::shared<double> data_;
      af::shared<double> sigmas_;
  };

  struct anomalous_combined_wrappers
  {
    typedef anomalous_combined w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("anomalous_combined", no_init)
        .def(init<cctbx::sgtbx::space_group const&,
                  af::const_ref<cctbx::miller::index<> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&>())
        .def("miller_indices", &w_t::miller_indices)
        .def("data", &w_t::data)
        .def("sigmas", &w_t::sigmas)
      ;
    }
  };

  void init_module()
  {
    anomalous_combined_wrappers::wrap();
  }

} // namespace <anonymous>

BOOST_PYTHON_MODULE(iotbx_dtrek_ext)
{
  init_module();
}
