#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cctbx/miller.h>
#include <cctbx/error.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <fstream>
#include <vector>

namespace {

  namespace miller = cctbx::miller;
  namespace af = scitbx::af;

  class read_formatted
  {
    public:
      read_formatted() {}

      void
      new_line(std::string const& line)
      {
        len_ = line.size();
        pos_ = 0;
        buffer_.reserve(len_+1);
        for(std::size_t i=0;i<len_;i++) buffer_[i] = line[i];
        buffer_[len_] = '\0';
      }

      const char*
      next_field(std::size_t w, char fill_char)
      {
        field_.reserve(w+1);
        std::size_t i;
        for(i=0; i < w && pos_ < len_; i++, pos_++) {
          field_[i] = buffer_[pos_];
        }
        for(; i < w; i++) {
          field_[i] = fill_char;
        }
        field_[w] = '\0';
        return &*field_.begin();
      }

      int
      next_int(std::size_t w)
      {
        return std::atoi(next_field(w, '0'));
      }

      double
      next_double(std::size_t w)
      {
        return std::atof(next_field(w, '0'));
      }

      miller::index<>
      next_miller_index(std::size_t w)
      {
        miller::index<> result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = next_int(w);
        }
        return result;
      }

    protected:
      std::size_t len_;
      std::size_t pos_;
      std::vector<char> buffer_;
      std::vector<char> field_;
  };

  class no_merge_original_index_arrays
  {
    public:
      no_merge_original_index_arrays(
        std::string const& file_name,
        std::size_t n_lines_skip)
      {
        std::ifstream cin(file_name.c_str());
        std::string line;
        for (std::size_t i=0;i<n_lines_skip;i++) {
          std::getline(cin, line);
          CCTBX_ASSERT(!cin.eof());
        }
        read_formatted rf;
        while (true) {
          std::getline(cin, line);
          if (cin.eof()) break;
          rf.new_line(line);
          original_indices_.push_back(rf.next_miller_index(4));
          unique_indices_.push_back(rf.next_miller_index(4));
          batch_numbers_.push_back(rf.next_int(6));
          centric_tags_.push_back(rf.next_int(2));
          spindle_flags_.push_back(rf.next_int(2));
          asymmetric_unit_indices_.push_back(rf.next_int(3));
          i_obs_.push_back(rf.next_double(8));
          sigmas_.push_back(rf.next_double(8));
        }
      }

      af::shared<miller::index<> >
      original_indices() const { return original_indices_; }

      af::shared<miller::index<> >
      unique_indices() const { return unique_indices_; }

      af::shared<int>
      batch_numbers() const { return batch_numbers_; }

      af::shared<int>
      centric_tags() const { return centric_tags_; }

      af::shared<int>
      asymmetric_unit_indices() const { return asymmetric_unit_indices_; }

      af::shared<int>
      spindle_flags() const { return spindle_flags_; }

      af::shared<double>
      i_obs() const { return i_obs_; }

      af::shared<double>
      sigmas() const { return sigmas_; }

    private:
      af::shared<miller::index<> > original_indices_;
      af::shared<miller::index<> > unique_indices_;
      af::shared<int> batch_numbers_;
      af::shared<int> centric_tags_;
      af::shared<int> asymmetric_unit_indices_;
      af::shared<int> spindle_flags_;
      af::shared<double> i_obs_;
      af::shared<double> sigmas_;
  };

  struct no_merge_original_index_arrays_wrappers
  {
    typedef no_merge_original_index_arrays w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("no_merge_original_index_arrays", no_init)
        .def(init<std::string const&, std::size_t>())
        .def("original_indices", &w_t::original_indices)
        .def("unique_indices", &w_t::unique_indices)
        .def("batch_numbers", &w_t::batch_numbers)
        .def("centric_tags", &w_t::centric_tags)
        .def("asymmetric_unit_indices", &w_t::asymmetric_unit_indices)
        .def("spindle_flags", &w_t::spindle_flags)
        .def("i_obs", &w_t::i_obs)
        .def("sigmas", &w_t::sigmas)
      ;
    }
  };

  void init_module()
  {
    no_merge_original_index_arrays_wrappers::wrap();
  }

} // namespace <anonymous>

BOOST_PYTHON_MODULE(iotbx_scalepack_ext)
{
  init_module();
}
