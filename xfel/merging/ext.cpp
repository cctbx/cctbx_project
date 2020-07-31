#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/miller.h>
#include <dials/array_family/reflection_table.h>
#include <algorithm>

typedef cctbx::miller::index<> miller_index_t;

namespace sx_merging {

  /**
   * This visitor solves 2 problems:
   * 1. Reflection tables have typed columns and the templated operator()
   * function allows us to modify the columns without knowing the type.
   * 2. Because of 1, we have to iterate over keys which means we
   * cannot use the reflection table operator [] because it checks
   * if the columns are all the same length. So here we extract the
   * columns first, then modify them.
   */
  struct hkl_splitter_visitor : public boost::static_visitor<void> {
    dials::af::reflection_table &src;
    typename dials::af::reflection_table::key_type key;
    const std::map<miller_index_t, size_t> chunk_lookup;
    const std::vector<dials::af::reflection_table> tables;
    dials::af::shared<miller_index_t> src_hkls;

    hkl_splitter_visitor(dials::af::reflection_table &src_,
                                   typename dials::af::reflection_table::key_type key_,
                                   const std::map<miller_index_t, size_t> &chunk_lookup_,
                                   const std::vector<dials::af::reflection_table> &tables_)

        : src(src_), key(key_), chunk_lookup(chunk_lookup_), tables(tables_) {
          src_hkls = src.get<miller_index_t>("miller_index_asymmetric");
        }

    template <typename U>
    void operator()(const dials::af::shared<U> &other_column) const {
      std::vector<dials::af::shared<U> > all_columns;
      for (size_t i = 0; i < tables.size(); i++) {
        dials::af::reflection_table::const_iterator it = tables[i].find(key);
        DIALS_ASSERT(it != tables[i].end());
        dials::af::shared<U> col = boost::get<dials::af::shared<U> >(it->second);
        all_columns.push_back(col);
      }

      for (size_t i = 0; i < src.size(); ++i) {
        miller_index_t hkl = src_hkls[i];
        std::map<miller_index_t, size_t>::const_iterator it = chunk_lookup.find(hkl);
        if (it == chunk_lookup.end())
          continue;
        size_t table_idx = it->second;
        all_columns[table_idx].push_back(other_column[i]);
      }
    }
  };


  void get_hkl_chunks_cpp(dials::af::reflection_table reflections,
                          const scitbx::af::shared<miller_index_t>& hkl_list,
                          const scitbx::af::shared<int>& chunk_id_list,
                          boost::python::list hkl_chunks_cpp){
    // set up a map hkl:chunk_id
    std::map<miller_index_t, size_t> chunk_lookup;
    for(size_t i = 0UL; i < hkl_list.size(); ++i){
      chunk_lookup[hkl_list[i]] = (size_t)chunk_id_list[i];
    }

    int n_chunks = boost::python::len(hkl_chunks_cpp);
    std::vector<dials::af::reflection_table> tables;
    for(size_t i=0UL; i < n_chunks; ++i){
      tables.push_back(boost::python::extract<dials::af::reflection_table>(hkl_chunks_cpp[i]));
    }

    // distribute reflections over chunks
    for (dials::af::reflection_table::const_iterator it = reflections.begin(); it != reflections.end(); ++it) {
      hkl_splitter_visitor visitor(reflections, it->first, chunk_lookup, tables);
      it->second.apply_visitor(visitor);
    }
  }

  /**
   * This visitor solves 2 problems:
   * 1. Reflection tables have typed columns and the templated operator()
   * function allows us to modify the columns without knowing the type.
   * 2. Because of 1, we have to iterate over keys which means we
   * cannot use the reflection table operator [] because it checks
   * if the columns are all the same length. So here we extract the
   * columns first, then modify them.
   */
  struct experiment_id_splitter_visitor : public boost::static_visitor<void> {
    dials::af::reflection_table &src;
    typename dials::af::reflection_table::key_type key;
    const std::map<std::string, size_t> chunk_lookup;
    const std::vector<dials::af::reflection_table> tables;
    dials::af::shared<std::string> expt_ids;

    experiment_id_splitter_visitor(dials::af::reflection_table &src_,
                                   typename dials::af::reflection_table::key_type key_,
                                   const std::map<std::string, size_t> &chunk_lookup_,
                                   const std::vector<dials::af::reflection_table> &tables_)

        : src(src_), key(key_), chunk_lookup(chunk_lookup_), tables(tables_) {
          expt_ids = src.get<std::string>("exp_id");
        }

    template <typename U>
    void operator()(const dials::af::shared<U> &other_column) const {
      std::vector<dials::af::shared<U> > all_columns;
      for (size_t i = 0; i < tables.size(); i++) {
        dials::af::reflection_table::const_iterator it = tables[i].find(key);
        DIALS_ASSERT(it != tables[i].end());
        dials::af::shared<U> col = boost::get<dials::af::shared<U> >(it->second);
        all_columns.push_back(col);
      }

      for (size_t i = 0; i < src.size(); ++i) {
        std::string experiment_id = expt_ids[i];
        std::map<std::string, size_t>::const_iterator it = chunk_lookup.find(experiment_id);
        DIALS_ASSERT(it != chunk_lookup.end());
        size_t table_idx = it->second;
        all_columns[table_idx].push_back(other_column[i]);
      }
    }
  };

  //template <typename T>
  void split_reflections_by_experiment_chunks_cpp(dials::af::reflection_table reflections,
                                                  const scitbx::af::shared<std::string>& exp_id_list,
                                                  const scitbx::af::shared<int>& chunk_id_list,
                                                  boost::python::list reflection_chunks){
    // set up a map exp_id:chunk_id
    std::map<std::string, size_t> chunk_lookup;
    for(size_t i = 0UL; i < exp_id_list.size(); ++i){
      chunk_lookup[exp_id_list[i]] = (size_t)chunk_id_list[i];
    }

    int n_chunks = boost::python::len(reflection_chunks);
    std::vector<dials::af::reflection_table> tables;
    for(size_t i=0UL; i < n_chunks; ++i){
      tables.push_back(boost::python::extract<dials::af::reflection_table>(reflection_chunks[i]));
    }

    // distribute reflections over chunks
    for (dials::af::reflection_table::const_iterator it = reflections.begin(); it != reflections.end(); ++it) {
      experiment_id_splitter_visitor visitor(reflections, it->first, chunk_lookup, tables);
      it->second.apply_visitor(visitor);
    }
  }

  dials::af::reflection_table isigi_dict_to_reflection_table(
                          dials::af::shared<cctbx::miller::index<> > miller_indices,
                          boost::python::dict ISIGI) {
    dials::af::reflection_table table;
    boost::python::list keys = ISIGI.keys();
    dials::af::shared<cctbx::miller::index<> > miller_index          = table.get<cctbx::miller::index<> >("miller_index");
    dials::af::shared<cctbx::miller::index<> > miller_index_original = table.get<cctbx::miller::index<> >("miller_index_original");
    dials::af::shared<double> scaled_intensity = table.get<double>("scaled_intensity");
    dials::af::shared<double> isigi            = table.get<double>("isigi");
    dials::af::shared<double> slope            = table.get<double>("slope");
    dials::af::shared<double> iobs             = table.get<double>("iobs");
    dials::af::shared<size_t> miller_id           = table.get<size_t>("miller_id");
    dials::af::shared<size_t> crystal_id          = table.get<size_t>("crystal_id");
    for (size_t i = 0; i < miller_indices.size(); i ++) {
      if (!ISIGI.has_key<cctbx::miller::index<> >(miller_indices[i])) continue;
      boost::python::list values = boost::python::extract<boost::python::list>(ISIGI[miller_indices[i]]);
      for (size_t j = 0; j < boost::python::len(values); j++) {
        miller_index.push_back(miller_indices[i]);
        boost::python::tuple item = boost::python::extract<boost::python::tuple>(values[j]);
        scaled_intensity.push_back(boost::python::extract<double>(item[0]));
        isigi.push_back(boost::python::extract<double>(item[1]));
        slope.push_back(boost::python::extract<double>(item[2]));
        miller_id.push_back(i);
        crystal_id.push_back(0);
        iobs.push_back(0);
        miller_index_original.push_back(cctbx::miller::index<int>(0,0,0));
      }
    }
    return table;
  }
}

using namespace boost::python;
namespace sx_merging{
namespace boost_python { namespace {

  void
  sx_merging_init_module() {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    def("get_hkl_chunks_cpp",&sx_merging::get_hkl_chunks_cpp);
    def("isigi_dict_to_reflection_table",&sx_merging::isigi_dict_to_reflection_table);
    def("split_reflections_by_experiment_chunks_cpp",&sx_merging::split_reflections_by_experiment_chunks_cpp);
  }

}
}} // namespace sx_merging::boost_python::<anonymous>

BOOST_PYTHON_MODULE(sx_merging_ext)
{
  sx_merging::boost_python::sx_merging_init_module();
}
