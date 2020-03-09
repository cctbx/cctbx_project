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
  void foo2(){
    std::cout <<"HELLO merging foo2"<< std::endl;
  }

  void get_hkl_chunks_cpp(dials::af::reflection_table reflections,
                          const scitbx::af::shared<miller_index_t>& hkl_list,
                          const scitbx::af::shared<int>& chunk_id_list,
                          boost::python::list hkl_chunks_cpp){
    // set up a map hkl:chunk_id
    std::map<miller_index_t, size_t> chunk_lookup;
    for(size_t i = 0UL; i < hkl_list.size(); ++i){
      chunk_lookup[hkl_list[i]] = (size_t)chunk_id_list[i];
    }

    SCITBX_ASSERT(reflections.contains("miller_index_asymmetric"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.value"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.variance"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.value.unmodified"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.variance.unmodified"));
    SCITBX_ASSERT(reflections.contains("exp_id"));

    scitbx::af::ref<miller_index_t> miller_index    = reflections["miller_index_asymmetric"];
    scitbx::af::ref<double> intensity               = reflections["intensity.sum.value"];
    scitbx::af::ref<double> variance                = reflections["intensity.sum.variance"];
    scitbx::af::ref<double> intensity_unmodified    = reflections["intensity.sum.value.unmodified"];
    scitbx::af::ref<double> variance_unmodified     = reflections["intensity.sum.variance.unmodified"];
    scitbx::af::ref<std::string> experiment_id      = reflections["exp_id"];

    miller_index_t* mi_ptr          = miller_index.begin();
    double* intensity_ptr           = intensity.begin();
    double* variance_ptr            = variance.begin();
    double* intensity_unmodified_ptr           = intensity_unmodified.begin();
    double* variance_unmodified_ptr            = variance_unmodified.begin();
    std::string* experiment_id_ptr  = experiment_id.begin();

    int n_chunks = boost::python::len(hkl_chunks_cpp);
    std::vector<dials::af::reflection_table> tables;
    for(size_t i=0UL; i < n_chunks; ++i){
      tables.push_back(boost::python::extract<dials::af::reflection_table>(hkl_chunks_cpp[i]));
      SCITBX_ASSERT(tables.back().contains("miller_index_asymmetric"));
    }

    // cache all columns for all hkl chunks
    std::vector<scitbx::af::shared<miller_index_t> >  mi_cols;
    std::vector<scitbx::af::shared<double> >          intensity_cols;
    std::vector<scitbx::af::shared<double> >          variance_cols;
    std::vector<scitbx::af::shared<double> >          intensity_unmodified_cols;
    std::vector<scitbx::af::shared<double> >          variance_unmodified_cols;
    std::vector<scitbx::af::shared<std::string> >     experiment_id_cols;

    for(size_t i=0UL; i < n_chunks; ++i){
      mi_cols.push_back(tables[i]["miller_index_asymmetric"]);
      intensity_cols.push_back(tables[i]["intensity.sum.value"]);
      variance_cols.push_back(tables[i]["intensity.sum.variance"]);
      intensity_unmodified_cols.push_back(tables[i]["intensity.sum.value.unmodified"]);
      variance_unmodified_cols.push_back(tables[i]["intensity.sum.variance.unmodified"]);
      experiment_id_cols.push_back(tables[i]["exp_id"]);
    }

    // distribute reflections over chunks
    for(size_t i=0UL; i < reflections.size(); ++i){
      miller_index_t  hkl       = *mi_ptr++;
      double          intensity = *intensity_ptr++;
      double          variance  = *variance_ptr++;
      double          intensity_unmodified = *intensity_unmodified_ptr++;
      double          variance_unmodified  = *variance_unmodified_ptr++;
      std::string     experiment_id = *experiment_id_ptr++;

      if( 0 != chunk_lookup.count(hkl) )
      {
          size_t chunk_id  = chunk_lookup[hkl];
          mi_cols[chunk_id].push_back(hkl);
          intensity_cols[chunk_id].push_back(intensity);
          variance_cols[chunk_id].push_back(variance);
          intensity_unmodified_cols[chunk_id].push_back(intensity_unmodified);
          variance_unmodified_cols[chunk_id].push_back(variance_unmodified);
          experiment_id_cols[chunk_id].push_back(experiment_id);
       }
    }
  }

  void split_reflections_by_experiment_chunks_cpp(dials::af::reflection_table reflections,
                                                  const scitbx::af::shared<std::string>& exp_id_list,
                                                  const scitbx::af::shared<int>& chunk_id_list,
                                                  boost::python::list reflection_chunks){
    // set up a map exp_id:chunk_id
    std::map<std::string, size_t> chunk_lookup;
    for(size_t i = 0UL; i < exp_id_list.size(); ++i){
      chunk_lookup[exp_id_list[i]] = (size_t)chunk_id_list[i];
    }

    SCITBX_ASSERT(reflections.contains("miller_index"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.value"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.variance"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.value.unmodified"));
    SCITBX_ASSERT(reflections.contains("intensity.sum.variance.unmodified"));
    SCITBX_ASSERT(reflections.contains("exp_id"));
    SCITBX_ASSERT(reflections.contains("s1"));

    scitbx::af::ref<miller_index_t> miller_index    = reflections["miller_index"];
    scitbx::af::ref<double> intensity               = reflections["intensity.sum.value"];
    scitbx::af::ref<double> variance                = reflections["intensity.sum.variance"];
    scitbx::af::ref<double> intensity_unmodified    = reflections["intensity.sum.value.unmodified"];
    scitbx::af::ref<double> variance_unmodified     = reflections["intensity.sum.variance.unmodified"];
    scitbx::af::ref<std::string> experiment_id      = reflections["exp_id"];
    scitbx::af::ref<scitbx::vec3<double> > s1       = reflections["s1"];

    miller_index_t* mi_ptr          = miller_index.begin();
    double* intensity_ptr           = intensity.begin();
    double* variance_ptr            = variance.begin();
    double* intensity_unmodified_ptr = intensity_unmodified.begin();
    double* variance_unmodified_ptr  = variance_unmodified.begin();
    std::string* experiment_id_ptr  = experiment_id.begin();
    scitbx::vec3<double>* s1_ptr    = s1.begin();

    int n_chunks = boost::python::len(reflection_chunks);
    std::vector<dials::af::reflection_table> tables;
    for(size_t i=0UL; i < n_chunks; ++i){
      tables.push_back(boost::python::extract<dials::af::reflection_table>(reflection_chunks[i]));
      SCITBX_ASSERT(tables.back().contains("miller_index"));
    }

    // cache all columns for all hkl chunks
    std::vector<scitbx::af::shared<miller_index_t> >  mi_cols;
    std::vector<scitbx::af::shared<double> >          intensity_cols;
    std::vector<scitbx::af::shared<double> >          variance_cols;
    std::vector<scitbx::af::shared<double> >          intensity_unmodified_cols;
    std::vector<scitbx::af::shared<double> >          variance_unmodified_cols;
    std::vector<scitbx::af::shared<std::string> >     experiment_id_cols;
    std::vector<scitbx::af::shared<scitbx::vec3<double> > >     s1_cols;

    for(size_t i=0UL; i < n_chunks; ++i){
      mi_cols.push_back(tables[i]["miller_index"]);
      intensity_cols.push_back(tables[i]["intensity.sum.value"]);
      variance_cols.push_back(tables[i]["intensity.sum.variance"]);
      intensity_unmodified_cols.push_back(tables[i]["intensity.sum.value.unmodified"]);
      variance_unmodified_cols.push_back(tables[i]["intensity.sum.variance.unmodified"]);
      experiment_id_cols.push_back(tables[i]["exp_id"]);
      s1_cols.push_back(tables[i]["s1"]);
    }

    // distribute reflections over chunks
    for(size_t i=0UL; i < reflections.size(); ++i){
      miller_index_t  hkl             = *mi_ptr++;
      double          intensity       = *intensity_ptr++;
      double          variance        = *variance_ptr++;
      double          intensity_unmodified = *intensity_unmodified_ptr++;
      double          variance_unmodified  = *variance_unmodified_ptr++;
      std::string     experiment_id   = *experiment_id_ptr++;
      scitbx::vec3<double> s1         = *s1_ptr++;

      if( 0 != chunk_lookup.count(experiment_id) )
      {
          size_t chunk_id = chunk_lookup[experiment_id];
          mi_cols[chunk_id].push_back(hkl);
          intensity_cols[chunk_id].push_back(intensity);
          variance_cols[chunk_id].push_back(variance);
          intensity_unmodified_cols[chunk_id].push_back(intensity_unmodified);
          variance_unmodified_cols[chunk_id].push_back(variance_unmodified);
          experiment_id_cols[chunk_id].push_back(experiment_id);
          s1_cols[chunk_id].push_back(s1);
       }
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

    def("foo2",&sx_merging::foo2);
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
