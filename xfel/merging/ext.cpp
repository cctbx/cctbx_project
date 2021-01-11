#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/tuple.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/miller.h>
#include <cctbx/miller/match_indices.h>
#include <cctbx/sgtbx/change_of_basis_op.h>
#include <dials/array_family/reflection_table.h>
#include <algorithm>
#include <tuple>
#include <string>
#include <scitbx/math/linear_correlation.h>

typedef cctbx::miller::index<> miller_index_t;
namespace af = scitbx::af;
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
    dials::af::reflection_table::key_type key;
    const std::map<miller_index_t, size_t> chunk_lookup;
    const std::vector<dials::af::reflection_table> tables;
    dials::af::shared<miller_index_t> src_hkls;

    hkl_splitter_visitor(dials::af::reflection_table &src_,
                                   dials::af::reflection_table::key_type key_,
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
    dials::af::reflection_table::key_type key;
    const std::map<std::string, size_t> chunk_lookup;
    const std::vector<dials::af::reflection_table> tables;
    dials::af::shared<std::string> expt_ids;

    experiment_id_splitter_visitor(dials::af::reflection_table &src_,
                                   dials::af::reflection_table::key_type key_,
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

  struct compute_rij_wij_detail {
    typedef cctbx::sgtbx::change_of_basis_op cbop_t;

    compute_rij_wij_detail(
        const af::shared<int> &lower_i,
        const af::shared<int> &upper_i,
        const af::shared<double> &data,
        const int &n_obs_min
        ) :
      lower_i_(lower_i),
      upper_i_(upper_i),
      data_(data),
      weights_("count"),
      n_obs_min_(n_obs_min) {
    }

    af::shared<cbop_t> cb_ops;
    af::shared<af::shared<cctbx::miller::index<> > > indices;
    af::shared<int> lower_i_;
    af::shared<int> upper_i_;
    af::shared<double> data_;
    int n_sym_ops;
    int n_obs_min_;
    std::string weights_;

    void set_indices(cbop_t cb_op,
                     af::shared<cctbx::miller::index<> > index_set){
      cb_ops.push_back(cb_op);
      indices.push_back(index_set);
      n_sym_ops = cb_ops.size();
    }

    /*
    void store_data(const af::shared<double> &data) {
      data_ = data; // consider if we need to copy it instead?
    }
    */



    boost::python::tuple
    compute_one_row(
        const int &n_lattices,
        const int &i_row) const
    {
      // port of _compute_rij_matrix_one_row_block from test.py
      typedef std::tuple<long, long, std::string> cache_key_t;
      typedef std::map<cache_key_t, std::tuple<double, int> > rij_cache_t;


      rij_cache_t rij_cache;
      cache_key_t cache_key;
      cctbx::sgtbx::rt_mx k_inv_kk;

      af::shared<int> rij_row, rij_col;
      af::shared<double> rij_data;
      af::shared<int> wij_row, wij_col;
      af::shared<double> wij_data;

      double corr_coeff;
      int n_obs;

      int i_lower = lower_i_[i_row];
      int i_upper = upper_i_[i_row];

      af::shared<double> intensities_i, intensities_j;

      af::shared<cctbx::miller::match_indices> matchers;
      for (int i=0; i<n_sym_ops; ++i) {
        af::shared<cctbx::miller::index<> > indices_i;
        for (int ii=i_lower; ii<i_upper; ++ii) {
          indices_i.push_back(indices[i][ii]);
        }
        matchers.push_back(cctbx::miller::match_indices(indices_i));
      }

      for (int j_col=i_row; j_col<n_lattices; ++j_col) {
        int j_lower = lower_i_[j_col], j_upper = upper_i_[j_col];

        af::shared<cctbx::miller::index<> > indices_j;

        for (int k=0; k<n_sym_ops; ++k) {
          cctbx::miller::match_indices &matcher = matchers[k];
          const cbop_t& cb_op_k = cb_ops[k];


          for (int kk=0; kk<n_sym_ops; ++kk) {
            const cbop_t& cb_op_kk = cb_ops[kk];
            k_inv_kk = cb_op_k.c_inv() * cb_op_kk.c();
            std::string k_inv_kk_str = k_inv_kk.as_xyz();
            cache_key = std::make_tuple(i_row, j_col, k_inv_kk_str);
            rij_cache_t::const_iterator found = rij_cache.find(cache_key);
            if (found != rij_cache.end()) {
              std::tuple<double, int> hit = found->second;
              corr_coeff = std::get<0>(hit);
              n_obs = std::get<1>(hit);
            }
            else {

              indices_j.clear();
              for (int i=j_lower; i<j_upper; ++i) {
                indices_j.push_back(indices[kk][i]);
              }

              matcher.match_cached(indices_j);

              af::shared<af::tiny<std::size_t, 2> > pairs = matcher.pairs();

              intensities_i.clear();
              intensities_j.clear();
              for (int i=0; i<pairs.size(); ++i) {
                intensities_i.push_back(data_[i_lower+pairs[i][0]]);
                intensities_j.push_back(data_[j_lower+pairs[i][1]]);
              }


              if (i_row==j_col && k==kk) continue;

              scitbx::math::linear_correlation<> corr(
                  af::make_const_ref(intensities_i),
                  af::make_const_ref(intensities_j)
                  );

              if (corr.is_well_defined()) {
                corr_coeff = corr.coefficient();
                n_obs = corr.n();
              }
              else {
                corr_coeff = -1.;
                n_obs = -1;
              }


              rij_cache[cache_key] = std::make_tuple(corr_coeff, n_obs);

            }

            if (n_obs==-1 || corr_coeff == -1. || n_obs < n_obs_min_)
              continue;
            else {
              int ik = i_row + (n_lattices * k);
              int jk = j_col + (n_lattices * kk);
              rij_row.push_back(ik);
              rij_col.push_back(jk);
              rij_data.push_back(corr_coeff);
              if (i_row != j_col) {
                rij_row.push_back(jk);
                rij_col.push_back(ik);
                rij_data.push_back(corr_coeff);
              }
              if (weights_ == "count") {
                wij_row.push_back(ik);
                wij_col.push_back(jk);
                wij_data.push_back(2*n_obs); //The factor of 2 is to match a bug
                                             //in the original implementation.
                if (i_row != j_col) {
                  wij_row.push_back(jk);
                  wij_col.push_back(ik);
                  wij_data.push_back(2*n_obs); //The factor of 2 is to match a bug
                                               //in the original implementation.
                }
              }
            }



          }
        }
      }
      return boost::python::make_tuple(rij_row, rij_col, rij_data, wij_row, wij_col, wij_data);
    }

    std::string foo(const int& n_lattices,
                  af::shared<int> lower_i,
                  af::shared<int> upper_i,
                  af::shared<double> data)
    {
    int NN = n_lattices * n_sym_ops;
    af::flex_double wij(af::flex_grid<>(NN,NN));
    double* wijbegin = wij.begin();
    af::flex_double rij(af::flex_grid<>(NN,NN));
    double* rijbegin = rij.begin();
    SCITBX_EXAMINE(n_sym_ops);
    for (int i = 0; i < n_lattices; i++){
        //def _compute_rij_matrix_one_row_block(i):
          //  print("compute_rij_matrix_one_row_block",i,"of",n_lattices)
          //  rij_cache = {}

          //  n_sym_ops = len(self._sym_ops)
          //  NN = n_lattices * n_sym_ops

          //  rij_row = []
          //  rij_col = []
          //  rij_data = []
          //  if self._weights is not None:
          //      wij_row = []
          //      wij_col = []
          //      wij_data = []
          //  else:
          //      wij = None
      SCITBX_EXAMINE(i);
      int i_lower = lower_i[i]; int i_upper = upper_i[i];
      double* ptr_intensities_i = &(data[i_lower]);
      for (int j = 0; j < n_lattices; j++){
        int j_lower = lower_i[j]; int j_upper = upper_i[j];
        double* ptr_intensities_j = &(data[j_lower]);

        typedef af::shared<cctbx::miller::index<> > indices_t;
        for (int k = 0; k < n_sym_ops; ++k){
          cbop_t& cb_op_k = cb_ops[k];
          indices_t indices_i;
          for (int idx = i_lower; idx < i_upper; idx++){
            indices_i.push_back(indices[k][idx]);
          } //only way I know of making a slice; giving the indices of one lattice reindexed by k

          for (int kk = 0; kk < n_sym_ops; ++kk){
            cbop_t& cb_op_kk = cb_ops[kk];
            if (i == j && k == kk){ continue; } //# don't include correlation of dataset with itself

            int ik = i + (n_lattices * k); // row
            int jk = j + (n_lattices * kk); // column

/*                        key = (i, j, str(cb_op_k.inverse() * cb_op_kk))
                        if use_cache and key in rij_cache:
                            cc, n = rij_cache[key]
                        else:
*/
            //this seems to be a big time sink and do we need to separately calculate indices_j,
            //if we store all instances of n_sym_ops x n_sym_ops of indices_i

            indices_t indices_j;
            for (int idx = j_lower; idx < j_upper; idx++){
              indices_j.push_back(indices[kk][idx]);
            }
/*

                            matches = miller.match_indices(indices_i, indices_j)
                            pairs = matches.pairs()
                            isel_i = pairs.column(0)
                            isel_j = pairs.column(1)
                            isel_i = isel_i.select(
                                self._patterson_group.epsilon(indices_i.select(isel_i))
                                == 1
                            )
                            isel_j = isel_j.select(
                                self._patterson_group.epsilon(indices_j.select(isel_j))
                                == 1
                            )
                            corr = flex.linear_correlation(
                                intensities_i.select(isel_i),
                                intensities_j.select(isel_j),
                            )

                            if corr.is_well_defined():
                                cc = corr.coefficient()
                                n = corr.n()
                                rij_cache[key] = (cc, n)
                            else:
                                cc = None
                                n = None
                        #print("n", n, "cc", cc, self._min_pairs)

                        if (
                            n is None
                            or cc is None
                            or (self._min_pairs is not None and n < self._min_pairs)
                        ):
                            continue

                        if self._weights == "count":
                            wij_row.extend([ik, jk])
                            wij_col.extend([jk, ik])
                            wij_data.extend([n, n])
                        elif self._weights == "standard_error":
                            assert n > 2
                            # http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
                            se = math.sqrt((1 - cc ** 2) / (n - 2))
                            wij = 1 / se
                            wij_row.extend([ik, jk])
                            wij_col.extend([jk, ik])
                            wij_data.extend([wij, wij])

                        rij_row.append(ik)
                        rij_col.append(jk)
                        rij_data.append(cc)
*/
          } //over kk
/*
            rij = sparse.coo_matrix((rij_data, (rij_row, rij_col)), shape=(NN, NN))
            if self._weights is not None:
                wij = sparse.coo_matrix((wij_data, (wij_row, wij_col)), shape=(NN, NN))

            return rij, wij
*/
        } // over k
      } //over j
    } // over i
    return ("Hello foo\n");
  } // foo
  };
}

using namespace boost::python;
namespace sx_merging{
namespace boost_python { namespace {

  void
  sx_merging_init_module() {
    using namespace boost::python;

    def("get_hkl_chunks_cpp",&sx_merging::get_hkl_chunks_cpp);
    def("isigi_dict_to_reflection_table",&sx_merging::isigi_dict_to_reflection_table);
    def("split_reflections_by_experiment_chunks_cpp",&sx_merging::split_reflections_by_experiment_chunks_cpp);

    class_<compute_rij_wij_detail>("compute_rij_wij_detail", no_init)
      .def(init<
        const af::shared<int>,
        const af::shared<int>,
        const af::shared<double>,
        const int >())
      .def("set_indices",&sx_merging::compute_rij_wij_detail::set_indices)
      .def("foo",&sx_merging::compute_rij_wij_detail::foo)
      .def("compute_one_row",&sx_merging::compute_rij_wij_detail::compute_one_row)
    ;


  }

}
}} // namespace sx_merging::boost_python::<anonymous>

BOOST_PYTHON_MODULE(sx_merging_ext)
{
  sx_merging::boost_python::sx_merging_init_module();
}
