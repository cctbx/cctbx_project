#ifndef SMTBX_STRUCTURE_FACTORS_DIRECT_TABLE_BASED_H
#define SMTBX_STRUCTURE_FACTORS_DIRECT_TABLE_BASED_H

#include <smtbx/error.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <smtbx/structure_factors/direct/standard_xray.h>
#include <cctbx/miller/lookup_utils.h>
#include <fstream>

namespace smtbx { namespace structure_factors { namespace table_based {
  template <typename FloatType>
  class table_based_anisotropic
    : public direct::one_scatterer_one_h::scatterer_contribution<FloatType>
  {
    typedef direct::one_scatterer_one_h::scatterer_contribution<FloatType>
      base_type;
  public:
    typedef FloatType float_type;
    typedef std::complex<float_type> complex_type;
  private:
    af::ref_owning_shared< xray::scatterer<float_type> > scatterers;
    miller::lookup_utils::lookup_tensor<float_type> mi_lookup;
    // rows in order of original hkl index -> scatterer contribution
    std::vector<std::vector<complex_type> > data;
    //af::shared<complex_type>
    bool use_ad;
  public:
    typedef FloatType float_t;
    // Copy constructor
    table_based_anisotropic(const table_based_anisotropic &tbsc)
      :
      scatterers(tbsc.scatterers),
      mi_lookup(tbsc.mi_lookup),
      data(tbsc.data),
      use_ad(tbsc.use_ad)
    {}

    table_based_anisotropic(
      af::shared< xray::scatterer<float_type> > const &scatterers)
      :
      scatterers(scatterers)
    {}

    void read_table(const std::string &file_name,
      sgtbx::space_group const &space_group,
      bool anomalous_flag)
    {
      using namespace std;
      ifstream in_file(file_name.c_str());
      string line;
      vector<std::string> toks;
      size_t lc = 0;
      vector<size_t> sc_indices(scatterers.size());
      af::shared<cctbx::miller::index<> > miller_indices;
      while (std::getline(in_file, line)) {
        lc++;
        boost::trim(line);
        if (line.empty()) {
          break;
        }
        toks.clear();
        // is header?
        if (lc <= 3) {
          boost::split(toks, line, boost::is_any_of(":"));
          SMTBX_ASSERT(toks.size() == 2);
          if (boost::iequals(toks[0], "scatterers")) {
            std::vector<std::string> stoks;
            boost::trim(toks[1]);
            boost::split(stoks, toks[1], boost::is_any_of(" "));
            SMTBX_ASSERT(stoks.size() == scatterers.size());
            map<string, size_t> sc_map;
            for (size_t sci = 0; sci < scatterers.size(); sci++) {
              sc_map[scatterers[sci].label] = sci;
            }
            for (size_t sci = 0; sci < scatterers.size(); sci++) {
              map<string, size_t>::iterator fsci = sc_map.find(stoks[sci]);
              SMTBX_ASSERT(fsci != sc_map.end());
              sc_indices[sci] = fsci->second;
            }
          }
          else if (boost::iequals(toks[0], "AD accounted")) {
            boost::trim(toks[1]);
            use_ad = boost::iequals(toks[1], "false");
          }
        }
        // data
        else {
          boost::split(toks, line, boost::is_any_of(" "));
          SMTBX_ASSERT(toks.size() == 3 + scatterers.size());
          cctbx::miller::index<> mi(
            boost::lexical_cast<double>(toks[0]),
            boost::lexical_cast<double>(toks[1]),
            boost::lexical_cast<double>(toks[2]));
          miller_indices.push_back(mi);
          vector<complex_type> row;
          row.reserve(scatterers.size());
          for (size_t sci = 3; sci < toks.size(); sci++) {
            size_t ci = toks[sci].find(',');
            if (ci != string::npos) {
              row.push_back(
                complex_type(
                  boost::lexical_cast<double>(toks[sci].substr(0, ci)),
                  boost::lexical_cast<double>(toks[sci].substr(ci+1))));
            }
            else {
              row.push_back(
                complex_type(
                  boost::lexical_cast<double>(toks[sci])));
            }
          }
          data.push_back(row);
        }
      }
      mi_lookup = miller::lookup_utils::lookup_tensor<float_type>(
        miller_indices.const_ref(),
        space_group,
        anomalous_flag);
    }

    virtual complex_type get(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      long h_idx = mi_lookup.find_hkl(h);
      SMTBX_ASSERT(h_idx >= 0);
      complex_type rv = data[static_cast<size_t>(h_idx)][scatterer_idx];
      if (use_ad) {
        xray::scatterer<> const &sc = scatterers[scatterer_idx];
        if (sc.flags.use_fp_fdp()) {
          return complex_type(rv.real() + sc.fp, rv.imag() + sc.fdp);
        }
        else {
          return rv;
        }
      }
      else {
        return rv;
      }
    }

    virtual base_type &at_d_star_sq(
      float_type d_star_sq)
    {
      return *this;
    }

    virtual base_type *raw_fork() const {
      return new table_based_anisotropic(*this);
    }

  };

}}} // smtbx::structure_factors::table_based

#endif // GUARD
