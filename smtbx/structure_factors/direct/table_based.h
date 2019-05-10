#ifndef SMTBX_STRUCTURE_FACTORS_DIRECT_TABLE_BASED_H
#define SMTBX_STRUCTURE_FACTORS_DIRECT_TABLE_BASED_H

#include <smtbx/error.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <smtbx/structure_factors/direct/standard_xray.h>
#include <cctbx/miller/lookup_utils.h>
#include <fstream>

namespace smtbx { namespace structure_factors { namespace table_based {

  template <typename FloatType>
  class table_data {
  public:
    typedef FloatType float_type;
    typedef std::complex<float_type> complex_type;
  protected:
    // rows in order of original hkl index -> scatterer contribution
    af::shared<std::vector<complex_type> > data_;
    af::shared<cctbx::miller::index<> > miller_indices_;
    bool use_ad;
    af::shared<sgtbx::rot_mx> rot_mxs_;
  public:
    af::shared<std::vector<complex_type> > const &data() const {
      return data_;
    }

    af::shared<sgtbx::rot_mx> const &rot_mxs() const {
      return rot_mxs_;
    }

    af::shared<cctbx::miller::index<> > const &miller_indices() const {
      return miller_indices_;
    }

    bool use_AD() const {
      return use_ad;
    }
  };

  template <typename FloatType>
  class table_reader : public table_data<FloatType> {
  public:
    typedef FloatType float_type;
    typedef std::complex<float_type> complex_type;
  private:
    typedef table_data<FloatType> parent_t;

    void read(af::shared<xray::scatterer<float_type> > const &scatterers,
      const std::string &file_name)
    {
      using namespace std;
      ifstream in_file(file_name.c_str());
      string line;
      vector<std::string> toks;
      size_t lc = 0;
      vector<size_t> sc_indices(scatterers.size());
      bool header_read = false;
      while (std::getline(in_file, line)) {
        lc++;
        boost::trim(line);
        if (line.empty()) {
          break;
        }
        toks.clear();
        // is header?
        if (!header_read) {
          boost::split(toks, line, boost::is_any_of(":"));
          SMTBX_ASSERT(toks.size() == 2);
          if (boost::iequals(toks[0], "scatterers")) {
            std::vector<std::string> stoks;
            boost::trim(toks[1]);
            boost::split(stoks, toks[1], boost::is_any_of(" "));
            SMTBX_ASSERT(stoks.size() == scatterers.size());
            map<string, size_t> sc_map;
            for (size_t sci = 0; sci < scatterers.size(); sci++) {
              sc_map[boost::to_upper_copy(scatterers[sci].label)] = sci;
            }
            for (size_t sci = 0; sci < scatterers.size(); sci++) {
              boost::to_upper(stoks[sci]);
              map<string, size_t>::iterator fsci = sc_map.find(stoks[sci]);
              SMTBX_ASSERT(fsci != sc_map.end());
              sc_indices[sci] = fsci->second;
            }
          }
          else if (boost::iequals(toks[0], "AD accounted")) {
            boost::trim(toks[1]);
            parent_t::use_ad = boost::iequals(toks[1], "false");
          }
          else if (boost::iequals(toks[0], "Symm")) {
            boost::trim(toks[1]);
            vector<std::string> symms_toks;
            boost::split(symms_toks, toks[1], boost::is_any_of(";"));
            for (size_t sti = 0; sti < symms_toks.size(); sti++) {
              boost::trim(symms_toks[sti]);
              if (symms_toks[sti].empty()) {
                break;
              }
              vector<std::string> symm_toks;
              boost::split(symm_toks, symms_toks[sti], boost::is_any_of(" "));
              SMTBX_ASSERT(symm_toks.size() == 9);
              sgtbx::rot_mx rmx;
              for (size_t mei = 0; mei < 9; mei++) {
                rmx[mei] = boost::lexical_cast<int>(symm_toks[mei]);
              }
              parent_t::rot_mxs_.push_back(rmx);
            }
          }
          else if (boost::iequals(toks[0], "data")) {
            header_read = true;
          }
        }
        // data
        else {
          boost::split(toks, line, boost::is_any_of(" "));
          SMTBX_ASSERT(toks.size() == 3 + scatterers.size());
          cctbx::miller::index<> mi(
            boost::lexical_cast<int>(toks[0]),
            boost::lexical_cast<int>(toks[1]),
            boost::lexical_cast<int>(toks[2]));
          parent_t::miller_indices_.push_back(mi);
          vector<complex_type> row;
          row.resize(scatterers.size());
          for (size_t sci = 3; sci < toks.size(); sci++) {
            size_t ci = toks[sci].find(',');
            if (ci != string::npos) {
              complex_type v(
                boost::lexical_cast<float_type>(toks[sci].substr(0, ci)),
                boost::lexical_cast<float_type>(toks[sci].substr(ci + 1)));
              row[sc_indices[sci - 3]] = v;
            }
            else {
              row[sc_indices[sci - 3]] = complex_type(
                boost::lexical_cast<float_type>(toks[sci]));
            }
          }
          parent_t::data_.push_back(row);
        }
      }
    }

  public:
    table_reader(af::shared<xray::scatterer<float_type> > const &scatterers,
      const std::string &file_name)
    {
      read(scatterers, file_name);
    }
  };

  template <typename FloatType>
  class table_based_isotropic
    : public direct::one_scatterer_one_h::scatterer_contribution<FloatType>
  {
    typedef direct::one_scatterer_one_h::scatterer_contribution<FloatType>
      base_type;
  public:
    typedef FloatType float_type;
    typedef std::complex<float_type> complex_type;
  private:
    miller::lookup_utils::lookup_tensor<float_type> mi_lookup;
    // hkl x scatterer x contribution
    af::shared <std::vector<complex_type> > data;
  public:
    // Copy constructor
    table_based_isotropic(const table_based_isotropic &tbsc)
      :
      mi_lookup(tbsc.mi_lookup),
      data(tbsc.data)
    {}

    table_based_isotropic(
      af::shared< xray::scatterer<float_type> > const &scatterers,
      table_reader<FloatType> const &data_,
      sgtbx::space_group const &space_group,
      bool anomalous_flag)
      :
      data(data_.miller_indices().size())
    {
      SMTBX_ASSERT(data_.rot_mxs().size() <= 1);
      for (size_t i = 0; i < data.size(); i++) {
        data[i].resize(scatterers.size());
        for (size_t j = 0; j < scatterers.size(); j++) {
          complex_type v = data_.data()[i][j];
          if (data_.use_AD()) {
            xray::scatterer<> const &sc = scatterers[j];
            if (sc.flags.use_fp_fdp()) {
              v = complex_type(v.real() + sc.fp, v.imag() + sc.fdp);
            }
          }
          data[i][j] = v;
        }
      }
      mi_lookup = miller::lookup_utils::lookup_tensor<float_type>(
        data_.miller_indices().const_ref(),
        space_group,
        anomalous_flag);
    }

    virtual complex_type get(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      long h_idx = mi_lookup.find_hkl(h);
      SMTBX_ASSERT(h_idx >= 0);
      return data[static_cast<size_t>(h_idx)][scatterer_idx];
    }

    virtual std::vector<complex_type> const &get_full(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      SMTBX_NOT_IMPLEMENTED();
      throw 1;
    }

    virtual base_type &at_d_star_sq(
      float_type d_star_sq)
    {
      return *this;
    }

    virtual bool is_spherical() const {
      return true;
    }

    virtual base_type *raw_fork() const {
      return new table_based_isotropic(*this);
    }
  };

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
    af::shared<cctbx::miller::index<> > lookup_indices;
    miller::lookup_utils::lookup_tensor<float_type> mi_lookup;
    // hkl x scatterer x hkl*r contribution
    af::shared <af::shared<std::vector<complex_type> > > data;
  public:
    // Copy constructor
    table_based_anisotropic(const table_based_anisotropic &tbsc)
      :
      scatterers(tbsc.scatterers),
      lookup_indices(tbsc.lookup_indices),
      mi_lookup(tbsc.mi_lookup),
      data(tbsc.data)
    {}

    table_based_anisotropic(
      af::shared< xray::scatterer<float_type> > const &scatterers,
      table_reader<FloatType> const &data_,
      sgtbx::space_group const &space_group,
      bool anomalous_flag)
      :
      scatterers(scatterers)
    {
      SMTBX_ASSERT(data_.rot_mxs().size() == space_group.n_smx());
      SMTBX_ASSERT((data_.data().size() % space_group.n_smx()) == 0);

      std::vector<size_t> r_map;
      r_map.resize(space_group.n_smx());
      for (std::size_t i = 0; i < space_group.n_smx(); i++) {
        sgtbx::rot_mx const& r = space_group.smx(i).r();
        bool found = false;
        for (size_t mi = 0; mi < space_group.n_smx(); mi++) {
          if (r == data_.rot_mxs()[mi]) {
            r_map[mi] = i;
            found = true;
            break;
          }
        }
        SMTBX_ASSERT(found);
      }
      data.resize(data_.data().size() / space_group.n_smx());
      lookup_indices.resize(data.size());
      for (size_t hi = 0; hi < data.size(); hi++) {
        af::shared<std::vector<complex_type> > row(scatterers.size());
        for (size_t sci = 0; sci < scatterers.size(); sci++) {
          std::vector<complex_type> h_row;
          h_row.resize(space_group.n_smx());
          for (size_t mi = 0; mi < space_group.n_smx(); mi++) {
            const size_t r_off = data.size() * mi;
            complex_type v = data_.data()[r_off + hi][sci];
            if (data_.use_AD()) {
              xray::scatterer<> const &sc = scatterers[sci];
              if (sc.flags.use_fp_fdp()) {
                v = complex_type(v.real() + sc.fp, v.imag() + sc.fdp);
              }
            }
            h_row[r_map[mi]] = v;
          }
          row[sci] = h_row;
        }
        data[hi] = row;
        lookup_indices[hi] = data_.miller_indices()[hi];
      }

      mi_lookup = miller::lookup_utils::lookup_tensor<float_type>(
        lookup_indices.const_ref(),
        space_group,
        anomalous_flag);
    }

    virtual complex_type get(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      SMTBX_NOT_IMPLEMENTED();
      throw 1;
    }

    virtual std::vector<complex_type> const &get_full(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      long h_idx = mi_lookup.find_hkl(h);
      SMTBX_ASSERT(h_idx >= 0);
      return data[static_cast<size_t>(h_idx)][scatterer_idx];
    }


    virtual base_type &at_d_star_sq(
      float_type d_star_sq)
    {
      return *this;
    }

    virtual bool is_spherical() const {
      return false;
    }

    virtual base_type *raw_fork() const {
      return new table_based_anisotropic(*this);
    }
  };

  template <typename FloatType>
  struct builder {
    static direct::one_scatterer_one_h::scatterer_contribution<FloatType> *
      build(
        af::shared< xray::scatterer<FloatType> > const &scatterers,
        std::string const &file_name,
        sgtbx::space_group const &space_group,
        bool anomalous_flag)
    {
      table_reader<FloatType> data(scatterers, file_name);
      if (data.rot_mxs().size() <= 1) {
        return new table_based_isotropic<FloatType>(
          scatterers,
          data,
          space_group,
          anomalous_flag);
      }
      return new table_based_anisotropic<FloatType>(
        scatterers,
        data,
        space_group,
        anomalous_flag);
    }

  };

}}} // smtbx::structure_factors::table_based

#endif // GUARD
