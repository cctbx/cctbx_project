#ifndef SMTBX_STRUCTURE_FACTORS_DIRECT_TABLE_BASED_H
#define SMTBX_STRUCTURE_FACTORS_DIRECT_TABLE_BASED_H

#include <smtbx/error.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <smtbx/structure_factors/direct/standard_xray.h>
#include <cctbx/miller/lookup_utils.h>
#include <cctbx/xray/scatterer_lookup.h>
#include <fstream>
#include <stdint.h>

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
    bool expanded;
    // one per scatterer of the structure: did the table have a column for it?
    // A table need not cover the whole structure, and the ones it misses are
    // given a spherical form factor rather than nothing.
    af::shared<bool> covered_;
  public:
    table_data()
    : expanded(false),
      use_ad(false)
    {}
    af::shared<std::vector<complex_type> > const &data() const {
      return data_;
    }

    af::shared<bool> const &covered() const {
      return covered_;
    }

    //! The scatterers the table had no column for, in order.
    af::shared<std::size_t> not_covered() const {
      af::shared<std::size_t> result;
      for (std::size_t i = 0; i < covered_.size(); i++) {
        if (!covered_[i]) {
          result.push_back(i);
        }
      }
      return result;
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

    bool is_expanded() const {
      return expanded;
    }
  };

  template <typename FloatType>
  class table_reader : public table_data<FloatType> {
  public:
    typedef FloatType float_type;
    typedef std::complex<float_type> complex_type;
  private:
    typedef table_data<FloatType> parent_t;
    const uctbx::unit_cell& u_cell;

    void read(af::shared<xray::scatterer<float_type> > const& scatterers,
      const std::string &file_name)
    {
      using namespace std;
    //   typedef cctbx::xray::scatterer_id_5<float_type, fractional<float_type>, 16> scatterer_id_t;
    typedef cctbx::xray::scatterer_id_big<float_type, fractional<float_type> > scatterer_id_t;
      ifstream in_file(file_name.c_str());
      string line;
      vector<std::string> toks;
      size_t lc = 0;
      // sized by the table's own column count once the header names them
      vector<size_t> sc_indices;
      parent_t::covered_ = af::shared<bool>(scatterers.size(), false);
      bool header_read = false, ids_read = false;
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
          if (toks.size() < 2) {
            continue;
          }
          if (boost::iequals(toks[0], "scatterers") && !ids_read) {
            std::vector<std::string> stoks;
            boost::trim(toks[1]);
            boost::split(stoks, toks[1], boost::is_any_of(" "));
            map<string, size_t> sc_map;
            for (size_t sci = 0; sci < scatterers.size(); sci++) {
              sc_map[boost::to_upper_copy(scatterers[sci].label)] = sci;
            }
            sc_indices.assign(stoks.size(), ~0);
            for (size_t sci = 0; sci < stoks.size(); sci++) {
              boost::to_upper(stoks[sci]);
              map<string, size_t>::iterator fsci = sc_map.find(stoks[sci]);
              if (fsci != sc_map.end()) {
                sc_indices[sci] = fsci->second;
                parent_t::covered_[fsci->second] = true;
              }
            }
          }
          else if (boost::iequals(toks[0], "scatterer_ids")) {
            std::vector<std::string> stoks;
            boost::trim(toks[1]);
            boost::split(stoks, toks[1], boost::is_any_of(" "));
            int nr_scat = scatterers.size();
            af::shared<int> data(nr_scat);
            for (size_t sci = 0; sci < nr_scat; sci++) {
                data[sci] = scatterers[sci].get_part();
            }

            cctbx::xray::scatterer_ID_lookup<FloatType> scatter_lookup(u_cell, scatterers, data);
            sc_indices.assign(stoks.size(), ~0);
            for (size_t sci = 0; sci < stoks.size(); sci++) {
              scatterer_id_t sc_id(stoks[sci]);
              size_t idx = scatter_lookup.find_index(sc_id);
              sc_indices[sci] = idx;
              if (idx != ~0) {
                parent_t::covered_[idx] = true;
              }
            }
            ids_read = true;
          }
          else if (boost::iequals(toks[0], "AD")) {
            boost::trim(toks[1]);
            parent_t::use_ad = boost::iequals(toks[1], "false");
          }
          else if (boost::iequals(toks[0], "Symm")) {
            boost::trim(toks[1]);
            if (boost::iequals(toks[1], "expanded")) {
              parent_t::expanded = true;
            }
            else {
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
          }
          else if (boost::iequals(toks[0], "data")) {
            header_read = true;
          }
        }
        // data
        else {
          boost::split(toks, line, boost::is_any_of(" "));
          // as many columns as the table names, which need not be as many
          // atoms as the structure has
          SMTBX_ASSERT(toks.size() == 3 + sc_indices.size());
          cctbx::miller::index<> mi(
            atoi(toks[0].c_str()),
            atoi(toks[1].c_str()),
            atoi(toks[2].c_str()));
          parent_t::miller_indices_.push_back(mi);
          vector<complex_type> row;
          // an atom the table misses keeps a zero and is served elsewhere
          row.assign(scatterers.size(), complex_type(0));
#pragma omp parallel for
          for (int sci = 3; sci < toks.size(); sci++) {
            if (sc_indices[sci - 3] == ~0) {
              continue;
            }
            size_t ci = toks[sci].find(',');
            if (ci != string::npos) {
              complex_type v(
                  atof(toks[sci].substr(0, ci).c_str()),
                  atof(toks[sci].substr(ci + 1).c_str()));
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

    void read_binary(af::shared<xray::scatterer<float_type> > const &scatterers,
      const std::string &file_name)
    {
      using namespace std;
    //   typedef cctbx::xray::scatterer_id_5<float_type, fractional<float_type>, 16> scatterer_id_t;
      typedef cctbx::xray::scatterer_id_big<float_type, fractional<float_type> > scatterer_id_t;
      ifstream tsc_file(file_name.c_str(), ios::binary);
      const size_t charsize = sizeof(char);
      const size_t intsize = sizeof(int);
      const size_t uint64size = sizeof(uint64_t);
      const size_t complex_doublesize = sizeof(complex<double>);
      const size_t complex_type_size = sizeof(complex_type);
      //If the size is not according to double type the binary will not be readable
      SMTBX_ASSERT(complex_doublesize == complex_type_size);
      int head = 0;
      tsc_file.read((char*)&head, intsize);
      const int nr_scat = scatterers.size();
      string header_str;
      if (head != 0) {
        vector<char> header(head);
        tsc_file.read(&header[0], head * charsize);
        header_str = string(header.begin(), header.end());
      }
      //read scatterer labels or ids and map onto scattterers list
      int sc_len = 0;
      tsc_file.read((char*)&sc_len, intsize);
      /* The table's column count and the structure's atom count are two
      different numbers. They usually agree, but a table need not cover the
      whole structure, so the file is read with its own width and each column
      is placed where it belongs -- or dropped, if the structure has no atom
      for it. Columns the structure has no atom for are read and discarded;
      atoms with no column are recorded in covered_ and later given a
      spherical form factor.
      */
      // sc_len is a count of ids in the one case and a count of bytes of
      // label text in the other, so the column count is settled per branch
      size_t n_columns = 0;
      vector<size_t> sc_indices;
      parent_t::covered_ = af::shared<bool>(nr_scat, false);
      if (boost::icontains(header_str, "SCATTERER_IDS")) {
        n_columns = sc_len;
        sc_indices.assign(n_columns, ~0);
        af::shared<int> data(nr_scat);
        for (size_t sci = 0; sci < nr_scat; sci++) {
            data[sci] = scatterers[sci].get_part();
        }

        cctbx::xray::scatterer_ID_lookup<FloatType> scatter_lookup(u_cell, scatterers, data);
        for (size_t sci = 0; sci < n_columns; sci++) {
          scatterer_id_t sc_id(tsc_file); //Read the raw bytes and convert to scatterer_id_t
          size_t idx = scatter_lookup.find_index(sc_id);
          sc_indices[sci] = idx;
          if (idx != ~0) {
            parent_t::covered_[idx] = true;
          }
        }
      }
      else {
        vector<char> scat_line(sc_len);
        tsc_file.read((char*)scat_line.data(), sc_len * charsize);
        string scat_str(scat_line.begin(),scat_line.end());
        vector<string> toks;
        boost::split(toks, scat_str, boost::is_any_of(" "));
        map<string, size_t> sc_map;
        for (size_t sci = 0; sci < nr_scat; sci++) {
          sc_map[boost::to_upper_copy(scatterers[sci].label)] = sci;
        }
        n_columns = toks.size();
        sc_indices.assign(n_columns, ~0);
        for (size_t sci = 0; sci < n_columns; sci++) {
          boost::to_upper(toks[sci]);
          map<string, size_t>::iterator fsci = sc_map.find(toks[sci]);
          if (fsci != sc_map.end()) {
            sc_indices[sci] = fsci->second;
            parent_t::covered_[fsci->second] = true;
          }
        }
      }
      //binary tsc files will only be written in expanded mode
      parent_t::expanded = true;
      //read number of indices in tscb file
      int nr_hkl[1] = { 0 };
      tsc_file.read((char*)&nr_hkl, intsize);
      //read indices and scattering factors row by row
      int index[3] = { 0,0,0 };
      // a row of the file is as wide as the table, a row of the model as wide
      // as the structure, and the two need not be the same
      vector<complex_type> buffer(n_columns);
      // one entry per reflection, not one per element of every reflection:
      // over-reserving here scales with the whole table
      parent_t::miller_indices_.reserve(nr_hkl[0]);
      // Size the table up front and fill it in place: pushing each row back
      // copies it, which is another pass over the whole table.
      parent_t::data_.resize(nr_hkl[0]);
      for (int run = 0; run < *nr_hkl; run++) {
        tsc_file.read((char*)&index, 3*intsize);
        cctbx::miller::index<> mi(index[0], index[1], index[2]);
        parent_t::miller_indices_.push_back(mi);
        // A read per scatterer spends more in stream bookkeeping than it
        // moves, and a table holds one row per reflection: take the row in one
        // read and put the scatterers in their places afterwards.
        tsc_file.read((char*)&buffer[0], n_columns * complex_doublesize);
        vector<complex_type> &target = parent_t::data_[run];
        // an atom the table misses keeps a zero here and is served elsewhere
        target.assign(nr_scat, complex_type(0));
        for (size_t i = 0; i < n_columns; i++) {
          if (sc_indices[i] != ~0) {
            target[sc_indices[i]] = buffer[i];
          }
        }
      }
      tsc_file.close();
      SMTBX_ASSERT(!tsc_file.bad());
    }

  public:
    table_reader(const uctbx::unit_cell& u_cell,
      af::shared<xray::scatterer<float_type> > const &scatterers,
      const std::string &file_name)
      : u_cell(u_cell)
    {
      if(file_name.find(".tscb")!=std::string::npos){
        read_binary(scatterers, file_name);
      }
      else{
        read(scatterers, file_name);
      }
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
      // shared, not copied: af::shared is reference counted
      data(data_.data())
    {
      SMTBX_ASSERT(data_.rot_mxs().size() <= 1);
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
    miller::lookup_utils::lookup_tensor<float_type> mi_lookup;
    // hkl x scatterer x hkl*r contribution
    af::shared <af::shared<std::vector<complex_type> > > data;
  public:
    // Copy constructor
    table_based_anisotropic(const table_based_anisotropic &tbsc)
      :
      mi_lookup(tbsc.mi_lookup),
      data(tbsc.data)
    {}

    table_based_anisotropic(
      af::shared< xray::scatterer<float_type> > const &scatterers,
      table_reader<FloatType> const &data_,
      sgtbx::space_group const &space_group,
      bool anomalous_flag)
    {
      SMTBX_ASSERT(data_.rot_mxs().size() == space_group.n_smx());
      SMTBX_ASSERT((data_.data().size() % space_group.n_smx()) == 0);

      std::vector<size_t> r_map;
      r_map.resize(space_group.n_smx());
      for (size_t i = 0; i < space_group.n_smx(); i++) {
        sgtbx::rot_mx const& r = data_.rot_mxs()[i];
        bool found = false;
        for (size_t mi = 0; mi < space_group.n_smx(); mi++) {
          if (r == space_group.smx(mi).r()) {
            r_map[i] = mi;
            found = true;
            break;
          }
        }
        SMTBX_ASSERT(found);
      }
      data.resize(data_.data().size() / space_group.n_smx());
      af::shared<cctbx::miller::index<> > lookup_indices(data.size());
      for (size_t hi = 0; hi < data.size(); hi++) {
        af::shared<std::vector<complex_type> > row(scatterers.size());
        for (size_t sci = 0; sci < scatterers.size(); sci++) {
          std::vector<complex_type> h_row(space_group.n_smx());
          for (size_t mi = 0; mi < space_group.n_smx(); mi++) {
            const size_t r_off = data.size() * mi;
            miller::index<> h =
              data_.miller_indices()[hi] * space_group.smx(r_map[mi]).r();
            SMTBX_ASSERT(h == data_.miller_indices()[r_off + hi]);
            h_row[r_map[mi]] = data_.data()[r_off + hi][sci];
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
  class lookup_based_anisotropic
    : public direct::one_scatterer_one_h::scatterer_contribution<FloatType>
  {
    typedef direct::one_scatterer_one_h::scatterer_contribution<FloatType>
      base_type;
  public:
    typedef FloatType float_type;
    typedef std::complex<float_type> complex_type;
  private:
    typedef std::map<cctbx::miller::index<>, std::size_t,
      cctbx::miller::fast_less_than<> > lookup_t;
    lookup_t mi_lookup;
    sgtbx::space_group const &space_group;
    af::shared<std::vector<complex_type> > data;
    bool anomalous_flag;
    mutable std::vector<complex_type> tmp;
  public:
    // Copy constructor
    lookup_based_anisotropic(const lookup_based_anisotropic &lbsc)
      :
      mi_lookup(lbsc.mi_lookup),
      space_group(lbsc.space_group),
      data(lbsc.data),
      anomalous_flag(lbsc.anomalous_flag),
      tmp(lbsc.tmp.size())
    {}

    lookup_based_anisotropic(
      af::shared< xray::scatterer<float_type> > const &scatterers,
      table_reader<FloatType> const &data_,
      sgtbx::space_group const &space_group,
      bool anomalous_flag)
      :
      space_group(space_group),
      // af::shared is reference counted, so this shares the table the reader
      // built rather than copying it. Copying it row by row doubled both the
      // time and the peak memory of loading a table.
      data(data_.data()),
      anomalous_flag(anomalous_flag),
      tmp(space_group.n_smx())
    {
      SMTBX_ASSERT(data_.rot_mxs().size() <= 1);
      SMTBX_ASSERT(data_.is_expanded());
      const af::shared<cctbx::miller::index<> > &indices = data_.miller_indices();
      for (size_t i = 0; i < data.size(); i++) {
        mi_lookup[indices[i]] = i;
      }
    }
    // for testing
    lookup_based_anisotropic(
      uctbx::unit_cell const &unit_cell,
      sgtbx::space_group const &space_group,
      af::shared<xray::scatterer<float_type> > const &scatterers,
      direct::one_scatterer_one_h::isotropic_scatterer_contribution<FloatType> &isc,
      af::shared<cctbx::miller::index<> > const &indices)
      :
      space_group(space_group),
      data(indices.size()*space_group.n_smx()),
      tmp(space_group.n_smx())
    {
      for (size_t i = 0; i < indices.size(); i++) {
        float_type d_star_sq = unit_cell.d_star_sq(indices[i]);
        direct::one_scatterer_one_h::scatterer_contribution<FloatType> const& sc =
          isc.at_d_star_sq(d_star_sq);
        for (size_t j = 0; j < space_group.n_smx(); j++) {
          size_t d_off = j * indices.size() + i;
          miller::index<> h = indices[i] * space_group.smx(j).r();
          mi_lookup[h] = d_off;
          data[d_off].resize(scatterers.size());
          for (size_t k = 0; k < scatterers.size(); k++) {
            complex_type v = sc.get(k, h);
            if (scatterers[k].flags.use_fp_fdp()) { // revert of applied...
              v = complex_type(v.real() - scatterers[k].fp, 0);
            }
            data[d_off][k] = v;
          }
        }
      }
    }

    virtual complex_type get(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      SMTBX_NOT_IMPLEMENTED();
      throw 1;
    }

    virtual std::vector<complex_type> const &get_full(std::size_t scatterer_idx,
      miller::index<> const &h_) const
    {
      for (std::size_t i = 0; i < space_group.n_smx(); i++) {
        miller::index<> h = h_ * space_group.smx(i).r();
        lookup_t::const_iterator l = mi_lookup.find(h);
        if (l == mi_lookup.end() && !anomalous_flag) {
          h *= -1;
          l = mi_lookup.find(h);
          SMTBX_ASSERT(l != mi_lookup.end())(h.as_string());
          tmp[i] = std::conj(data[l->second][scatterer_idx]);
        }
        else{
          SMTBX_ASSERT(l != mi_lookup.end())(h.as_string());
          tmp[i] = data[l->second][scatterer_idx];
        }
      }
      return tmp;
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
      return new lookup_based_anisotropic(*this);
    }
  };

  /* A table that covers only part of the structure, spherical for the rest.

  A tabulated table is made for a particular set of atoms, and the structure it
  is used with may have gained an atom since, or never have been fully covered.
  Rather than refuse the table altogether, the atoms it names are taken from it
  and the others get the spherical form factor they would have had with no
  table at all. The mixture is not free of consequence -- part of the structure
  is then modelled less well than the rest -- so which atoms fell back is
  reported, for whoever has to say so.
  */
  template <typename FloatType>
  class partially_tabulated
    : public direct::one_scatterer_one_h::scatterer_contribution<FloatType>
  {
    typedef direct::one_scatterer_one_h::scatterer_contribution<FloatType>
      base_type;
  public:
    typedef FloatType float_type;
    typedef std::complex<float_type> complex_type;
  private:
    /* isotropic_scatterer_contribution keeps the registry by reference, so
    something has to own the one it is built against for as long as it lives.
    The caller's cannot be relied on: from Python the registry is a temporary
    that dies at the end of the call. Hence a copy here -- it is refcounted
    inside, so this is cheap -- declared before spherical_ so that it outlives
    it, and rebuilt rather than forked on copy for the same reason.
    */
    af::shared<xray::scatterer<float_type> > scatterers_;
    xray::scattering_type_registry registry_;
    boost::shared_ptr<base_type> tabulated_, spherical_;
    af::shared<bool> covered_;
    std::size_t n_smx_;
    mutable std::vector<complex_type> tmp_;
  public:
    partially_tabulated(const partially_tabulated &other)
      : scatterers_(other.scatterers_),
        registry_(other.registry_),
        tabulated_(other.tabulated_->raw_fork()),
        spherical_(new direct::one_scatterer_one_h::
          isotropic_scatterer_contribution<float_type>(scatterers_, registry_)),
        covered_(other.covered_),
        n_smx_(other.n_smx_),
        tmp_(other.n_smx_)
    {}

    partially_tabulated(base_type *tabulated,
      af::shared<xray::scatterer<float_type> > const &scatterers,
      xray::scattering_type_registry const &registry,
      af::shared<bool> const &covered, std::size_t n_smx)
      : scatterers_(scatterers),
        registry_(registry),
        tabulated_(tabulated),
        spherical_(new direct::one_scatterer_one_h::
          isotropic_scatterer_contribution<float_type>(scatterers_, registry_)),
        covered_(covered),
        n_smx_(n_smx),
        tmp_(n_smx)
    {}

    virtual complex_type get(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      return covered_[scatterer_idx] ? tabulated_->get(scatterer_idx, h)
                                     : spherical_->get(scatterer_idx, h);
    }

    /* A spherical form factor depends on the reflection only through its
    resolution, which every symmetry equivalent shares, so the entries this
    returns for an atom off the table are all the same number.
    */
    virtual std::vector<complex_type> const &get_full(std::size_t scatterer_idx,
      miller::index<> const &h) const
    {
      if (covered_[scatterer_idx]) {
        return tabulated_->get_full(scatterer_idx, h);
      }
      tmp_.assign(n_smx_, spherical_->get(scatterer_idx, h));
      return tmp_;
    }

    virtual base_type &at_d_star_sq(float_type d_star_sq) {
      tabulated_->at_d_star_sq(d_star_sq);
      spherical_->at_d_star_sq(d_star_sq);
      return *this;
    }

    //! Mixed, so not the uniform spherical case callers may optimise for.
    virtual bool is_spherical() const {
      return false;
    }

    virtual af::shared<std::size_t> scatterers_not_in_table() const {
      af::shared<std::size_t> result;
      for (std::size_t i = 0; i < covered_.size(); i++) {
        if (!covered_[i]) {
          result.push_back(i);
        }
      }
      return result;
    }

    virtual base_type *raw_fork() const {
      return new partially_tabulated(*this);
    }
  };

  template <typename FloatType>
  struct builder {
    /* Without a registry a table which misses an atom is an error, as before.
    build_with_fallback supplies one and gets the spherical fallback instead.
    */
    static direct::one_scatterer_one_h::scatterer_contribution<FloatType> *
      build(
        const uctbx::unit_cell& u_cell,
        af::shared< xray::scatterer<FloatType> > const &scatterers,
        std::string const &file_name,
        sgtbx::space_group const &space_group,
        bool anomalous_flag)
    {
      return build_impl(u_cell, scatterers, file_name, space_group,
        anomalous_flag, 0);
    }

    static direct::one_scatterer_one_h::scatterer_contribution<FloatType> *
      build_with_fallback(
        const uctbx::unit_cell& u_cell,
        af::shared< xray::scatterer<FloatType> > const &scatterers,
        std::string const &file_name,
        sgtbx::space_group const &space_group,
        bool anomalous_flag,
        xray::scattering_type_registry const &scattering_type_registry)
    {
      return build_impl(u_cell, scatterers, file_name, space_group,
        anomalous_flag, &scattering_type_registry);
    }

    static direct::one_scatterer_one_h::scatterer_contribution<FloatType> *
      build_impl(
        const uctbx::unit_cell& u_cell,
        af::shared< xray::scatterer<FloatType> > const &scatterers,
        std::string const &file_name,
        sgtbx::space_group const &space_group,
        bool anomalous_flag,
        xray::scattering_type_registry const *scattering_type_registry)
    {
      table_reader<FloatType> data(u_cell, scatterers, file_name);
      if(!data.use_AD()){
        anomalous_flag = false;
      }
      direct::one_scatterer_one_h::scatterer_contribution<FloatType> *result;
      if (data.rot_mxs().size() <= 1) {
        if (data.is_expanded()) {
          result = new lookup_based_anisotropic<FloatType>(
            scatterers,
            data,
            space_group,
            anomalous_flag);
        }
        else {
          result = new table_based_isotropic<FloatType>(
            scatterers,
            data,
            space_group,
            anomalous_flag);
        }
      }
      else {
        result = new table_based_anisotropic<FloatType>(
          scatterers,
          data,
          space_group,
          anomalous_flag);
      }
      /* A table which covers every atom is used as it stands. One which does
      not is paired with spherical form factors for the rest, so the refinement
      can go on -- but only if the caller supplied the registry to compute them
      from. Without it there is nothing to fall back on, and a table that does
      not fit the structure is an error as it was before.
      */
      af::shared<std::size_t> missing = data.not_covered();
      if (missing.size() != 0) {
        if (scattering_type_registry == 0) {
          delete result;
          throw SMTBX_ERROR("the table does not cover the whole structure and"
            " no scattering type registry was given to fall back on");
        }
        result = new partially_tabulated<FloatType>(
          result,
          scatterers,
          *scattering_type_registry,
          data.covered(),
          space_group.n_smx());
      }
      return result;
    }

    static direct::one_scatterer_one_h::scatterer_contribution<FloatType> *
      build_lookup_based_for_tests(
        uctbx::unit_cell const &unit_cell,
        sgtbx::space_group const &space_group,
        af::shared<xray::scatterer<FloatType> > const &scatterers,
        xray::scattering_type_registry const &scattering_type_registry,
        af::shared<cctbx::miller::index<> > const &indices)
    {
      direct::one_scatterer_one_h::isotropic_scatterer_contribution<FloatType>
        isc(scatterers, scattering_type_registry);
      return new lookup_based_anisotropic<FloatType>(
        unit_cell,
        space_group,
        scatterers,
        isc,
        indices);
    }
  };

}}} // smtbx::structure_factors::table_based

#endif // GUARD
