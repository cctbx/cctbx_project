#include <iotbx/mtz/batch.h>
#include <iotbx/mtz/column.h>

namespace iotbx { namespace mtz {

  af::shared<std::size_t>
  cmtz_struct_sizes()
  {
    af::shared<std::size_t> result;
    result.push_back(sizeof(CMtz::MTZCOL));
    result.push_back(sizeof(CMtz::MTZSET));
    result.push_back(sizeof(CMtz::MTZXTAL));
    result.push_back(sizeof(CMtz::MTZBAT));
    result.push_back(sizeof(CMtz::SYMGRP));
    result.push_back(sizeof(CMtz::MTZ));
    return result;
  }

  object::object()
  :
    ptr_(CMtz::MtzMalloc(0, 0), ptr_deleter)
  {
    if (ptr_.get() == 0) throw cctbx::error("MtzMalloc failed.");
    ptr_->refs_in_memory = true;
  }

  object::object(const char* file_name)
  {
    CCTBX_ASSERT(file_name != 0);
    ptr_ = boost::shared_ptr<CMtz::MTZ>(
      CMtz::MtzGet(file_name, true), ptr_deleter);
    if (ptr_.get() == 0) {
      throw cctbx::error(std::string("MTZ file read error: ") + file_name);
    }
  }

  std::string
  object::title() const
  {
    char result[sizeof(ptr()->title)];
    int title_length = CMtz::ccp4_lrtitl(ptr(), result);
    return std::string(result, title_length);
  }

  object&
  object::set_title(const char* title, bool append)
  {
    CCTBX_ASSERT(title != 0);
    int set_title_success = CMtz::ccp4_lwtitl(ptr(), title, append);
    CCTBX_ASSERT(set_title_success);
    return *this;
  }

  af::shared<std::string>
  object::history() const
  {
    CMtz::MTZ* p = ptr();
    af::shared<std::string> result((af::reserve(p->histlines)));
    for(int i=0;i<p->histlines;i++) {
      const char* line = p->hist+MTZRECORDLENGTH*i;
      int j = 0;
      for(;j<MTZRECORDLENGTH;j++) {
        if (line[j] == '\0') break;
      }
      result.push_back(std::string(line, j));
    }
    return result;
  }

  object&
  object::add_history(af::const_ref<std::string> const& lines)
  {
    boost::shared_ptr<char> buffer(
      CMtz::MtzCallocHist(lines.size()), CMtz::MtzFreeHist);
    for(std::size_t i=0;i<lines.size();i++) {
      strncpy(
        buffer.get()+i*MTZRECORDLENGTH,
        lines[i].c_str(),
        std::min(
          static_cast<std::size_t>(MTZRECORDLENGTH),
          lines[i].size()));
    }
    int add_history_success = CMtz::MtzAddHistory(
      ptr(),
      reinterpret_cast<char (*)[MTZRECORDLENGTH]>(buffer.get()),
      lines.size());
    CCTBX_ASSERT(add_history_success);
    return *this;
  }

  object&
  object::add_history(const char* line)
  {
    CCTBX_ASSERT(line != 0);
    return add_history(af::tiny<std::string, 1>(line).const_ref());
  }

  object&
  object::set_space_group_name(const char* name)
  {
    CCTBX_ASSERT(name != 0);
    char* target = ptr()->mtzsymm.spcgrpname;
    const unsigned target_size = sizeof(ptr()->mtzsymm.spcgrpname);
    strncpy(target, name, target_size-1);
    target[target_size-1] = '\0';
    return *this;
  }

  object&
  object::set_point_group_name(const char* name)
  {
    CCTBX_ASSERT(name != 0);
    char* target = ptr()->mtzsymm.pgname;
    const unsigned target_size = sizeof(ptr()->mtzsymm.pgname);
    strncpy(target, name, target_size-1);
    target[target_size-1] = '\0';
    return *this;
  }

  cctbx::sgtbx::space_group
  object::space_group() const
  {
    CMtz::MTZ* p = ptr();
    cctbx::sgtbx::space_group result;
    scitbx::mat3<double> rm;
    scitbx::vec3<double> tv;
    for(int im=0;im<p->mtzsymm.nsym;im++) {
      for (int ir=0;ir<3;ir++) {
        for (int ic=0;ic<3;ic++) {
          rm(ir,ic) = p->mtzsymm.sym[im][ir][ic];
        }
        tv[ir] = p->mtzsymm.sym[im][ir][3];
      }
      result.expand_smx(cctbx::sgtbx::rt_mx(rm, tv));
    }
    return result;
  }

  object&
  object::set_space_group(cctbx::sgtbx::space_group const& space_group)
  {
    CMtz::MTZ* p = ptr();
    CCTBX_ASSERT(sizeof(p->mtzsymm.sym) / sizeof(*p->mtzsymm.sym)
              >= space_group.order_z());
    p->mtzsymm.nsymp = static_cast<int>(space_group.order_p());
    p->mtzsymm.nsym = static_cast<int>(space_group.order_z());
    for (int im=0;im<p->mtzsymm.nsym;im++) {
      cctbx::sgtbx::rt_mx sm = space_group(im).mod_positive();
      cctbx::sgtbx::rot_mx rm = sm.r();
      cctbx::sgtbx::tr_vec tv = sm.t();
      scitbx::mat3<int> rm_num = rm.num();
      float rm_den = rm.den();
      scitbx::vec3<int> tv_num = tv.num();
      float tv_den = tv.den();
      for (int ir=0;ir<3;ir++) {
        for (int ic=0;ic<3;ic++) {
          p->mtzsymm.sym[im][ir][ic] = rm_num(ir,ic)/rm_den;
        }
        p->mtzsymm.sym[im][ir][3] = tv_num[ir]/tv_den;
      }
      for (int ic=0;ic<3;ic++) {
        p->mtzsymm.sym[im][3][ic] = 0.;
      }
      p->mtzsymm.sym[im][3][3] = 1.;
    }
    return *this;
  }

  void
  object::reserve(int capacity)
  {
    CMtz::MTZ* p = ptr();
    if (!p->refs_in_memory) return;
    for(int i=0;i<p->nxtal;i++) {
      for(int j=0;j<p->xtal[i]->nset;j++) {
        for(int k=0;k<p->xtal[i]->set[j]->ncol;k++) {
          capacity = std::max(
            capacity,
            column_array_capacity(p->xtal[i]->set[j]->col[k]));
        }
      }
    }
    for(int i=0;i<p->nxtal;i++) {
      for(int j=0;j<p->xtal[i]->nset;j++) {
        for(int k=0;k<p->xtal[i]->set[j]->ncol;k++) {
          ccp4array_reserve(p->xtal[i]->set[j]->col[k]->ref, capacity);
        }
      }
    }
  }

  void
  object::adjust_column_array_sizes(int new_nref)
  {
    CMtz::MTZ* p = ptr();
    if (!p->refs_in_memory) return;
    if (new_nref > p->nref) {
      reserve(new_nref);
      for(int i=0;i<p->nxtal;i++) {
        for(int j=0;j<p->xtal[i]->nset;j++) {
          for(int k=0;k<p->xtal[i]->set[j]->ncol;k++) {
            CMtz::MTZCOL* col_k = p->xtal[i]->set[j]->col[k];
            int old_size = column_array_size(col_k);
            if (new_nref > old_size) {
              ccp4array_resize(col_k->ref, new_nref);
              for(int iref=old_size;iref<new_nref;iref++) {
                *(reinterpret_cast<union float_uint_uchar*>(
                  &col_k->ref[iref])) = CCP4::ccp4_nan();
              }
            }
          }
        }
      }
    }
  }

  af::shared<batch>
  object::batches() const
  {
    af::shared<batch> result((af::reserve(n_batches())));
    for(int i_batch=0;i_batch<n_batches();i_batch++) {
      result.push_back(batch(*this, i_batch));
    }
    return result;
  }

  batch
  object::add_batch()
  {
    CMtz::MTZBAT* p = ptr()->batch;
    CMtz::MTZBAT* p_tail = p;
    int max_batch_number = 0;
    int i_batch = 0;
    for(;;i_batch++) {
      if (p == 0) break;
      max_batch_number = std::max(max_batch_number, p->num);
      p_tail = p;
      p = p->next;
    }
    std::vector<float> buf(
      static_cast<std::size_t>(NBATCHINTEGERS + NBATCHREALS),
      static_cast<float>(0));
#if defined(__DECCXX_VER)
    BOOST_STATIC_ASSERT(sizeof(float) == sizeof(int));
#else
    CCTBX_ASSERT(sizeof(float) == sizeof(int));
#endif
    CCTBX_ASSERT(CMtz::ccp4_lwbat(
      ptr(), 0, max_batch_number+1, &*buf.begin(), "") == 1);
    p = (p_tail == 0 ? ptr()->batch : p_tail->next);
    CCTBX_ASSERT(p != 0);
    CCTBX_ASSERT(p->next == 0);
    CCTBX_ASSERT(p->num == max_batch_number+1);
    return batch(*this, i_batch);
  }

  af::shared<crystal>
  object::crystals() const
  {
    af::shared<crystal> result((af::reserve(n_crystals())));
    for(int i_crystal=0;i_crystal<n_crystals();i_crystal++) {
      result.push_back(crystal(*this, i_crystal));
    }
    return result;
  }

  crystal
  object::add_crystal(
    const char* name,
    const char* project_name,
    af::double6 const& unit_cell_parameters)
  {
    CCTBX_ASSERT(name != 0);
    CCTBX_ASSERT(project_name != 0);
    float uc_params[6];
    for(int i=0;i<6;i++) uc_params[i] = unit_cell_parameters[i];
    int i_crystal = n_crystals();
    CMtz::MTZXTAL* crystal_ptr = CMtz::MtzAddXtal(
      ptr(), name, project_name, uc_params);
    CCTBX_ASSERT(crystal_ptr != 0);
    CCTBX_ASSERT(n_crystals() == i_crystal+1);
    crystal result(*this, i_crystal);
    CCTBX_ASSERT(result.ptr() == crystal_ptr);
    return result;
  }

  crystal
  object::add_crystal(
    const char* name,
    const char* project_name,
    cctbx::uctbx::unit_cell const& unit_cell)
  {
    return add_crystal(name, project_name, unit_cell.parameters());
  }

  bool
  object::has_column(const char* label) const
  {
    CCTBX_ASSERT(label != 0);
    for(int i_crystal=0;i_crystal<n_crystals();i_crystal++) {
      crystal x(*this, i_crystal);
      for(int i_dataset=0;i_dataset<x.n_datasets();i_dataset++) {
        dataset s(x, i_dataset);
        for(int i_column=0;i_column<s.n_columns();i_column++) {
          column c(s, i_column);
          if (CMtz::MtzPathMatch(label, c.label())) {
            return true;
          }
        }
      }
    }
    return false;
  }

  column
  object::get_column(const char* label) const
  {
    CCTBX_ASSERT(label != 0);
    for(int i_crystal=0;i_crystal<n_crystals();i_crystal++) {
      crystal x(*this, i_crystal);
      for(int i_dataset=0;i_dataset<x.n_datasets();i_dataset++) {
        dataset s(x, i_dataset);
        for(int i_column=0;i_column<s.n_columns();i_column++) {
          column c(s, i_column);
          if (CMtz::MtzPathMatch(label, c.label())) {
            return c;
          }
        }
      }
    }
    throw cctbx::error(std::string("Unknown MTZ column label: ") + label);
  }

  hkl_columns
  object::get_hkl_columns() const
  {
    return hkl_columns(
      get_column("H"),
      get_column("K"),
      get_column("L"));
  }

  af::double2
  object::max_min_resolution() const
  {
    double d_star_sq_min = -1;
    double d_star_sq_max = -1;
    int n_refl = n_reflections();
    int n_crys = n_crystals();
    if (n_refl > 0 && n_crys > 0) {
      hkl_columns hkl = get_hkl_columns();
      std::vector<cctbx::uctbx::unit_cell> unit_cells;
      for(int i_crystal=0;i_crystal<n_crys;i_crystal++) {
        unit_cells.push_back(crystal(*this, i_crystal).unit_cell());
      }
      for(int i_refl=0;i_refl<n_refl;i_refl++) {
        for(int i_crystal=0;i_crystal<n_crys;i_crystal++) {
          double d_star_sq = unit_cells[i_crystal].d_star_sq(
            hkl.get_miller_index(i_refl));
          if (d_star_sq_min > d_star_sq || d_star_sq_min < 0) {
              d_star_sq_min = d_star_sq;
          }
          if (d_star_sq_max < d_star_sq) {
              d_star_sq_max = d_star_sq;
          }
        }
      }
    }
    return af::double2(
      d_star_sq_min <= 0 ? -1 : 1/std::sqrt(d_star_sq_min),
      d_star_sq_max <= 0 ? -1 : 1/std::sqrt(d_star_sq_max));
  }

  namespace {
    inline
    std::complex<double>
    polar_deg(double ampl, double phi)
    {
      return std::polar(ampl, phi*scitbx::constants::pi_180);
    }
  }

  af::shared<cctbx::miller::index<> >
  object::extract_miller_indices() const
  {
    int n_refl = n_reflections();
    af::shared<cctbx::miller::index<> > result((af::reserve(n_refl)));
    hkl_columns hkl = get_hkl_columns();
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      result.push_back(hkl.get_miller_index(i_refl));
    }
    return result;
  }

  void
  object::replace_miller_indices(
    af::const_ref<cctbx::miller::index<> > const& miller_indices)
  {
    CCTBX_ASSERT(miller_indices.size() == n_reflections());
    hkl_columns hkl = get_hkl_columns();
    for(int i_refl=0;i_refl<miller_indices.size();i_refl++) {
      hkl.replace_miller_index(i_refl, miller_indices[i_refl]);
    }
  }

  integer_group
  object::extract_integers(
    const char* column_label) const
  {
    int n_refl = n_reflections();
    integer_group result(false, n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data(get_column(column_label));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (!data.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(data.int_datum(i_refl));
      }
    }
    return result;
  }

  af::shared<int>
  object::extract_integers(
    af::const_ref<int> const& mtz_reflection_indices,
    const char* column_label) const
  {
    int n_refl = n_reflections();
    af::shared<int> result((af::reserve(mtz_reflection_indices.size())));
    column data(get_column(column_label));
    for(int i=0;i<mtz_reflection_indices.size();i++) {
      int i_refl = mtz_reflection_indices[i];
      CCTBX_ASSERT(i_refl >= 0 && i_refl < n_refl);
      CCTBX_ASSERT(!data.is_ccp4_nan(i_refl));
      result.push_back(data.int_datum(i_refl));
    }
    return result;
  }

  integer_group
  object::extract_integers_anomalous(
    const char* column_label_plus,
    const char* column_label_minus) const
  {
    int n_refl = n_reflections();
    integer_group result(true, 2*n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data_plus(get_column(column_label_plus));
    column data_minus(get_column(column_label_minus));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (!data_plus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(data_plus.int_datum(i_refl));
      }
      if (!data_minus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(-hkl.get_miller_index(i_refl));
        result.data.push_back(data_minus.int_datum(i_refl));
      }
    }
    return result;
  }

  real_group
  object::extract_reals(
    const char* column_label) const
  {
    int n_refl = n_reflections();
    real_group result(false, n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data(get_column(column_label));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (!data.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(data.float_datum(i_refl));
      }
    }
    return result;
  }

  af::shared<double>
  object::extract_reals(
    af::const_ref<int> const& mtz_reflection_indices,
    const char* column_label) const
  {
    int n_refl = n_reflections();
    af::shared<double> result((af::reserve(mtz_reflection_indices.size())));
    column data(get_column(column_label));
    for(int i=0;i<mtz_reflection_indices.size();i++) {
      int i_refl = mtz_reflection_indices[i];
      CCTBX_ASSERT(i_refl >= 0 && i_refl < n_refl);
      CCTBX_ASSERT(!data.is_ccp4_nan(i_refl));
      result.push_back(data.float_datum(i_refl));
    }
    return result;
  }

  real_group
  object::extract_reals_anomalous(
    const char* column_label_plus,
    const char* column_label_minus) const
  {
    int n_refl = n_reflections();
    real_group result(true, 2*n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data_plus(get_column(column_label_plus));
    column data_minus(get_column(column_label_minus));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (!data_plus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(data_plus.float_datum(i_refl));
      }
      if (!data_minus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(-hkl.get_miller_index(i_refl));
        result.data.push_back(data_minus.float_datum(i_refl));
      }
    }
    return result;
  }

  hl_group
  object::extract_hendrickson_lattman(
    const char* column_label_a,
    const char* column_label_b,
    const char* column_label_c,
    const char* column_label_d) const
  {
    int n_refl = n_reflections();
    hl_group result(false, n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data_a(get_column(column_label_a));
    column data_b(get_column(column_label_b));
    column data_c(get_column(column_label_c));
    column data_d(get_column(column_label_d));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (   data_b.is_ccp4_nan(i_refl) != data_a.is_ccp4_nan(i_refl)
          || data_c.is_ccp4_nan(i_refl) != data_a.is_ccp4_nan(i_refl)
          || data_d.is_ccp4_nan(i_refl) != data_a.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting Hendrickson-Lattman array"
          " from columns: ")
          + column_label_a + ", "
          + column_label_b + ", "
          + column_label_c + ", "
          + column_label_d);
      }
      if (!data_a.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(cctbx::hendrickson_lattman<>(
          data_a.float_datum(i_refl),
          data_b.float_datum(i_refl),
          data_c.float_datum(i_refl),
          data_d.float_datum(i_refl)));
      }
    }
    return result;
  }

  hl_group
  object::extract_hendrickson_lattman_anomalous(
    const char* column_label_a_plus,
    const char* column_label_b_plus,
    const char* column_label_c_plus,
    const char* column_label_d_plus,
    const char* column_label_a_minus,
    const char* column_label_b_minus,
    const char* column_label_c_minus,
    const char* column_label_d_minus) const
  {
    int n_refl = n_reflections();
    hl_group result(true, n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data_ap(get_column(column_label_a_plus));
    column data_bp(get_column(column_label_b_plus));
    column data_cp(get_column(column_label_c_plus));
    column data_dp(get_column(column_label_d_plus));
    column data_am(get_column(column_label_a_minus));
    column data_bm(get_column(column_label_b_minus));
    column data_cm(get_column(column_label_c_minus));
    column data_dm(get_column(column_label_d_minus));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (   data_bp.is_ccp4_nan(i_refl) != data_ap.is_ccp4_nan(i_refl)
          || data_cp.is_ccp4_nan(i_refl) != data_ap.is_ccp4_nan(i_refl)
          || data_dp.is_ccp4_nan(i_refl) != data_ap.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting Hendrickson-Lattman array"
          " from columns: ")
          + column_label_a_plus + ", "
          + column_label_b_plus + ", "
          + column_label_c_plus + ", "
          + column_label_d_plus);
      }
      if (!data_ap.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(cctbx::hendrickson_lattman<>(
          data_ap.float_datum(i_refl),
          data_bp.float_datum(i_refl),
          data_cp.float_datum(i_refl),
          data_dp.float_datum(i_refl)));
      }
      if (   data_bm.is_ccp4_nan(i_refl) != data_am.is_ccp4_nan(i_refl)
          || data_cm.is_ccp4_nan(i_refl) != data_am.is_ccp4_nan(i_refl)
          || data_dm.is_ccp4_nan(i_refl) != data_am.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting Hendrickson-Lattman array"
          " from columns: ")
          + column_label_a_minus + ", "
          + column_label_b_minus + ", "
          + column_label_c_minus + ", "
          + column_label_d_minus);
      }
      if (!data_am.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(-hkl.get_miller_index(i_refl));
        result.data.push_back(cctbx::hendrickson_lattman<>(
          data_am.float_datum(i_refl),
          data_bm.float_datum(i_refl),
          data_cm.float_datum(i_refl),
          data_dm.float_datum(i_refl)));
      }
    }
    return result;
  }

  observations_group
  object::extract_observations(
    const char* column_label_data,
    const char* column_label_sigmas) const
  {
    int n_refl = n_reflections();
    observations_group result(false, n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data(get_column(column_label_data));
    column sigmas(get_column(column_label_sigmas));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (data.is_ccp4_nan(i_refl) != sigmas.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting observation array from columns: ")
          + column_label_data + ", " + column_label_sigmas);
      }
      if (!data.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(data.float_datum(i_refl));
        result.sigmas.push_back(sigmas.float_datum(i_refl));
      }
    }
    return result;
  }

  observations_group
  object::extract_observations_anomalous(
    const char* column_label_data_plus,
    const char* column_label_sigmas_plus,
    const char* column_label_data_minus,
    const char* column_label_sigmas_minus) const
  {
    int n_refl = n_reflections();
    observations_group result(true, 2*n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data_plus(get_column(column_label_data_plus));
    column sigmas_plus(get_column(column_label_sigmas_plus));
    column data_minus(get_column(column_label_data_minus));
    column sigmas_minus(get_column(column_label_sigmas_minus));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (data_plus.is_ccp4_nan(i_refl) != sigmas_plus.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting observation array from columns: ")
          + column_label_data_plus + ", " + column_label_sigmas_plus);
      }
      if (!data_plus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(data_plus.float_datum(i_refl));
        result.sigmas.push_back(sigmas_plus.float_datum(i_refl));
      }
      if (data_minus.is_ccp4_nan(i_refl) != sigmas_minus.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting observation array from columns: ")
          + column_label_data_minus + ", " + column_label_sigmas_minus);
      }
      if (!data_minus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(-hkl.get_miller_index(i_refl));
        result.data.push_back(data_minus.float_datum(i_refl));
        result.sigmas.push_back(sigmas_minus.float_datum(i_refl));
      }
    }
    return result;
  }

  observations_group
  object::extract_delta_anomalous(
    const char* column_label_f_data,
    const char* column_label_f_sigmas,
    const char* column_label_d_data,
    const char* column_label_d_sigmas) const
  {
    int n_refl = n_reflections();
    observations_group result(true, 2*n_refl);
    hkl_columns hkl = get_hkl_columns();
    column f_data(get_column(column_label_f_data));
    column f_sigmas(get_column(column_label_f_sigmas));
    column d_data(get_column(column_label_d_data));
    column d_sigmas(get_column(column_label_d_sigmas));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (   f_sigmas.is_ccp4_nan(i_refl) != f_data.is_ccp4_nan(i_refl)
          || d_sigmas.is_ccp4_nan(i_refl) != d_data.is_ccp4_nan(i_refl)
          || (f_data.is_ccp4_nan(i_refl) && !d_data.is_ccp4_nan(i_refl))) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_f_data   + ", "
          + column_label_f_sigmas + ", "
          + column_label_d_data   + ", "
          + column_label_d_sigmas);
      }
      if (!f_data.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        double fd = f_data.float_datum(i_refl);
        double fs = f_sigmas.float_datum(i_refl);
        if (d_data.is_ccp4_nan(i_refl)) {
          result.data.push_back(fd);
          result.sigmas.push_back(fs);
        }
        else {
          double ddh = d_data.float_datum(i_refl) * .5;
          double dsh = d_sigmas.float_datum(i_refl) * .5;
          double s = std::sqrt(fs*fs + dsh*dsh);
          result.mtz_reflection_indices.push_back(i_refl);
          result.indices.push_back(-hkl.get_miller_index(i_refl));
          result.data.push_back(fd + ddh);
          result.data.push_back(fd - ddh);
          result.sigmas.push_back(std::sqrt(s));
          result.sigmas.push_back(std::sqrt(s));
        }
      }
    }
    return result;
  }

  complex_group
  object::extract_complex(
    const char* column_label_ampl,
    const char* column_label_phi) const
  {
    int n_refl = n_reflections();
    complex_group result(false, n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data_ampl(get_column(column_label_ampl));
    column data_phi(get_column(column_label_phi));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (data_ampl.is_ccp4_nan(i_refl) != data_phi.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_ampl + ", " + column_label_phi);
      }
      if (!data_ampl.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(polar_deg(
          data_ampl.float_datum(i_refl), data_phi.float_datum(i_refl)));
      }
    }
    return result;
  }

  complex_group
  object::extract_complex_anomalous(
    const char* column_label_ampl_plus,
    const char* column_label_phi_plus,
    const char* column_label_ampl_minus,
    const char* column_label_phi_minus) const
  {
    int n_refl = n_reflections();
    complex_group result(true, n_refl);
    hkl_columns hkl = get_hkl_columns();
    column data_ampl_plus(get_column(column_label_ampl_plus));
    column data_phi_plus(get_column(column_label_phi_plus));
    column data_ampl_minus(get_column(column_label_ampl_minus));
    column data_phi_minus(get_column(column_label_phi_minus));
    for(int i_refl=0;i_refl<n_refl;i_refl++) {
      if (   data_ampl_plus.is_ccp4_nan(i_refl)
          != data_phi_plus.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_ampl_plus + ", " + column_label_phi_plus);
      }
      if (!data_ampl_plus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(hkl.get_miller_index(i_refl));
        result.data.push_back(polar_deg(
          data_ampl_plus.float_datum(i_refl),
          data_phi_plus.float_datum(i_refl)));
      }
      if (   data_ampl_minus.is_ccp4_nan(i_refl)
          != data_phi_minus.is_ccp4_nan(i_refl)) {
        throw cctbx::error(std::string(
          "Unexpected NAN while extracting complex array from columns: ")
          + column_label_ampl_minus + ", " + column_label_phi_minus);
      }
      if (!data_ampl_minus.is_ccp4_nan(i_refl)) {
        result.mtz_reflection_indices.push_back(i_refl);
        result.indices.push_back(-hkl.get_miller_index(i_refl));
        result.data.push_back(polar_deg(
          data_ampl_minus.float_datum(i_refl),
          data_phi_minus.float_datum(i_refl)));
      }
    }
    return result;
  }

  void
  object::write(const char* file_name)
  {
    CCTBX_ASSERT(file_name != 0);
    if (!CMtz::MtzPut(ptr(), file_name)) {
      throw cctbx::error("MTZ write failed.");
    }
  }

  void
  object::ptr_deleter(CMtz::MTZ* ptr)
  {
    if (ptr != 0) {
      if (ptr->batch != 0 && ptr->n_orig_bat == 0) {
        ptr->n_orig_bat = 1; // force MtzFreeBatch
      }
      CCTBX_ASSERT(CMtz::MtzFree(ptr));
    }
  }

}} // namespace iotbx::mtz
