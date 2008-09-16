#include <iotbx/mtz/dataset.h>
#include <iotbx/error.h>

namespace iotbx { namespace mtz {

  crystal::crystal(object const& mtz_object, int i_crystal)
  :
    mtz_object_(mtz_object),
    i_crystal_(i_crystal)
  {
    IOTBX_ASSERT(i_crystal >= 0);
    IOTBX_ASSERT(i_crystal < mtz_object.n_crystals());
  }

  CMtz::MTZXTAL*
  crystal::ptr() const
  {
    IOTBX_ASSERT(mtz_object_.n_crystals() > i_crystal_);
    return CMtz::MtzIxtal(mtz_object_.ptr(), i_crystal_);
  }

  crystal&
  crystal::set_id(int id)
  {
    if (ptr()->xtalid != id) {
      CMtz::MTZ* p = mtz_object().ptr();
      IOTBX_ASSERT(p->refs_in_memory);
      for(int i=0;i<p->nxtal;i++) {
        IOTBX_ASSERT(p->xtal[i]->xtalid != id);
      }
      ptr()->xtalid = id;
    }
    return *this;
  }

  crystal&
  crystal::set_name(const char* new_name)
  {
    IOTBX_ASSERT(new_name != 0);
    IOTBX_ASSERT(std::strlen(new_name) < sizeof(ptr()->xname));
    if (std::strcmp(ptr()->xname, new_name) == 0) return *this;
    if (mtz_object().has_crystal(new_name)) {
      throw std::runtime_error(
          std::string("mtz::crystal::set_name(new_name=\"")
        + new_name
        + "\"): new_name is used already for another crystal.");
    }
    std::strcpy(ptr()->xname, new_name);
    return *this;
  }

  crystal&
  crystal::set_project_name(const char* new_project_name)
  {
    IOTBX_ASSERT(new_project_name != 0);
    IOTBX_ASSERT(std::strlen(new_project_name) < sizeof(ptr()->pname));
    std::strcpy(ptr()->pname, new_project_name);
    return *this;
  }

  af::small<double, 6>
  crystal::unit_cell_parameters() const
  {
    af::small<double, 6> result;
    float* cell = ptr()->cell;
    for(std::size_t i=0;i<6;i++) result.push_back(cell[i]);
    return result;
  }

  crystal&
  crystal::set_unit_cell_parameters(af::small<double, 6> const& parameters)
  {
    float* cell = ptr()->cell;
    for(std::size_t i=0;i<6;i++) cell[i] = static_cast<float>(parameters[i]);
    return *this;
  }

  af::shared<dataset>
  crystal::datasets() const
  {
    af::shared<dataset> result((af::reserve(n_datasets())));
    for(int i_dataset=0;i_dataset<n_datasets();i_dataset++) {
      result.push_back(dataset(*this, i_dataset));
    }
    return result;
  }

  dataset
  crystal::add_dataset(
    const char *name,
    double wavelength)
  {
    int i_dataset = n_datasets();
    CMtz::MTZSET* dataset_ptr = CMtz::MtzAddDataset(
      mtz_object().ptr(), ptr(), name, static_cast<float>(wavelength));
    IOTBX_ASSERT(dataset_ptr != 0);
    IOTBX_ASSERT(n_datasets() == i_dataset+1);
    dataset result(*this, i_dataset);
    IOTBX_ASSERT(result.ptr() == dataset_ptr);
    return result;
  }

  bool
  crystal::has_dataset(const char* name) const
  {
    IOTBX_ASSERT(name != 0);
    for(int i_dataset=0;i_dataset<n_datasets();i_dataset++) {
      dataset s(*this, i_dataset);
      if (std::strcmp(s.name(), name) == 0) {
        return true;
      }
    }
    return false;
  }

}} // namespace iotbx::mtz
