#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/mat3.h>
#include <iostream>
#include <iotbx/mtzwriter.h>
#include <string.h>

namespace iotbx { namespace mtz {

MtzWriter::MtzWriter()
:
  mtz(CMtz::MtzMalloc(0,0))
{
  CCTBX_ASSERT(mtz != NULL);
}

MtzWriter::~MtzWriter() { CMtz::MtzFree(mtz); }

void
MtzWriter::setTitle(const std::string& title)
{
  const std::size_t max_chars = sizeof(mtz->title);
  strncpy(mtz->title, title.c_str(), max_chars-1);
  mtz->title[max_chars-1] = '\0';
}

void
MtzWriter::setSpaceGroup(
  const cctbx::sgtbx::space_group& sg,
  const std::string& symbol)
{
  int val_nsymx = sg.order_z();
  int val_nsympx= sg.order_p();
  float rsymx[192][4][4];

  CCTBX_ASSERT(sg.order_z() <= 192);
  for(std::size_t i=0;i<sg.order_z();i++) {
    cctbx::sgtbx::rt_mx s = sg(i);
    cctbx::sgtbx::rot_mx r = s.r();
    cctbx::sgtbx::tr_vec t = s.t();
    scitbx::mat3<int> r_num = r.num();
    float r_den = r.den();
    scitbx::vec3<int> t_num = t.num();
    float t_den = t.den();
    for (std::size_t p=0;p<3;p++) {
      for (std::size_t q=0;q<3;q++) {
        rsymx[i][p][q] = r_num(p,q)/r_den;
      }
      rsymx[i][p][3] = t_num[p]/t_den;
      rsymx[i][3][p] = 0.0;
    }
    rsymx[i][3][3] = 1.0;
  }

  CCTBX_ASSERT(sg.conventional_centring_type_symbol() != '\0');
  char val_ltypex = sg.conventional_centring_type_symbol();

  cctbx::sgtbx::space_group_type sgtype(sg.type());
  int val_nspgrx = sgtype.number();

  char sgtype_s[11];
  CCTBX_ASSERT(strlen(symbol.c_str()) <= 10);
  sprintf(sgtype_s,"%s",symbol.c_str());

  cctbx::sgtbx::matrix_group::code mgcode = sg.point_group_type();
  char pgname[11];
  CCTBX_ASSERT(strlen(mgcode.label()) <= 10);
  sprintf(pgname,"%s",mgcode.label());

  CMtz::ccp4_lwsymm(mtz, val_nsymx, val_nsympx, rsymx, &val_ltypex,
                    val_nspgrx, sgtype_s, pgname);
}

void
MtzWriter::oneCrystal(
  const std::string& crystal,
  const std::string& project,
  const cctbx::uctbx::unit_cell& uc)
{
  std::size_t ind;
  std::string strip_crystal = crystal;
  std::string strip_project = project;
  while ((ind=strip_crystal.find(' '))!=std::string::npos) {
    strip_crystal.erase(ind,1);
  }
  while ((ind=strip_project.find(' '))!=std::string::npos) {
    strip_project.erase(ind,1);
  }
  float cell[6];
  for (ind=0;ind<6;++ind) {cell[ind]=uc.parameters()[ind];}
  onextal = CMtz::MtzAddXtal(mtz, strip_crystal.c_str(), strip_project.c_str(),
            cell);
}

void
MtzWriter::oneDataset(
  const std::string& dataset,
  const double& w)
{
  float wavelength=w;
  std::size_t ind;
  std::string strip_dataset = dataset;
  while ((ind=strip_dataset.find(' '))!=std::string::npos) {
    strip_dataset.erase(ind,1);
  }
  oneset = CMtz::MtzAddDataset(mtz,onextal,strip_dataset.c_str(),wavelength);
  CMtz::MtzAddColumn(mtz,oneset,"H","H");
  CMtz::MtzAddColumn(mtz,oneset,"K","H");
  CMtz::MtzAddColumn(mtz,oneset,"L","H");
}

void
MtzWriter::safe_ccp4_lwrefl(
  const float adata[],
  CMtz::MTZCOL *lookup[],
  const int ncol,
  const int iref)
{
  int i,j,k;
  /* if this is extra reflection, check memory */
  if (mtz->refs_in_memory && iref > mtz->nref) {
   /* Loop over crystals */
    for (i = 0; i < mtz->nxtal; ++i) {
   /* Loop over datasets for each crystal */
     for (j = 0; j < mtz->xtal[i]->nset; ++j) {
    /* Loop over columns for each dataset */
      for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
       if (iref > ccp4array_size(mtz->xtal[i]->set[j]->col[k]->ref)){
         //Always reallocate in increments of 8K
         int newsize = 8192 * (1 + iref/8192);
         ccp4array_resize(mtz->xtal[i]->set[j]->col[k]->ref, newsize);
       }
      }
     }
    }
  }

  if (mtz->refs_in_memory && iref > mtz->nref) {
  /* Loop over crystals */
   for (i = 0; i < mtz->nxtal; ++i) {
  /* Loop over datasets for each crystal */
    for (j = 0; j < mtz->xtal[i]->nset; ++j) {
   /* Loop over columns for each dataset */
     for (k = 0; k < mtz->xtal[i]->set[j]->ncol; ++k) {
      mtz->xtal[i]->set[j]->col[k]->ref[iref-1] = CCP4::ccp4_nan().f;
     }
    }
   }
  }

  for (i = 0; i < ncol; ++i) {
    if (lookup[i]) {
      /* update reflection in memory or add to refldata array. */
      if (mtz->refs_in_memory) {
        lookup[i]->ref[iref-1] = adata[i];
      }
      /* update column ranges */
      if (!CMtz::ccp4_ismnf(mtz, adata[i])) {
        if (adata[i] < lookup[i]->min) lookup[i]->min = adata[i];
        if (adata[i] > lookup[i]->max) lookup[i]->max = adata[i];
      }
    }
  }

  /* increment nref if we are adding new reflections */
  if (iref > mtz->nref)
    mtz->nref = iref;
}

void
MtzWriter::addColumn(
  const std::string& name,
  char type_code,
  af::const_ref<cctbx::miller::index<> > const& miller_indices,
  af::const_ref<double> const& data)
{
  using namespace cctbx;
  CMtz::MTZCOL* write_columns[4];
  write_columns[0]=CMtz::MtzColLookup(mtz,"H");
  write_columns[1]=CMtz::MtzColLookup(mtz,"K");
  write_columns[2]=CMtz::MtzColLookup(mtz,"L");
  if (CMtz::MtzColLookup(mtz,name.c_str())!=NULL)
    throw Error("Attempt to overwrite existing column "+name);
  write_columns[3]=CMtz::MtzAddColumn(mtz,oneset,name.c_str(),&type_code);
  CCTBX_ASSERT(miller_indices.size() == data.size());
  typedef std::map<miller::index<>, std::size_t> lookup_dict_type;
  lookup_dict_type lookup_dict;
  for( std::size_t i=0; i<mtz->nref; i++ ) {
    miller::index<> M((int)write_columns[0]->ref[i],
                      (int)write_columns[1]->ref[i],
                      (int)write_columns[2]->ref[i]);
    lookup_dict[M] = i;
  }
  for(std::size_t i=0; i<miller_indices.size(); i++) {
    std::size_t i_mtz;
    lookup_dict_type::const_iterator
      ld_pos = lookup_dict.find(miller_indices[i]);
    if (ld_pos != lookup_dict.end()) {
      i_mtz = ld_pos->second;
      write_columns[3]->ref[i_mtz] = (float)data[i];
    }
    else {
      i_mtz = mtz->nref; //NUMBER_OF_EXISTING_INDICES
      i_mtz += 1;
      float adata[4] = {miller_indices[i][0],
                        miller_indices[i][1],
                        miller_indices[i][2],
                        data[i]}; //ADD miller_indices[i] TO THE MTZ OBJECT
      //ccp4_lwrefl(mtz, adata, write_columns, 4, i_mtz);
      safe_ccp4_lwrefl(adata, write_columns, 4, i_mtz);

      //FILL IN MISSING VALUES FOR ALL OTHER COLUMNS
    }
  }
}

void
MtzWriter::write(const std::string& filename)
{
  CMtz::MtzPut(mtz, const_cast<char*>(filename.c_str()));
}

}} // namespace iotbx::mtz
