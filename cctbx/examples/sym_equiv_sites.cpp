#include <cctbx/sgtbx/sym_equiv_sites.h>
#include <scitbx/array_family/simple_io.h>

int main()
{
  using namespace cctbx;
  using std::cout;
  using std::endl;

  uctbx::unit_cell unit_cell(scitbx::af::double6(11,12,13,90,100,90));
  cout << "unit cell: " << unit_cell.parameters().const_ref() << endl;

  sgtbx::space_group_type space_group_type("C 2");
  cout << "space group: " << space_group_type.lookup_symbol() << endl;
  cout << endl;

  sgtbx::space_group const& space_group = space_group_type.group();
  CCTBX_ASSERT(space_group.is_compatible_unit_cell(unit_cell));

  cout << "symmetry operations: " << endl;
  for (std::size_t i = 0; i < space_group.order_z(); i++) {
    cout << "  " << space_group(i).as_xyz() << endl;
  }
  cout << endl;

  for (double z=0.51; z < 0.54; z += 0.02) {

    fractional<> site(0.0, 0.13, z);

    sgtbx::site_symmetry site_symmetry(
      unit_cell,
      space_group,
      site,
      /*min_distance_sym_equiv*/ 0.5,
      /*assert_min_distance_sym_equiv*/ true);

    sgtbx::sym_equiv_sites<> sym_equiv_sites(site_symmetry);

    cout << "original_site:"
         << sym_equiv_sites.original_site().const_ref() << endl;
    cout << "special_op: "
         << sym_equiv_sites.special_op().as_xyz() << endl;
    cout << "is_special_position: "
         << sym_equiv_sites.is_special_position() << endl;

    af::const_ref<scitbx::vec3<double> >
      coordinates = sym_equiv_sites.coordinates().const_ref();
    for (std::size_t i = 0; i < coordinates.size(); i++) {
      cout << "coordinates[" << i << "] = "
           << coordinates[i].const_ref() << endl;
    }
    cout << endl;

  }

  return 0;
}
