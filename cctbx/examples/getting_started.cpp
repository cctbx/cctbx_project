#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/space_group_type.h>
#include <scitbx/array_family/simple_io.h>

int main()
{
  using cctbx::uctbx::unit_cell;
  using cctbx::sgtbx::space_group_type;
  using cctbx::sgtbx::space_group;
  using std::cout;
  using std::endl;

  unit_cell ucell(scitbx::af::double6(11,12,13,90,100,90));
  cout << "unit cell: " << ucell.parameters().const_ref() << endl;

  space_group_type sg_type("C 2");
  cout << "space group: " << sg_type.lookup_symbol() << endl;

  space_group const& sg = sg_type.group();
  CCTBX_ASSERT(sg.is_compatible_unit_cell(ucell));

  cout << "symmetry operations: " << endl;
  for (std::size_t i = 0; i < sg.order_z(); i++) {
    cout << "  " << sg(i).as_xyz() << endl;
  }

  return 0;
}
