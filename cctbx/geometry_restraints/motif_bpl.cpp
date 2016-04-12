#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_arg.hpp>
#include <cctbx/geometry_restraints/motif.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct motif_atom_wrappers : boost::python::pickle_suite
  {
    typedef motif::atom w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.name,
        self.scattering_type,
        self.nonbonded_type,
        self.partial_charge
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("motif_atom", no_init)
        .def(init<
          const char*,
            optional<const char*, const char*, double> >((
              arg("name"),
              arg("scattering_type")="",
              arg("nonbonded_type")="",
              arg("partial_charge")=0)))
        .def_readwrite("name", &w_t::name)
        .def_readwrite("scattering_type", &w_t::scattering_type)
        .def_readwrite("nonbonded_type", &w_t::nonbonded_type)
        .def_readwrite("partial_charge", &w_t::partial_charge)
        .def_pickle(motif_atom_wrappers())
        ;
    }
  };

  struct motif_bond_wrappers : boost::python::pickle_suite
  {
    typedef motif::bond w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.atom_names,
        self.type,
        self.distance_ideal,
        self.weight,
        self.id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif_bond", no_init)
        .def(init<
          af::tiny<std::string, 2>,
            optional<const char*, double, double, const char*> >((
              arg("atom_names"),
              arg("type")="",
              arg("distance_ideal")=0,
              arg("weight")=0,
              arg("id")="")))
        .add_property("atom_names",
          make_getter(&w_t::atom_names, rbv()),
          make_setter(&w_t::atom_names, dcp()))
        .def_readwrite("type", &w_t::type)
        .def_readwrite("distance_ideal", &w_t::distance_ideal)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("id", &w_t::id)
        .def_pickle(motif_atom_wrappers())
        ;
    }
  };

  struct motif_angle_wrappers : boost::python::pickle_suite
  {
    typedef motif::angle w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.atom_names,
        self.angle_ideal,
        self.weight,
        self.id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif_angle", no_init)
        .def(init<
          af::tiny<std::string, 3>,
            optional<double, double, const char*> >((
              arg("atom_names"),
              arg("angle_ideal")=0,
              arg("weight")=0,
              arg("id")="")))
        .add_property("atom_names",
          make_getter(&w_t::atom_names, rbv()),
          make_setter(&w_t::atom_names, dcp()))
        .def_readwrite("angle_ideal", &w_t::angle_ideal)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("id", &w_t::id)
        .def_pickle(motif_angle_wrappers())
        ;
    }
  };

  struct motif_dihedral_wrappers : boost::python::pickle_suite
  {
    typedef motif::dihedral w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.atom_names,
        self.angle_ideal,
        self.weight,
        self.periodicity,
        self.id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif_dihedral", no_init)
        .def(init<
          af::tiny<std::string, 4>,
            optional<double, double, int, const char*> >((
              arg("atom_names"),
              arg("angle_ideal")=0,
              arg("weight")=0,
              arg("periodicity")=0,
              arg("id")="")))
        .add_property("atom_names",
          make_getter(&w_t::atom_names, rbv()),
          make_setter(&w_t::atom_names, dcp()))
        .def_readwrite("angle_ideal", &w_t::angle_ideal)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("periodicity", &w_t::periodicity)
        .def_readwrite("id", &w_t::id)
        .def_pickle(motif_dihedral_wrappers())
        ;
    }
  };

  struct motif_chirality_wrappers : boost::python::pickle_suite
  {
    typedef motif::chirality w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.atom_names,
        self.volume_sign,
        self.both_signs,
        self.volume_ideal,
        self.weight,
        self.id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif_chirality", no_init)
        .def(init<
          af::tiny<std::string, 4>,
            optional<const char*, bool, double, double, const char*> >((
              arg("atom_names"),
              arg("volume_sign")="",
              arg("both_signs")=false,
              arg("volume_ideal")=0,
              arg("weight")=0,
              arg("id")="")))
        .add_property("atom_names",
          make_getter(&w_t::atom_names, rbv()),
          make_setter(&w_t::atom_names, dcp()))
        .def_readwrite("volume_sign", &w_t::volume_sign)
        .def_readwrite("both_signs", &w_t::both_signs)
        .def_readwrite("volume_ideal", &w_t::volume_ideal)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("id", &w_t::id)
        .def_pickle(motif_chirality_wrappers())
        ;
    }
  };

  struct motif_planarity_wrappers : boost::python::pickle_suite
  {
    typedef motif::planarity w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.atom_names,
        self.weights,
        self.id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif_planarity", no_init)
        .def(init<>())
        .def(init<
          af::shared<std::string> const&,
          af::shared<double> const&,
            optional<const char*> >((
              arg("atom_names"), arg("weights"), arg("id")="")))
        .add_property("atom_names",
          make_getter(&w_t::atom_names, rbv()),
          make_setter(&w_t::atom_names, dcp()))
        .add_property("weights",
          make_getter(&w_t::weights, rbv()),
          make_setter(&w_t::weights, dcp()))
        .def_readwrite("id", &w_t::id)
        .def_pickle(motif_planarity_wrappers())
        ;
    }
  };

  // XXX move to scitbx::boost_python, also create+move array_from_list
  template <typename ArrayType>
  boost::python::list
  array_as_list(ArrayType const& array)
  {
    boost::python::list result;
    typedef typename ArrayType::const_iterator aci;
    aci ae = array.end();
    for(aci ai=array.begin();ai!=ae;ai++) {
      result.append(*ai);
    }
    return result;
  }

  struct motif_wrappers : boost::python::pickle_suite
  {
    typedef motif w_t;

#define CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(type, member) \
    static boost::python::list \
    member##_as_list(w_t const& self) \
    { \
      return array_as_list(self.member.const_ref()); \
    } \
\
    static void \
    set_##member(w_t& self, boost::python::object const& sequence) \
    { \
      unsigned size = boost::python::len(sequence); \
      self.member = af::shared<motif::type>(af::reserve(size)); \
      for(unsigned i=0;i<size;i++) { \
        boost::python::extract<motif::type> proxy(sequence[i]); \
        self.member.push_back(proxy()); \
      } \
    }

    CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(atom, atoms)
    CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(bond, bonds)
    CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(angle, angles)
    CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(dihedral, dihedrals)
    CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(chirality, chiralities)
    CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(planarity, planarities)

    static boost::python::tuple
    getstate(w_t const& self)
    {
      return boost::python::make_tuple(
         self.id,
         self.description,
         self.info,
         self.manipulation_ids,
         self.atoms,
         self.bonds,
         self.angles,
         self.dihedrals,
         self.chiralities,
         self.planarities
        );
    }

    static void
    setstate(w_t& self, boost::python::tuple state)
    {
      self.id =                boost::python::extract< std::string >                   (state[0]);
      self.description =       boost::python::extract<std::string>                     (state[1]);
      self.info =              boost::python::extract< af::shared<std::string> >       (state[2]);
      self.manipulation_ids =  boost::python::extract< af::shared<std::string> >       (state[3]);
      self.atoms =             boost::python::extract< af::shared<motif::atom> >       (state[4]);
      self.bonds =             boost::python::extract< af::shared<motif::bond> >       (state[5]);
      self.angles =            boost::python::extract< af::shared<motif::angle> >      (state[6]);
      self.dihedrals =         boost::python::extract< af::shared<motif::dihedral> >   (state[7]);
      self.chiralities =       boost::python::extract< af::shared<motif::chirality> >  (state[8]);
      self.planarities =       boost::python::extract< af::shared<motif::planarity> >  (state[9]);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif", no_init)
        .def(init<>())
        .def_readwrite("id", &w_t::id)
        .def_readwrite("description", &w_t::description)
        .add_property("info",
          make_getter(&w_t::info, rbv()),
          make_setter(&w_t::info, dcp()))
        .add_property("manipulation_ids",
          make_getter(&w_t::manipulation_ids, rbv()),
          make_setter(&w_t::manipulation_ids, dcp()))
        .def("atoms_as_list", atoms_as_list)
        .def("set_atoms", set_atoms)
        .def("bonds_as_list", bonds_as_list)
        .def("set_bonds", set_bonds)
        .def("angles_as_list", angles_as_list)
        .def("set_angles", set_angles)
        .def("dihedrals_as_list", dihedrals_as_list)
        .def("set_dihedrals", set_dihedrals)
        .def("chiralities_as_list", chiralities_as_list)
        .def("set_chiralities", set_chiralities)
        .def("planarities_as_list", planarities_as_list)
        .def("set_planarities", set_planarities)
        .def_pickle(motif_wrappers())
        ;
    }
  };

  struct motif_alteration_wrappers : boost::python::pickle_suite
  {
    typedef motif::alteration w_t;

    static std::string
    get_action(w_t const& self) { return self.action.description(); }

    static void
    set_action(w_t& self, std::string const& description)
    {
      self.action = w_t::action_type(description);
    }

    static std::string
    get_operand(w_t const& self) { return self.operand.description(); }

    static void
    set_operand(w_t& self, std::string const& description)
    {
      self.operand = w_t::operand_type(description);
    }

    static boost::python::list
    planarity_atom_actions_as_list(w_t const& self)
    {
      boost::python::list result;
      typedef w_t::action_type const* aci;
      aci ae = self.planarity_atom_actions.end();
      for(aci ai=self.planarity_atom_actions.begin();ai!=ae;ai++) {
        result.append(ai->description());
      }
      return result;
    }

    static void
    set_planarity_atom_actions(
      w_t& self, boost::python::object const& sequence)
    {
      self.planarity_atom_actions.clear(); // for the case of an exception
      unsigned size = boost::python::len(sequence);
      af::shared<w_t::action_type> new_actions((af::reserve(size)));
      for(unsigned i=0;i<size;i++) {
        boost::python::extract<const char*> proxy(sequence[i]);
        new_actions.push_back(w_t::action_type(proxy()));
      }
      self.planarity_atom_actions = new_actions;
    }

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.action, self.operand);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif_alteration", no_init)
        .def(init<optional<std::string const&, std::string const&> >((
          arg("action"), arg("operand"))))
        .add_property("action",
          make_function(get_action),
          make_function(set_action))
        .add_property("operand",
          make_function(get_operand),
          make_function(set_operand))
        .add_property("motif_ids",
          make_getter(&w_t::motif_ids, rbv()),
          make_setter(&w_t::motif_ids, dcp()))
        .def_readwrite("atom", &w_t::atom)
        .def_readwrite("motif_atom_name", &w_t::motif_atom_name)
        .def_readwrite("bond", &w_t::bond)
        .def_readwrite("angle", &w_t::angle)
        .def_readwrite("dihedral", &w_t::dihedral)
        .def_readwrite("chirality", &w_t::chirality)
        .def_readwrite("planarity", &w_t::planarity)
        .def_readwrite("planarity_motif_id", &w_t::planarity_motif_id)
        .def("planarity_atom_actions_as_list", planarity_atom_actions_as_list)
        .def("set_planarity_atom_actions", set_planarity_atom_actions)
        .def("change_partial_charge", &w_t::change_partial_charge)
        .def("set_change_partial_charge", &w_t::set_change_partial_charge,
          (arg("state")), return_self<>())
        .def("change_distance_ideal", &w_t::change_distance_ideal)
        .def("set_change_distance_ideal", &w_t::set_change_distance_ideal,
          (arg("state")), return_self<>())
        .def("change_weight", &w_t::change_weight)
        .def("set_change_weight", &w_t::set_change_weight,
          (arg("state")), return_self<>())
        .def("change_angle_ideal", &w_t::change_angle_ideal)
        .def("set_change_angle_ideal", &w_t::set_change_angle_ideal,
          (arg("state")), return_self<>())
        .def("change_periodicity", &w_t::change_periodicity)
        .def("set_change_periodicity", &w_t::set_change_periodicity,
          (arg("state")), return_self<>())
        .def("change_both_signs", &w_t::change_both_signs)
        .def("set_change_both_signs", &w_t::set_change_both_signs,
          (arg("state")), return_self<>())
        .def("change_volume_ideal", &w_t::change_volume_ideal)
        .def("set_change_volume_ideal", &w_t::set_change_volume_ideal,
          (arg("state")), return_self<>())
        .def_pickle(motif_alteration_wrappers())
        ;
    }
  };

  struct motif_manipulation_wrappers : boost::python::pickle_suite
  {
    typedef motif::manipulation w_t;

    CCTBX_GEOMETRY_RESTRAINTS_MOTIF_LIST_MEMBER(alteration, alterations)

      static boost::python::tuple
      getstate(w_t const& self)
    {
      return boost::python::make_tuple(
        self.id,
        self.description,
        self.info,
        self.alterations
        );
    }

    static void
      setstate(w_t& self, boost::python::tuple state)
    {
      self.id =          boost::python::extract< std::string >            (state[0]);
      self.description = boost::python::extract< std::string >            (state[1]);
      self.info =        boost::python::extract< af::shared<std::string> >(state[2]);
      self.alterations = boost::python::extract< af::shared<motif::alteration> >(state[3]);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      class_<w_t>("motif_manipulation", no_init)
        .def(init<>())
        .def_readwrite("id", &w_t::id)
        .def_readwrite("description", &w_t::description)
        .add_property("info",
          make_getter(&w_t::info, rbv()),
          make_setter(&w_t::info, dcp()))
        .def("alterations_as_list", alterations_as_list)
        .def("set_alterations", set_alterations)
        .def_pickle(motif_wrappers())
        ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    motif_atom_wrappers::wrap();
    motif_bond_wrappers::wrap();
    motif_angle_wrappers::wrap();
    motif_dihedral_wrappers::wrap();
    motif_chirality_wrappers::wrap();
    motif_planarity_wrappers::wrap();
    motif_wrappers::wrap();
    motif_alteration_wrappers::wrap();
    motif_manipulation_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_motif() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
