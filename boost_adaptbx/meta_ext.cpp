/* The only purpose of this module is to provide access to the
   Boost.Python metaclass via meta_ext.empty.__class__.
   See also:
     boost/libs/python/doc/tutorial/doc/quickstart.txt, keyword injector
     boost.python.injector (boost/python.py)
 */

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>

namespace { struct empty {}; }

BOOST_PYTHON_MODULE(boost_python_meta_ext)
{
  boost::python::class_<empty>("empty", boost::python::no_init);
}
