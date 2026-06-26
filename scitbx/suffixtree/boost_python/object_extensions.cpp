#include <boost/python/object.hpp>
#include <boost/python/extract.hpp>

#include <boost/functional/hash.hpp>

#include <scitbx/suffixtree/boost_python/object_extensions.hpp>

namespace boost
{
namespace python
{

std::size_t
hash_value(const object& obj)
{
#if PY_MAJOR_VERSION >= 3
  typedef Py_hash_t hash_t;
#else // Py2
  typedef long hash_t;
#endif
  hash_t raw = extract<hash_t>( obj.attr( "__hash__" )() );
  return boost::hash_value( raw );
}

std::ostream&
operator <<(std::ostream& s, const object& obj)
{
  return s << ( extract< const char* >( obj.attr( "__str__" )() ) );
}

} // namespace python
} // namespace boost
