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
  return boost::hash_value( PyObject_Hash( obj.ptr() ) );
}

std::ostream&
operator <<(std::ostream& s, const object& obj)
{
  return s << ( extract< const char* >( obj.attr( "__str__" )() ) );
}

} // namespace python
} // namespace boost
