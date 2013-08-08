#include <ostream>
#include <boost/python/object_fwd.hpp>

namespace boost
{
namespace python
{

std::size_t hash_value(const object& obj);

std::ostream& operator <<(std::ostream& s, const object& obj);

} // namespace python
} // namespace boost
