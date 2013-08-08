#ifndef SUFFIXTREE_EXCEPTION_HPP_
#define SUFFIXTREE_EXCEPTION_HPP_

#include <exception>

namespace scitbx
{

namespace suffixtree
{

class bad_edge_type : public std::exception
{
  const char* what() const throw()
  {
    return "Invalid edge type for operation";
  }
};


class bad_state : public std::exception
{
  const char* what() const throw()
  {
    return "Incorrect state for operation";
  }
};


class bad_tree : public std::exception
{
  const char* what() const throw()
  {
    return "Invalid tree";
  }
};

} // namespace suffixtree
} // namespace scitbx

#endif // SUFFIXTREE_EXCEPTION_HPP_
