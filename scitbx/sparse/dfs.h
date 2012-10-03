#ifndef SCITBX_SPARSE_DFS_H
#define SCITBX_SPARSE_DFS_H

#include <boost/tuple/tuple.hpp>
#include <stack>

namespace scitbx { namespace sparse {

/** A DFS through the graph associated to a sparse matrix, starting
from the nonzero elements of a vector.
*/
template<class Matrix>
class depth_first_search
{
public:
  typedef typename Matrix::value_type value_type;
  typedef typename Matrix::column_type column_type;
  typedef typename Matrix::index_type index_type;
  typedef typename Matrix::const_row_iterator const_row_iterator;

private:
  enum colour_type {white, gray, black};
  std::vector<colour_type> colour;

  typedef typename std::vector<index_type>::iterator row_idx_iter;

public:
    // Construct a DFS to visit the elements of a m x n matrix
    depth_first_search(index_type m, index_type n)
    : colour(std::max(m,n), white)
  {}

  // Perform the DFS through the graph of M, starting from the nonzeros of d
  template<class Visitor>
  void operator()(const Matrix& M, const column_type &d, Visitor& vis);
};

template<class Matrix>
template<class Visitor>
void depth_first_search<Matrix>::operator()(const Matrix& M,
                                            const column_type &d,
                                            Visitor& vis)
{
  std::vector<index_type> marked; // grayed and blackened vertices

  vis.dfs_started();

  for (const_row_iterator i_d=d.begin(); i_d != d.end(); i_d++) {
    /* DFS from nonzero of d */
    std::stack<boost::tuple<index_type,
                            const_row_iterator,
                            const_row_iterator> > stack;
    index_type k = vis.permute_rhs( i_d.index() );
    if (colour[k] != white) continue;
    vis.dfs_started_from_vertex(k);
    if (vis.dfs_shall_cut_tree_rooted_at(k)) continue;
    const_row_iterator col_iter = M.col(k).begin(),
                       col_end  = M.col(k).end();
    stack.push(boost::make_tuple(k, col_iter, col_end));
    while (!stack.empty()) {
      index_type l;
      boost::tuples::tie(l, col_iter, col_end) = stack.top();
      stack.pop();
      while (col_iter != col_end) {
        index_type k = vis.permute( col_iter.index() );
        if (colour[k] == white) {
          colour[k] = gray; marked.push_back(k);
          vis.dfs_found_tree_edge(l, k);
          if (vis.dfs_shall_cut_tree_edge(l, k)) continue;
          stack.push(boost::make_tuple(l, ++col_iter, col_end));
          l = k;
          col_iter = M.col(l).begin();
          col_end  = M.col(l).end();
        }
        else {
          ++col_iter;
        }
      }
      colour[l] = black; marked.push_back(l);
      vis.dfs_finished_vertex(l);
    }
  }
  /* Reset colour array for the next depth-first search.
    That's an important trick to keep the total number of ops
    proportional to the number of non-zeroes: if we were to reset
    the entire vector colour, we would run in O(n)
  */
  for (row_idx_iter p = marked.begin(); p != marked.end(); p++) {
    colour[*p] = white;
  }
}

}} // scitbx::sparse

#endif
