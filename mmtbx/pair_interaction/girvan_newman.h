#ifndef MMTBX_PAIR_INTERACTION_GIRVAN_NEWMAN_H
#define MMTBX_PAIR_INTERACTION_GIRVAN_NEWMAN_H

#include <scitbx/array_family/shared.h>
#include <scitbx/error.h>
#include <vector>
#include <set>
#include <stack>
#include <queue>
#include <utility>
#include <algorithm>

namespace mmtbx { namespace pair_interaction {

namespace af = scitbx::af;

/*
  Girvan-Newman edge removal, as boost::betweenness_centrality_clustering
  with bc_clustering_threshold(threshold, g, normalize=false) on an
  undirected adjacency_list<setS, vecS>: while the largest edge betweenness
  (Brandes, halved for the undirected graph) is >= threshold, remove the
  first edge in insertion order that has it. Boost recomputes the
  betweenness of the whole graph after every removal; the betweenness of an
  edge depends only on its connected component, so this implementation
  recomputes only the component(s) of the removed edge, and never touches a
  component whose largest betweenness is already below the threshold. The
  removal sequence and the values are the same: within a component the
  sources are visited in increasing vertex order and the neighbours of a
  vertex in increasing vertex order (setS order), as boost does over the
  whole graph, and the other components contribute nothing.
*/
class girvan_newman
{
public:
  girvan_newman(
    int n_vertices,
    af::shared<int> const& edge_list,
    double threshold)
  :
    n(n_vertices), threshold(threshold)
  {
    SCITBX_ASSERT(n_vertices >= 0);
    SCITBX_ASSERT(edge_list.size() % 2 == 0);
    // deduplicate: setS ignores a repeated or reversed pair; a self-pair is
    // kept as an edge (never removed: it lies on no shortest path, its
    // betweenness is 0) but is not traversed
    std::set<std::pair<int,int> > seen;
    for(std::size_t k = 0; k < edge_list.size(); k += 2) {
      int i = edge_list[k], j = edge_list[k+1];
      SCITBX_ASSERT(i >= 0 && i < n && j >= 0 && j < n);
      std::pair<int,int> key(std::min(i,j), std::max(i,j));
      if(seen.find(key) != seen.end()) continue;
      seen.insert(key);
      edges.push_back(std::make_pair(i, j));
    }
    m = edges.size();
    alive.assign(m, true);
    bc.assign(m, 0.0);
    adj.assign(n, std::vector<std::pair<int,int> >());
    for(int e = 0; e < m; e++) {
      if(edges[e].first == edges[e].second) continue;
      adj[edges[e].first ].push_back(std::make_pair(edges[e].second, e));
      adj[edges[e].second].push_back(std::make_pair(edges[e].first,  e));
    }
    for(int v = 0; v < n; v++) std::sort(adj[v].begin(), adj[v].end());
    comp.assign(n, -1);
    for(int v = 0; v < n; v++) {
      if(comp[v] >= 0) continue;
      int c = new_component();
      label(v, c);
      compute(c);
    }
    run();
  }

  af::shared<int> remaining_edges() const
  {
    af::shared<int> result;
    for(int e = 0; e < m; e++) {
      if(!alive[e]) continue;
      result.push_back(edges[e].first);
      result.push_back(edges[e].second);
    }
    return result;
  }

private:
  int n, m;
  double threshold;
  std::vector<std::pair<int,int> > edges;
  std::vector<bool> alive;
  std::vector<double> bc;
  std::vector<std::vector<std::pair<int,int> > > adj; // (neighbour, edge)
  std::vector<int> comp;
  std::vector<std::vector<int> > comp_vertices;
  std::vector<bool> comp_final;
  std::vector<int> comp_max_edge;
  std::vector<double> comp_max_val;

  int new_component()
  {
    comp_vertices.push_back(std::vector<int>());
    comp_final.push_back(false);
    comp_max_edge.push_back(-1);
    comp_max_val.push_back(0.0);
    return comp_vertices.size() - 1;
  }

  // Label the component of v (over the alive edges) with c.
  void label(int v, int c)
  {
    comp_vertices[c].clear();
    std::vector<int> stack;
    stack.push_back(v);
    comp[v] = c;
    while(!stack.empty()) {
      int w = stack.back(); stack.pop_back();
      comp_vertices[c].push_back(w);
      for(std::size_t k = 0; k < adj[w].size(); k++) {
        if(!alive[adj[w][k].second]) continue;
        int u = adj[w][k].first;
        if(comp[u] == c) continue;
        comp[u] = c;
        stack.push_back(u);
      }
    }
    std::sort(comp_vertices[c].begin(), comp_vertices[c].end());
  }

  // Brandes edge betweenness of component c (halved: undirected), its
  // largest value and the first edge carrying it, and whether c is final.
  void compute(int c)
  {
    std::vector<int> const& vs = comp_vertices[c];
    for(std::size_t a = 0; a < vs.size(); a++) {
      int v = vs[a];
      for(std::size_t k = 0; k < adj[v].size(); k++) {
        if(alive[adj[v][k].second]) bc[adj[v][k].second] = 0.0;
      }
    }
    std::vector<std::size_t> path_count(n, 0); // boost: degree_size_type
    std::vector<int> distance(n, -1);
    std::vector<double> dependency(n, 0.0);
    std::vector<std::vector<std::pair<int,int> > > incoming(n); // (v, edge)
    for(std::size_t a = 0; a < vs.size(); a++) {
      int s = vs[a];
      for(std::size_t b = 0; b < vs.size(); b++) {
        int w = vs[b];
        incoming[w].clear();
        path_count[w] = 0;
        distance[w] = -1;
        dependency[w] = 0.0;
      }
      path_count[s] = 1;
      distance[s] = 0;
      std::stack<int> ordered;
      std::queue<int> q;
      q.push(s);
      while(!q.empty()) {
        int v = q.front(); q.pop();
        ordered.push(v);
        for(std::size_t k = 0; k < adj[v].size(); k++) {
          int e = adj[v][k].second;
          if(!alive[e]) continue;
          int w = adj[v][k].first;
          if(distance[w] < 0) {            // tree edge
            distance[w] = distance[v] + 1;
            path_count[w] = path_count[v];
            incoming[w].push_back(std::make_pair(v, e));
            q.push(w);
          }
          else if(distance[w] == distance[v] + 1) { // non-tree edge
            path_count[w] += path_count[v];
            incoming[w].push_back(std::make_pair(v, e));
          }
        }
      }
      while(!ordered.empty()) {
        int w = ordered.top(); ordered.pop();
        for(std::size_t k = 0; k < incoming[w].size(); k++) {
          int v = incoming[w][k].first, e = incoming[w][k].second;
          double factor = double(path_count[v]) / double(path_count[w]);
          factor *= (1.0 + dependency[w]);
          dependency[v] += factor;
          bc[e] += factor;
        }
      }
    }
    int best = -1;
    double best_val = 0.0;
    for(int e = 0; e < m; e++) {
      if(!alive[e] || comp[edges[e].first] != c) continue;
      bc[e] = bc[e] / 2.0;
      if(best < 0 || bc[e] > best_val) { best = e; best_val = bc[e]; }
    }
    comp_max_edge[c] = best;
    comp_max_val[c] = best_val;
    comp_final[c] = (best < 0) || (best_val < threshold);
  }

  void run()
  {
    while(true) {
      int best_c = -1;
      for(std::size_t c = 0; c < comp_vertices.size(); c++) {
        if(comp_final[c]) continue;
        if(best_c < 0 ||
           comp_max_val[c] > comp_max_val[best_c] ||
           (comp_max_val[c] == comp_max_val[best_c] &&
            comp_max_edge[c] < comp_max_edge[best_c])) {
          best_c = c;
        }
      }
      if(best_c < 0) break;
      int e = comp_max_edge[best_c];
      alive[e] = false;
      int i = edges[e].first, j = edges[e].second;
      // unlabel the old component, so that the piece of j is told apart
      // from the piece of i (both carried best_c)
      std::vector<int> const& old = comp_vertices[best_c];
      for(std::size_t k = 0; k < old.size(); k++) comp[old[k]] = -1;
      label(i, best_c);
      int c2 = -1;
      if(comp[j] < 0) {
        c2 = new_component();
        label(j, c2);
      }
      compute(best_c);
      if(c2 >= 0) compute(c2);
    }
  }
};

inline af::shared<int> girvan_newman_remaining_edges(
  int n_vertices,
  af::shared<int> const& edges,
  double threshold)
{
  girvan_newman gn(n_vertices, edges, threshold);
  return gn.remaining_edges();
}

}} // namespace mmtbx::pair_interaction

#endif // MMTBX_PAIR_INTERACTION_GIRVAN_NEWMAN_H
