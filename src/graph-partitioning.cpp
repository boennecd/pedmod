#include <vector>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <list>
#include <limits>
#include <string>
#include <iostream>

namespace {

class vertex;

/// a weighted directed edge used to keep track of neighbors
struct weighted_edge {
  vertex const *v;
  double weight;

  weighted_edge(vertex const *v, double weight): v(v), weight(weight) { }

  inline operator vertex const*() const noexcept {
    return v;
  }
};

/// simple vertex class
class vertex {
private:
  using const_iterator = std::vector<weighted_edge>::const_iterator;

  /// vector of edges
  std::vector<weighted_edge> edges;

public:
  unsigned const id;
  double const weight;

  vertex(unsigned const id, double const weight = 1): id(id), weight(weight) { }

  inline const_iterator begin() const noexcept {
    return edges.begin();
  }

  inline const_iterator end() const noexcept {
    return edges.end();
  }

  inline size_t n_edges() const noexcept {
    return edges.size();
  }

  inline void add_edge(vertex *other, double const weight){
    if(other != this and other){
      edges       .emplace_back(other, weight);
      other->edges.emplace_back(this , weight);
    }
  }
};

/// simple edge class
struct edge {
  vertex const * x,
               * y;

  edge(): x(nullptr), y(nullptr) { }
  edge(vertex const * x, vertex const * y):
    x(x), y(y) { }
};

/// check if x is almost equal y
inline bool is_almost_equal(double const x, double const y){
  constexpr double const eps = 1000 * std::numeric_limits<double>::epsilon();
  return std::abs(x - y) < (std::abs(y) + eps) * eps;
}

/***
 * creates a vector with vertices from an edge list. The ids should be in the
 * range of 0,...,n-1 where n is the largest id. The latter three arguments can
 * be used to provide weights for the vertices which have weights that are not
 * equal to one.
 */
std::vector<vertex> create_vertices
  (int const *from, int const *to, int const n_edges,
   int const *id, double const *vertex_weight, int const n_weights,
   double const *edge_weight){
  // find the number greatest id
  int const * f_max = std::max_element(from, from + n_edges),
            * t_max = std::max_element(to  , to   + n_edges);
  if(!f_max or *f_max < 0L)
    throw std::invalid_argument("invalid from list");
  if(!t_max or *t_max < 0L)
    throw std::invalid_argument("invalid to list");
  int const n_obs = std::max(*f_max, *t_max) + 1L;

  // create map with the weights
  std::unordered_map<int, double> vertex_weight_map;
  for(int i = 0; i < n_weights; ++i, ++id, ++vertex_weight)
    vertex_weight_map.insert({ *id, *vertex_weight });

  // it is important first to create the vertices and then the edges as the
  // edges use pointers and thus will not be valid if the std::vector for what
  // ever reason resizes
  std::vector<vertex> out;
  out.reserve(n_obs);
  for(int i = 0; i < n_obs; ++i){
    auto vertex_weight_map_entry = vertex_weight_map.find(i);
    double const w =
      vertex_weight_map_entry != vertex_weight_map.end() ?
      vertex_weight_map_entry->second : 1;
    out.emplace_back(i, w);
  }

  // add the edges
  for(int i = 0; i < n_edges; ++i, ++from, ++to, ++edge_weight){
    // checks
    if(!from or *from > n_obs or *from < 0)
      throw std::invalid_argument("invalid from");
    if(!to or *to > n_obs or *to < 0)
      throw std::invalid_argument("invalid to");

    if(*from < n_obs and *to < n_obs)
      out[*from].add_edge(&out[*to], *edge_weight);
  }

  return out;
}

/// represents a block (2-connected component / biconnected component)
struct block {
  // perhaps we should not use unordered_set here
  /// the vertices in the block
  std::unordered_set<vertex const*> vertices;
  /// the cut vertices in the block
  std::unordered_set<vertex const*> cut_verices;

  inline bool has_vertex(vertex const *vertex){
    return vertices.count(vertex);
  }
};

/// block cut tree like structure with blocks and the edges that connect them
class block_cut_tree {
public:
  struct block_edge {
    block const * v1,
                * v2;
  };

  /// the blocks that defines the tree
  std::vector<block> const &blocks;
  /// the edges between the blocks
  std::vector<block_edge> block_edges;
  /// the map from a cut point to the block_edges
  std::unordered_map<vertex const *, std::vector<size_t> >
    cut_point_edges;

  block_cut_tree(std::vector<block> const &blocks):
    blocks(blocks) {
    //
    // create a map from cut vertices to the blocks and create a map from the
    // blocks to their indices
    std::unordered_map<vertex const *, std::deque<block const *> >
      cut_vertex_to_blocks;
    for(block const &b : blocks){
      for(vertex const *v : b.cut_verices){
        auto ele = cut_vertex_to_blocks.insert({ v, { } });
        ele.first->second.emplace_back(&b);
      }
    }

    // add the edges
    for(auto &ele : cut_vertex_to_blocks){
      std::deque<block const *> &bs = ele.second;
      size_t const n_blocks = bs.size();

      for(size_t i = 0; i < n_blocks; ++i)
        for(size_t j = i + 1; j < n_blocks; ++j){
          block_edges.push_back(block_edge { bs[i], bs[j] } );
          auto new_cut_point_edges = cut_point_edges.insert({ ele.first, { } });
          new_cut_point_edges.first->second.emplace_back(
              block_edges.size() - 1L);
        }
    }
  }
};

/**
 * struct to find the biconnected components using the method suggest
 * by Hopcroft and Tarjan (1973). See
 *   https://en.wikipedia.org/wiki/Biconnected_component
 */
class biconnected_components {
  /// class to hold extra information for vertices during the optimization
  class vertex_w_info  {
    using const_iterator = typename std::vector<vertex_w_info*>::const_iterator;

    /// vector of edges
    std::vector<vertex_w_info*> edges;

  public:
    /// the original vertex
    vertex const * const org_vertex;

    vertex_w_info(vertex const * const org_vertex):
      org_vertex(org_vertex) {
      if(org_vertex)
        edges.reserve(org_vertex->n_edges());
    }

    inline const_iterator begin() const noexcept {
      return edges.begin();
    }

    inline const_iterator end() const noexcept {
      return edges.end();
    }

    inline void add_edge(vertex_w_info *other){
      if(other != this and other){
        edges       .push_back(other);
        other->edges.push_back(this);
      }
    }

    bool visited = false;
    unsigned depth = 0,
             low   = 0;
    vertex_w_info const * parent = nullptr;
  };

  // perhaps we should not use unordered_map here
  /// maps from the old ids to 0,...,n-1
  std::unordered_map<vertex const*, unsigned> id_map;
  /// vertices with extra info needed for the method
  std::vector<vertex_w_info> vertex_info;
  /// queue used during the algorithm to keep track of the used edges
  std::deque<edge> edges;

  void check_if_cut_point(
      vertex_w_info &x, unsigned depth,
      std::vector<block> &bicon_comps,
      std::deque<edge> &edges){
    x.visited = true;
    x.depth = depth;
    x.low = depth;
    unsigned child_count(0);

    for(auto neighbor : x){
      if(!neighbor->visited){
        neighbor->parent = &x;
        edges.emplace_back(x.org_vertex, neighbor->org_vertex);
        check_if_cut_point(*neighbor, depth + 1L, bicon_comps,
                           edges);
        ++child_count;
        x.low = std::min(x.low, neighbor->low);

        bool is_cut_point =
          ( x.parent and neighbor->low >= x.depth) or
          (!x.parent and child_count > 1);
        if(!is_cut_point)
          continue;

        bicon_comps.emplace_back();
        auto &new_component = bicon_comps.back();
        while(edges.size() > 0 and
                (edges.back().x != x.org_vertex or
                   edges.back().y != neighbor->org_vertex)){
          auto &last_edge = edges.back();
          new_component.vertices.insert(last_edge.x);
          new_component.vertices.insert(last_edge.y);
          edges.pop_back();
        }
        if(edges.size() < 1)
          continue;
        auto &last_edge = edges.back();
        new_component.vertices.insert(last_edge.x);
        new_component.vertices.insert(last_edge.y);
        new_component.cut_verices.insert(x.org_vertex);
        edges.pop_back();

      } else if(neighbor != x.parent){
        x.low = std::min(x.low, neighbor->depth);
        if(neighbor->depth < x.depth)
          edges.emplace_back(x.org_vertex, neighbor->org_vertex);
      }
    }
  }

  // same as the previous method but this only returns the cuts points
  void check_if_cut_point_only_points(
      vertex_w_info &x, unsigned depth,
      std::unordered_set<vertex const *> &cut_points){
    x.visited = true;
    x.depth = depth;
    x.low = depth;
    unsigned child_count(0);

    for(auto neighbor : x){
      if(!neighbor->visited){
        neighbor->parent = &x;
        check_if_cut_point_only_points(
          *neighbor, depth + 1L, cut_points);
        ++child_count;
        x.low = std::min(x.low, neighbor->low);

        bool is_cut_point =
          ( x.parent and neighbor->low >= x.depth) or
          (!x.parent and child_count > 1);
        if(is_cut_point)
          cut_points.emplace(x.org_vertex);

      } else if(neighbor != x.parent)
        x.low = std::min(x.low, neighbor->depth);
    }
  }

public:
  // set the vertices but only add the edge if do_add_edge returns TRUE
  template<typename Tcontainer, typename Func>
  void set_vertices_if(Tcontainer const &x, Func do_add_edge){
    id_map.clear();
    vertex_info.clear();
    vertex_info.reserve(x.size());
    {
      unsigned id(0L);
      for(vertex const *xi : x){
        id_map.emplace(std::make_pair(xi, id++));
        vertex_info.emplace_back(xi);
      }
    }

    unsigned id(0L);
    auto vi = vertex_info.begin();
    for(vertex const *xi : x){
      for(weighted_edge const &neighbor : *xi){
        if(!do_add_edge(*neighbor.v, *xi))
          continue;

        unsigned const id_neighbor = id_map.at(neighbor);
        if(id_neighbor > id)
          vi->add_edge(&vertex_info[id_neighbor]);
      }
      ++id;
      ++vi;
    }
  }

  template<typename Tcontainer>
  void set_vertices(Tcontainer const &x){
    set_vertices_if(x, [](vertex const&,vertex const&) { return true; });
  }

  biconnected_components() = default;
  biconnected_components(std::vector<vertex> const &x){
    std::vector<vertex const *> tmp;
    tmp.reserve(x.size());
    for(vertex const &xi : x)
      tmp.emplace_back(&xi);
    set_vertices(tmp);
  }

  std::vector<block> get(){
    if(vertex_info.size() < 1)
      return { };

    std::vector<block> out;
    for(auto &v : vertex_info)
      if(!v.visited)
        check_if_cut_point(v, 0, out, edges);

    if(edges.size() > 0){
      out.emplace_back();
      auto &new_component = out.back().vertices;
      for(auto &e : edges){
        new_component.insert(e.x);
        new_component.insert(e.y);
      }
      edges.clear();
    }

    // have to add cut_vertices from other blocks. The algorithm only adds them
    // to one of the blocks
    std::unordered_set<vertex const *> cut_vertices;
    for(block const &b : out)
      cut_vertices.insert(b.cut_verices.begin(), b.cut_verices.end());
    for(block &b : out)
      for(vertex const * v : b.vertices)
        if(cut_vertices.count(v))
          b.cut_verices.insert(v);

    return out;
  }

  void get_cut_points(std::unordered_set<vertex const *> &cut_points){
    cut_points.clear();
    if(vertex_info.size() < 1)
      return;

    for(auto &v : vertex_info)
      if(!v.visited)
        check_if_cut_point_only_points(v, 0, cut_points);
  }
};

/**
 * constructs an approximate maximally balanced connected partition using the
 * method suggested by Chlebíková (1996). All member functions assume that the
 * graph is connected to start with.
 */
template<typename Tstream>
class max_balanced_partition {
  /**
   * struct to hold the weights needed for the algorithm. The class computes
   * the weights once and avoids repeated computation of the same quantities.
   */
  struct block_info {
    /// holds the total weight in some direction from a cut point
    struct directed_edge_w_weight {
      /// the block that is pointed to
      block_info const *other_block;
      /// weight of vertices in direction other_block less the common cut point
      double const weight;
    };

    /// the block which this object holds information for
    block const &b;
    /// the total weight of all vertices in the block including cut points
    double const inner_weight;

    /// map from cut points to directed_edge_w_weight objects
    std::unordered_map<vertex const *, std::vector<directed_edge_w_weight> >
      edges_to_other_blocks;

    block_info(block const &b):
    b(b),
    inner_weight(([&]() -> double {
      double out(0.);
      for(vertex const * v : b.vertices)
        out += v->weight;
      return out;
    })()) { }

    /**
     * sets the edges_to_other_blocks. This needs to be called for each
     * block_info object. The last argument is to avoid and inf recursion.
     */
    void set_edges_to_other_blocks(
        block_cut_tree const &bct,
        std::unordered_map<block const *, block_info> &b_info,
        std::unordered_set<vertex const *> &visited_cut_points){
      for(vertex const *v : b.cut_verices){
        // get or create the vector of directed edges with weights
        auto vec_of_w_edges_pair = edges_to_other_blocks.insert( { v, { } } );
        std::vector<directed_edge_w_weight> &vec_of_w_edges =
          vec_of_w_edges_pair.first->second;

        if(!visited_cut_points.insert(v).second)
          // the cut point is already being processed
          continue;

        // loop over the edges of this cut point. Add needed edges if we need to
        std::vector<size_t> const &edge_indices = bct.cut_point_edges.at(v);
        for(size_t edge_idx : edge_indices){
          // get the pointer to the other block
          block_cut_tree::block_edge const &edge = bct.block_edges[edge_idx];
          block const *other = edge.v1 != &b ? edge.v1 : edge.v2;

          // only step further if the edge does not exists
          {
            // TODO: possibly done many times. Perhaps use a set instead?
            auto ptr_match = std::find_if(
              vec_of_w_edges.begin(), vec_of_w_edges.end(),
              [&](directed_edge_w_weight const &x) -> bool {
                return &x.other_block->b == other;
              });
            if(ptr_match != vec_of_w_edges.end())
              continue;
          }

          // add the edges to other block if needed
          block_info &other_info = b_info.at(other);
          other_info.set_edges_to_other_blocks(bct, b_info, visited_cut_points);

          // compute the weight and add the edge
          double weight(other_info.inner_weight - v->weight);
          for(auto const &other_edges : other_info.edges_to_other_blocks){
            if(other_edges.first != v)
              for(auto const &w_e : other_edges.second)
                weight += w_e.weight;
          }

          // add the new edge
          vec_of_w_edges.push_back({ &other_info, weight });
        }
      }
    }
  };

  block_cut_tree const &bct;
  std::unordered_map<block const *, block_info> b_info;
  /**
   * as a sanity check, we compute the total weight which we can tests against
   * later. If we fail to get this value in any block then something is wrong
   * (e.g. the graph is not connected).
   */
  double const total_weight;
  /// holds whether to print log info during the computation
  unsigned const trace;
  /// the stream used for tracing
  Tstream &Tout;

  /**
   * recursively add all the block_infos from the tree if everything is
   * connected.
   */
  static void add_block(block const &b, block_cut_tree const &tree,
                        std::unordered_map<block const *, block_info> &out){
    auto res = out.insert({ &b, block_info(b) });
    if(!res.second)
      // the element was already there
      return;

    // add connected blocks
    for(vertex const * v : b.cut_verices){
      std::vector<size_t> const &edges = tree.cut_point_edges.at(v);
      for(size_t edge_idx : edges){
        block_cut_tree::block_edge const &b_edge = tree.block_edges[edge_idx];
        block const * other_block = b_edge.v1 != &b ? b_edge.v1 : b_edge.v2;
        add_block(*other_block, tree, out);
      }
    }
  }

  biconnected_components bicon_comp;

  /// keeps track of the best found solution of the connected mcb problem
  struct solution {
    /// the set we grow
    std::unordered_set<vertex const *> v1;
    /// the value of the weights in the smallest set
    double criterion;
    /// the block which contains the partition
    block_info const * block_with_solution;
  };

  /// keeps track of vertices that we add or remove
  struct working_set {
    /// the set of vertices that are in the set
    std::unordered_set<vertex const *> set;
    /// weights of the vertices
    std::unordered_map<vertex const *, double> &vertex_weight;
    /// vertices in the block that are connected to the set but not in the set
    std::unordered_set<vertex const *> connected_vertices;
    /// sum of the weights
    double weight = 0;
    /// the block we are working with
    block const * b = nullptr;

    working_set(std::unordered_map<vertex const *, double> &vertex_weight):
      vertex_weight(vertex_weight) { }

    /// adds a vertex
    void add_vertex(vertex const *new_vertex){
      if(!set.insert(new_vertex).second)
        // we did not add the vertex as it was already there
        return;

      // update the weight
      weight += vertex_weight[new_vertex];

      // remove the vertex from the connected vertices
      {
        auto con_match = connected_vertices.find(new_vertex);
        if(con_match != connected_vertices.end())
          connected_vertices.erase(con_match);
      }

      // update the connected vertices
      for(weighted_edge const &v : *new_vertex)
        if(!set.count(v) and b->vertices.count(v))
          connected_vertices.insert(v);
    }

    /**
     * removes a vertex. It is assumed that the vertex being removed is still
     * a connected to some vertex in the set.
     */
    void remove_vertex(vertex const *old_vertex){
      auto ptr = set.find(old_vertex);
      if(ptr == set.end())
        // was not in the set
        return;
      set.erase(ptr);

      // update the weight
      weight -= vertex_weight[old_vertex];

      // add the vertex the connected vertices
      connected_vertices.insert(old_vertex);

      // update the connected vertices
      for(vertex const * v : *old_vertex){
        auto ptr = connected_vertices.find(v);
        if(ptr == connected_vertices.end())
          // must already be in the set
          continue;

        // the question now is whether the vertex is still connected through
        // some other edge
        bool is_connected(false);
        for(vertex const *o : *v)
          if(set.count(o)){
            is_connected = true;
            break;
          }

        if(!is_connected)
          connected_vertices.erase(ptr);
      }
    }

    /// resets the set and changes the block
    void reset(block const * new_b){
      set.clear();
      connected_vertices.clear();
      weight = 0;
      b = new_b;
    }

    /***
     * resets the set and changes the block and set the vertices to the
     * complement of the passed set
     */
    void reset(block const * new_b,
               std::unordered_set<vertex const *> const &other_set){
      set.clear();
      connected_vertices.clear();
      weight = 0;
      b = new_b;
      for(vertex const * v : b->vertices)
        if(!other_set.count(v))
          add_vertex(v);
    }
  };

  /// sets the vertex weights and returns the vertex with the largest weight
  vertex const * set_vertex_weights
    (block_info const &b_info,
     std::unordered_map<vertex const *, double> &vertex_weight){
    vertex_weight.clear();
    block const &cur_block = b_info.b;

    // we also have to keep track of the largest element
    vertex const * max_weight_vertex = nullptr;
    double max_weight(-std::numeric_limits<double>::max());
    auto update_max_weight_vertex =
      [&](vertex const *new_ver, double const weight) -> void {
        if(!max_weight_vertex or weight > max_weight){
          max_weight = weight;
          max_weight_vertex = new_ver;
        }
      };

    double weight_sum(0.); // we compute the sum as a sanity check
    for(vertex const * v : cur_block.vertices){
      if(!cur_block.cut_verices.count(v)){
        // not a cut point. Simply use the weight of the vertex
        vertex_weight.insert({ v, v->weight });
        weight_sum += v->weight;
        update_max_weight_vertex(v, v->weight);
        continue;
      }

      // it is a cut point we have to sum up the cost
      double new_weight(v->weight);
      std::vector<typename  block_info::directed_edge_w_weight> const
        &edges_out = b_info.edges_to_other_blocks.at(v);
      for(auto const &e : edges_out)
        new_weight += e.weight;
      vertex_weight.insert({ v, new_weight });
      update_max_weight_vertex(v, new_weight);
      weight_sum += new_weight;
    }

    // perform a sanity check on the the sum of the weights
    if(!is_almost_equal(weight_sum, total_weight))
      throw std::runtime_error(
          "sum of weights in block is not equal to the total weight (" +
            std::to_string(weight_sum) + " and " +
            std::to_string(total_weight) + ")");

    return max_weight_vertex;
  }

  /// runs the algorithm suggested by Chlebíková
  solution get_connected_mcb_solution(){
    // declare the objects we will use. Start with the solution object
    solution res;
    res.criterion = -std::numeric_limits<double>::max();
    // the working v2 set which we remove elements from
    std::unordered_set<vertex const *> other_set;
    // the map from a vertex to its weight
    std::unordered_map<vertex const *, double> vertex_weight;
    // cut points used during the algorithm
    std::unordered_set<vertex const *> cut_points;
    // work set for the one we add edges to
    working_set v1_set(vertex_weight);

    auto add_solution_if_better =
      [&res, this](working_set &working_set,
                   block_info const * block_with_solution) -> void {
       double const criterion = std::min(working_set.weight,
                                         total_weight - working_set.weight);
       if(criterion > res.criterion){
         std::swap(res.v1, working_set.set);
         res.criterion = criterion;
         res.block_with_solution = block_with_solution;
       }
     };

    // used to report if trace is high enough
    auto report =
      [this](std::string const &prefix, double const set_size) -> void {
        Tout << prefix << ". Balance criterion is "
             << std::min(total_weight - set_size, set_size) << '\n';
    };

    // find the best solution
    for(auto const &b : b_info){
      // clean up from the previous iterations
      block const &cur_block = b.second.b;
      v1_set.reset(&cur_block);
      other_set.clear();
      vertex_weight.clear();

      // set the vertex weights and setup the set which we add to
      {
        vertex const * max_weight_vertex =
          set_vertex_weights(b.second, vertex_weight);

        // add the vertex with largest weight to the new set and check if we can
        // exit early
        v1_set.add_vertex(max_weight_vertex);
        if(v1_set.weight >= total_weight / 2 or cur_block.vertices.size() < 3){
          add_solution_if_better(v1_set, &b.second);
          if(trace > 1)
            report("Exited early", v1_set.weight);
          continue;
        }

        // create the list with the other set
        other_set.insert(cur_block.vertices.begin(), cur_block.vertices.end());
        other_set.erase(other_set.find(max_weight_vertex));
      }

      for(unsigned i = 0; i < cur_block.vertices.size() - 2L; ++i){
        // find the cut points. For this, we should only include "internal"
        // edges in the set where we remove edges from
        auto is_internal_edge = [&](vertex const &x, vertex const &y) -> bool {
          return other_set.count(&x) and other_set.count(&y);
        };
        bicon_comp.set_vertices_if(other_set, is_internal_edge);
        bicon_comp.get_cut_points(cut_points);

        // find the best solution among the intersection minus the cut points
        vertex const * min_weight_vertex = nullptr;
        double min_weight = std::numeric_limits<double>::max();
        for(vertex const * v : v1_set.connected_vertices)
          if(!cut_points.count(v)){
            double const w_weight = vertex_weight.at(v);
            if(
              // no vertex was set before
              !min_weight_vertex or
              // this one has a smaller weight
              w_weight < min_weight or
              // the weight is the same but to make sure that we get the same
              // on all platforms, then we take the one with the smallest id.
              // We are not guaranteed to get this otherwise because of the
              // arbitrary order of unordered_set
              (is_almost_equal(w_weight, min_weight) and
                 v->id < min_weight_vertex->id)){
              min_weight_vertex = v;
              min_weight = w_weight;
            }
          }

        if(!min_weight_vertex){
          // no vertex satisfied the criteria
          if(trace > 1)
            Tout << "Found no vertex to move\n";
          break;
        }

        if(min_weight >= total_weight - 2 * v1_set.weight){
          if(trace > 1)
            Tout << "Cannot improve balance criterion further\n";
          break;
        }

        // add the new vertex
        v1_set.add_vertex(min_weight_vertex);
        other_set.erase(other_set.find(min_weight_vertex));
      }

      if(trace > 1)
        report("Found split", v1_set.weight);
      add_solution_if_better(v1_set, &b.second);

      // check if we can exit early
      if(is_almost_equal(2 * res.criterion, total_weight)){
        if(trace > 0)
          report("Found perfect balanced connected partition early",
                 v1_set.weight);
        break;
      }
    }

    if(trace > 0)
      report("Found approximately balanced connected partition", res.criterion);

    return res;
  }

public:
  max_balanced_partition(block_cut_tree const &bct,
                         unsigned const trace, Tstream &Tout):
  bct(bct),
  b_info(([&]() -> std::unordered_map<block const *, block_info> {
    std::unordered_map<block const *, block_info> out;
    if(bct.blocks.size() > 0)
      add_block(bct.blocks[0], bct, out);
    return out;
  })()),
  total_weight(([&]() -> double {
    double out(0.);
    std::unordered_set<vertex const *> visited;
    for(auto const &b_i : b_info){
      for(vertex const * v : b_i.first->vertices)
        if(visited.insert(v).second)
          out += v->weight;
    }

    return out;
  })()),
  trace(trace), Tout(Tout) {
    std::unordered_set<vertex const *> visited_cut_points;
    for(auto &b_i : b_info){
      visited_cut_points.clear();
      b_i.second.set_edges_to_other_blocks(bct, b_info, visited_cut_points);
    }
  }

  /**
   * holds information of the solution of maximally balanced connected partition
   * problem.
   */
  struct mbcp_result {
    /// the criterion value
    double balance_criterion;
    /// vector with the removed edges
    std::vector<edge> removed_edges;
    /// vertices in the two sets in the block where the partition is
    std::vector<vertex const *> s1, s2;
  };

  /**
   * runs the approximate maximally balanced connected partition algorithm.
   * See
   *    https://cstheory.stackexchange.com/q/48864/62145
   */
  mbcp_result get(double const slack, unsigned const max_kl_it_inner,
                  unsigned const max_kl_it){
    // computes the object to return
    auto get_return_object = [](solution &res) -> mbcp_result {
      // we need to format the solution and return
      mbcp_result out;
      out.balance_criterion = res.criterion;
      if(res.v1.size() > 0){
        // we loop through the vertices in the block and add the edges in the cut
        block const &block_w_solution = res.block_with_solution->b;
        for(vertex const * v : block_w_solution.vertices){
          if(res.v1.count(v)){
            // the vertex is in the other set so we continue (do not add edges
            // twice)
            out.s1.emplace_back(v);
            continue;
          }
          out.s2.emplace_back(v);

          for(weighted_edge const &e : *v)
            if(res.v1.count(e))
              out.removed_edges.emplace_back(v, e);
        }
      }

      return out;
    };

    // possibly improve the result by finding a cut which is better
    solution res = get_connected_mcb_solution();
    if(slack <= 0)
      return get_return_object(res);

    if(!res.block_with_solution)
      throw std::runtime_error("block_with_solution not set");

    // the map from a vertex to its weight
    std::unordered_map<vertex const *, double> vertex_weight;
    set_vertex_weights(*res.block_with_solution, vertex_weight);

    // the block with the solution
    block const &b = res.block_with_solution->b;

    // the two sets we work with
    working_set v1_set(vertex_weight),
                v2_set(vertex_weight);
    v2_set.reset(&b, res.v1);
    v1_set.reset(&b, v2_set.set);

    // create a map with from vertices to the gain of moving them
    std::unordered_map<vertex const *, double> edge_gain;
    // updates the cost of a given vertex. The update_adjacent treats the
    // vertex as if it has just been moved and updates its neighbors
    // accordingly
    auto update_egde_cost =
      [&v1_set, &v2_set, &edge_gain]
      (vertex const * x, bool const update_adjacent) -> void {
        bool const is_in_v1 = v1_set.set.count(x);
        std::unordered_set<vertex const *>
          &current_set = is_in_v1 ? v1_set.set : v2_set.set,
          &other_set   = is_in_v1 ? v2_set.set : v1_set.set;

        double out(0);
        for(weighted_edge const &v : *x)
          if       (current_set.count(v)){
            // they are in the same set so moving the vertex is going to cost
            out -= v.weight;
            if(update_adjacent)
              // have to switch from plus to minus
              edge_gain[v.v] -= 2 * v.weight;
          } else if(other_set.count(v)){
            // they are not in the same set so moving the vertex is going to
            // yield a plus
            out += v.weight;
            if(update_adjacent)
              // have to switch from minus to plus
              edge_gain[v.v] += 2 * v.weight;
          }

        edge_gain[x] = out;
    };

    // insert the first values
    for(vertex const * v : b.vertices)
      update_egde_cost(v, false);

    // keeps track of the vertex, the gain, where the vertices is moved from
    // and to, and the balance criterion
    struct gain_score {
      vertex const *x = nullptr;
      double gain = 0;
      double balance_criterion = 0;
      working_set * moved_from = nullptr,
                  * moved_to   = nullptr;

      gain_score() = default;
      gain_score(vertex const *x, double const gain,
                 double const balance_criterion,
                 working_set * moved_from,
                 working_set * moved_to):
        x(x), gain(gain), balance_criterion(balance_criterion),
        moved_from(moved_from), moved_to(moved_to) { }
    };

    // queue with the moves we made
    std::deque<gain_score> scores;
    // set with the vertices we have used
    std::unordered_set<vertex const *> used_vertices;
    // cut points used during the algorithm
    std::unordered_set<vertex const *> cut_points;
    // minimum value for the balance criterion
    double const min_balance_criterion(total_weight / 2 - slack * total_weight);
    // the number of iterations to make in the inner part of the algorithm
    unsigned const n_it_take =
      std::min<unsigned>(max_kl_it_inner, b.vertices.size());
    // the balance criterion of the current solution
    double max_balance_crit = res.criterion;

    if(trace > 0)
      Tout << "Starting to reduce the cost of the cut in a block of size "
           << b.vertices.size() << ". The partition in the block consists of two sets of size "
           << v1_set.set.size() << ' ' << v2_set.set.size() << '\n';

    // perform the main loop
    for(unsigned i = 0; i < max_kl_it; ++i){
      scores.clear();
      used_vertices.clear();
      if(trace > 0)
        Tout << "Starting iteration " << i + 1 << '\n';

      for(unsigned j = 0; j < n_it_take; ++j){
        gain_score best_gain;

        auto find_best_move =
          [&](working_set &from, working_set &to) -> void {
            auto is_internal_edge =
              [&](vertex const &x, vertex const &y) -> bool {
                return from.set.count(&x) and from.set.count(&y);
            };
            bicon_comp.set_vertices_if(from.set, is_internal_edge);
            bicon_comp.get_cut_points(cut_points);

            for(vertex const * v : to.connected_vertices){
              if(used_vertices.count(v))
                continue;
              if(cut_points.count(v))
                continue;

              double const weight = vertex_weight.at(v),
                     balance_crit = std::min(from.weight - weight,
                                             to.weight   + weight);
              if(balance_crit < min_balance_criterion)
                continue;

              double const gain_v = edge_gain.at(v);
              if(
                // there was no solution before
                !best_gain.x or
                // the solution is better
                best_gain.gain < gain_v or
                // the solution is not better but the balance criterion is
                // better
                (is_almost_equal(best_gain.gain, gain_v) and
                      best_gain.balance_criterion < balance_crit) or
                // everything is the same except the id which is smaller. This
                // yields reproducible results despite using unordered_set
                (is_almost_equal(best_gain.gain, gain_v) and
                   is_almost_equal(best_gain.balance_criterion,
                                   balance_crit) and
                   v->id < best_gain.x->id))
                best_gain = { v, gain_v, balance_crit, &from, &to };
            }
          };

        // find the best move
        find_best_move(v1_set, v2_set);
        find_best_move(v2_set, v1_set);

        // perform the move if any
        if(!best_gain.x)
          // no move to perform
          break;

        best_gain.moved_to  ->add_vertex   (best_gain.x);
        best_gain.moved_from->remove_vertex(best_gain.x);
        update_egde_cost(best_gain.x, true);

        used_vertices.insert(best_gain.x);
        scores.emplace_back(std::move(best_gain));
      }

      if(trace > 0)
        Tout << "Found " << scores.size() << " vertices to move\n";
      if(scores.size() < 1)
        // found no swap
        break;

      // find best number of moves to keep
      unsigned max_idx(0);
      double  max_gain(0),
              gain_sum(0.);
      unsigned idx(0);
      for(gain_score const &s : scores){
        ++idx;
        gain_sum += s.gain;
        if(gain_sum > max_gain or
             (gain_sum == max_gain and
                s.balance_criterion > max_balance_crit)){
          max_gain = gain_sum;
          max_idx = idx;
          max_balance_crit = s.balance_criterion;
        }
      }

      if(trace > 0){
        if(max_idx > 0){
          Tout << "Keept " << max_idx << " moves with a gain of "
               << max_gain << " and a balance criterion of "
               << max_balance_crit << '\n';
        } else
          Tout << "Keept " << max_idx << " moves\n";
      }

      // undo moves
      unsigned const n_scores = scores.size();
      for(unsigned idx = max_idx; idx < n_scores; ++idx){
        gain_score const &s = scores.back();
        s.moved_to  ->remove_vertex(s.x);
        s.moved_from->add_vertex   (s.x);
        update_egde_cost(s.x, true);
        scores.pop_back();
      }

      if(max_idx < 1L)
        // no solution found
        break;
    }

    std::swap(res.v1, v1_set.set);
    res.criterion = std::min(v1_set.weight, v2_set.weight);

    return get_return_object(res);
  }
};

} // namespace

#include <Rcpp.h>

Rcpp::List biconnected_components_to_rcpp_list
  (std::vector<block> const &bicon_comps){

  Rcpp::List out(bicon_comps.size());
  for(size_t i = 0; i < bicon_comps.size(); ++i){
    auto &comp_i = bicon_comps[i];
    Rcpp::IntegerVector ele(comp_i.vertices.size());
    {
      size_t j(0);
      for(vertex const *comp_ele : comp_i.vertices)
        ele[j++] = comp_ele->id + 1L;
    }

    Rcpp::IntegerVector cut_verices(comp_i.cut_verices.size());
    {
      size_t j(0);
      for(vertex const *comp_ele : comp_i.cut_verices)
        cut_verices[j++] = comp_ele->id + 1L;
    }

    ele.attr("cut_verices") = cut_verices;
    out[i] = ele;
  }

  return out;
}

// [[Rcpp::export(.get_biconnected_components, rng = false)]]
Rcpp::List get_biconnected_components(
    Rcpp::IntegerVector const from, Rcpp::IntegerVector const to,
    Rcpp::IntegerVector const weights_ids, Rcpp::NumericVector const weights,
    Rcpp::NumericVector const edge_weights){
  if(from.size() != to.size())
    throw std::invalid_argument("size of from does not match size of to");
  if(edge_weights.size() != to.size())
    throw std::invalid_argument("size of edge_weights does not match size of to");
  if(weights_ids.size() != weights.size())
    throw std::invalid_argument("size of weights_ids does not match size of weights");

  std::vector<vertex> vertices = create_vertices(
    &from[0], &to[0], to.size(), &weights_ids[0], &weights[0],
    weights_ids.size(), &edge_weights[0]);
  auto bicon_comps = biconnected_components(vertices).get();
  return biconnected_components_to_rcpp_list(bicon_comps);
}

Rcpp::List block_cut_tree_to_rcpp_list
  (block const &b, block_cut_tree const &bct,
   std::unordered_set<vertex const*> &used_cut_points){
  Rcpp::IntegerVector ele(b.vertices.size());
  {
    size_t j(0);
    for(vertex const *comp_ele : b.vertices)
      ele[j++] = comp_ele->id + 1L;
  }

  Rcpp::IntegerVector cut_verices(b.cut_verices.size());
  {
    size_t j(0);
    for(vertex const *comp_ele : b.cut_verices)
      cut_verices[j++] = comp_ele->id + 1L;
  }

  Rcpp::List leafs;
  {
    for(vertex const *comp_ele : b.cut_verices)
      if(!used_cut_points.count(comp_ele)){
        used_cut_points.insert(comp_ele);
        auto &edges = bct.cut_point_edges.at(comp_ele);
        for(size_t edge_idx : edges){
          auto &edge = bct.block_edges[edge_idx];
          auto other = edge.v1 != &b ? edge.v1 : edge.v2;
          leafs.push_back(block_cut_tree_to_rcpp_list(
              *other, bct, used_cut_points));
        }
      }
  }

  return Rcpp::List::create(
    Rcpp::Named("vertices") = ele,
    Rcpp::Named("cut_vertices") = cut_verices,
    Rcpp::Named("leafs") = leafs);
}

// [[Rcpp::export(.get_block_cut_tree, rng = false)]]
Rcpp::List get_block_cut_tree(
    Rcpp::IntegerVector const from, Rcpp::IntegerVector const to,
    Rcpp::IntegerVector const weights_ids, Rcpp::NumericVector const weights,
    Rcpp::NumericVector const edge_weights){
  if(from.size() != to.size())
    throw std::invalid_argument("size of from does not match size of to");
  if(edge_weights.size() != to.size())
    throw std::invalid_argument("size of edge_weights does not match size of to");
  if(weights_ids.size() != weights.size())
    throw std::invalid_argument("size of weights_ids does not match size of weights");

  std::vector<vertex> vertices = create_vertices(
    &from[0], &to[0], to.size(), &weights_ids[0], &weights[0],
    weights_ids.size(), &edge_weights[0]);
  auto bicon_comps = biconnected_components(vertices).get();
  block_cut_tree bct(bicon_comps);

  std::unordered_set<vertex const*> used_cut_points;
  return block_cut_tree_to_rcpp_list(bicon_comps[0], bct, used_cut_points);
}

using max_balanced_partition_rcpp =
  max_balanced_partition<decltype(Rcpp::Rcout)>;

Rcpp::List mbcp_result_to_rcpp_list
  (max_balanced_partition_rcpp::mbcp_result const &res){
  Rcpp::NumericVector balance_criterion = { res.balance_criterion };
  Rcpp::IntegerMatrix removed_edges(res.removed_edges.size(), 2L);
  for(unsigned i = 0; i < res.removed_edges.size(); ++i){
    edge const &e = res.removed_edges[i];
    removed_edges(i, 0) = e.x->id + 1L;
    removed_edges(i, 1) = e.y->id + 1L;
  }

  auto get_set =
    [](std::vector<vertex const *> const &s) -> Rcpp::IntegerVector{
      Rcpp::IntegerVector out(s.size());
      for(unsigned i = 0; i < s.size(); ++i)
        out[i] = s[i]->id + 1L;
      return out;
  };

  return Rcpp::List::create(
    Rcpp::Named("balance_criterion") = balance_criterion,
    Rcpp::Named("removed_edges")     = removed_edges,
    Rcpp::Named("set_1")             = get_set(res.s1),
    Rcpp::Named("set_2")             = get_set(res.s2));
}

// [[Rcpp::export(.get_max_balanced_partition, rng = false)]]
Rcpp::List get_max_balanced_partition(
    Rcpp::IntegerVector const from, Rcpp::IntegerVector const to,
    Rcpp::IntegerVector const weights_ids, Rcpp::NumericVector const weights,
    Rcpp::NumericVector const edge_weights,
    double const slack, unsigned const max_kl_it_inner,
    unsigned const max_kl_it, unsigned const trace){
  if(from.size() != to.size())
    throw std::invalid_argument("size of from does not match size of to");
  if(edge_weights.size() != to.size())
    throw std::invalid_argument("size of edge_weights does not match size of to");
  if(weights_ids.size() != weights.size())
    throw std::invalid_argument("size of weights_ids does not match size of weights");
  if(slack >= .5)
    throw std::invalid_argument("invalid slack value");

  std::vector<vertex> vertices = create_vertices(
    &from[0], &to[0], to.size(), &weights_ids[0], &weights[0],
    weights_ids.size(), &edge_weights[0]);
  auto bicon_comps = biconnected_components(vertices).get();
  block_cut_tree bct(bicon_comps);
  max_balanced_partition_rcpp max_part(bct, trace, Rcpp::Rcout);
  max_balanced_partition_rcpp::mbcp_result const res =
    max_part.get(slack, max_kl_it_inner, max_kl_it);

  return mbcp_result_to_rcpp_list(res);
}
