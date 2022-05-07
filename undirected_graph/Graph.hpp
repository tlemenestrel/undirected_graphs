#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include "thrust/iterator/counting_iterator.h"
#include "thrust/iterator/transform_iterator.h"

#include <algorithm>
#include <vector>
#include <cassert>
#include <string>
#include <unordered_map>
#include <map>


#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

  // ===========================================================================
  // GRAPH
  // ===========================================================================

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 *
 * This code uses the proxy design pattern: the Node and Edge classes are 
 * used as proxies, while the internal_node and internal_edge structs contain
 * the underlyind data of those structures, allowing better data encapsulation.
 * 
 * The code has two functions, fetch_node() and fetch_edge(), that allow to 
 * access the data in the internal structs when we have a Node or Edge object.
 */

template <typename V, typename E>
class Graph: private totally_ordered <Graph<V, E>> {

  // Predeclare the internal struct for node and eges
  struct internal_edge; 
  struct internal_node;

public:

  // ===========================================================================
  // PUBLIC TYPE DEFINITIONS
  // ===========================================================================

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  struct NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  using node_value_type = V;
  using edge_value_type = E;

  // ===========================================================================
  // CONSTRUCTORS AND DESTRUCTOR
  // ===========================================================================

  /** Constructs an empty graph. */
  Graph() : nodes_(), edges_(), internal_nodes_(), internal_edges_(), 
  nodes_to_edge_(), next_node_id_(0), next_edge_id_(0), node_i2u_(), edge_i2u_() {};

  /** Default destructor */
  ~Graph() = default;

  // ===========================================================================
  // NODES
  // ===========================================================================

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered <Node> {
   
   public:
    /** Constructs an invalid node. */
    Node() {}

    /** Returns this node's graph pointer (cannot be modified). */
    const Graph* graph() const 
    {
      return graph_;
    }

    /** Returns this node's position (cannot be modified). */
    const Point& position() const 
    {
      return fetch_node().pt_;
    }

    /** Returns this node's position (can be modified). */
    Point& position() 
    {
      return fetch_node().pt_;
    }

    /** Returns this node's position by value */
    Point position_val() 
    {
      return fetch_node().pt_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const 
    {
      return fetch_node().get_node_idx_val();
    }

    IncidentIterator edge_begin() const
    {
      return IncidentIterator(graph_, get_incident_edges(), 0, node_id_);
    }

    IncidentIterator edge_end() const
    {
      return IncidentIterator(graph_, get_incident_edges(), degree(), node_id_);
    }

    /** Return the number of incident edges of a node. **/
    size_type degree() const 
    {
      std::vector<size_type> incident_edges = get_incident_edges();
      return incident_edges.size();
    }

    /** Return the vector containing the indexes of the incident edges of this node **/
    std::vector<size_type> get_incident_edges() const 
    {
      return this->graph_->adj_edges_.at(get_node_id());
    }

    /** Return this node's value */
    node_value_type& value()
    {
      return fetch_node().get_node_value();
    }

    /** Return this node's value (cannot be modified) */
    const node_value_type& value() const
    {
      return fetch_node().get_node_value();
    }

    /** Return this node's id */
    size_type get_node_id() const
    {
      return fetch_node().get_node_id();
    }

    /** Return this node's index */
    size_type& get_node_idx()
    {
      return fetch_node().get_node_idx();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const 
    {
      if(n.get_node_id() == get_node_id() && n.graph_ == graph_)
        return true;
      return false;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const 
    {
      if(n.get_node_id() > get_node_id() && n.graph_ == graph_)
        return true;
      else if (n.graph_ > graph_)
        return true;
      return false;
    }

    /** Helper method to return the appropriate element.
     * This loops over the elements until it finds the element with the
     * correct uid.
     */
    internal_node& fetch_node() const 
    {
        return this->graph_->internal_nodes_.at(node_id_);
    }

    private:
      friend class Graph;
      Graph* graph_; // Pointer to the graph of the node
      size_type node_id_; // Index of the node

      /** Construct a valid node using initalizer list. */
      Node(const Graph* graph, size_type id):
        graph_(const_cast<Graph*>(graph)), node_id_(id){}
  };

  // ===========================================================================
  // GRAPH - NODES METHOD
  // ===========================================================================

  /**
   * @brief Given a node iterator @a n, remove the node it points to from the graph
   * @param[in]  a Node Iterator pointing to the Node to be deleted
   * @returns return NodeIterator that points to the beginning of iterator 
   *
   * Complexity: O(num_nodes())
   */
  node_iterator remove_node(node_iterator n_it) 
  {
    Node n = (*n_it);
    remove_node(n);
    return node_begin();
  }

  /**
   * @brief Given a node @a n, remove it from the graph
   * @param[in]  a Node to be deleted
   * @returns return 0 if not deleted; 1 f deleted;
   *
   * Complexity: O(num_nodes())
   */
  size_type remove_node(const Node& n) 
  {
    if(!has_node(n)) 
    {
      return 0;
    }

    size_type n_id  = n.get_node_id();
    size_type n_idx = n.index();
    IncidentIterator it = n.edge_begin();

    // Remove adjacent edges

    while(n.degree() != 0) 
    {
    for (;it != n.edge_end();) 
    {
      Edge e = (*it);
      remove_edge(e);
      it = n.edge_begin();
    }
  }

    // Remove key from adjacent edges
    adj_edges_.erase(n_id);

    if(n_idx == (num_nodes()-1)) {
      node_i2u_.pop_back();
      return 1;
    }

    // Perform the swap
    std::iter_swap(node_i2u_.begin() + n_idx, node_i2u_.end() - 1);

    // Update the index of the node that was swapped
    auto n_old = nodes_.at((*(node_i2u_.begin() + n_idx)));
    n_old.get_node_idx() = n_idx;

    // Remove the old edge
    node_i2u_.pop_back();

    return 1;
  } 

  void default_initialize_nodes(int a) const
  {
    // Declare the iterators
    node_iterator niter_first = node_begin();

    // Default initialize values of all nodes to -1.
    for(; niter_first != node_end(); ++niter_first)
      (*niter_first).value() = a;
  }

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const 
  {
      return node_i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const 
  {
    return size();
  }

  Node add_node(const Point& position, const node_value_type& a = node_value_type()) 
  {
    // Declare a new node and add it to the vector of nodes
    internal_node new_internal_node = internal_node(position, a, num_nodes(), next_node_id_); 

    // Create the pair for the map of internal nodes
    std::pair<size_type, internal_node> key_val_node (next_node_id_, new_internal_node);
    internal_nodes_.insert(key_val_node);

    // Create the vector to store the future adjacent edges
    std::vector<size_type> adj_edges;
    std::pair<size_type, std::vector<size_type>> key_val_edge (next_node_id_, adj_edges);
    adj_edges_.insert(key_val_edge);

    // Declare new node and push it to vector
    node_type new_node = Node(this, next_node_id_); 
    nodes_.push_back(new_node);
    node_i2u_.push_back(next_node_id_);
    next_node_id_++;

    // Return the new node using the current graph and index
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const 
  {  
    if(this == n.graph_ && n.index() < size())
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const 
  {
    size_type n_id = node_i2u_.at(i);
    return Node(this, n_id);
  }

  std::string create_node_id(const Node&a, const Node&b) const
  {
    // Get the nodes ids
    std::string a_id = std::to_string(a.get_node_id());
    std::string b_id = std::to_string(b.get_node_id());

    // Create key
    std::string key = a_id + "-" + b_id;

    return key;
  }

  // ===========================================================================
  // EDGES
  // ===========================================================================
  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:

    /** Construct an invalid Edge. */
    Edge() {}

    /** Return this edge's value */
    edge_value_type& value()
    {
      return fetch_edge().get_edge_value();
    }
    
    /** Return this edge's value (cannot be modified) */
    const edge_value_type& value() const
    {
      return fetch_edge().get_edge_value();
    }
    
    /** Returns the first node of this Edge */
    Node node1() const 
    {
      return Node(this->graph_, fetch_edge().get_edge_n1());      
    }

    /** Return the second node of this Edge */
    Node node2() const 
    {
      return Node(this->graph_, fetch_edge().get_edge_n2());      
    }

    /** Returns the first node index of this Edge */
    size_type& node1_idx() const 
    {
      return fetch_edge().first_node_id_;      
    }

    /** Returns the first node index of this Edge */
    size_type node1_idx_val() const 
    {
      return fetch_edge().first_node_id_;      
    }
    /** Returns the second node index of this Edge */
    size_type& node2_idx() const 
    {
      return fetch_edge().second_node_id_;      
    }

    /** Returns the second node index of this Edge */
    size_type node2_idx_val() const 
    {
      return fetch_edge().second_node_id_;      
    }

    /** Returns the id of this Edge */
    size_type get_edge_id() const 
    {
      return fetch_edge().get_edge_id();
    }

    /** Returns the index of this Edge */
    size_type& get_edge_idx() const 
    {
      return fetch_edge().get_edge_idx();
    }

    /** Returns the index of this Edge */
    size_type get_edge_idx_val() const 
    {
      return fetch_edge().get_edge_idx_val();
    }

    /** Returns the length between the two nodes of this Edge */
    double length() const 
    {
      // Get the nodes coordinates
      Point p1 = node1().fetch_node().pt_;
      Point p2 = node2().fetch_node().pt_;

      // Compute distance
      float distance = sqrt(
        pow((p1.x - p2.x), 2) + 
        pow((p1.y - p2.y), 2) +
        pow((p1.z - p2.z), 2));

      return distance;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * First, check if both edges are from same graph. Then, check if 
     * their nodes match (in one direction or the other as graph is undirected)
     */
    bool operator==(const Edge& e) const 
    {
      if((this->graph_ == e.graph_) && ((fetch_edge().first_node_id_== 
        e.fetch_edge().first_node_id_ && fetch_edge().second_node_id_ == 
        e.fetch_edge().second_node_id_) || (fetch_edge().first_node_id_ == 
        e.fetch_edge().second_node_id_ && fetch_edge().second_node_id_
        == e.fetch_edge().second_node_id_)))
        return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const 
    {
      assert(e.graph_ != NULL and this->graph_ != NULL);

      // Get the min and max of nodes indexes for the two edges 
      size_type e_min = std::min(e.fetch_edge().first_node_id_, 
        e.fetch_edge().second_node_id_);
      size_type e_max = std::max(e.fetch_edge().first_node_id_, 
        e.fetch_edge().second_node_id_);

      size_type this_min = std::min(fetch_edge().first_node_id_, 
        fetch_edge().second_node_id_);
      size_type this_max = std::max(fetch_edge().first_node_id_, 
        fetch_edge().second_node_id_);

      if ((graph_ < e.graph_) ||((graph_ == e.graph_) && ((this_min < e_min) ||
       ((this_min == e_min) && (this_max < e_max))))) 
      {
        return true;
    }
    return false;
    }

    // Fetching function to get the underlying data of an Edge
    internal_edge& fetch_edge() const 
    {
      return this->graph_->internal_edges_.at(edge_id_);
    }

    private:
      // Allow Graph to access Edge's private member data and functions.
      friend class Graph;

      // Edge's private data attributes
      Graph* graph_;
      size_type edge_id_;

      //Construct a valid Edge using initalizer list.
      Edge(const Graph* graph, size_type edge_id):
        graph_(const_cast<Graph*>(graph)),
        edge_id_(edge_id){};
  };

  // ===========================================================================
  // GRAPH - EDGES METHOD
  // ===========================================================================

  /**
   * @brief Given an edge iterator, remove the edge that it points to
   * @param[in,out] e_it   An edge iterator pointing to the edge to be deleted
   * @returns return EdgeIterator object to the new @a graph_.edge_begin()
   *
   * @post If (has_edge(@a a, @a b) == true) {old @a num_edges_ == new @a num_edges_ + 1}
   *  else { old @a num_edges_ == new @a num_edges_ }
   *
   * Complexity: O(num nodes() + num edges())
   *
   */
  edge_iterator remove_edge(edge_iterator e_it) 
  {
    Edge e = (*e_it);
    remove_edge(e);
    return edge_begin();
  }

  /**
   * @brief Given two nodes, remove the connecting edge between them from a graph
   * @param[in,out] a   A node of the edge to be deleted
   * @param[in,out] b   Another node of the edge to be deleted
   * @returns return 0 if not deleted; 1 if successfully deleted
   *
   * @post If (has_edge(@a a, @a b) == true) {old @a num_edges_ == new @a num_edges_ + 1}
   *  else { old @a num_edges_ == new @a num_edges_ }
   *
   * Complexity: O(num nodes() + num edges())
   *
   */
  size_type remove_edge(const Node& a, const Node& b) 
  {
    if(has_edge(a,b)) 
    {
      Edge e = find_edge(a, b);
      remove_edge(e);
      return 1;
    }
    return 0;
  }

  /**
   * @brief Given an edge, remove the connecting edge between them from a graph
   * @param[in,out] e   An edge to be deleted
   * @returns return 0 if not deleted; 1 if successfully deleted
   *
   * @post If (has_edge(@a a, @a b) == true) {old @a num_edges_ == new @a num_edges_ + 1}
   *  else { old @a num_edges_ == new @a num_edges_ }
   *
   * Complexity: O(num nodes() + num edges())
   *
   */
  size_type remove_edge(const Edge& e) 
  {
    // Get the data of the edge
    Node n1          = e.node1();
    Node n2          = e.node2();
    size_type n1_id  = n1.get_node_id();
    size_type n2_id  = n2.get_node_id();
    size_type e_id   = e.get_edge_id();
    size_type e_idx   = e.get_edge_idx();

    // If the edges does not exist
    if (!has_edge(n1, n2)) {
      return 0;
    }

    // If the edge exists
    else 
    {
    // 1) Remove edge from nodes to edge  
    std::string key1 = create_node_id(n1, n2);
    std::string key2 = create_node_id(n2, n1);
    nodes_to_edge_.erase(key1);
    nodes_to_edge_.erase(key2);
    // 3) Remove the edge from adjacent edges for each node
    // Node 1
    std::vector<size_type>& vec1 = adj_edges_.at(n1_id);
    std::vector<size_type>::iterator position_1 = std::find(
      vec1.begin(), 
      vec1.end(), 
      e_id);
    if (position_1 != vec1.end()) { 
      vec1.erase(position_1);
    }

    // Node 2
    std::vector<size_type>& vec2 = adj_edges_.at(n2_id);
    std::vector<size_type>::iterator position_2 = std::find(
      vec2.begin(), 
      vec2.end(), 
      e_id);
    if (position_2 != vec2.end()) {
      vec2.erase(position_2);
    }

    if(e_idx == (num_edges()-1)) {
      edge_i2u_.pop_back();
      return 1;
    }

    // Perform the swap
    std::iter_swap(edge_i2u_.begin() + e_idx, edge_i2u_.end() - 1);
    // Update the index of the node that was swapped
    auto e_old = edges_.at((*(edge_i2u_.begin() + e_idx)));
    e_old.get_edge_idx() = e_idx;
    // Remove the old edge
    edge_i2u_.pop_back();
    return 1;
  }
}

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const 
  {
    return edge_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const 
  {
    size_type e_id = edge_i2u_.at(i);
    return Edge(this, e_id);
  }

  Edge find_edge(const Node& a, const Node& b) 
  {
    std::string key1 = create_node_id(a, b);
    std::string key2 = create_node_id(b, a);

    if(nodes_to_edge_.find(key1) == nodes_to_edge_.end()) 
    {
      size_type edge_id = nodes_to_edge_.find(key2)->second;
      return Edge(this, edge_id);
    }
    else if(nodes_to_edge_.find(key2) == nodes_to_edge_.end()) 
    {
      size_type edge_id = nodes_to_edge_.find(key1)->second;
      return Edge(this, edge_id);
    }
    assert(false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) 
  {
    std::string key1 = create_node_id(a, b);
    std::string key2 = create_node_id(b, a);
    if(nodes_to_edge_.find(key1) == nodes_to_edge_.end() &&
      nodes_to_edge_.find(key2) == nodes_to_edge_.end()) {
        return false;
      }
    else
      return true;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& c = edge_value_type()) 
  {
    // Check if the edge already exists. If yes, return it
    if(has_edge(a, b)) 
      return find_edge(a, b);

    // Create a key, val for the unordered map
    std::string key = create_node_id(a, b);
    std::pair<std::string, size_type> key_val (key, next_edge_id_);
    nodes_to_edge_.insert(key_val);

    // Create a new edge to be added and the pair for the map of internal edges
    internal_edge new_internal_edge = internal_edge(a.node_id_, b.node_id_, c, num_edges(), next_edge_id_);
    std::pair<size_type, internal_edge> key_val_2 (next_edge_id_, new_internal_edge);
    internal_edges_.insert(key_val_2);

    // Declare new edge 
    edge_type new_edge = Edge(this, next_edge_id_); 

    // Add edges to adjacency list of each nodes
    adj_edges_.at(a.get_node_id()).push_back(next_edge_id_);
    adj_edges_.at(b.get_node_id()).push_back(next_edge_id_);

    // Push it to vector of edges
    edges_.push_back(new_edge);
    edge_i2u_.push_back(next_edge_id_);

    next_edge_id_++;

    // Return the new Edge
    return new_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() 
  {
    nodes_.clear(); edges_.clear(); internal_nodes_.clear();internal_edges_.clear(); 
    nodes_to_edge_.clear(); node_i2u_.clear(); edge_i2u_.clear();
    next_edge_id_ = 0;
    next_node_id_ = 0;
  }

  // ===========================================================================
  // NODE ITERATOR
  // ===========================================================================

  /* Helper struct to get the node at an index */
  struct IndexToNode 
  {
    IndexToNode(const Graph* graph) : graph_(const_cast<Graph*>(graph)){};
    Node operator() (size_type node_idx){return graph_->node(node_idx);}
   private:
    Graph* graph_;
  };

  struct NodeIterator:thrust::transform_iterator<IndexToNode,thrust::counting_iterator<size_type>,node_type> 
  {
  using super_t = thrust::transform_iterator<IndexToNode, thrust::counting_iterator<size_type>, node_type>;

  // Default (invalid) constructor.
  NodeIterator(){};

  // KLUDGE conversion constructor
  NodeIterator(const super_t& ti) : super_t{ti} {};

  NodeIterator(const Graph* graph, const size_type nodeiter_idx = 0)
    : super_t{
      thrust::counting_iterator<size_type>(nodeiter_idx),
      IndexToNode(graph)} {};
  };

  // ===========================================================================
  // GRAPH - NODE ITERATOR METHODS
  // ===========================================================================

    NodeIterator node_begin() const
    {
      return NodeIterator(this, 0);
    }

    NodeIterator node_end() const
    {
      return NodeIterator(this, num_nodes());
   }

  // ===========================================================================
  // INCIDENT ITERATOR
  // ===========================================================================

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {};

    Edge operator*() const 
    {
      Edge e = Edge(graph_ptr_inc_itr_, incident_edges_.at(edge_iter_incident_idx_));

      if(e.node1_idx() == this->node_idx_) 
      {
        return e;
      }
      else
      {
        size_type e1_idx = e.node1_idx_val();
        size_type e2_idx = e.node2_idx_val();

        // Perform the swap
        e.node1_idx() = e2_idx;
        e.node2_idx() = e1_idx;
        Edge e = Edge(graph_ptr_inc_itr_, incident_edges_.at(edge_iter_incident_idx_));
        return e;
      }
      assert(false);
    }

    IncidentIterator& operator++() 
    {
      edge_iter_incident_idx_++;
      return (*this);
    }

    bool operator==(const IncidentIterator& i) const 
    {
      if(this->graph_ptr_inc_itr_ == i.graph_ptr_inc_itr_ 
        && this->incident_edges_ == i.incident_edges_
        && this->edge_iter_incident_idx_ == i.edge_iter_incident_idx_)
        return true;
      return false;
    }

    //Defines inequality between two iterators
    bool operator!=(const IncidentIterator& i) const
    {
        // Reuse the == operator
        return !(i == (*this));
    }

    // Get incident edges vector of an IncidentIterator
    std::vector<size_type>& get_incident_edges_vec() 
    {
      return incident_edges_;
    }

   private:
    friend class Graph;
    friend class Node;
    Graph* graph_ptr_inc_itr_;
    std::vector<size_type> incident_edges_;
    size_type edge_iter_incident_idx_;
    size_type node_idx_;

    IncidentIterator(
      const Graph* graph, 
      std::vector<size_type> incident_edges,
      size_type edge_iter_incident_idx,
      size_type node_idx)
      :
      graph_ptr_inc_itr_(const_cast<Graph*>(graph)),
      incident_edges_(incident_edges),
      edge_iter_incident_idx_(edge_iter_incident_idx),
      node_idx_(node_idx)
      {};
  };

  // ===========================================================================
  // EDGE ITERATOR
  // ===========================================================================

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    EdgeIterator& operator++() 
    {
      edge_iter_idx_++;
      return (*this);
    }

    bool operator==(const EdgeIterator& n) const
    {
      if(this->edge_iter_idx_ == n.edge_iter_idx_ 
        && this->graph_iter_ == n.graph_iter_)
        return true;
      return false;
    }

    bool operator!=(const EdgeIterator& e) const
    {
        // Reuse the == operator
        return !(e == (*this));
    }

    Edge operator*() const
    {
      size_type edge_id = graph_iter_->edge_i2u_.at(edge_iter_idx_);
      return Edge(graph_iter_, edge_id);
    }

   private:
    friend class Graph;
    size_type edge_iter_idx_;
    Graph* graph_iter_; // Pointer to the graph of the edge

    EdgeIterator(Graph* graph, size_type edge_iter_idx):
    edge_iter_idx_(edge_iter_idx),
    graph_iter_(graph) {};

  };

  // ===========================================================================
  // GRAPH - EDGE ITERATOR METHODS
  // ===========================================================================

    EdgeIterator edge_begin()
    {
      return EdgeIterator(this, 0);
    }

    EdgeIterator edge_end()
    {
      return EdgeIterator(this, num_edges());
    }

private:

  // ===========================================================================
  // INTERNAL NODE
  // ===========================================================================
  
    struct internal_node 
    {
    // Data attributes of the internal_node struct: a Point object and unique id
    Point pt_;
    node_value_type value_;
    size_type node_idx_;
    size_type node_id_;

    /** Return this node's value */
    node_value_type& get_node_value()
    {
      return value_;
    }

    /** Return this node's index */
    size_type& get_node_idx() 
    {
      return node_idx_;
    }

    /** Return this node's index */
    size_type get_node_idx_val()
    {
      return node_idx_;
    }

    /** Return this node's index */
    size_type get_node_id()
    {
      return node_id_;
    }
    // Default constructor
    internal_node(){};

    // Parameterized Constructor (with value)
    internal_node(Point position, node_value_type value, size_type node_idx, size_type node_id): 
    pt_(position), 
    value_(value),
    node_idx_(node_idx),
    node_id_(node_id)
    {};
  };

  // ===========================================================================
  // INTERNAL EDGE
  // ===========================================================================

    struct internal_edge 
    {
    // Data attributes of the internal_edge struct: a unique edge idx
    // the ids of the nodes it is connecting
    size_type first_node_id_;
    size_type second_node_id_;
    edge_value_type value_;
    size_type edge_idx_;
    size_type edge_id_;

    /** Return this edge's value */
    edge_value_type& get_edge_value()
    {
      return value_;
    }

    /** Return this edge's first node index */
    size_type& get_edge_n1()
    {
      return first_node_id_;
    }

    /** Return this edge's second node index */
    size_type& get_edge_n2()
    {
      return second_node_id_;
    }

    /** Return this edge's index */
    size_type& get_edge_idx()
    {
      return edge_idx_;
    }

    /** Return this edge's index */
    size_type get_edge_idx_val()
    {
      return edge_idx_;
    }

    /** Return this edge's id */
    size_type get_edge_id()
    {
      return edge_id_;
    }

    // Default constructor
    internal_edge() {};

    // Parameterized Constructor
    internal_edge(size_type first_id, size_type second_id, edge_value_type value, size_type edge_idx, size_type edge_id): 
    first_node_id_(first_id), 
    second_node_id_(second_id),
    value_(value),
    edge_idx_(edge_idx),
    edge_id_(edge_id)
    {};
  };

  // ===========================================================================
  // GRAPH DATA ATTRIBUTES
  // ===========================================================================

  // Vectors to store the Nodes and Edges in the Graph
  std::vector<node_type> nodes_;
  std::vector<edge_type> edges_;

  // Maps to store the Internal Nodes and Edges in the Graph
  std::map<size_type, internal_node> internal_nodes_;
  std::map<size_type, internal_edge> internal_edges_;

  // Map to store the adjacent edges of nodes
  std::map<size_type, std::vector<size_type>> adj_edges_;

  // Unordered map to store the Internal Nodes and Edges in the Graph
  std::unordered_map<std::string, size_type> nodes_to_edge_;

  size_type next_node_id_;
  size_type next_edge_id_;

  std::vector<size_type> node_i2u_;
  std::vector<size_type> edge_i2u_;
};

#endif // CME212_GRAPH_HPP
