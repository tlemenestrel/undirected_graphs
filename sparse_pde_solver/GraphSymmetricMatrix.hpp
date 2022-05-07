/**
 * @file GraphSymmetricMatrix.hpp
 * Implimentation file for treating the Graph as a MTL Matrix
 */

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include <fstream>
#include "CME212/Color.hpp"
#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"

using GraphType = Graph<char,char>;  //<  DUMMY Placeholder
using NodeType  = typename GraphType::node_type;

/* Check if a node is a boundary node. */
bool boundary(const NodeType& n) 
{
  // Get the position of the node
  const Point p = n.position();

  // 1st condition
  if (norm_inf(p) == 1)
    return true;

  Point p1 = p - Point( 0.6,  0.6, 0);
  Point p2 = p - Point(-0.6,  0.6, 0);
  Point p3 = p - Point( 0.6, -0.6, 0);
  Point p4 = p - Point(-0.6, -0.6, 0);

  // 2nd condition
  if(norm_inf(p1) < 0.2 or norm_inf(p2) < 0.2 or norm_inf(p3) < 0.2 or norm_inf(p4) < 0.2)
    return true;

  Point p1_bb = Point(-0.6, -0.2, -1);
  Point p2_bb = Point( 0.6,  0.2,  1);
  Box3D bb = Box3D(p1_bb, p2_bb);

  // 3rd condition
  if(bb.contains(p))
    return true;
  return false;
} 

/** A color functor. */
struct NodeColor
{
  NodeColor() {};
  CME212::Color operator()(const NodeType& n)
  {
    double nm = norm(n.position()); double x = n.position().z;
    return CME212::Color::make_heat(std::abs(x) / nm);
  }
};

/* A position functor */
struct NodePosition
{
  NodePosition(mtl::vec::dense_vector<double> x) : x_(x){}

  Point& operator()(NodeType& n)
  {
    size_t i = n.index();
    Point& p = n.position();
    p.z = x_[i];
    return p;
  }
  mtl::vec::dense_vector<double> x_;
};

/** Remove all the nodes in graph @a g whose posiiton is within Box3D @a bb.
 * @param[in,out] g  The Graph to remove nodes from
 * @param[in]    bb  The BoundingBox, all nodes inside this box will be removed
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const Box3D& bb) 
{
  for(auto i = g.node_begin(); i != g.node_end(); ++i)
  {
    auto n = *i;
    if(bb.contains(n.position()))
      g.remove_node(*i);
  }
}

/** Function g **/
double g(const Point& p) 
{
  // Define the points for the 2nd condition
  Point p1 = p - Point( 0.6,  0.6, 0);
  Point p2 = p - Point(-0.6,  0.6, 0);
  Point p3 = p - Point( 0.6, -0.6, 0);
  Point p4 = p - Point(-0.6, -0.6, 0);

  // Define the points and bounding box for the 3rd condition
  Point p1_bb = Point(-0.6, -0.2, -1);
  Point p2_bb = Point( 0.6,  0.2,  1);
  Box3D bb = Box3D(p1_bb, p2_bb);

  if (norm_inf(p) == 1)
    return 0.0;
  else if(norm_inf(p1) < 0.2 or norm_inf(p2) < 0.2 or norm_inf(p3) < 0.2 or norm_inf(p4) < 0.2)
    return -0.2;
  else if(bb.contains(p))
    return 1.0;
  return(1.0);
}

/* Function f. */
double f(Point& p)
{
  return 5 * cos(norm_1(p));
}

  // ===========================================================================
  // GRAPH SYMMETRIC MATRIC
  // ===========================================================================

class GraphSymmetricMatrix 
{
	public:
		GraphSymmetricMatrix(GraphType* graph):graph_(graph){};

  /* Get an element of the matrix. */
  double element(size_t i, size_t j) 
  {
    auto ni = graph_->node(i);
    auto nj = graph_->node(j);

    // Check 1st condition
    if(i == j and boundary(ni))
      return 1.0;

    // Check 2nd condition
    else if(i != j and (boundary(ni) or boundary(nj)))
      return 0.0;

    // Check 3rd condition
    else {
      if(i == j)
        return -double(ni.degree());
      if(graph_->has_edge(ni, nj))
        return 1.0;
      else
        return 0.0;
    }
  }

  /** Helper function to perform multiplication. Allows for delayed
  * evaluation of results.
  * Assign::apply(a, b) resolves to an assignment operation such as
  a*  += b, a -= b, or a = b. @pre @a size(v) == size(w) */
  template <typename VectorIn, typename VectorOut, typename Assign> 
  void mult(const VectorIn& v, VectorOut& w, Assign) const 
  {
    size_t s = graph_->size();

    for(size_t i = 0; i < s; ++i)
    {
      double curr = 0;
      for (size_t j = indptr_[i]; j < indptr_[i+1]; ++j)
      {
        curr += v[indices_[j]] * elements_[j];
      }
      Assign::apply(w[i], curr);
    }
  }

  /* * Matvec forwards to MTL ’s lazy m at _ c ve c _ mu l t ip l i er operator */
  template <typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
  operator *(const Vector& v) const{
    return {*this, v};
  }

  /** Make sparse matrix. **/
  void make_sparse_mat() 
  {
    indptr_.push_back(0);
    for(size_t i = 0; i < num_rows(); ++i)
    {
      for (size_t j = 0; j < num_cols(); ++j) 
      {
        if(element(i, j) != 0) 
        {
          elements_.push_back(element(i, j));
          indices_.push_back(j);
        }
      }
      indptr_.push_back(indices_.size());
    }
  }

  // Getter methods for the attributes of the matrix
  size_t num_rows() const {return graph_->size();}
  size_t num_cols() const {return graph_->size();}
  size_t m_size() const {return graph_->size() * graph_->size();}

  private:
    GraphType* graph_;
    std::vector<size_t> indptr_;
    std::vector<size_t> indices_;
    std::vector<double> elements_;
};

/** The number of rows in the matrix . */
inline std::size_t num_rows(const GraphSymmetricMatrix& A){
  return A.num_rows();
}

/** The number of columns in the matrix . */
inline std::size_t num_cols(const GraphSymmetricMatrix& A){
  return A.num_cols();
}

/** The number of elements in the matrix . */
inline std::size_t size(const GraphSymmetricMatrix& A){
  return A.m_size();
}
