/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <fstream>
#include <chrono>
#include <thread>

#include "thrust/for_each.h"
#include "thrust/iterator/transform_iterator.h"
#include "thrust/system/omp/execution_policy.h"

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "CME212/BoundingBox.hpp"

#include "Graph.hpp"
#include "SpaceSearcher.hpp"

// Gravity in meters/sec^2; initialize K and damping const c
static constexpr double grav = 9.81;
static constexpr double K = 100;
static constexpr double c = 1;

/** Custom structure of data to store with Nodes */
struct NodeData 
{
  Point vel;       
  double mass;     
  Point position;  
  NodeData() : vel(0), mass(1), position(0) {};
};

/** Custom structure of data to store with Edges */
struct EdgeData
{
  double K;
  double L;
  EdgeData() : K(100.0), L(0.0) {};
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = unsigned;

/* Functor to update nodes positions */
struct FunctorNodePosition  
{
    FunctorNodePosition(double t_dif) : d_(t_dif){};
    void operator()(Node n){ n.position() += n.value().vel * d_; }
    double d_;
}; 

/* Functor to update nodes velocity */
template <typename F>
struct FunctorNodeVelocity  
{
  FunctorNodeVelocity(double t, double t_dif, F& force)
    : f_(force), t_(t), d_(t_dif){};
    void operator()(Node n){n.value().vel+=f_(n,t_) * (d_/n.value().mass);}
    F f_; double t_; double d_; 
}; 

template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) 
{
  // Update the position of the nodes
  thrust::for_each(
    thrust::omp::par, 
    g.node_begin(), 
    g.node_end(),
    FunctorNodePosition(dt));
  constraint(g,t);
  // Update the velocity of the nodes
  thrust::for_each(
    thrust::omp::par, 
    g.node_begin(), 
    g.node_end(), 
    FunctorNodeVelocity<F>(t,dt,force));
  return t + dt;
}

class Force 
{
  public:
    Force() {}
    virtual Point operator()(Node n, double t) const 
    {
      (void) n;
      (void) t;
      return Point(0, 0, 0);
    }
    virtual ~Force() {}
};

class GravityForce : public Force 
{
  public:
    GravityForce() {};
    Point operator()(Node n, double t) const 
    {
      (void) t;
      Point grav_force = Point(0,0,-n.value().mass*grav);
      return grav_force;
    }
};

class MassSpringForce : public Force 
{
  public:
    MassSpringForce() {};
    Point operator()(Node n, double t) const 
    {
      (void) t;
      Point f_spring(0,0,0);
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) 
      {
        auto e1 = (*it);
        f_spring += -e1.value().K*(n.position() - e1.node2().position()) * (e1.length()- e1.value().L) / e1.length();
      }
      return f_spring;
    }
};

class DampingForce : public Force 
{
  public:
    DampingForce() {};
    Point operator()(Node n, double t) const 
    {
      (void) t;
      Point f_damp = -c * n.value().vel;
      return f_damp;
    }
};

struct CombinedForce 
{
  std::vector<const Force*> force_vec_;
  CombinedForce(std::vector<const Force*> force_vec): force_vec_(force_vec) {};
  Point operator()(Node n, double t) const {
    (void) t;
    Point f_comb(0,0,0);
    for (unsigned int i = 0; i < force_vec_.size(); i++) 
      f_comb += (*(force_vec_[i]))(n,t);
    return f_comb;
  }
};

CombinedForce make_combined_force(
  const Force& force_1, 
  const Force& force_2, 
  const Force& force_3 = Force()) 
{
  std::vector<const Force*> force_vec;

  force_vec.push_back(&force_1);
  force_vec.push_back(&force_2);
  force_vec.push_back(&force_3);

  return CombinedForce(force_vec);
}

/** Force function object for HW2 #1. */
struct Problem1Force 
{
  template <typename NODE>
  Point operator()(NODE n, double t) 
  {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) 
      return Point(0,0,0);
    Point f_spring(0,0,0);
    Point grav_force = Point(0,0,-n.value().mass*grav);

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = (*it);
      f_spring += -e.value().K * (n.position()-e.node2().position()) *
                                (e.length()-e.value().L)/e.length();
    }
    (void) t;
    return (grav_force + f_spring);
  }
};

class Constraint 
{
  public:
    Constraint() {};
    virtual void operator()(GraphType& g, double t) const 
    {
      (void) g;
      (void) t;
    }
    virtual ~Constraint() {}
};

class PinConstraint : public Constraint 
{
  public:
    PinConstraint() {};
    void operator()(GraphType& g, double t) const 
    {
      (void) t;
      // Fix points
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        if ((*it).value().position == Point(0,0,0) || (*it).value().position == Point(1,0,0)) 
        {
          (*it).value().vel = Point(0,0,0);
          (*it).position() = (*it).value().position;
        }
      }
    }
};

class PlaneConstraint : public Constraint 
{
  public:
    PlaneConstraint() {};
    void operator()(GraphType& g, double t) const {
      (void) t;
      for (auto it = g.node_begin(); it != g.node_end(); ++it) 
      {
        if (dot((*it).position(), Point(0,0,1)) < -0.75) 
        {
          (*it).position().elem[2]=-0.75;
          (*it).value().vel.elem[2]= 0;
        }
      }
    }
};

class SphereConstraint : public Constraint 
{
  public:
    SphereConstraint() {};
    void operator()(GraphType& g, double t) const 
    {
      (void) t;
      Point center(0.5,0.5,-0.5);
      double radius = 0.15;
      for (auto it = g.node_begin(); it != g.node_end(); ++it) 
      {
        if (norm((*it).position() - center) < radius) 
        {
          Point diff_vec = (*it).position() - center;
          double dist_c = norm(diff_vec);

          (*it).position() = center + diff_vec * radius / dist_c;
          (*it).value().vel -= dot((*it).value().vel, diff_vec/dist_c) * diff_vec / dist_c;
        }
      }
    }
};

class TearConstraint : public Constraint 
{
  public:
    TearConstraint() {};
    void operator()(GraphType& g, double t) const 
    {
      (void) t;
      Point center(0.5,0.5,-0.5);
      double radius = 0.15;
      // Remove nodes under constraint
      auto it = g.node_begin();
      while (it != g.node_end()) {
        if (norm((*it).position() - center) < radius) 
          g.remove_node((*it));
        else 
          ++it;
      }
    }
};

struct CombinedConstraints 
{
  std::vector<const Constraint*> c_vec_;
  CombinedConstraints(std::vector<const Constraint*> constr_vec): c_vec_(constr_vec) {}
  void operator()(GraphType& graph, double t) const 
  {
    (void) t;
    for (unsigned int i = 0; i < c_vec_.size(); i++) 
      (*(c_vec_[i]))(graph, t);
  }
};

CombinedConstraints make_combined_constraints(
  const Constraint& constr_1, 
  const Constraint& constr_2,    
  const Constraint& constr_3 = Constraint()) 
{
  std::vector<const Constraint*> c_vec;
  c_vec.push_back(&constr_1);
  c_vec.push_back(&constr_2);
  c_vec.push_back(&constr_3);
  return CombinedConstraints(c_vec);
}

struct ReduceVel 
{
    ReduceVel(Node& n, double radius): n_(n), radius_(radius){};
    void operator()(Node n)
    {
        Point point_r = n_.position() - n.position();
        double norm_l2 = normSq(point_r);
        if (n_ != n && norm_l2 < radius_) 
            n_.value().vel -= (dot(point_r, n_.value().vel) / norm_l2) * point_r;
    }
    Node n_;
    double radius_;
};

Box3D InterBBx(Box3D& box_tiny, Box3D& box_large)
{
    Point tiny_max  = box_tiny.max();
    Point tiny_min  = box_tiny.min();

    Point large_max = box_large.max();
    Point large_min = box_large.min();

    Point max_result = tiny_max;
    Point min_result = tiny_min;

    for(Point::size_type i = 0; i < tiny_min.size(); ++i)
    {
        if(tiny_max[i] > large_max[i])
            max_result[i] = large_max[i];
        if(tiny_min[i] < large_min[i])
            min_result[i] = large_min[i];
    }

    for(Point::size_type i = 0; i < tiny_min.size(); ++i)
    {
        if(max_result[i] < large_min[i])
            max_result[i] = large_min[i];
        if(min_result[i] > large_max[i])
            min_result[i] = large_max[i];
    }

    return Box3D(min_result, max_result);
}

struct GetInf 
{
    GetInf(SpaceSearcher<Node>& ss):ss_(ss){};
    void operator()(Node node1)
    {
        const Point& c = node1.position();
        double r = std::numeric_limits<double>::max();

        for (auto eit = node1.edge_begin(); eit != node1.edge_end(); ++eit)
            r = std::min(r, normSq((*eit).node2().position() - c));

        r *= 0.9;
        Point up  = c - sqrt(r);
        Point low = c + sqrt(r);

        Box3D tiny_bb(low, up);
        Box3D large_bb = ss_.bounding_box();
        Box3D bb_final = InterBBx(tiny_bb, large_bb);

        thrust::for_each(ss_.begin(bb_final), ss_.end(bb_final), ReduceVel(node1, r));
    }
    private:
      SpaceSearcher<Node>& ss_;
};

struct SelfCollisionConstraint : public Constraint 
{
    void operator()(GraphType& graph, double t) 
    {
        // Silence compiler
        (void) t;
        auto node_to_position = [](const Node& n) { return n.position();};

        Box3D large_bb = Box3D(
          thrust::make_transform_iterator(
            graph.node_begin(), 
            node_to_position),
          thrust::make_transform_iterator(
            graph.node_end(), 
            node_to_position));

        Point low_extend = large_bb.min();
        Point up_extend  = large_bb.max();
        for(size_type i = 0; i < low_extend.size(); ++i)
        {
            low_extend[i] = -abs(low_extend[i]) * 1.5;
            up_extend[i]  =  abs(up_extend[i])  * 1.5;
        }
        large_bb = Box3D(low_extend, up_extend);

        SpaceSearcher<Node> search(large_bb, graph.node_begin(), graph.node_end(), node_to_position);
        thrust::for_each(graph.node_begin(), graph.node_end(), GetInf(search));
    }
};
