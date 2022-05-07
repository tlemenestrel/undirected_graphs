#include <gtest/gtest.h>
#include <fstream>

#include "CME212/Util.hpp"
#include "Graph.hpp"
#include "shortest_path.hpp"
//#include "subgraph.hpp"

class GraphPointFixture : public ::testing::Test {
 protected:
   //Define types
  using GraphType    = Graph<int, int>;
  using NodeType     = typename GraphType::node_type;
  using NodeIter     = typename GraphType::node_iterator;
  using EdgeType     = typename GraphType::edge_type;
  using EdgeIter     = typename GraphType::edge_iterator;
  using IncidentIter = typename GraphType::incident_iterator;

  //Set up Graph and Points
  GraphType graph;
  std::vector<Point> points;
  virtual void SetUp() {
    for(int i = 0; i < 10; i++)
      points.push_back(Point(i));
  }
};

// Test node iterator
TEST_F(GraphPointFixture, NodeIter)
{
  NodeType n1 = graph.add_node(points[0]);
  NodeType n2 = graph.add_node(points[1]);
  NodeType n3 = graph.add_node(points[2]);
  
  int iter = 0;
  NodeIter ni = graph.node_begin();

  EXPECT_EQ((*ni), n1) << " error in node iteration ";
  ++ni;
  EXPECT_EQ((*ni), n2) << " error in node iteration ";
  ++ni;
  EXPECT_EQ((*ni), n3) << " error in node iteration ";

  NodeIter ni2 = graph.node_begin();
  for(;ni2 != graph.node_end(); ++ni2)
    iter++;
  EXPECT_EQ(iter, 3)    << " error in node iteration ";
}

// Test edge iterator
TEST_F(GraphPointFixture, EdgeIter)
{
  NodeType n1 = graph.add_node(points[0]);
  NodeType n2 = graph.add_node(points[1]);
  NodeType n3 = graph.add_node(points[2]);

  EdgeType e1 = graph.add_edge(n1, n2);
  EdgeType e2 = graph.add_edge(n1, n3);
  EdgeType e3 = graph.add_edge(n2, n3);

  EdgeIter ei = graph.edge_begin();
  EXPECT_EQ(e1, (*ei)) << " error in edge iteration ";
  ASSERT_NE(e2, (*ei))  << " error in edge iteration ";
  ASSERT_NE(e3, (*ei))  << " error in edge iteration ";
  ++ei;
  EXPECT_EQ(e2, (*ei)) << " error in edge iteration ";
  ASSERT_NE(e1, (*ei))  << " error in edge iteration ";
  ASSERT_NE(e3, (*ei))  << " error in edge iteration ";
  ++ei;
  EXPECT_EQ(e3, (*ei)) << " error in edge iteration ";
  ASSERT_NE(e1, (*ei))  << " error in edge iteration ";
  ASSERT_NE(e2, (*ei))  << " error in edge iteration ";
}

// Test degree function
TEST_F(GraphPointFixture, Degree)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);

  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n0, n3);

  EXPECT_EQ(0, n4.degree()) << "n4 degree is 0";
  EXPECT_EQ(3, n0.degree()) << "n1 degree is 3";
}

// Test get incident edges
TEST_F(GraphPointFixture, GetIncidentEdges)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);
  NodeType n5 = graph.add_node(points[5]);
  NodeType n6 = graph.add_node(points[6]);

  EdgeType e0 = graph.add_edge(n0, n1);
  graph.add_edge(n1, n3);
  graph.add_edge(n1, n4);

  EdgeType e4 = graph.add_edge(n0, n2);
  graph.add_edge(n2, n5);
  graph.add_edge(n2, n6);

  std::vector<size_type> vec  = n0.get_incident_edges();
  std::vector<size_type> vec2 = n1.get_incident_edges();
  std::vector<size_type> vec3 = n5.get_incident_edges();

  std::vector<size_type>vec4;
  vec4.push_back(e0.get_edge_id());
  vec4.push_back(e4.get_edge_id());

  EXPECT_EQ(2, vec.size()) << "error in get incident edges";
  EXPECT_EQ(3, vec2.size()) << "error in get incident edges";
  EXPECT_EQ(1, vec3.size()) << "error in get incident edges";

  EXPECT_EQ(vec, vec4)     << "error in get incident edges";
  EXPECT_EQ(vec.at(0), vec4.at(0))     << "error in get incident edges";
  EXPECT_EQ(vec.at(1), vec4.at(1))     << "error in get incident edges";
  ASSERT_NE(vec.at(0), vec4.at(1))     << "error in get incident edges";
  ASSERT_NE(vec.at(1), vec4.at(0))     << "error in get incident edges";
}

// Test Incident iterator
TEST_F(GraphPointFixture, IncidentIterator)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);
  NodeType n5 = graph.add_node(points[5]);
  NodeType n6 = graph.add_node(points[6]);

  EdgeType e0 = graph.add_edge(n0, n1);
  EdgeType e1 = graph.add_edge(n1, n3);
  EdgeType e2 = graph.add_edge(n1, n4);

  graph.add_edge(n0, n2);
  graph.add_edge(n2, n5);
  graph.add_edge(n2, n6);

  IncidentIter it = n1.edge_begin();
  
  EXPECT_EQ(e0, (*it))  << " error in incident iteration ";
  ASSERT_NE(e1, (*it))  << " error in incident iteration ";
  ASSERT_NE(e2, (*it))  << " error in incident iteration ";
  ++it;

  EXPECT_EQ(e1, (*it))  << " error in incident iteration ";
  ASSERT_NE(e0, (*it))  << " error in incident iteration ";
  ASSERT_NE(e2, (*it))  << " error in incident iteration ";
  ++it;

  EXPECT_EQ(e2, (*it))  << " error in incident iteration ";
  ASSERT_NE(e0, (*it))  << " error in incident iteration ";
  ASSERT_NE(e1, (*it))  << " error in incident iteration ";

  // Iterate over its ajdacent edges
  IncidentIter r_begin = n1.edge_begin();
  IncidentIter r_end   = n1.edge_end();

  int iter = 0;

  for(;r_begin != r_end; ++r_begin) 
    iter++;

  EXPECT_EQ(iter, 3)    << " error in node iteration ";
}

// Test Incident iterator Nodes
TEST_F(GraphPointFixture, IncidentIteratorNodes)
{
  NodeType n1 = graph.add_node(points[1]);
  graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);

  //EdgeType e0 = graph.add_edge(n0, n1);
  EdgeType e1 = graph.add_edge(n1, n3);
  EdgeType e2 = graph.add_edge(n1, n4);

  IncidentIter it = n1.edge_begin();
  
  EXPECT_EQ(e1, (*it))  << " error in incident iteration ";
  EXPECT_EQ(e1.node1(), (*it).node1())  << " error in incident iteration ";
  EXPECT_EQ(e1.node2(), (*it).node2())  << " error in incident iteration ";
  ASSERT_NE(e1.node1(), (*it).node2())  << " error in incident iteration ";
  ASSERT_NE(e1.node2(), (*it).node1())  << " error in incident iteration ";
  ASSERT_NE(e2, (*it))  << " error in incident iteration ";
}

// Test Functor Node Comparison
TEST_F(GraphPointFixture, NodeComparison)
{
  NodeType n1 = graph.add_node(Point(0));
  NodeType n2 = graph.add_node(Point(1));

  NodeType n3 = graph.add_node(Point(3));
  NodeType n4 = graph.add_node(Point(8));

  // Declare the functor 
  NodeComparison node_comp;

  Point p1 = Point(0.8, 0.8, 0.8);
  Point p2 = Point(6, 5, 7);

  EXPECT_EQ(false, node_comp(n1, n2, p1)) << " error in nearest node" ;
  EXPECT_EQ(true , node_comp(n4, n3, p2)) << " error in nearest node" ;
}

// Test Nearest Node
TEST_F(GraphPointFixture, NearestNode)
{
  NodeType n1 = graph.add_node(points[0]);
  NodeType n2 = graph.add_node(points[1]);
  graph.add_node(points[2]);
  graph.add_node(points[3]);

  Point p1 = Point(0.8, 0.8, 0.8);

  NodeIter nearest = nearest_node(graph, p1);
  NodeType node_comp = *nearest;

  EXPECT_EQ(n2, node_comp) << " error in nearest node" ;
  ASSERT_NE(n1, node_comp) << " error in nearest node" ;
}

// Test Default Initialize nodes
TEST_F(GraphPointFixture, DefaultInitializeNodes) 
{
  NodeType n0 = graph.add_node(points[0], 0);
  NodeType n1 = graph.add_node(points[1], 1);
  NodeType n2 = graph.add_node(points[2], 2);

  graph.default_initialize_nodes(-1);

  EXPECT_EQ(-1, n0.value()) << "error in default initialize nodes";
  EXPECT_EQ(-1, n1.value()) << "error in default initialize nodes";
  EXPECT_EQ(-1, n2.value()) << "error in default initialize nodes";
  ASSERT_NE(0,  n2.value()) << "error in default initialize nodes";
}

// Test Set node val
TEST_F(GraphPointFixture, SetNodeVal) 
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[0], 1);

  n0.value() = 3;
  n1.value() = 2;

  EXPECT_EQ(3, n0.value()) << "error in Set node val";
  EXPECT_EQ(2, n1.value()) << "error in Set node val";
  ASSERT_NE(0, n0.value()) << "error in Set node val";
  ASSERT_NE(1, n1.value()) << "error in Set node val";
}

// Test Shortest path length
TEST_F(GraphPointFixture, ShortestPathLength)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);
  NodeType n4 = graph.add_node(points[4]);
  NodeType n5 = graph.add_node(points[5]);
  NodeType n6 = graph.add_node(points[6]);
  NodeType n7 = graph.add_node(points[7]);

  graph.add_edge(n0, n1);
  graph.add_edge(n1, n3);
  graph.add_edge(n1, n4);

  graph.add_edge(n0, n2);
  graph.add_edge(n2, n5);
  graph.add_edge(n2, n6);
  graph.add_edge(n5, n4);
  graph.add_edge(n5, n7);
  graph.add_edge(n4, n7);

  Point p1 = Point(0.1, 0.1, 0.1);
  NodeIter nearest = nearest_node(graph, p1);
  NodeType root = *nearest;

  int path = shortest_path_lengths(graph, root);

  EXPECT_EQ(n0, root) << " error in shortest path length";
  ASSERT_NE(n1, root) << " error in shortest path length";

  EXPECT_EQ(0,  n0.value()) << " error in shortest path length";
  EXPECT_EQ(1,  n1.value()) << " error in shortest path length";
  EXPECT_EQ(1,  n2.value()) << " error in shortest path length";
  EXPECT_EQ(2,  n3.value()) << " error in shortest path length";
  EXPECT_EQ(2,  n4.value()) << " error in shortest path length";
  EXPECT_EQ(2,  n5.value()) << " error in shortest path length";
  EXPECT_EQ(2,  n6.value()) << " error in shortest path length";

  ASSERT_NE(-1, n0.value()) << " error in shortest path length";
  ASSERT_NE(-1, n1.value()) << " error in shortest path length";
  ASSERT_NE(-1, n2.value()) << " error in shortest path length";
  ASSERT_NE(-1, n3.value()) << " error in shortest path length";
  ASSERT_NE(-1, n4.value()) << " error in shortest path length";
  ASSERT_NE(-1, n5.value()) << " error in shortest path length";
  ASSERT_NE(-1, n6.value()) << " error in shortest path length";

  EXPECT_EQ(path, 3)  << " error in shortest path length";
}
