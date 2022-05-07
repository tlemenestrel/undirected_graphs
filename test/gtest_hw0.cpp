#include <gtest/gtest.h>
#include <fstream>
#include "CME212/Util.hpp"
#include "Graph.hpp"

class GraphPointFixture : public ::testing::Test 
{
 protected:
   //Define types
  using GraphType   = Graph<int, int>;
  using NodeType    = typename GraphType::node_type;
  using EdgeType    = typename GraphType::edge_type;
  using size_type = unsigned;

  //Set up Graph and Points
  GraphType graph;
  std::vector<Point> points;
  virtual void SetUp() {
    for(int i = 0; i < 10; i++)
      points.push_back(Point(i));
  }
};

// Add a node
TEST_F(GraphPointFixture, AddNode)
{
  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  graph.add_node(points[2]);
  graph.add_node(points[2]);

  EXPECT_EQ(graph.num_nodes(), 5) << " error when adding nodes " ;
}

// Test has edge
TEST_F(GraphPointFixture, HasEdge)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);

  graph.add_edge(n0, n1);
  
  EXPECT_EQ(true, graph.has_edge(n0, n1)) << " add edge creates duplicate edges " ;
}

// Add existing edge doesn't create another edge
TEST_F(GraphPointFixture, AddEdge)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);

  graph.add_edge(n0, n1);
  graph.add_edge(n0, n1);
  graph.add_edge(n1, n0);

  EXPECT_EQ(graph.num_edges(), 1) << " add edge creates duplicate edges " ;
}

// Test has_node function
TEST_F(GraphPointFixture, HasNode)
{
  GraphType::node_type n0 = graph.add_node(points[0]);
  EXPECT_TRUE( graph.has_node(n0) ) << "has_node did not find n0";
}

// Test num nodes/size functions
TEST_F(GraphPointFixture, Size)
{
  EXPECT_EQ(graph.num_nodes(),graph.size()) << "num_nodes and size are different"  ;
  EXPECT_EQ(graph.size(), 0) << "starting size is not 0" ;

  graph.add_node(points[0]);
  graph.add_node(points[1]);
  graph.add_node(points[2]);
  graph.add_node(points[3]);

  EXPECT_EQ(graph.num_nodes(), graph.size()) << "num_nodes and size are different";
  EXPECT_EQ(graph.num_nodes(), 4) << "size is incorrect";
}

// Make sure the number of edges are counted correctly
// Should be 6 and not 12 as duplicated edges are not allowed
// Cannot have (a, b) and (b, a)

TEST_F(GraphPointFixture, NumEdges)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  NodeType n3 = graph.add_node(points[3]);

  graph.add_edge(n0, n1);
  graph.add_edge(n0, n2);
  graph.add_edge(n0, n3);
  graph.add_edge(n1, n0);
  graph.add_edge(n1, n2);
  graph.add_edge(n1, n3);
  graph.add_edge(n2, n0);
  graph.add_edge(n2, n1);
  graph.add_edge(n2, n3);
  graph.add_edge(n3, n0);
  graph.add_edge(n3, n1);
  graph.add_edge(n3, n2);

  EXPECT_EQ(graph.num_edges(), 6) << " num_edges() doesn't find the right num of edges " ;
}

// Verify only one of e0 < e1 or e1 < e0 is true
TEST_F(GraphPointFixture, Tricotomy)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  NodeType n2 = graph.add_node(points[2]);
  
  EdgeType e0 = graph.add_edge(n0, n1);
  EdgeType e1 = graph.add_edge(n1, n2);

  EXPECT_TRUE( (e0 < e1) ^ (e1 < e0) ) << "error in edge comparison";
}

// Test edge function
TEST_F(GraphPointFixture, Edge)
{
  NodeType n0 = graph.add_node(points[0]);
  NodeType n1 = graph.add_node(points[1]);
  EdgeType e0 = graph.add_edge(n0, n1);
  
  EXPECT_EQ(e0, graph.edge(0)) << "error in edge retrieval"  ;
}

// Test adding node with default value
TEST_F(GraphPointFixture, DefaultNodeVal)
{
  graph.add_node(points[0]);
  EXPECT_EQ( graph.node(0).value(), 0 ) << "add_node does not initalize node vale with a default 0 value";
}

// Test fetch node
TEST_F(GraphPointFixture, FetchNode)
{
  NodeType n0 = graph.add_node(points[0], 3);
  NodeType n1 = graph.add_node(points[1], 4);
  NodeType n2 = graph.add_node(points[2], 5);

  EXPECT_EQ(3, n0.fetch_node().get_node_value()) << "error in fetch node";
  EXPECT_EQ(4, n1.fetch_node().get_node_value()) << "error in fetch node";
  EXPECT_EQ(5, n2.fetch_node().get_node_value()) << "error in fetch node";
  ASSERT_NE(4, n2.fetch_node().get_node_value()) << "error in fetch node";
}

// Test fetch edge
TEST_F(GraphPointFixture, FetchEdge)
{
  NodeType n0 = graph.add_node(points[0], 3);
  NodeType n1 = graph.add_node(points[1], 4);
  NodeType n2 = graph.add_node(points[2], 5);

  EdgeType e0 = graph.add_edge(n0, n1);
  EdgeType e1 = graph.add_edge(n1, n2);

  EXPECT_EQ(e0.fetch_edge().get_edge_n1(), n0.index()) << "error in fetch edge";
  EXPECT_EQ(e1.fetch_edge().get_edge_n2(), n2.index()) << "error in fetch edge";
  ASSERT_NE(e1.fetch_edge().get_edge_n2(), n1.index()) << "error in fetch edge";
}

// Test find edge
TEST_F(GraphPointFixture, FindEdge)
{
  NodeType n0 = graph.add_node(points[0], 3);
  NodeType n1 = graph.add_node(points[1], 4);
  NodeType n2 = graph.add_node(points[2], 5);

  EdgeType e0 = graph.add_edge(n0, n1);
  EdgeType e1 = graph.add_edge(n0, n2);

  EXPECT_EQ(e0, graph.find_edge(n0, n1)) << "error in find edge";
  EXPECT_EQ(e0, graph.find_edge(n1, n0)) << "error in find edge";
  ASSERT_NE(e1, graph.find_edge(n0, n1)) << "error in fetch edge";
  ASSERT_NE(e0, graph.find_edge(n0, n2)) << "error in fetch edge";
}

// Change Node Index
TEST_F(GraphPointFixture, ChangeNodeIndex)
{
  NodeType n0 = graph.add_node(points[0], 3);
  NodeType n1 = graph.add_node(points[1], 4);
  NodeType n2 = graph.add_node(points[2], 5);

  EXPECT_EQ(0, n0.index()) << "error in node swap";
  EXPECT_EQ(1, n1.index()) << "error in node swap";
  EXPECT_EQ(2, n2.index()) << "error in node swap";

  n2.get_node_idx() = n0.index();
  EXPECT_EQ(0, n2.index()) << "error in node swap";
  EXPECT_EQ(0, n0.index()) << "error in node swap";
}
