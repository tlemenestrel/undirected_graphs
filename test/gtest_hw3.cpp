#include <gtest/gtest.h>
#include <fstream>
#include "CME212/Util.hpp"
#include "Graph.hpp"
#include <vector>
#include "GraphSymmetricMatrix.hpp"

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

TEST_F(GraphPointFixture, ChangeNodeCoord) 
{
  NodeType n1 = graph.add_node(Point(0, 0, 10));
  EXPECT_EQ(10, n1.position().z) << "Error: changing point position";
  EXPECT_EQ(0,  n1.position().x) << "Error: changing point position";

}