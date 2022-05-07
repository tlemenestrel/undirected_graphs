/**
 * @file shortest_path.cpp
 * Implimentation file for using our templated Graph to determine shortest paths.
 */

#include <fstream>
#include <math.h>
#include <queue>
#include <vector>

#include "CME212/SFML_Viewer.hpp"
#include "CME212/Util.hpp"
#include "CME212/Color.hpp"

#include "Graph.hpp"

// Define our types
using GraphType = Graph<int, int>;
using NodeType  = typename GraphType::node_type;
using NodeIter  = typename GraphType::node_iterator;
using IncidentIter  = typename GraphType::incident_iterator;
using size_type = unsigned;

// Functor for String comparison
struct NodeComparison
{
    bool operator() (const NodeType& n1, const NodeType& n2, const Point& pt_comp) 
    {
    	// Get the points coordinates of each point
        Point p1 = n1.fetch_node().pt_; 
        Point p2 = n2.fetch_node().pt_;

        // Compute distances
        float d1 = sqrt(
        	pow((p1.x - pt_comp.x), 2) + 
        	pow((p1.y - pt_comp.y), 2) +
        	pow((p1.z - pt_comp.z), 2));

        float d2 = sqrt(
        	pow((p2.x - pt_comp.x), 2) + 
        	pow((p2.y - pt_comp.y), 2) +
        	pow((p2.z - pt_comp.z), 2));

        if(d1 < d2)
        	return true;
        return false;
    }
};

struct ColorPathFunctor
{
  int maxpath_;

  ColorPathFunctor(const int p) : maxpath_(p) {};

  CME212::Color operator()(const Graph<int, int>::node_type n)
  {
  	float c = 1. - (1.+(float)(n.value()))/(1.+(float)(maxpath_));
    return CME212::Color::make_heat(1.0 - c);
  }
};

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Euclidean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */
NodeIter nearest_node(const GraphType& g, const Point& point)
{
	// Start the iterator
    NodeIter it = g.node_begin();

    // Declare the functor 
    NodeComparison node_comp;

    // Initiliaze the closest point as first point iterated over in the graph
    NodeType closest_point = *it;

    // Move to the next point
    ++it;

    // While not done with the search
    while(it != g.node_end())
    {
        if(node_comp((*it), closest_point, point))
        {
            closest_point = *it;
        }
        ++it;
    }
    const GraphType* graph = &g;

  	return NodeIter(graph, closest_point.index());
}

/** Update a graph with the shortest path lengths from a root node.
 * @param[in,out] g     Input graph
 * @param[in,out] root  Root node to start the search.
 * @return The maximum path length found.
 *
 * @post root.value() == 0
 * @post Graph has modified node values indicating the minimum path length
 *           to the root.
 * @post Graph nodes that are unreachable from the root have value() == -1.
 *
 * This sets all nodes' value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(const GraphType& g, NodeType& root)
{	
	// Default initialize all nodes' values to -1
	g.default_initialize_nodes(-1);

	// Mark root node as visited.
	root.value() = 0;

	// Declare the queue of nodes waiting to be evaluated and push the root.
	std::queue<NodeType> waiting;
	waiting.push(root);

 	int max = 0;

	// While the BFS is still running
	while (!waiting.empty())
	{
		NodeType r = waiting.front();
		int cur = r.value();

		// Iterate over its ajdacent edges
		IncidentIter r_begin = r.edge_begin();
		IncidentIter r_end   = r.edge_end();

		for(;r_begin != r_end; ++r_begin) 
		{
			if((*r_begin).node2().value() == -1) 
			{
				(*r_begin).node2().value() = cur + 1;
	        	if((*r_begin).node2().value() > max) 
	        		max = (*r_begin).node2().value();
				waiting.push((*r_begin).node2());
			}
			if((*r_begin).node2().value() > cur + 1) 
			{ 
				(*r_begin).node2().value() = cur + 1;
				if((*r_begin).node2().value() > max) 
					max = (*r_begin).node2().value();
			}
		}
		waiting.pop();
	}
	return max;
}
