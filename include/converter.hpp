#include <vector>
#include "edge.hpp"
/*
* alias Graph_ is shortcut for graph representation in algorithms and Converter class
*/
using Graph_ = std::vector<std::vector<std::pair<int, long long>>>; // first = dest, second = weight

struct Converter {
    /*
    * Converter transforms input graph representations to our uniform format
    * Possible input formats: adjacency list, adjacency matrix, list of edges
    * Output format: adjacency list (using Graph_)
    */
    Graph_ graph; // Graph_ graph is a result of convertions

    Graph_ convert(long long** matrix, int size1d);
    
    Graph_
    convert(const std::vector<Edge>& edge_list, bool is_directed);

    // Do nothing, input format = output format
    Converter(const Graph_& graph);

    // Adjancency matrix
    Converter(long long** matrix, int size1d);

    // Edge list
    Converter(const std::vector<Edge>& edge_list, bool is_directed);
};