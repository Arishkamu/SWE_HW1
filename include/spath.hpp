#ifndef SPATH_HPP_
#define SPATH_HPP_
#include "converter.hpp"
#include "edge.hpp"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

/*
 * alias Graph_ is shortcut for graph representation in algorithms and Converter
 * class
 */
using Graph_ =
    std::vector<std::vector<std::pair<int, long long>>>; // first = dest, second
                                                         // = weight

/*
 * Edge structure describes an edge in a graph.
 * source - vertex where an edge starts from
 * dest - vertex where an edge points to
 * weight - number. Our graphs are weighted
 * INF is used in all algorithms as maximum possible value of weight in terms of
 * limits of long long
 */

std::pair<long long, std::vector<int>>
dijkstra_high_density(int start, int finish, Converter c);

std::pair<long long, std::vector<int>>
dijkstra_low_density(int start, int finish, Converter c);

std::pair<std::vector<long long>, std::vector<int>> bellman_ford(int start,
                                                                 Converter c);

std::vector<long long> johnson(Graph_ &g);

std::pair<long long, std::vector<int>> astar(int start, int finish, Converter c,
                                             long long (*heuristic)(int goal,
                                                                    int other));

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<long long>>>
floyd_warshall(Converter c);

std::vector<int> path_from_next(const std::vector<std::vector<int>> &next,
                                int start, int finish);

std::vector<int> path_from_parent(const std::vector<int> &parent, int start,
                                  int finish);

std::pair<long long, std::vector<int>>
bellman_for_two_vertices(int start, int finish, Converter c);

std::pair<long long, std::vector<int>> bfs(int start, int finish, Converter c);

std::vector<std::pair<int, int>> lee_algorithm(int start, int finish,
                                               Converter converter);

std::pair<long long, std::vector<int>>
naive_shortest_path(int start, int finish, Converter c);

std::pair<double, std::vector<int>>
karp_algorithm(Converter c);

#endif
