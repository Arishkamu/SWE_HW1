#include "spath.hpp"

std::vector<int> path_from_parent(const std::vector<int> &parent, int start,
                                  int finish) {
  std::vector<int> path{finish};
  int cur = finish;
  while (parent[cur] != start) {
    path.push_back(parent[cur]);
    cur = parent[cur];
  }
  path.push_back(start);
  std::reverse(path.begin(), path.end());
  return path;
}

std::pair<long long, std::vector<int>>
dijkstra_high_density(int start, int finish, Converter converter) {
  /*
   * Optimal for graphs with high density
   * Complexity: O(n^2 + m)
   */
  const Graph_ &graph = converter.graph;
  int n = graph.size();
  // start - start vertex, finish - finish vertex
  assert(start >= 0 && start < n);
  assert(finish >= 0 && finish < n);

  // safe, we do only comparisons, no arithmetic operations with INF
  long long INF = std::numeric_limits<long long>::max();
  std::vector<long long> dist(n, INF);
  std::vector<int> parent(n, -1);
  std::vector<bool> used(n, false);

  dist[start] = 0;
  for (int i = 0; i < n; i++) {
    int vertex = -1;
    for (int j = 0; j < n; j++) {
      if (!used[j] && (vertex == -1 || dist[j] < dist[vertex]))
        vertex = j;
    }
    if (dist[vertex] == INF)
      break;
    used[vertex] = true;

    for (int j = 0; j < graph[vertex].size(); j++) {
      int u = graph[vertex][j].first; // u is adjacent vertex to vertex
      long long weight = graph[vertex][j].second;

      if (dist[u] > dist[vertex] + weight) {
        dist[u] = dist[vertex] + weight;
        parent[u] = vertex;
      }
    }
  }
  if (dist[finish] == INF) {
    return {-1, std::vector<int>()};
  }
  auto path = path_from_parent(parent, start, finish);
  return {dist[finish], path};
}

std::pair<long long, std::vector<int>>
dijkstra_low_density(int start, int finish, Converter converter) {
  /*
   * Optimal for graphs with low density
   * Complexity: O(m * log(n))
   */
  const Graph_ &graph = converter.graph;
  int n = graph.size();
  // start - start vertex, finish - finish vertex
  assert(start >= 0 && start < n);
  assert(finish >= 0 && finish < n);

  // safe, we do only comparisons, no arithmetic operations with INF
  long long INF = std::numeric_limits<long long>::max();
  std::vector<long long> dist(n, INF);
  std::vector<int> parent(n, -1);
  dist[start] = 0;

  std::set<std::pair<long long, int>> estimates;
  estimates.insert({0, start});

  while (!estimates.empty()) {
    auto [est, v] = *estimates.begin();
    estimates.erase(estimates.begin());
    for (auto &item : graph[v]) {
      int to = item.first;
      long long weight = item.second;
      if (dist[to] > dist[v] + weight) {
        estimates.erase({dist[to], to});
        dist[to] = dist[v] + weight;
        parent[to] = v;
        estimates.insert({dist[to], to});
      }
    }
  }

  if (dist[finish] == INF) {
    return {-1, std::vector<int>()};
  }
  auto path = path_from_parent(parent, start, finish);
  return {dist[finish], path};
}

std::pair<std::vector<long long>, std::vector<int>>
bellman_ford(int start, Converter converter) {
  /*
   * Weight of every cycle should be non-negative
   * Complexity: O(mn)
   */
  const Graph_ &graph = converter.graph;
  int n = graph.size();
  // start - start vertex
  assert(start >= 0 && start < n);

  // safe, we do only comparisons, no arithmetic operations with INF
  long long INF = std::numeric_limits<long long>::max();
  // O(n) memory
  std::vector<int> parent(n, -1);
  std::vector<long long> dist(n, INF);
  dist[start] = 0;

  for (int k = 1; k < n; k++) {
    bool changed = false; // Flag for early break

    for (int v = 0; v < n; v++) {
      for (auto &item : graph[v]) {
        int u = item.first;
        long long weight = item.second;
        if (dist[u] != INF && dist[u] + weight < dist[v]) {
          dist[v] = dist[u] + weight;
          parent[v] = u;
          changed = true;
        }
      }
    }

    if (!changed)
      break;
  }
  return {dist, parent};
}

std::pair<long long, std::vector<int>>
bellman_for_two_vertices(int start, int finish, Converter converter) {
  auto [dist, parent] = bellman_ford(start, converter);

  // safe, we do only comparisons, no arithmetic operations with INF
  long long INF = std::numeric_limits<long long>::max();
  if (dist[finish] == INF) {
    return {-1, std::vector<int>()};
  }
  std::vector<int> path = path_from_parent(parent, start, finish);
  return std::pair<long long, std::vector<int>>{dist[finish], path};
}

std::vector<long long> johnson(Graph_ &graph) {
  /*
   * Complexity: O(nm)
   * returns potentials!!!
   * After that you can apply more fast Dijkstra to find needed shortest paths.
   * It is well when you need to find multiple times the shortest path on a
   * given graph
   */

  int n = graph.size();

  /*
   * fakeVertex is used for creating one temp vertex which is needed
   * for calculating potentials
   * After processing temp vertex will be removed
   */
  std::vector<std::pair<int, long long>> fakeVertex(n);
  for (int i = 0; i < n; i++)
    fakeVertex[i] = {i, 0};

  graph.push_back(std::move(fakeVertex)); // adding temp vertex

  // Find potentials with Ford-Bellman
  std::vector<long long> potentials = bellman_ford(n, graph).first;
  potentials.pop_back(); // deleting potential for temp vertex

  graph.pop_back(); // deleting temp vertex
  return potentials;
}

std::pair<long long, std::vector<int>>
astar(int start, int finish, Converter converter,
      long long (*heuristic)(int goal, int other)) {
  /*
   * Complexity: O(m*log(n)*complexity(heuristic))
   * modified Dijkstra's algorithm for searching path to a one specific point
   * by using some "good" estimator of distance to it for giving priority
   * instead of pure distance
   */
  const Graph_ &graph = converter.graph;
  int n = graph.size();
  // start - start vertex, finish - finish vertex
  assert(start >= 0 && start < n);
  assert(finish >= 0 && finish < n);

  // safe, we do only comparisons, no arithmetic operations with INF
  long long INF = std::numeric_limits<long long>::max();
  std::vector<long long> dist(n, INF);
  std::vector<int> parent(n, -1);
  dist[start] = 0;
  std::set<std::pair<long long, long long>> estimates;
  estimates.insert({heuristic(finish, start), start});
  while (!estimates.empty()) {
    auto [est, v] = *estimates.begin();
    if (v == finish)
      break;
    estimates.erase(estimates.begin());
    for (auto &[neigh, weight] : graph[v]) {
      if (dist[neigh] > dist[v] + weight) {
        dist[neigh] = dist[v] + weight;
        estimates.insert({dist[neigh] + heuristic(finish, neigh), neigh});
        parent[neigh] = v;
      }
    }
  }
  if (dist[finish] == INF) {
    return {-1, std::vector<int>()};
  }
  auto path = path_from_parent(parent, start, finish);
  return {dist[finish], path};
}

std::pair<std::vector<std::vector<int>>, std::vector<std::vector<long long>>>
floyd_warshall(Converter converter) {
  /*
   * Complexity: O(n^3)
   * computes distance from all to all and returns matrix next[u][v] which holds
   * value for the next vertice on the path from u to v, which can be then used
   * in function path_from_next to restore any path.
   */
  const Graph_ &graph = converter.graph;

  // safe, we do only comparisons, no arithmetic operations with INF
  long long INF = std::numeric_limits<long long>::max();
  std::vector<std::vector<long long>> dist(
      graph.size(), std::vector<long long>(graph.size(), INF));
  std::vector<std::vector<int>> next(graph.size(),
                                     std::vector<int>(graph.size(), -1));
  for (int i = 0; i < graph.size(); i++) {
    dist[i][i] = 0;
    next[i][i] = i;
  }
  for (int i = 0; i < graph.size(); i++) {
    for (const auto &[v, weight] : graph[i]) {
      dist[i][v] = weight;
      next[i][v] = v;
    }
  }
  for (int i = 0; i < graph.size(); i++) {
    for (int u = 0; u < graph.size(); u++) {
      for (int v = 0; v < graph.size(); v++) {
        if (dist[u][i] != INF && dist[i][v] != INF) {
          if (dist[u][i] + dist[i][v] < dist[u][v]) {
            dist[u][v] = dist[u][i] + dist[i][v];
            next[u][v] = next[u][i];
          }
        }
      }
    }
  }
  return {next, dist};
}

std::pair<long long, std::vector<int>> bfs(int start, int finish,
                                           Converter converter) {
  /*
   * BFS algorithm is used to find the shortest path between start and finish in
   * an unweighted graph Time complexity: O(V + E) Space complexity: O(V) V -
   * number of vertices E - number of edges
   */
  const Graph_ &graph = converter.graph;
  std::size_t n = graph.size();

  std::queue<int> queue;
  queue.push(start);

  long long INF = std::numeric_limits<long long>::max();

  // Distance from the start to all vertices
  std::vector<long long> dist(n, INF);
  dist[start] = 0;

  // Vertex parents (for path trace)
  std::vector<int> parent(n, -1);

  // Loop for BFS
  while (!queue.empty()) {
    int from = queue.front();
    queue.pop();
    for (auto u : graph[from]) {
      int to = u.first;
      // A shorter route has been found
      if (dist[to] > dist[from] + 1) {
        dist[to] = dist[from] + 1;
        parent[to] = from;
        queue.push(to);
      }
    }
  }

  // The path to the final is not found
  if (dist[finish] == INF) {
    return {-1, std::vector<int>()};
  }

  // Trace the path
  auto path = path_from_parent(parent, start, finish);
  return {dist[finish], path};
}

std::pair<long long, std::vector<int>>
naive_shortest_path(int start, int finish, Converter converter) {
  /*
   * Naive algorithm is used to find the shortest path between start and finish
   * in a weighted graph. It works by generating all possible paths from start
   * to finish and selecting the one with the minimum total weight.
   *
   * Time complexity: O(V! * E)
   * Space complexity: O(V)
   * V - number of vertices
   * E - number of edges
   */
  const Graph_ &graph = converter.graph;
  std::size_t n = graph.size();

  long long INF = std::numeric_limits<long long>::max();

  // List of vertices to be permuted (excluding start and finish)
  std::vector<int> vertices;
  for (int i = 0; i < n; ++i) {
    if (i != start && i != finish) {
      vertices.push_back(i);
    }
  }

  // Initialize minimum distance and best path
  long long min_distance = INF;
  std::vector<int> best_path;

  // Check direct edge from start to finish
  for (const auto &edge : graph[start]) {
    if (edge.first == finish) {
      min_distance = edge.second;
      best_path = {start, finish};
      break;
    }
  }

  // Generate all possible subsets of intermediate vertices
  for (int mask = 0; mask < (1 << vertices.size()); ++mask) {
    std::vector<int> subset;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
      if (mask & (1 << i)) {
        subset.push_back(vertices[i]);
      }
    }

    // Generate all permutations of the subset
    do {
      long long current_distance = 0;
      int current_vertex = start;
      bool path_exists = true;

      // Check if the current permutation forms a valid path
      for (int next_vertex : subset) {
        bool edge_found = false;
        for (const auto &edge : graph[current_vertex]) {
          if (edge.first == next_vertex) {
            current_distance += edge.second;
            current_vertex = next_vertex;
            edge_found = true;
            break;
          }
        }
        if (!edge_found) {
          path_exists = false;
          break;
        }
      }

      // If the path is valid, check the edge to the finish vertex
      if (path_exists) {
        bool edge_found = false;
        for (const auto &edge : graph[current_vertex]) {
          if (edge.first == finish) {
            current_distance += edge.second;
            edge_found = true;
            break;
          }
        }
        // If a shorter path is found, update the minimum distance and best path
        if (edge_found && current_distance < min_distance) {
          min_distance = current_distance;
          best_path = subset;
          best_path.insert(best_path.begin(), start);
          best_path.push_back(finish);
        }
      }
    } while (std::next_permutation(subset.begin(), subset.end()));
  }

  // If no path is found, return {-1, empty vector}
  if (min_distance == INF) {
    return {-1, std::vector<int>()};
  } else {
    // Return the minimum distance and the best path
    return {min_distance, best_path};
  }
}
std::vector<int> path_from_next(const std::vector<std::vector<int>> &next,
                                int start, int finish) {
  /*
   * It is applied to 'next' from floyd_warshall algorithm to reconstruct the
   * path from start vertex (start) to finish vertex (finish)
   */
  int cur = start;
  std::vector<int> path;
  while (cur != finish) {
    path.push_back(cur);
    cur = next[cur][finish];
  }

  path.push_back(finish);
  return path;
}
