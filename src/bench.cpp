#include <benchmark/benchmark.h>
#include "graph.h"
#include "../include/spath.hpp"
#include <vector>
#include <random>
#include <memory>
#include <limits>
#include <atomic>

// Global counters for mem statistic
std::atomic<int64_t> g_allocated_memory(0); // number of allocated bytes
std::atomic<int64_t> g_num_allocs(0);  // number of allocations 

/*
* Redefinition of operator new since we should intercept all mem allocations
*/
void* operator new(size_t size) {
    g_allocated_memory += size;
    g_num_allocs++;
    return malloc(size);
}


/*
* This class is responsible for calculating consumed memory during bench-tests
*/
class CustomMemoryManager: public benchmark::MemoryManager {
public:
    // executed before each test
    void Start() BENCHMARK_OVERRIDE {
        g_num_allocs = 0;
        g_allocated_memory = 0;
    }

    // executed after each test
    void Stop(Result& result) BENCHMARK_OVERRIDE {
        result.num_allocs = g_num_allocs;
        result.max_bytes_used = g_allocated_memory;
    }
};
    
std::unique_ptr<CustomMemoryManager> mm(new CustomMemoryManager());

Graph_ 
generate_random_graph(int num_nodes, int num_edges,
                 int min_weight, int max_weight, bool is_undirected = true) {
    Graph_ graph(num_nodes);

    // graph generation 
    auto generated_proxy_graph = Graph::random(num_nodes, num_edges).allowLoops().connected();
    if (!is_undirected) {
        generated_proxy_graph = generated_proxy_graph.directed();
    }

    // Graph_ as edge list
    std::vector<std::pair<int, int>> generated_graph = generated_proxy_graph.g().shuffled().edges();

    // weight generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> weight_dist(min_weight, max_weight);

    for (int i = 0; i < generated_graph.size(); ++i) {
        int weight = weight_dist(gen);
        auto pp = generated_graph[i];
        graph[pp.first].push_back({pp.second, weight});
        if (is_undirected) 
            graph[pp.second].push_back({pp.first, weight});
    }

    return graph;
}

static void 
BM_Dijkstra_high_density(benchmark::State& state) {
    int num_nodes = state.range(0);
    int num_edges = state.range(1);
    int min_weight = state.range(2);
    int max_weight = state.range(3);

    Graph_ graph = generate_random_graph(num_nodes, num_edges, min_weight, max_weight);
    Converter converter(graph);
    int start_node = 0;
    int end_node = num_nodes - 1;

    for (auto _ : state) {
        benchmark::DoNotOptimize(dijkstra_high_density(start_node, end_node, converter));
    }

    // For report
    state.counters["Nodes"] = num_nodes;
    state.counters["Edges"] = num_edges;
}

static void 
BM_Dijkstra_low_density(benchmark::State& state) {
    int num_nodes = state.range(0);
    int num_edges = state.range(1);
    int min_weight = state.range(2);
    int max_weight = state.range(3);

    Graph_ graph = generate_random_graph(num_nodes, num_edges, min_weight, max_weight);
    Converter converter(graph);
    int start_node = 0;
    int end_node = num_nodes - 1;

    for (auto _ : state) {
        benchmark::DoNotOptimize(dijkstra_low_density(start_node, end_node, converter));
    }

    // For report
    state.counters["Nodes"] = num_nodes;
    state.counters["Edges"] = num_edges;
}

static void 
BM_bellman_ford(benchmark::State& state) {
    int num_nodes = state.range(0);
    int num_edges = state.range(1);
    int min_weight = state.range(2);
    int max_weight = state.range(3);

    Graph_ graph = generate_random_graph(num_nodes, num_edges, min_weight, max_weight);
    Converter converter(graph);
    int start_node = 0;

    for (auto _ : state) {
        benchmark::DoNotOptimize(bellman_ford(start_node, converter));
    }

    // For report
    state.counters["Nodes"] = num_nodes;
    state.counters["Edges"] = num_edges;
}

static void 
BM_floyd_warshall(benchmark::State& state) {
    int num_nodes = state.range(0);
    int num_edges = state.range(1);
    int min_weight = state.range(2);
    int max_weight = state.range(3);

    Graph_ graph = generate_random_graph(num_nodes, num_edges, min_weight, max_weight);
    Converter converter(graph);

    for (auto _ : state) {
        benchmark::DoNotOptimize(floyd_warshall(converter));
    }

    // For report
    state.counters["Nodes"] = num_nodes;
    state.counters["Edges"] = num_edges;
}


BENCHMARK(BM_Dijkstra_high_density)
    ->Args({100, 500, 1, 100})    // 100 vertexes, 500 edges, weight up to 100
    ->Args({1000, 10000, 1, 100}) // 1000 vertexes, 10K edges, weight up to 100
    ->Args({5000, 50000, 1, 100}) // 5K vertexes, 50K edges, weight up to 100
    ->Args({10000, 400000, 1, 100}) // 10K vertexes, 400K edges, weight up to 100 
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_Dijkstra_low_density)
    ->Args({100, 500, 1, 100})    // 100 vertexes, 500 edges, weight up to 100
    ->Args({1000, 10000, 1, 100}) // 1000 vertexes, 10K edges, weight up to 100
    ->Args({5000, 50000, 1, 100}) // 5K vertexes, 50K edges, weight up to 100
    ->Args({10000, 400000, 1, 100}) // 10K vertexes, 400K edges, weight up to 100
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_bellman_ford)
    ->Args({100, 500, -150, 100})    // 100 vertexes, 500 edges, weight up to 100
    ->Args({1000, 10000, -150, 100}) // 1000 vertexes, 10K edges, weight up to 100
    ->Args({5000, 50000, -150, 100}) // 5K vertexes, 50K edges, weight up to 100
    ->Args({5000, 400000, -150, 100}) // 5K vertexes, 400K edges, weight up to 100
    ->Unit(benchmark::kMillisecond);

BENCHMARK(BM_floyd_warshall)
    ->Args({100, 500, -150, 100})    // 100 vertexes, 500 edges, weight up to 100
    ->Args({1000, 1000, -150, 100}) // 1000 vertexes, 1K edges, weight up to 100
    ->Unit(benchmark::kMillisecond);

int main(int argc, char** argv)
{
    config.generateLargeObjects = true;
    ::benchmark::RegisterMemoryManager(mm.get());
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}
