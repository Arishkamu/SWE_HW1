#include "../include/spath.hpp"
#include <iostream>
#include <catch2/catch_test_macros.hpp>

void test_shortest_path(Converter adjacency_list_converter, Converter adjacency_matrix_converter, Converter edge_list_converter, 
int start, int finish, std::pair<long long, std::vector<int>> expected, 
std::pair<long long, std::vector<int>> (*algorithm)(int start, int finish , Converter converter)) {
    auto res1 = algorithm(start,finish,adjacency_list_converter);
    auto res2 = algorithm(start,finish,adjacency_matrix_converter);
    auto res3 = algorithm(start,finish,edge_list_converter);
    REQUIRE( res1 == expected );
    REQUIRE( res2 == expected );
    REQUIRE( res3 == expected );
}

void test_floyd_warshall(Converter adjacency_list_converter, Converter adjacency_matrix_converter, Converter edge_list_converter, 
    std::vector<std::vector<std::pair<long long, std::vector<int>>>>&expected_mtx) {
        auto [next1, dist1] = floyd_warshall(adjacency_list_converter);
        auto [next2, dist2] = floyd_warshall(adjacency_matrix_converter);
        auto [next3, dist3] = floyd_warshall(edge_list_converter);
        REQUIRE(dist1 == dist2);
        REQUIRE(dist2 == dist3);
        REQUIRE(next1 == next2);
        REQUIRE(next2 == next3);
        int n = adjacency_list_converter.graph.size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                long long INF = std::numeric_limits<long long>::max();
                std::vector<int>path1;
                long long d1;
                if (dist1[i][j] == INF) {
                    path1 = {};
                    d1 = -1;
                } else {
                    path1 = path_from_next(next1, i, j);
                    d1 = dist1[i][j];
                }
                auto ans1 = std::make_pair(d1, path1);
                REQUIRE( expected_mtx[i][j] == ans1 );
            }
        }
    }


TEST_CASE("Converter", "[Converter]") {
    SECTION("adjacency list") {
        int n = 5;
        Graph_ adjacency_list(n);
        adjacency_list[0] = {{1, 10}, {2, 20}, {4, 50}};
        adjacency_list[1] = {{0, 10}, {3, 30}};
        adjacency_list[2] = {{0, 20}, {3, 40}, {4, 60}};
        adjacency_list[3] = {{1, 30}, {2, 40}, {4, 70}};
        adjacency_list[4] = {{0, 50}, {2, 60}, {3, 70}};

        Converter converter(adjacency_list);
        REQUIRE(converter.graph == adjacency_list);
    }

    SECTION("adjacency matrix") {
        int n = 5;
        long long** adjacency_matrix = new long long*[n];
        for (int i = 0; i < n; i++) {
            adjacency_matrix[i] = new long long[n]{};
        }

        adjacency_matrix[0][1] = 10;
        adjacency_matrix[0][2] = 20;
        adjacency_matrix[0][4] = 50;
        adjacency_matrix[1][3] = 30;
        adjacency_matrix[2][3] = 40;
        adjacency_matrix[2][4] = 60;
        adjacency_matrix[3][4] = 70;
        adjacency_matrix[1][0] = 10;
        adjacency_matrix[2][0] = 20;
        adjacency_matrix[3][1] = 30;
        adjacency_matrix[3][2] = 40;
        adjacency_matrix[4][0] = 50;
        adjacency_matrix[4][2] = 60;
        adjacency_matrix[4][3] = 70;

        Converter converter(adjacency_matrix, n);
        Graph_ expected = {
            {{1, 10}, {2, 20}, {4, 50}},
            {{0, 10}, {3, 30}},
            {{0, 20}, {3, 40}, {4, 60}},
            {{1, 30}, {2, 40}, {4, 70}},
            {{0, 50}, {2, 60}, {3, 70}}
        };

        REQUIRE(converter.graph == expected);

        for (int i = 0; i < n; i++) {
            delete[] adjacency_matrix[i];
        }
        delete[] adjacency_matrix;
    }

    SECTION("edge list (undirected)") {
        std::vector<Edge> edge_list = {
            {0, 1, 10}, {0, 2, 20}, {0, 4, 50},
            {1, 3, 30}, {2, 3, 40}, {2, 4, 60}, {3, 4, 70}
        };

        Converter converter(edge_list, false);
        Graph_ expected = {
            {{1, 10}, {2, 20}, {4, 50}},
            {{0, 10}, {3, 30}},
            {{0, 20}, {3, 40}, {4, 60}},
            {{1, 30}, {2, 40}, {4, 70}},
            {{0, 50}, {2, 60}, {3, 70}}
        };

        REQUIRE(converter.graph == expected);
    }

    SECTION("edge list (directed)") {
        std::vector<Edge> edge_list = {
            {0, 1, 10}, {0, 2, 20}, {0, 4, 50},
            {1, 3, 30}, {2, 3, 40}, {2, 4, 60}, {3, 4, 70}
        };

        Converter converter(edge_list, true);
        Graph_ expected = {
            {{1, 10}, {2, 20}, {4, 50}},
            {{3, 30}},
            {{3, 40}, {4, 60}},
            {{4, 70}},
            {}
        };

        REQUIRE(converter.graph == expected);
    }

    SECTION("empty edge list") {
        std::vector<Edge> edge_list = {};
        Converter converter(edge_list, false);
        REQUIRE(converter.graph.empty());
    }

    SECTION("empty adjacency list") {
        Graph_ empty_list = {};
        Converter converter(empty_list);
        REQUIRE(converter.graph.empty());
    }

    SECTION("empty adjacency matrix") {
        int n = 5;
        long long** adjacency_matrix = new long long*[n];
        for (int i = 0; i < n; i++) {
            adjacency_matrix[i] = new long long[n]{};
        }

        Converter converter(adjacency_matrix, n);
        REQUIRE(converter.graph == Graph_(n));

        for (int i = 0; i < n; i++) {
            delete[] adjacency_matrix[i];
        }
        delete[] adjacency_matrix;
    }
}



TEST_CASE( "Correctness", "[Requirement 1]" ) {
    SECTION("empty graph") {
        int n = 5;
        std::vector<std::vector<std::pair<int, long long>>>adjacency_list(n);
        long long** adjacency_matrix = new long long*[n];
        for (int i = 0; i < n; i++) {
            adjacency_matrix[i] = new long long[n];
            for (int j = 0; j < n; j++) {
                adjacency_matrix[i][j] = 0;
            }
        }
        std::vector<Edge>edge_list;

        Converter adjacency_list_converter = Converter(adjacency_list);
        Converter adjacency_matrix_converter = Converter(adjacency_matrix, 5);
        Converter edge_list_converter = Converter(edge_list, false);
        assert(edge_list_converter.graph.empty());
        std::pair<long long, std::vector<int>> expected = {-1,std::vector<int>()};

        // Edge list is empty, so inner assertion fails and crushes the program. Can't be "checked" normally.
        SECTION("Dijkstra for high density") {
            auto res1 = dijkstra_high_density(0,2,adjacency_list_converter);
            auto res2 = dijkstra_high_density(0,2,adjacency_matrix_converter);

            REQUIRE( res1 == expected );
            REQUIRE( res2 == expected );
        }

        SECTION("Dijkstra for low density") {
            auto res1 = dijkstra_low_density(0,2,adjacency_list_converter);
            auto res2 = dijkstra_low_density(0,2,adjacency_matrix_converter);

            REQUIRE( res1 == expected );
            REQUIRE( res2 == expected );
        }

        SECTION("Bellman-Ford") {
            auto res1 = bellman_for_two_vertices(0,2,adjacency_list_converter);
            auto res2 = bellman_for_two_vertices(0,2,adjacency_matrix_converter);

            REQUIRE( res1 == expected );
            REQUIRE( res2 == expected );
        }
        SECTION("Floyd-Warshall") {
            auto [next1, dist1] = floyd_warshall(adjacency_list_converter);
            auto [next2, dist2] = floyd_warshall(adjacency_matrix_converter);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    long long INF = std::numeric_limits<long long>::max();
                    std::vector<int>path1,path2;
                    long long d1,d2;
                    if (dist1[i][j] == INF ) {
                        path1 = path2 = {};
                        d1 = d2 = -1;
                    } else {
                        path1 = path_from_next(next1, i, j);
                        path2 = path_from_next(next2, i, j);
                        d1 = dist1[i][j];
                        d2 = dist2[i][j];
                    }
                    auto ans1 = std::make_pair(d1, path1);
                    auto ans2 = std::make_pair(d2, path2);
                    if (i == j) {

                        REQUIRE( ans1 == std::pair<long long, std::vector<int>>{0,{i}} );
                        REQUIRE( ans2 == std::pair<long long, std::vector<int>>{0,{i}} );
                    } else {
                        REQUIRE( ans1 == expected );
                        REQUIRE( ans2 == expected );
                    }
                }
            }
        }


        for (int i = 0; i < n; i++) {
            delete adjacency_matrix[i];
        }
        delete adjacency_matrix;
    }

    SECTION("full graph") {
        int n = 5;
        std::vector<std::vector<std::pair<int, long long>>>adjacency_list(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                adjacency_list[i].emplace_back(j, 1);
            }
        }

        long long** adjacency_matrix = new long long*[n];
        for (int i = 0; i < n; i++) {
            adjacency_matrix[i] = new long long[n];
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                adjacency_matrix[i][j] = 1;
            }
        }
        std::vector<Edge>edge_list;

        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                if (i == j) {
                    continue;
                }
                Edge e = {i, j, 1};
                edge_list.push_back(e);
            }
        }

        Converter adjacency_list_converter = Converter(adjacency_list);
        Converter adjacency_matrix_converter = Converter(adjacency_matrix, 5);
        Converter edge_list_converter = Converter(edge_list, false);
        std::pair<long long, std::vector<int>> expected = {1,{0,3}};

        SECTION("Dijkstra for high density") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 0, 3, expected, dijkstra_high_density);
        }
        SECTION("Dijkstra for low density") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 0, 3, expected, dijkstra_low_density);
        }
        SECTION("Bellman-Ford") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 0, 3, expected, bellman_for_two_vertices);
        }
        SECTION("Floyd-Warshall") {
            std::vector<std::vector<std::pair<long long, std::vector<int>>>>expected_mtx(n);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i == j) {
                        expected_mtx[i].emplace_back(0, std::vector<int>{i});
                    } else {
                        expected_mtx[i].emplace_back(1, std::vector<int>{i,j});
                    }
                }
            }
            test_floyd_warshall(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, expected_mtx);
        }

        for (int i = 0; i < n; i++) {
            delete adjacency_matrix[i];
        }
        delete adjacency_matrix;

    }
    SECTION("complex graph") {
        int n = 5;
        std::vector<std::vector<std::pair<int, long long>>>adjacency_list(n);
        adjacency_list[0] = {{1,100}, {4, 10}, {3, 20}};
        adjacency_list[1] = {{2,7}, {0,100}};
        adjacency_list[2] = {{1,7}, {3,20}};
        adjacency_list[3] = {{4,15},{0,20},{2,20}};
        adjacency_list[4] = {{0,10},{3,15}};

        long long** adjacency_matrix = new long long*[n];
        for (int i = 0; i < n; i++) {
            adjacency_matrix[i] = new long long[n];
            for (int j = 0; j < n; j++) {
                adjacency_matrix[i][j] = 0;
            }
        }
        adjacency_matrix[0][1] = adjacency_matrix[1][0] = 100;
        adjacency_matrix[0][4] = adjacency_matrix[4][0] = 10;
        adjacency_matrix[0][3] = adjacency_matrix[3][0] = 20;
        adjacency_matrix[1][2] = adjacency_matrix[2][1] = 7;
        adjacency_matrix[2][3] = adjacency_matrix[3][2] = 20;
        adjacency_matrix[3][4] = adjacency_matrix[4][3] = 15;

        std::vector<Edge>edge_list = {{0,1,100},{0,4,10},{0,3,20},{1,2,7},{2,3,20},{3,4,15}};

        Converter adjacency_list_converter = Converter(adjacency_list);
        Converter adjacency_matrix_converter = Converter(adjacency_matrix, 5);
        Converter edge_list_converter = Converter(edge_list, false);

        std::pair<long long, std::vector<int>> expected = {47,{0,3,2,1}};

        SECTION("Dijkstra for high density") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 0, 1, expected, dijkstra_high_density);
        }
        SECTION("Dijkstra for low density") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 0, 1, expected, dijkstra_low_density);
        }
        SECTION("Bellman-Ford") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 0, 1, expected, bellman_for_two_vertices);
        }
        SECTION("Floyd-Warshall") {
            std::vector<std::vector<std::pair<long long, std::vector<int>>>>expected_mtx(n);
            expected_mtx[0] = {{0, {0}}, {47, {0,3,2,1}}, {40, {0,3,2}}, {20, {0,3}}, {10, {0,4}}};
            expected_mtx[1] = {{47, {1,2,3,0}}, {0,{1}}, {7,{1,2}}, {27, {1,2,3}}, {42, {1,2,3,4}}};
            expected_mtx[2] = {{40, {2,3,0}}, {7,{2,1}}, {0,{2}}, {20,{2,3}}, {35,{2,3,4}}};
            expected_mtx[3] = {{20,{3,0}}, {27, {3,2,1}}, {20, {3,2}}, {0,{3}}, {15, {3,4}}};
            expected_mtx[4] = {{10, {4,0}}, {42, {4,3,2,1}}, {35,{4,3,2}}, {15,{4,3}}, {0, {4}}};
            test_floyd_warshall(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, expected_mtx);
        }


        for (int i = 0; i < n; i++) {
            delete adjacency_matrix[i];
        }
        delete adjacency_matrix;


    }

    SECTION("two components") {
        int n = 5;
        std::vector<std::vector<std::pair<int, long long>>>adjacency_list(n);
        adjacency_list[0] = {{4,10}, {3,20}};
        adjacency_list[1] = {{2,7}};
        adjacency_list[2] = {{1,7}};
        adjacency_list[3] = {{0,20},{4,15}};
        adjacency_list[4] = {{0,10}, {3,15}};

        long long** adjacency_matrix = new long long*[n];
        for (int i = 0; i < n; i++) {
            adjacency_matrix[i] = new long long[n];
            for (int j = 0; j < n; j++) {
                adjacency_matrix[i][j] = 0;
            }
        }
        adjacency_matrix[0][3] = adjacency_matrix[3][0] = 20;
        adjacency_matrix[0][4] = adjacency_matrix[4][0] = 10;
        adjacency_matrix[1][2] = adjacency_matrix[2][1] = 7;
        adjacency_matrix[3][4] = adjacency_matrix[4][3] = 15;

        std::vector<Edge>edge_list = {{0,3,20},{0,4,10},{1,2,7},{3,4,15}};
        Converter adjacency_list_converter = Converter(adjacency_list);
        Converter adjacency_matrix_converter = Converter(adjacency_matrix, 5);
        Converter edge_list_converter = Converter(edge_list, false);
        std::pair<long long, std::vector<int>> expected = {-1,std::vector<int>()};

        SECTION("Dijkstra for high density") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 2, 4, expected, dijkstra_high_density);
        }
        SECTION("Dijkstra for low density") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 2, 4, expected, dijkstra_low_density);
        }
        SECTION("Bellman-Ford") {
            test_shortest_path(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, 2, 4, expected, bellman_for_two_vertices);
        }
        SECTION("Floyd-Warshall") {
            std::vector<std::vector<std::pair<long long, std::vector<int>>>>expected_mtx(n);
            expected_mtx[0] = {{0, {0}}, {-1, {}}, {-1, {}}, {20, {0,3}}, {10, {0,4}}};
            expected_mtx[1] = {{-1, {}}, {0,{1}}, {7,{1,2}}, {-1, {}}, {-1, {}}};
            expected_mtx[2] = {{-1, {}}, {7,{2,1}}, {0,{2}}, {-1,{}}, {-1,{}}};
            expected_mtx[3] = {{20,{3,0}}, {-1, {}}, {-1, {}}, {0,{3}}, {15, {3,4}}};
            expected_mtx[4] = {{10, {4,0}}, {-1, {}}, {-1,{}}, {15,{4,3}}, {0, {4}}};
            test_floyd_warshall(adjacency_list_converter, adjacency_matrix_converter, edge_list_converter, expected_mtx);
        }

        for (int i = 0; i < n; i++) {
            delete adjacency_matrix[i];
        }
        delete adjacency_matrix;
    }
}

