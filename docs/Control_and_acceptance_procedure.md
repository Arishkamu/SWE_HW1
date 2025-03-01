## Control and acceptance procedure

The library should meet all requirements described in Requirements.md
Tests to meet the requirements:
- #### Test group 1 (Requirements 1,2) 
    Check if all algorithms give the expected answer on 4 types of graph with $\geq$ 5 vertices: empty, full, a complex connected graph with cycles, a graph with $\geq$ 2 non-trivial components (and start and finish are in different ones). All of the graphs should be examined in all of the three representations and give the same answer. 
    A* should be tested with Euclidean distance as heuristic.
- #### Script (Requirement 3)
    There should be a correct shell script/sequence of commands of $\leq$ 8 lines of code that will allow using library. A video of installation and usage in Ubuntu Linux and an arbitrary IDE with successful compilation should be recorded.
- #### Test group 2 (Requirement 4)
    Test generates 100 random graphs for n = 10000 in adjacency list representation and checks if time for computing Floyd's and Johnson's (asymptotically the heaviest algorithms with V^3 and VElogV respectively) are both less than 2 minutes.
