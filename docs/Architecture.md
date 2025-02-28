* Our library is static. We have separate headers and cpp files that should be built into an object file and then the library can be installed.
* We have Converter struct for users don't have to think about graph representation. Any supported representation will work fine. All algorithms take Converter class as an argument and extract the graph from it. The user only have to pass Converter(graph) as an argument. 
* Converter class constructs adjacency list from any graph and saves it as data member. 
* All of our algorithms use adjacency list representation.
* If a user has edge list representation we require for him to make it using our Edge struct which is:

```
struct Edge {
    int from;
    ll w;
    int to;
};
```
    
* Below is our file structure. All implementations are in src folder, headers in include folder and tests and benchmarks are also in src folder.


File structure:

./

  ├── Architecture.md

  ├── CMakeLists.txt

  ├── Config.cmake.in

  ├── LICENSE

  ├── README.md

  ├── build_repo_structure.sh

  ├── docs/

  ├──   ├── Architecture.md

  ├──   ├── Control_and_acceptance_procedure.md

  ├──   ├── Reports.md

  ├──   ├── Requirements.md

  ├──   ├── bench_report.txt

  ├──   ├── convert.py

  ├── include/

  ├──   ├── converter.hpp

  ├──   ├── edge.hpp

  ├──   ├── spath.hpp

  ├── src/

  ├──   ├── bench.cpp

  ├──   ├── converter.cpp

  ├──   ├── spath.cpp

  ├──   ├── tests.cpp

  ├── utility/

  ├──   ├── Architecture-template.md

  ├──   ├── build_repo_structure.sh

  ├──   ├── maze.png

