# Benchmark reports on different OS
These file contains different bench reports ran on 3 OS (Linux distributions).
File `docs/bench_report.txt` reflects benchmark results from CI. 

Platform parameters:
```
Model CPU name: AMD Ryzen 7 4800H with Radeon Graphics
Caches (sum of all):      
  L1d:                    256 KiB (8 instances)
  L1i:                    256 KiB (8 instances)
  L2:                     4 MiB (8 instances)
  L3:                     8 MiB (2 instances)
CPU max MHz:          2900,0000
CPU min MHz:          1400,0000
Memory:
    Type: DDR4
    Speed: 3200 MT/s
    RAM size: 16GB
```
## Preparation
To execute tests we used OS on a virtual machine (full virtualization). We built the object file, restarted the system, turned off battery safe mode, turned off auto updates:
```
sudo systemctl stop apt-daily.service
sudo systemctl stop apt-daily-upgrade.service
```
Ran with taskset -c 0 to attach our executable to a single core and avoid context switch:
```
taskset -c 0 ./Bench --benchmark_format=json > bench_result.json
``` 
## Ubuntu jammy
```
Benchmark Report
==================================================

Algorithm: BM_Dijkstra_high_density/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.02 ms
Real Time: 0.02 ms
Allocations per Iteration: 220.06
Max Used (Bytes / MBytes): 489109 / 0.489109
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 1.47 ms
Real Time: 1.47 ms
Allocations per Iteration: 2436.75
Max Used (Bytes / MBytes): 8768109 / 8.768109
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 34.69 ms
Real Time: 34.69 ms
Allocations per Iteration: 12112.56
Max Used (Bytes / MBytes): 44436117 / 44.436117
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 143.7 ms
Real Time: 143.7 ms
Allocations per Iteration: 119801.0
Max Used (Bytes / MBytes): 185397605 / 185.397605
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.02 ms
Real Time: 0.02 ms
Allocations per Iteration: 411.75
Max Used (Bytes / MBytes): 639408 / 0.639408
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 0.54 ms
Real Time: 0.54 ms
Allocations per Iteration: 4949.25
Max Used (Bytes / MBytes): 10723984 / 10.723984
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 3.5 ms
Real Time: 3.5 ms
Allocations per Iteration: 24793.12
Max Used (Bytes / MBytes): 54287752 / 54.287752
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 14.01 ms
Real Time: 14.01 ms
Allocations per Iteration: 80589.44
Max Used (Bytes / MBytes): 358462468 / 358.462468
--------------------------------------------------

Algorithm: BM_bellman_ford/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.11 ms
Real Time: 0.11 ms
Allocations per Iteration: 216.0
Max Used (Bytes / MBytes): 506948 / 0.506948
--------------------------------------------------

Algorithm: BM_bellman_ford/1000/10000/-150/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 26.85 ms
Real Time: 26.85 ms
Allocations per Iteration: 2433.31
Max Used (Bytes / MBytes): 8963504 / 8.963504
--------------------------------------------------

Algorithm: BM_bellman_ford/5000/50000/-150/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 1107.41 ms
Real Time: 1107.46 ms
Allocations per Iteration: 118686.0
Max Used (Bytes / MBytes): 17957568 / 17.957568
--------------------------------------------------

Algorithm: BM_bellman_ford/5000/400000/-150/100
Nodes: 5000.0
Edges: 400000.0
CPU Time: 4098.71 ms
Real Time: 4098.93 ms
Allocations per Iteration: 484724.0
Max Used (Bytes / MBytes): 131323712 / 131.323712
--------------------------------------------------

Algorithm: BM_floyd_warshall/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.79 ms
Real Time: 0.79 ms
Allocations per Iteration: 618.75
Max Used (Bytes / MBytes): 4480726 / 4.480726
--------------------------------------------------

Algorithm: BM_floyd_warshall/1000/1000/-150/100
Nodes: 1000.0
Edges: 1000.0
CPU Time: 44.14 ms
Real Time: 44.14 ms
Allocations per Iteration: 5838.08
Max Used (Bytes / MBytes): 290611022 / 290.611022
--------------------------------------------------

Algorithm: BM_bfs/100/500/1/1
Nodes: 100.0
Edges: 500.0
CPU Time: 0.0 ms
Real Time: 0.0 ms
Allocations per Iteration: 220.12
Max Used (Bytes / MBytes): 499208 / 0.499208
--------------------------------------------------

Algorithm: BM_bfs/1000/10000/1/1
Nodes: 1000.0
Edges: 10000.0
CPU Time: 0.08 ms
Real Time: 0.08 ms
Allocations per Iteration: 2444.5
Max Used (Bytes / MBytes): 8814152 / 8.814152
--------------------------------------------------

Algorithm: BM_bfs/5000/50000/1/1
Nodes: 5000.0
Edges: 50000.0
CPU Time: 0.49 ms
Real Time: 0.49 ms
Allocations per Iteration: 12156.19
Max Used (Bytes / MBytes): 44891040 / 44.89104
--------------------------------------------------

Algorithm: BM_bfs/10000/400000/1/1
Nodes: 10000.0
Edges: 400000.0
CPU Time: 2.65 ms
Real Time: 2.65 ms
Allocations per Iteration: 44403.19
Max Used (Bytes / MBytes): 331313840 / 331.31384
--------------------------------------------------

Algorithm: BM_lee/100/500/1/1
Nodes: 100.0
Edges: 500.0
CPU Time: 0.0 ms
Real Time: 0.0 ms
Allocations per Iteration: 219.31
Max Used (Bytes / MBytes): 491192 / 0.491192
--------------------------------------------------

Algorithm: BM_lee/1000/10000/1/1
Nodes: 1000.0
Edges: 10000.0
CPU Time: 0.07 ms
Real Time: 0.07 ms
Allocations per Iteration: 2443.12
Max Used (Bytes / MBytes): 8792368 / 8.792368
--------------------------------------------------

Algorithm: BM_lee/5000/50000/1/1
Nodes: 5000.0
Edges: 50000.0
CPU Time: 0.36 ms
Real Time: 0.36 ms
Allocations per Iteration: 12157.0
Max Used (Bytes / MBytes): 44487976 / 44.487976
--------------------------------------------------

Algorithm: BM_lee/10000/400000/1/1
Nodes: 10000.0
Edges: 400000.0
CPU Time: 1.57 ms
Real Time: 1.57 ms
Allocations per Iteration: 44397.44
Max Used (Bytes / MBytes): 330465184 / 330.465184
--------------------------------------------------

Algorithm: BM_naive_shortest_path/5/10/1/10
Nodes: 5.0
Edges: 10.0
CPU Time: 0.0 ms
Real Time: 0.0 ms
Allocations per Iteration: 34.06
Max Used (Bytes / MBytes): 14547 / 0.014547
--------------------------------------------------

Algorithm: BM_naive_shortest_path/7/15/1/10
Nodes: 7.0
Edges: 15.0
CPU Time: 0.0 ms
Real Time: 0.0 ms
Allocations per Iteration: 99.12
Max Used (Bytes / MBytes): 29079 / 0.029079
--------------------------------------------------

Algorithm: BM_naive_shortest_path/10/20/1/10
Nodes: 10.0
Edges: 20.0
CPU Time: 1.13 ms
Real Time: 1.13 ms
Allocations per Iteration: 845.19
Max Used (Bytes / MBytes): 177495 / 0.177495
--------------------------------------------------

Algorithm: BM_dag_shortest_paths/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.0 ms
Real Time: 0.0 ms
Allocations per Iteration: 118.38
Max Used (Bytes / MBytes): 191110 / 0.19111
--------------------------------------------------

Algorithm: BM_dag_shortest_paths/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 0.05 ms
Real Time: 0.05 ms
Allocations per Iteration: 1168.94
Max Used (Bytes / MBytes): 2635830 / 2.63583
--------------------------------------------------

Algorithm: BM_dag_shortest_paths/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 0.41 ms
Real Time: 0.41 ms
Allocations per Iteration: 5821.62
Max Used (Bytes / MBytes): 13557174 / 13.557174
--------------------------------------------------

Algorithm: BM_dag_shortest_paths/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 1.79 ms
Real Time: 1.79 ms
Allocations per Iteration: 13733.12
Max Used (Bytes / MBytes): 75073446 / 75.073446
--------------------------------------------------

Algorithm: BM_dag_shortest_path/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.01 ms
Real Time: 0.01 ms
Allocations per Iteration: 196.06
Max Used (Bytes / MBytes): 298453 / 0.298453
--------------------------------------------------

Algorithm: BM_dag_shortest_path/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 0.1 ms
Real Time: 0.1 ms
Allocations per Iteration: 2095.5
Max Used (Bytes / MBytes): 4324581 / 4.324581
--------------------------------------------------

Algorithm: BM_dag_shortest_path/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 0.66 ms
Real Time: 0.66 ms
Allocations per Iteration: 10269.25
Max Used (Bytes / MBytes): 21884549 / 21.884549
--------------------------------------------------

Algorithm: BM_dag_shortest_path/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 2.25 ms
Real Time: 2.25 ms
Allocations per Iteration: 23517.88
Max Used (Bytes / MBytes): 129828661 / 129.828661
--------------------------------------------------

Algorithm: BM_karp_algorithm/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.18 ms
Real Time: 0.18 ms
Allocations per Iteration: 430.38
Max Used (Bytes / MBytes): 3034590 / 3.03459
--------------------------------------------------

Algorithm: BM_karp_algorithm/1000/10000/-150/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 46.16 ms
Real Time: 46.16 ms
Allocations per Iteration: 4547.67
Max Used (Bytes / MBytes): 205087766 / 205.087766
--------------------------------------------------

Algorithm: BM_karp_algorithm/5000/50000/-150/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 4143.86 ms
Real Time: 4144.08 ms
Allocations per Iteration: 128714.0
Max Used (Bytes / MBytes): 322257754 / 322.257754
--------------------------------------------------
```


