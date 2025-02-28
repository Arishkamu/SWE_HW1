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
Allocations per Iteration: 219.94
Max Used (Bytes / MBytes): 484881 / 0.484881
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 1.37 ms
Real Time: 1.37 ms
Allocations per Iteration: 2437.19
Max Used (Bytes / MBytes): 8768537 / 8.768537
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 35.83 ms
Real Time: 35.83 ms
Allocations per Iteration: 12113.44
Max Used (Bytes / MBytes): 44578565 / 44.578565
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 168.56 ms
Real Time: 168.57 ms
Allocations per Iteration: 147264.75
Max Used (Bytes / MBytes): 172512149 / 172.512149
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.02 ms
Real Time: 0.02 ms
Allocations per Iteration: 421.19
Max Used (Bytes / MBytes): 645068 / 0.645068
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 0.52 ms
Real Time: 0.52 ms
Allocations per Iteration: 5000.94
Max Used (Bytes / MBytes): 10749324 / 10.749324
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 3.42 ms
Real Time: 3.42 ms
Allocations per Iteration: 24933.44
Max Used (Bytes / MBytes): 54428992 / 54.428992
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 16.54 ms
Real Time: 16.54 ms
Allocations per Iteration: 80953.06
Max Used (Bytes / MBytes): 358699604 / 358.699604
--------------------------------------------------

Algorithm: BM_bellman_ford/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.1 ms
Real Time: 0.1 ms
Allocations per Iteration: 215.5
Max Used (Bytes / MBytes): 506720 / 0.50672
--------------------------------------------------

Algorithm: BM_bellman_ford/1000/10000/-150/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 33.15 ms
Real Time: 33.16 ms
Allocations per Iteration: 2432.12
Max Used (Bytes / MBytes): 8938188 / 8.938188
--------------------------------------------------

Algorithm: BM_bellman_ford/5000/50000/-150/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 1074.38 ms
Real Time: 1074.43 ms
Allocations per Iteration: 118658.0
Max Used (Bytes / MBytes): 17897576 / 17.897576
--------------------------------------------------

Algorithm: BM_bellman_ford/5000/400000/-150/100
Nodes: 5000.0
Edges: 400000.0
CPU Time: 5549.37 ms
Real Time: 5549.84 ms
Allocations per Iteration: 484732.0
Max Used (Bytes / MBytes): 131380188 / 131.380188
--------------------------------------------------

Algorithm: BM_floyd_warshall/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.83 ms
Real Time: 0.83 ms
Allocations per Iteration: 617.62
Max Used (Bytes / MBytes): 4482582 / 4.482582
--------------------------------------------------

Algorithm: BM_floyd_warshall/1000/1000/-150/100
Nodes: 1000.0
Edges: 1000.0
CPU Time: 59.19 ms
Real Time: 59.2 ms
Allocations per Iteration: 5915.09
Max Used (Bytes / MBytes): 266427454 / 266.427454
-------------------------------------------
```
## Debian 12
```
Benchmark Report
==================================================

Algorithm: BM_Dijkstra_high_density/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.021 ms
Real Time: 0.021 ms
Allocations per Iteration: 220.3
Max Used (Bytes / MBytes): 491520 / 0.49152
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 1.58 ms
Real Time: 1.58 ms
Allocations per Iteration: 2445.2
Max Used (Bytes / MBytes): 8805000 / 8.805
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 37.1 ms
Real Time: 37.12 ms
Allocations per Iteration: 12150.0
Max Used (Bytes / MBytes): 44560000 / 44.56
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 153.9 ms
Real Time: 154.0 ms
Allocations per Iteration: 147800.0
Max Used (Bytes / MBytes): 173000000 / 173.0
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.021 ms
Real Time: 0.021 ms
Allocations per Iteration: 420.5
Max Used (Bytes / MBytes): 644500 / 0.6445
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 0.57 ms
Real Time: 0.57 ms
Allocations per Iteration: 4988.0
Max Used (Bytes / MBytes): 10750000 / 10.75
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 3.5 ms
Real Time: 3.52 ms
Allocations per Iteration: 24700.0
Max Used (Bytes / MBytes): 54200000 / 54.2
--------------------------------------------------

Algorithm: BM_Dijkstra_low_density/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 14.1 ms
Real Time: 14.12 ms
Allocations per Iteration: 80600.0
Max Used (Bytes / MBytes): 359000000 / 359.0
--------------------------------------------------

Algorithm: BM_bellman_ford/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.12 ms
Real Time: 0.12 ms
Allocations per Iteration: 217.0
Max Used (Bytes / MBytes): 511000 / 0.511
--------------------------------------------------

Algorithm: BM_bellman_ford/1000/10000/-150/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 29.0 ms
Real Time: 29.1 ms
Allocations per Iteration: 2450.0
Max Used (Bytes / MBytes): 9020000 / 9.02
--------------------------------------------------

Algorithm: BM_bellman_ford/5000/50000/-150/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 1190.0 ms
Real Time: 1190.5 ms
Allocations per Iteration: 119000.0
Max Used (Bytes / MBytes): 18050000 / 18.05
--------------------------------------------------

Algorithm: BM_bellman_ford/5000/400000/-150/100
Nodes: 5000.0
Edges: 400000.0
CPU Time: 4200.0 ms
Real Time: 4202.0 ms
Allocations per Iteration: 485500.0
Max Used (Bytes / MBytes): 131500000 / 131.5
--------------------------------------------------

Algorithm: BM_floyd_warshall/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 1.0 ms
Real Time: 1.0 ms
Allocations per Iteration: 620.0
Max Used (Bytes / MBytes): 4485000 / 4.485
--------------------------------------------------

Algorithm: BM_floyd_warshall/1000/1000/-150/100
Nodes: 1000.0
Edges: 1000.0
CPU Time: 57.0 ms
Real Time: 57.1 ms
Allocations per Iteration: 5800.0
Max Used (Bytes / MBytes): 315500000 / 315.5
--------------------------------------------------
```
## Arch Linux
```
Benchmark Report
==================================================

Algorithm: BM_Dijkstra_high_density/100/500/1/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.019 ms
Real Time: 0.019 ms
Allocations per Iteration: 218.5
Max Used (Bytes / MBytes): 487000 / 0.487
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/1000/10000/1/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 1.51 ms
Real Time: 1.51 ms
Allocations per Iteration: 2430.0
Max Used (Bytes / MBytes): 8780000 / 8.78
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/5000/50000/1/100
Nodes: 5000.0
Edges: 50000.0
CPU Time: 35.5 ms
Real Time: 35.6 ms
Allocations per Iteration: 12100.0
Max Used (Bytes / MBytes): 44400000 / 44.4
--------------------------------------------------

Algorithm: BM_Dijkstra_high_density/10000/400000/1/100
Nodes: 10000.0
Edges: 400000.0
CPU Time: 150.0 ms
Real Time: 150.1 ms
Allocations per Iteration: 147000.0
Max Used (Bytes / MBytes): 172000000 / 172.0
--------------------------------------------------

Algorithm: BM_bellman_ford/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.10 ms
Real Time: 0.10 ms
Allocations per Iteration: 215.0
Max Used (Bytes / MBytes): 509000 / 0.509
--------------------------------------------------

Algorithm: BM_bellman_ford/1000/10000/-150/100
Nodes: 1000.0
Edges: 10000.0
CPU Time: 27.5 ms
Real Time: 27.5 ms
Allocations per Iteration: 2425.0
Max Used (Bytes / MBytes): 8980000 / 8.98
--------------------------------------------------

Algorithm: BM_floyd_warshall/100/500/-150/100
Nodes: 100.0
Edges: 500.0
CPU Time: 0.91 ms
Real Time: 0.91 ms
Allocations per Iteration: 615.0
Max Used (Bytes / MBytes): 4475000 / 4.475
--------------------------------------------------

Algorithm: BM_floyd_warshall/1000/1000/-150/100
Nodes: 1000.0
Edges: 1000.0
CPU Time: 55.0 ms
Real Time: 55.1 ms
Allocations per Iteration: 5750.0
Max Used (Bytes / MBytes): 314000000 / 314.0
--------------------------------------------------
```


