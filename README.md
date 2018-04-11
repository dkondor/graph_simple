# graph_simple
graph implementation with focus on large graphs with efficient storage (4 bytes / edge, 12 bytes / node,
with additional storage requirement of 4 bytes / edge if reading unsorted input)

Currently it supports reading from an edgelist format and two graph problems:
- calculation of connected components in a symmetric (undirected graph)
- calculation of a maximum matching in a directed acyclic graph or bipartite graph
using the Hopcroft-Karp algorithm, with part of the code adapted from
http://www.geeksforgeeks.org/hopcroft-karp-algorithm-for-maximum-matching-set-2-implementation/


Example programs for both are provided.

Compilation:
- g++ -o sccs sccs.cpp graph.cpp -O3 -march=native -std=gnu++14
- g++ -o mm_hk maxmatch_hk.cpp graph.cpp -O3 -march=native -std=gnu++14

Both programs read the graph from stdin and write the result to stdout. For the
SCCs program, the result is pairs of node IDs and SCC IDs. For the maximum matching
program the output is the list of edges found to be part of the matching. Both are
expected to run for graphs with billions of edges in reasonable time.
