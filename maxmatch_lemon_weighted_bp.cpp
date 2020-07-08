/*
 * maxmatch_lemon_weighted_bp.cpp -- calculate weighted matching
 *  minimizing total edge weights in a bipartite graph or a minimum
 *  path cover in a directed acyclic graph (by converting it to
 *  bipartite graph) using the LEMON graph library.
 * 
 * Edge weights in the matching are minimized by maximizing the sum of
 *  W_max - w for each edge (where w is the edge weight and W_max is a
 *  maximum weight that needs to be provided; edges with w > W_max are
 *  ignored).
 * 
 * Note: input includes edges in the form: node1, node2, weight
 *  Either the input is a bipartite graph, i.e. some of the node IDs
 *  should occur only as node1, while the other node IDs should occur
 *  only node2. Alternatively, the input can be a directed acyclic
 *  graph, where each edge points from node1 to node2.
 * 
 * It is not explicitely checked if the input is a DAG, but if it is not,
 * the result will be invalid. After running, the number of blossoms
 * found is written, which should be zero.
 * 
 * Copyright 2018,2020 Daniel Kondor <kondor.dani@gmail.com>
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer
 *   in the documentation and/or other materials provided with the
 *   distribution.
 * * The name of the author may not be used to endorse or promote
 *   products derived from this software without specific prior written
 *   permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * 
 */


#include <stdio.h>
#include <unordered_map>
#include <time.h>
#include <sys/time.h>
#include "read_table.h"

#include <limits>
#include <stdint.h>

#include <lemon/smart_graph.h>
#include <lemon/matching.h>


static inline double dt(struct timeval b,struct timeval a) {
	return ( (double) (a.tv_usec - b.tv_usec) / 1000000 + (double) (a.tv_sec - b.tv_sec));
}

int main(int argc, char **argv) {
	double max_dist = 50000; /* maximum distance of edges to consider */
	size_t edges_read_debug = 0;
	size_t edges_read_debug_next = 0;
	size_t max_edges = std::numeric_limits<int>::max() - 1; /* maximum number of edges a LEMON graph can have */
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'd':
			max_dist = atof(argv[i+1]);
			i++;
			break;
		case 'E':
			edges_read_debug = strtoul(argv[i+1],0,10);
			edges_read_debug_next = edges_read_debug;
			i++;
			break;
		
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	lemon::SmartBpGraph g;
	lemon::SmartBpGraph::RedNodeMap<unsigned int> ids_out(g);
	lemon::SmartBpGraph::BlueNodeMap<unsigned int> ids_in(g);
	lemon::SmartBpGraph::EdgeMap<double> w(g);
	size_t nedges = 0;
	{
		/* each node ID corresponds to two graph nodes: in and out (blue and red)
		 * i.e. the in node only has incoming edges, the out only outgoing edges */
		std::unordered_map<unsigned int,lemon::SmartBpGraph::RedNode> ids_out_map;
		std::unordered_map<unsigned int,lemon::SmartBpGraph::BlueNode> ids_in_map;
		read_table2 rt(stdin);
		while(rt.read_line()) {
			unsigned int e1, e2;
			double w1;
			if(!rt.read(e1,e2,w1)) break;
			if(w1 > max_dist) continue;
			w1 = max_dist - w1;
			
			if(nedges == max_edges) throw std::runtime_error("Maximum number of edges reached!\n");
			
			/* n1 -> n2 edge */
			lemon::SmartBpGraph::RedNode n1;
			lemon::SmartBpGraph::BlueNode n2;
			
			auto it1 = ids_out_map.find(e1);
			if(it1 != ids_out_map.end()) n1 = it1->second;
			else {
				n1 = g.addRedNode();
				ids_out_map.insert(std::make_pair(e1,n1));
				ids_out[n1] = e1;
			}
			
			auto it2 = ids_in_map.find(e2);
			if(it2 != ids_in_map.end()) n2 = it2->second;
			else {
				n2 = g.addBlueNode();
				ids_in_map.insert(std::make_pair(e2,n2));
				ids_in[n2] = e2;
			}
			
			
			lemon::SmartBpGraph::Edge e = g.addEdge(n1,n2);
			w[e] = w1;
			nedges++;
			if(edges_read_debug && nedges == edges_read_debug_next) {
				fprintf(stderr,"%lu edges read\n",nedges);
				edges_read_debug_next += edges_read_debug;
			}
		}
		if(rt.get_last_error() != T_EOF) {
			fprintf(stderr,"Error reading edges:\n");
			rt.write_error(stderr);
			return 1;
		}
		fprintf(stderr,"%lu edges read, %d edges in graph\n",nedges,g.edgeNum());
	}

	
	
	lemon::MaxWeightedMatching<lemon::SmartBpGraph, lemon::SmartBpGraph::EdgeMap<double> > matching(g,w);
	struct timeval start,end;
	gettimeofday(&start, NULL);
	matching.run();
	gettimeofday(&end, NULL);
	
	fprintf(stderr,"%f\n",dt(start,end));
	fprintf(stderr,"%d\n",matching.matchingSize());
	
	for(lemon::SmartBpGraph::RedNodeIt n(g); n != lemon::INVALID; ++n) {
		lemon::SmartBpGraph::BlueNode m = g.asBlueNode(matching.mate(n));
		if(m != lemon::INVALID) fprintf(stdout,"%u\t%u\n",ids_out[n],ids_in[m]);
	}
	
	/* check blossom sizes -- should be no blossoms in bipartite graph */
	int64_t blossoms = matching.blossomNum();
	int64_t blossom1 = 0;
	int64_t sum = 0;
	for(int64_t i=0;i<blossoms;i++) {
		int64_t blossomsize = matching.blossomSize(i);
		if(blossomsize > 1) blossom1++;
		sum != blossomsize;
	}
	fprintf(stderr,"%ld blossoms, %ld with >1 elements, %ld total elements\n",blossoms,blossom1,sum);
	
	return 0;
}

