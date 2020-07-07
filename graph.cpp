/*
 * graph.cpp -- simple graph implementation storing it as a list of edges
 * 	includes implementation for finding connected components of a symmetric graph
 * 	and an implementation of finding a maximum matching for a bipartite graph
 * 		or DAG using the Hopcroft-Karp algorithm
 * 
 * Copyright 2018 Daniel Kondor <kondor.dani@gmail.com>
 * 
 * Hopcroft-Karp algorithm adapted from
 * http://www.geeksforgeeks.org/hopcroft-karp-algorithm-for-maximum-matching-set-2-implementation/
 * (no license provided for it)
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following disclaimer
 *   in the documentation and/or other materials provided with the
 *   distribution.
 * * Neither the name of the  nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
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
 */



#include "graph.h"

/* 3-4. common interface for creating the edges after copying the out-edges to the edges array and using a separate in-edges array */
int graph::create_graph(unsigned int* in_edges, std::unordered_map<unsigned int,unsigned int>* ids_map) {
	/* note: clear has been called before at this point and out-edges are stored in the edges array */
	/* 1. sort edges using zip iterators */
	zi::zip_it<unsigned int*,unsigned int*> e = zi::make_zip_it(in_edges,edges);
	zi::zip_it<unsigned int*,unsigned int*> end = zi::make_zip_it(in_edges+nedges,edges+nedges);
	std::sort(e,end,[](const auto& a, const auto& b) {
			return a.first < b.first;
		});
	/* 2. call the general interface creating the node objects */
	return create_graph_sorted(e,end,ids_map,false);
}
		
		
/* 4. create graph of edges supplied as two C-style arrays
 * edges point as e1[i] -> e2[i]
 * the e2 array is "taken over" by this instance, it should not be freed by the caller later
 * both arrays are modified by sorting them */
int graph::create_graph_arrays(unsigned int* e1, unsigned int* e2, size_t e_size, std::unordered_map<unsigned int, unsigned int>* ids_map) {
	clear();
	if( ! (e1 && e2 && e_size) ) return 1;
	edges = e2;
	nedges = e_size;
	edges_size = e_size;
	return create_graph(e1,ids_map);
}
/* 4. create graph of edges supplied as two std::vector<unsiged int>
 * edges point as e1[i] -> e2[i]
 * the e2 vector is "taken over", i.e. it is swapped by a local vector */
int graph::create_graph_vectors(std::vector<unsigned int>& e1, std::vector<unsigned int>& e2, std::unordered_map<unsigned int, unsigned int>* ids_map) {
	clear();
	if(e1.size() != e2.size() || e1.size() == 0) return 0;
	edges_vect.swap(e2);
	edges_owned = false;
	nedges = edges_vect.size();
	edges_size = edges_vect.capacity();
	edges = edges_vect.data();
	return create_graph(e1.data(),ids_map);
}

/* Helper function for the above to handle reading the graph from a file.
 * If partitioned == true, it expects already partitioned input; otherwise
 * it sorts the input.
 * If ids != null, it replaces ids in the edges to form a continuous
 * range 0...N-1, and the real ids from the file are stored in the supplied
 * vector */
int graph::read_graph(FILE* f, bool partitioned, std::vector<unsigned int>* ids) {
	if(!f) return 1;
	if(ids) ids->clear();
	std::unordered_map<unsigned int,unsigned int> ids_map;
	std::unordered_map<unsigned int,unsigned int>* p_ids_map = 0;
	if(ids) p_ids_map = &ids_map;
	
	tsv_iterator itf(f,0,0,true);
	
	int ret;
	if(partitioned) ret = create_graph_partitioned(itf,tsv_iterator_sentinel(),p_ids_map);
	else ret = create_graph_copy_sort(itf,tsv_iterator_sentinel(),p_ids_map);
	if(ret || itf.is_error()) {/* check and signal parse errors */
		clear();
		fprintf(stderr,"graph::read_graph(): error reading the graph!\n");
		return 1;
	}
	
	if(ids) {
		ids->resize(nnodes,0);
		for(auto it = ids_map.begin();it != ids_map.end();++it) (*ids)[it->second] = it->first;
	}
	return 0;
}



/* find connected components in a symmetric graph */
unsigned int graph::find_sccs(std::vector<unsigned int>& sccids) const {
	sccids.clear();
	sccids.resize(nnodes,nnodes);
	
	unsigned int sccid = 0;
	std::vector<unsigned int> search_next;
	std::vector<unsigned int> search_cur;
	for(unsigned int i=0;i<nnodes;i++) if(sccids[i] == nnodes) {
		/* perform search from node i */
		sccids[i] = sccid;
		search_cur.push_back(i);
		do {
			for(unsigned int j : search_cur) {
				/* go through the neighbors of j and add them to the current
				 * SCC and the list to be searched next if not there already */
				for(edges_iterator it(this,j,true);it != edges_end();++it) {
					unsigned int k = it->second;
					if(sccids[k] != sccid) {
						if(sccids[k] != nnodes) {
							fprintf(stderr,"Inconsistent graph while searching for connected components!\nNote: symmetric graph required\n");
							sccids.clear();
							return 0;
						}
						sccids[k] = sccid;
						search_next.push_back(k);
					}
				}
			}
			search_cur.clear();
			search_cur.swap(search_next);
		} while(search_cur.size() > 0);
		sccid++;
	}
	return sccid;
}


		
/* calculate maximum matching using the Hopcroft-Karp algorithm
 * store the edges in the maximum matching in the provided vector (res)
 * if use_r == false, use a version of dfs without recursion
 * (using recursion might result in stack overflow even for
 * moderate size graphs as well)
 * 
 * code adapted from
 * http://www.geeksforgeeks.org/hopcroft-karp-algorithm-for-maximum-matching-set-2-implementation
 */
int graph::maxmatch_hk(std::vector<std::pair<unsigned int, unsigned int> >& res, bool use_r) const {
	// Initialize result
    res.clear();
    if(!nedges) return 0;
    if(nedges == 1) { res.push_back(*edges_begin()); return 0; }
	
	// pairU[u] stores pair of u in matching where u
    // is a vertex on left side of Bipartite Graph.
    // If u doesn't have any pair, then pairU[u] is NIL
    unsigned int* pairU = (unsigned int*)malloc(sizeof(unsigned int)*(3*nnodes+nnodes_v+4));
    if(!pairU) return 1;
    uint64_t* dfsQ = 0;
    if(!use_r) {
		dfsQ = (uint64_t*)malloc(sizeof(uint64_t)*(nnodes+1));
		if(!dfsQ) {
			free(pairU);
			return 1;
		}
	}
    
 
    // pairV[v] stores pair of v in matching. If v
    // doesn't have any pair, then pairU[v] is NIL
    unsigned int* pairV = pairU + nnodes + 1;
 
    // dist[u] stores distance of left side vertices
    // dist[u] is one more than dist[u'] if u is next
    // to u'in augmenting path
    unsigned int* dist = pairV + nnodes_v + 1;
    
    // working space for bfs and non-recursive dfs as the queue to use
    // as each node can be added maximum once to the queue, it's size is maximum nnodes+1
    unsigned int* Q = dist + nnodes + 1;
 
    // Initialize NIL as pair of all vertices
    for(unsigned int u=0; u<=nnodes; u++) pairU[u] = NIL;
    for(unsigned int v=0; v<=nnodes_v; v++) pairV[v] = NIL;
 
 
    // Keep updating the result while there is an
    // augmenting path.
    while (bfs(pairU,pairV,dist,Q))
        // Find a free vertex
        for (unsigned int u=0; u<nnodes; u++)
            // If current vertex is free and there is
            // an augmenting path from current vertex
            if (pairU[u]==NIL) {
				if(use_r) dfs(u,pairU,pairV,dist);
				else dfs_nor(u,pairU,pairV,dist,Q,dfsQ);
			}
    
    
    // save result
    for(unsigned int u=0;u<nnodes;u++) if(pairU[u] != NIL) res.push_back(std::make_pair(u,pairU[u]));
    free(pairU);
    return 0;
}




	
/* bfs for maximum matching
 * Returns true if there is an augmenting path, else returns false */
bool graph::bfs(unsigned int* pairU, unsigned int* pairV, unsigned int* dist, unsigned int* Q) const {
	const unsigned int INF = NIL+1; // note: path lengths will be < nnodes -> < NIL
	unsigned int ql = 0; // length of the queue (i.e. u nodes to consider)
 
    // First layer of vertices (set distance as 0)
    for (unsigned int u=0; u<nnodes; u++)
    {
        // If this is a free vertex, add it to queue
        if (pairU[u]==NIL)
        {
            // u is not matched
            dist[u] = 0;
            if(ql >= nnodes+1) handle_error("graph::bfs(): maximum queue length reached!\n"); // this should not happen, just to make sure
            Q[ql] = u;
            ql++;
        }
 
        // Else set distance as infinite so that this vertex
        // is considered next time
        else dist[u] = INF;
    }
 
    // Initialize distance to NIL as infinite
    dist[NIL] = INF;
 
    // Q is going to contain vertices of left side only. 
    for(unsigned int i=0;i<ql;i++) {
        // Dequeue a vertex
        int u = Q[i];
 
        // If this node is not NIL and can provide a shorter path to NIL
        if (dist[u] < dist[NIL])
        {
            // Get all adjacent vertices of the dequeued vertex u
            for (edges_iterator it = edges_begin_n(u); it != edges_end(); ++it)
            {
                unsigned int v = it->second;
 
                // If pair of v is not considered so far
                // (v, pairV[V]) is not yet explored edge.
                if (dist[pairV[v]] == INF)
                {
                    // Consider the pair and add it to queue
                    dist[pairV[v]] = dist[u] + 1;
                    if(ql >= nnodes+1) handle_error("graph::bfs(): maximum queue length reached!\n"); // this should not happen, just to make sure
		            Q[ql] = pairV[v];
		            ql++;
                }
            }
        }
    }
 
    // If we could come back to NIL using alternating path of distinct
    // vertices then there is an augmenting path
    return (dist[NIL] != INF);
}

/* dfs for maximum matching
 * Returns true if there is an augmenting path beginning with vertex u (recursive call) */
bool graph::dfs(unsigned int u, unsigned int* pairU, unsigned int* pairV, unsigned int* dist) const {
    if (u != NIL) {
		const unsigned int INF = NIL+1; // note: path lengths will be < nnodes -> < NIL
		for (edges_iterator it = edges_begin_n(u); it != edges_end(); ++it) {
            // Adjacent to u
            unsigned int v = it->second;
 
            // Follow the distances set by BFS
            if (dist[pairV[v]] == dist[u]+1)
            {
                /* If dfs for pair of v also returns true
                 * note: recursion can continue for nnodes depth, this can result in stack overflow */
                if (dfs(pairV[v], pairU, pairV, dist) == true)
                {
                    pairV[v] = u;
                    pairU[u] = v;
                    return true;
                }
            }
        }
 
        // If there is no augmenting path beginning with u.
        dist[u] = INF;
        return false;
    }
    return true;
}

/* dfs for maximum matching
 * Returns true if there is an augmenting path beginning with vertex u (non-recursive version) */
void graph::dfs_nor(unsigned int u, unsigned int* pairU, unsigned int* pairV, unsigned int* dist, unsigned int* path, uint64_t* ix) const {

	const unsigned int INF = NIL+1; // note: path lengths will be < nnodes -> < NIL
	//~ unsigned int* ix = path + nnodes + 1; -- note: ix provided as separate array as it needs to be 64-bit
	
	unsigned int l = 0;
	uint64_t j = idx[u];
	while(1) {
		// look at edges starting from u
		bool v_found = false;
		bool nil_found = false;
		unsigned int v;
		for(;j<nedges && j<idx[u] + outdeg[u];j++) {
            // Adjacent to u
            v = edges[j];
            // Follow the distances set by BFS
            if (dist[pairV[v]] == dist[u]+1)
            {
				if(pairV[v] == NIL) {
					nil_found = true;
					break;
				}
				
                path[l] = u;
                ix[l] = j;
                l++;
                u = pairV[v];
                j = idx[u];
                v_found = true;
                break;
            }
        }
        
        if(v_found) continue; // continue the recursion
        
        if(nil_found) {
			// found a good path, go back and update
			while(1) {
				pairV[v] = u;
				pairU[u] = v;
				if(l == 0) break;
				l--;
				u = path[l];
				j = ix[l];
				v = edges[j];
			}
			break; // and end the run
		}
        
		// did not find anything for u, step back and continue the loop
		dist[u] = INF;
		if(l == 0) break; // in this case, did not find anything in the search
		l--;
		u = path[l];
		j = ix[l] + 1;
	}
}




unsigned int graph::real_deg(int n) const {
	uint64_t deg;
	if(n+1 < nnodes) deg = idx[n+1] - idx[n];
	else deg = nedges - idx[n];
	return (unsigned int)deg; /* note: assuming degrees are always < 2^32 */
}

/* make the graph symmetric by simply growing the edges structure */
int graph::make_symmetric() {
	/* 1. check if the graph is symmetric, increase node degrees */
	size_t extra_size = 0;
	for(int i=0;i<nnodes;i++) {
		/* note: degrees should not be relied during running, they are being updated
		 * calculate the current real degree */
		unsigned int deg = real_deg(i);
		for(unsigned int k=0;k<deg;k++) {
			unsigned int j = edges[idx[i] + k]; /* i -- j edge */
			/* check if there is a j -- i edge as well */
			unsigned int degj = real_deg(j);
			unsigned int* it = std::lower_bound(edges + idx[j], edges + idx[j] + degj, i);
			if(it >= edges + idx[j] + degj || *it != i) {
				/* not found, will need to add j -- i edge
				 * at this point only increase the degree of j */
				outdeg[j]++;
				extra_size++;
			}
		}
	}
	
	if(extra_size == 0) return 0; /* no extra edges need to be added */
	/* 2. realloc the edges array, alloc temporary space for old node degrees */
	size_t old_size = nedges;
	if(edges_owned) {
		if(nedges + extra_size > edges_size) {
			unsigned int* tmp1 = (unsigned int*)realloc(edges,(nedges + extra_size)*sizeof(unsigned int));
			if(!tmp1) return 1;
			edges = tmp1;
			nedges += extra_size;
			edges_size = nedges;
		}
		else nedges += extra_size;
	}
	else {
		edges_vect.resize(nedges + extra_size);
		edges = edges_vect.data();
		nedges += extra_size;
		edges_size = nedges;
	}
	
	unsigned int* tmp = (unsigned int*)malloc(sizeof(unsigned int)*nnodes);
	if(!tmp) return 1;
	
	/* 3. move all edges according to the new degrees,
	 * store old degrees to be able to add the new edges
	 * need to iterate in reverse
	 * 	-- note: i needs to be signed for the following loop to work! */
	size_t newidx = nedges;
	for(unsigned int i1 = nnodes; i1 > 0; i1--) {
		unsigned int i = i1-1;
		/* calculate the old degree of i */
		unsigned int deg = old_size - idx[i]; /* again, assuming that all degrees < 2^32 */
		tmp[i] = deg;
		if(outdeg[i] > newidx) {
			fprintf(stderr,"graph::make_symmetric(): inconsistency found while adding edges!\n");
			free(tmp);
			return 1;
		}
		newidx -= outdeg[i]; /* new start for edges */
		if(newidx != idx[i]) for(unsigned int j=deg;j>0;j--) edges[newidx+j-1] = edges[idx[i]+j-1];
		outdeg[i] = deg; /* reset the old degree so we know how to grow */
		old_size = idx[i];
		idx[i] = newidx; /* save the new start of edges */
	}
	
	/* 4. add missing edges */
	size_t edges_added = 0;
	for(unsigned int i = 0; i < nnodes; i++) {
		unsigned int deg = tmp[i]; /* degrees saved previously */
		for(unsigned int k = 0; k < deg; k++) {
			unsigned int j = edges[idx[i] + k];
			size_t degj = tmp[j];
			unsigned int* it = std::lower_bound(edges + idx[j], edges + idx[j] + degj, i);
			if(it >= edges + idx[j] + degj || *it != i) {
				/* add j -- i edge */
				unsigned int max_deg = real_deg(j);
				if(outdeg[j] >= max_deg) {
					fprintf(stderr,"graph::make_symmetric(): inconsistency found while adding edges!\n");
					free(tmp);
					return 1;
				}
				edges[idx[j] + outdeg[j]] = i;
				outdeg[j]++;
				edges_added++;
			}
		}
	}
	
	free(tmp);
	if(edges_added != extra_size) {
		fprintf(stderr,"graph::make_symmetric(): inconsistency found while adding edges!\n");
		return 1;
	}
	return 0;
}


