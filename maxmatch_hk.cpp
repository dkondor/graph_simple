/*
 * maxmatch_hk.cpp -- calculate maximum matching in a directed acyclic graph
 * 	or in a bipartite graph
 * 
 * Copyright 2018 Daniel Kondor <kondor.dani@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <stdio.h>
#include <time.h>
#include "graph.h"

int main(int argc, char **argv)
{
	bool use_r = true;
	bool partitioned = false;
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'r':
			use_r = false;
			break;
		case 'p':
			partitioned = true;
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	std::vector<unsigned int> ids;
	graph g;
	g.read_graph(stdin,partitioned,&ids);
	fprintf(stderr,"%lu edges read, %lu total memory required\n",g.num_edges(),g.get_memory_hk());
	
	struct timespec t1;
	clock_gettime(CLOCK_MONOTONIC,&t1);
	std::vector<std::pair<unsigned int,unsigned int> > res;
	g.maxmatch_hk(res,use_r);
	struct timespec t2;
	clock_gettime(CLOCK_MONOTONIC,&t2);
	
	for(auto it = res.begin(); it != res.end(); ++it) {
		fprintf(stdout,"%u\t%u\n",ids[it->first],ids[it->second]);
	}
	
	double runtime = t2.tv_sec - t1.tv_sec + (t2.tv_nsec - t1.tv_nsec) / 1e9;
	fprintf(stderr,"%lu matches found, %f s total runtime\n",res.size(),runtime);
	
	return 0;
}

