/*
 * sccs.cpp -- find connected components in an undirected graph using
 * 	depth-first search
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
#include "graph.h"

int main(int argc, char** argv) {
	bool read_part = false; /* read already partitioned graph */
	for(int i=1;i<argc;i++) if(argv[i][0] == '-') switch(argv[i][1]) {
		case 'p':
			read_part = true;
			break;
		default:
			fprintf(stderr,"Unknown parameter: %s!\n",argv[i]);
			break;
	}
	
	graph g;
	std::vector<unsigned int> ids;
	std::vector<unsigned int> sccids;
	g.read_graph(stdin,read_part,&ids);
	g.make_symmetric();
	
	unsigned int nsccs = g.find_sccs(sccids);
	fprintf(stderr,"%u connected components found\n",nsccs);
	if(sccids.size() != ids.size()) fprintf(stderr,"Inconsistent results!\n");
	else for(size_t i=0;i<ids.size();i++) fprintf(stdout,"%u\t%u\n",ids[i],sccids[i]);
	return 0;
}


