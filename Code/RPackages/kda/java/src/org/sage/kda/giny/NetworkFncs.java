/*
  Copyright (c) 2010, Sage Bionetworks

  This library is free software; you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published
  by the Free Software Foundation; either version 2.1 of the License, or
  any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  The software and
  documentation provided hereunder is on an "as is" basis, and the
  Institute for Systems Biology and the Whitehead Institute
  have no obligations to provide maintenance, support,
  updates, enhancements or modifications.  In no event shall the
  Institute for Systems Biology and the Whitehead Institute
  be liable to any party for direct, indirect, special,
  incidental or consequential damages, including lost profits, arising
  out of the use of this software and its documentation, even if the
  Institute for Systems Biology and the Whitehead Institute
  have been advised of the possibility of such damage.  See
  the GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this library; if not, write to the Free Software Foundation,
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/
package org.sage.kda.giny;

import static org.sage.kda.giny.SetUtils.diff;
import giny.filter.Filter;
import giny.model.Edge;
import giny.model.GraphPerspective;
import giny.model.Node;

import java.util.HashSet;
import java.util.Set;

public class NetworkFncs {
	
	public static GraphPerspective findNLayerNeighbors_GraphPerspective(final GraphPerspective network, final Set<Node> subnetNodes, int nlayers, boolean directed){
		final Set<Node> subnetwork = findNLayerNeighbors(network, subnetNodes, new HashSet<Node>(), nlayers, directed);
		Filter f = new Filter() {
			public boolean passesFilter(Object obj) {
				if(obj instanceof Node){
					return isGood((Node)obj);
				}else if(obj instanceof Edge){
					return isGood(((Edge)obj).getSource()) && isGood(((Edge)obj).getTarget());
				}
				throw new IllegalArgumentException("Unexpected obj in filter" + obj.getClass());
			}
			
			private boolean isGood(Node n){
				return subnetwork.contains(n) || subnetNodes.contains(n);
			}
		};
		return network.createGraphPerspective(f);
	}
	

	
	public static Set<Node> findNLayerNeighbors(GraphPerspective network, Set<Node> subnetNodes, Set<Node> all, int nlayers, boolean directed){
		
		all.addAll(subnetNodes);
		
		//base case
		if(nlayers < 0)
			throw new IllegalArgumentException("nlayers argument must be 1 or greater");
		else if(nlayers == 0)
			return all;
		
		Node[] subnetArr = subnetNodes.toArray(new Node[subnetNodes.size()]);
		HashSet<Node> neighbors = new HashSet<Node>();
		for(int i = 0; i < subnetArr.length; ++i){
			int[] adjacentEdges = network.getAdjacentEdgeIndicesArray(
					subnetArr[i].getRootGraphIndex(),
					!directed,
					!directed,
					true);
			for(int j = 0; j < adjacentEdges.length; ++j){
				neighbors.add(network.getNode(network.getEdgeTargetIndex(adjacentEdges[j])));
				if(!directed)
					neighbors.add(network.getNode(network.getEdgeSourceIndex(adjacentEdges[j])));
			}
		}
		
		Set<Node> diff = diff(all, neighbors);
		
		if(neighbors.size() == 0 || diff.size() == 0)
			return all;
		else
			return findNLayerNeighbors(network, diff, all, nlayers-1, directed);
		
	}
	
	public static Set<Node> downStreamGenes(GraphPerspective network, Set<Node> seedNodes, int nLayers, boolean directed){
		Set<Node> prenodes = seedNodes;
		Set<Node> curnodes = new HashSet<Node>(seedNodes);
		int cnt = nLayers;
		
		while(cnt > 0){
			Set<Node> neighbors = findNLayerNeighbors(network, prenodes, new HashSet<Node>(), 1, directed);
		
			curnodes.addAll(neighbors);
			
			Set<Node> pcdiff = diff(curnodes, prenodes);
			
			if(pcdiff.size() == 0){
				break;
			}
			
			prenodes = new HashSet<Node>(pcdiff);
			
			cnt--;
		}
		
		return curnodes;
	}
}
