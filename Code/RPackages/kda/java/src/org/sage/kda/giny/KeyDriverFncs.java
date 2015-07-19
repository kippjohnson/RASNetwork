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

import static org.sage.kda.giny.NetworkFncs.downStreamGenes;
import static org.sage.kda.giny.SetUtils.diff;
import static org.sage.kda.giny.SetUtils.intersect;
import giny.model.GraphPerspective;
import giny.model.Node;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class KeyDriverFncs {
	

	public static List<KeyDriver> keyDriverInSubnetwork(GraphPerspective network, 
			Set<Node> signature, int nLayers, boolean directed)
	{
		
		HashSet<Node> allNodes = new HashSet<Node>();
		for(java.util.Iterator<Node> iter = network.nodesIterator(); iter.hasNext(); )
			allNodes.add(iter.next());
		
		int[] nodeIdxs = network.getNodeIndicesArray();
		int noSubnetsize = nodeIdxs.length;
		
		boolean network_as_signature = network.getNodeCount() == signature.size() &&
			signature.containsAll(allNodes);
		
	
		List<Set<Node>> downStreamGenes = new ArrayList<Set<Node>>(noSubnetsize);
		
		for(int i = 0; i < nodeIdxs.length; ++i){
			Node n = network.getNode(nodeIdxs[i]);
			HashSet<Node> seed = new HashSet<Node>(1);
			seed.add(n);
			Set<Node> idn = downStreamGenes(network, seed, nLayers, directed);
			downStreamGenes.add(diff(idn,seed));
		}
		
//		int mymincut = min_downstreamnodes;
//		if(mymincut <= 0){
//			int[] dssizes = new int[downStreamGenes.size()];
//			for(int i = 0; i < dssizes.length; ++i)
//				dssizes[i] = downStreamGenes.get(i).size();
//			
//			mymincut = (int)(MathUtils.mean(dssizes) + MathUtils.sd(dssizes));
//		}
		
		ArrayList<KeyDriver> results = new ArrayList<KeyDriver>();
		for(int i = 0; i < noSubnetsize; ++i){
			
				Set<Node> hits = intersect(downStreamGenes.get(i), signature);
//				double pv = new HypergeometricDistributionImpl(
//						noSubnetsize,
//						downStreamGenes.get(i).size(), 
//						signature.size()).upperCumulativeProbability(hits.size());
				
				double pv = MathUtils.hyperUpperCumulative(hits.size(), signature.size(), downStreamGenes.get(i).size(), noSubnetsize);
				KeyDriver r = new KeyDriver(
						network,
						network.getNode(nodeIdxs[i]), 
						hits.size(), 
						signature.size(), 
						noSubnetsize, 
						pv, 
						downStreamGenes.get(i));
				//if(r.pv < pvThreshold){
					results.add(r);
				//}
			//}
		}
		

		localAnalysis(results, network);
		boostHubs(results, network, directed);
		
		return results;
	}
	
	static void localAnalysis(List<KeyDriver> keyDrivers, GraphPerspective network){
		
		int[] idns = new int[keyDrivers.size()];
		for(int i = 0; i < idns.length; ++i){
			idns[i] = keyDrivers.get(i).downstreamGenes.size();
		}
		
		for(int i= 0; i < keyDrivers.size(); ++i){
			ArrayList<Set<Node>> check = new ArrayList<Set<Node>>();
			KeyDriver r = keyDrivers.get(i);
			for(int j = 0; j < idns.length; ++j){
				if(idns[j] > idns[i]){
					check.add(keyDrivers.get(j).downstreamGenes);
				}
			}
			r.islocal = (check.size() > 0)
				? !SetUtils.isSetInSets(r.downstreamGenes, check)
				: true;
		}
	}
	
	static void boostHubs(List<KeyDriver> keyDrivers, GraphPerspective network, boolean directed){
		
		int[] nodes = network.getNodeIndicesArray();
		int[] degrees = new int[nodes.length];
		for(int i= 0; i < nodes.length; ++i){
			degrees[i] = directed
					? network.getOutDegree(nodes[i])
					: network.getDegree(nodes[i]);
		}
		
		HashSet<Node> sel = new HashSet<Node>();
		double cutoff = MathUtils.mean(degrees) + 2 * MathUtils.sd(degrees);
		
		for(int i = 0; i < nodes.length; ++i){
			int v = directed ? network.getOutDegree(nodes[i]) : network.getDegree(nodes[i]);
			if(v > cutoff)
				sel.add(network.getNode(nodes[i]));
		}
		
		for(KeyDriver kd : keyDrivers)
			kd.isboosted = sel.contains(kd.node); 
		
	}
	
	static String printNodes(Collection c){
		StringBuffer sb = new StringBuffer();
		for(Object o : c){
			if(o instanceof Node){
				sb.append(((Node)o).getIdentifier());
			}else{
				sb.append(o.toString());
			}
			sb.append(",");
		}
		return sb.toString();
	}
	
}
