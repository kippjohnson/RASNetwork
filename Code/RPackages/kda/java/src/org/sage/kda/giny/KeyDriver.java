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

import giny.model.GraphPerspective;
import giny.model.Node;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

public class KeyDriver {

	public final Node node;
	public final int no_hits;
	public final int signature_length;
	public final int subnet_size;
	public final double pv;
	public final double pv_bonferroni;
	public final Set<Node> downstreamGenes;
	public final GraphPerspective graph;
	// annotations
	public boolean islocal;
	public boolean isboosted;
	public double qvalue; // compute only when permutations are run
	
	KeyDriver(GraphPerspective graph, Node node, int no_hits, int signature_length, int subnet_size, double pv, Set<Node> downStreamGenes){
		this.graph = graph;
		this.node = node;
		this.no_hits = no_hits;
		this.signature_length = signature_length;
		this.subnet_size = subnet_size;
		this.pv = pv;
		this.pv_bonferroni = Math.min(pv * subnet_size, 1.0); // bonferroni corrected
		this.downstreamGenes = downStreamGenes;
	}
	
	static class PValueComparator implements Comparator<KeyDriver>{

		public int compare(KeyDriver o1, KeyDriver o2) {
			return Double.compare(o1.pv, o2.pv);
		}

	} 
	
	public static List<KeyDriver> pvalueFilter(List<KeyDriver> kd, double pvalue){
		ArrayList<KeyDriver> nkd = new ArrayList<KeyDriver>();
		for(KeyDriver k : kd){
			if(k.pv <= pvalue)
				nkd.add(k);
		}
		return nkd;
	}
}
