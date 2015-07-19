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

import fing.model.FingRootGraphFactory;
import giny.model.GraphPerspective;
import giny.model.Node;
import giny.model.RootGraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ImportExport {
	
	public static void output(Writer w, Map<String,List<KeyDriver>> results) throws IOException{
		DecimalFormat formatter = new DecimalFormat("0.###E0");
		// write header
		w.write("module\tkeydriver\tpvalue\tcorrecte-pvalue\tqvalue\thits\tdownstream\ttype\n");
		List<String> moduleNames = new ArrayList<String>(results.keySet());
		Collections.sort(moduleNames);
		for(String module : moduleNames){
			for(KeyDriver kd : results.get(module)){
				w.write(module);
				w.write("\t");
				w.write(kd.node.getIdentifier());
				w.write("\t");
				w.write(formatter.format(kd.pv));
				w.write("\t");
				w.write(formatter.format(kd.pv_bonferroni));
				w.write("\t");
				w.write(formatter.format(kd.qvalue));
				w.write("\t");
				w.write(Integer.toString(kd.no_hits));
				w.write("\t");
				w.write(Integer.toString(kd.downstreamGenes.size()));
				w.write("\t");
				w.write(kd.islocal ? "local" : "global");
				w.write("\n");
			}
		}
	}
	
	public static Map<String,Set<Node>> loadModules(String f, Network network)
		throws IOException 
	{
		HashMap<String,Set<Node>> modules = new HashMap<String,Set<Node>>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line;
		while((line = br.readLine()) != null){
			line = line.trim();
			if(line.equals(""))
				continue;
			String[] arr = line.split("\t");
			Set<Node> set = modules.get(arr[1]);
			if(set == null){
				set = new HashSet<Node>();
				modules.put(arr[1], set);
			}
			Node node = network.nodeMap.get(arr[0]);
			if(node == null){
				System.err.println("Warning: geneset node '" + arr[0] + "' not in network.");
			}else{
				set.add(node);
			}
		}
		return modules;
	}

	public static Network loadNetwork(String f, boolean directed) 
		throws IOException
	{
		Network network = new Network();
		RootGraph graph = FingRootGraphFactory.instantiateRootGraph();
		
		BufferedReader br =null;
		try{
			br = new BufferedReader(
				new InputStreamReader(ImportExport.class.getResourceAsStream(f)));
		}catch(NullPointerException e){
			br = new BufferedReader(
					new FileReader(f));
		}
		String line = null;
		
		HashMap<String,Node> nodeMap = new HashMap<String,Node>();
		
		while((line = br.readLine()) != null){
			line = line.trim();
			if(line.equals(""))
				continue;
			String[] arr = line.split("\t");
			Node n1 = nodeMap.get(arr[0]);
			Node n2 = nodeMap.get(arr[1]);
			if(n1 == null){
				n1 = graph.getNode(graph.createNode());
				n1.setIdentifier(arr[0]);
				nodeMap.put(arr[0], n1);
			}
			if(n2 == null){
				n2 = graph.getNode(graph.createNode());
				n2.setIdentifier(arr[1]);
				nodeMap.put(arr[1], n2);
			}
			graph.createEdge(n1.getRootGraphIndex(), n2.getRootGraphIndex());
		}
		network.nodeMap = nodeMap;
		network.graph = graph.createGraphPerspective(graph.getNodeIndicesArray(), graph.getEdgeIndicesArray());
		return network;
	}
	
	static class Network {
		HashMap<String,Node> nodeMap;
		GraphPerspective graph;
	}
}
