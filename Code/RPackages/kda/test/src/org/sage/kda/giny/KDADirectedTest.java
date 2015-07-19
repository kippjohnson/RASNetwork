package org.sage.kda.giny;

import giny.model.GraphPerspective;
import giny.model.Node;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class KDADirectedTest {

	private ImportExport.Network _network;
	
	@Before
	public void setup() throws Exception {
		_network = ImportExport.loadNetwork("./data/directedNetwork.txt",true);
	}
	
	@Test
	public void testNeighbors(){
		String[] nodes = {"a","b","c","d","e","f","g","h","i"};
		
		for(int i = 0; i < nodes.length; ++i){
			HashSet<Node> module = new HashSet<Node>();
			module.add(_network.nodeMap.get(nodes[i]));
			for(int j = 0; j < nodes.length-1; ++j){
				Set<Node> neighbors = NetworkFncs.findNLayerNeighbors(_network.graph, module, new HashSet<Node>(), j+1, true);
				//System.out.println(nodes[i] + " " + i + " " + (j+1)  + " " + neighbors.size());
				assertEquals(neighbors.size(), Math.min(j+2, nodes.length -i));
			}
		}
	}
	

	@Test
	public void testNeighbors_Graph(){
		String[] nodes = {"a","b","c","d","e","f","g","h","i"};
		
		for(int i = 0; i < nodes.length; ++i){
			HashSet<Node> module = new HashSet<Node>();
			module.add(_network.nodeMap.get(nodes[i]));
			for(int j = 0; j < nodes.length-1; ++j){
				GraphPerspective graph = NetworkFncs.findNLayerNeighbors_GraphPerspective(
						_network.graph, module, j+1, true);
				//System.out.println(nodes[i] + " " + i + " " + (j+1)  + " " + neighbors.size());
				
				assertEquals(graph.getNodeIndicesArray().length, Math.min(j+2, nodes.length - i));
			}
		}
	}
	
	@Test
	public void testDownstream(){
		String[] nodes = {"a","b","c","d","e","f","g","h","i"};
		
		for(int i = 0; i < nodes.length; ++i){
			HashSet<Node> module = new HashSet<Node>();
			module.add(_network.nodeMap.get(nodes[i]));
			for(int j = 0; j < nodes.length-1; ++j){
				Set<Node> neighbors = NetworkFncs.downStreamGenes(_network.graph, module, j+1, true);
				//System.out.println(nodes[i] + " " + i + " " + (j+1)  + " " + neighbors.size());
				assertEquals(neighbors.size(), Math.min(j+2, nodes.length - i));
			}
		}
	}
	
//	@Test
//	public void testKDA(){
//		String[] nodes = {"b","c","d"};
//
//		HashSet<Node> module = new HashSet<Node>();
//		for(int i = 0; i < nodes.length; ++i){
//			module.add(_network.nodeMap.get(nodes[i]));
//		}
//			
//		
//		List<KeyDriver> kds = KeyDriverAnalysis.kda(_network.graph, module, 1.0, 5, 1, true);
//		for(KeyDriver kd : kds){
//			System.out.println(kd.node.getIdentifier() + "\t" + kd.no_hits + "\t" + kd.downstreamGenes.size() + "\t" + kd.pv);
//		}
//	}
}
