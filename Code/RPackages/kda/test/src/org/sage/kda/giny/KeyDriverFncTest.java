package org.sage.kda.giny;


import giny.model.Node;
import giny.model.RootGraph;

import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class KeyDriverFncTest {

	@Test
	public void downStreamGenesTestDirected_1() throws Exception {
		ImportExport.Network n = ImportExport.loadNetwork("network2.txt", true);
		Set<Node> dsg = NetworkFncs.downStreamGenes(n.graph, makeSet(n, new String[] {"a"}), 1, true);
		assertTrue(dsg.containsAll(makeSet(n, new String[] {"b","c" })));
	}
	
	@Test
	public void downStreamGenesTestDirected_2() throws Exception {
		ImportExport.Network n = ImportExport.loadNetwork("network2.txt", true);
		Set<Node> dsg = NetworkFncs.downStreamGenes(n.graph, makeSet(n, new String[] {"a"}), 2, true);
		assertTrue(dsg.containsAll(makeSet(n, new String[] {"b","c","d" })));
	}
	
	@Test
	public void downStreamGenesTestUndirected_1() throws Exception {
		ImportExport.Network n  = ImportExport.loadNetwork("network2.txt", false);
		Set<Node> dsg = NetworkFncs.downStreamGenes(n.graph, makeSet(n, new String[] {"a"}), 1, false);
		assertTrue(dsg.containsAll(makeSet(n, new String[] {"b","c","e" })));
	}
	
	@Test
	public void downStreamGenesTestUndirected_2() throws Exception {
		ImportExport.Network n = ImportExport.loadNetwork("network2.txt", false);
		Set<Node> dsg = NetworkFncs.downStreamGenes(n.graph, makeSet(n, new String[] {"a"}), 2, false);
		assertTrue(dsg.containsAll(makeSet(n, new String[] {"b","c","e","d","f" })));
	}
	
	@Test
	public void nLayerNeighborTestDirected() throws Exception {
		ImportExport.Network n = ImportExport.loadNetwork("network2.txt", true);
		Set<Node> genes = NetworkFncs.findNLayerNeighbors(
				n.graph, makeSet(n, new String[] {"a"}), new HashSet<Node>(), 1, true);
		assertEquals(3, genes.size()); 
		genes = NetworkFncs.findNLayerNeighbors(n.graph, makeSet(n, new String[] {"a"}),new HashSet<Node>(), 2, true);
		assertEquals(4, genes.size());
	}
	
	@Test
	public void nLayerNeighborTestDirected_multi() throws Exception {
		ImportExport.Network n = ImportExport.loadNetwork("network2.txt", true);
		Set<Node> genes = NetworkFncs.findNLayerNeighbors(n.graph, 
				makeSet(n, new String[] {"a","b"}),new HashSet<Node>(), 1,true);
		assertEquals(4, genes.size());
		
		genes = NetworkFncs.findNLayerNeighbors(n.graph, 
				makeSet(n, new String[] {"a","b"}),new HashSet<Node>(), 2, true);
		assertEquals(5, genes.size());
	}
	
	@Test
	public void nLayerNeighborTestUnDirected1() throws Exception {
		ImportExport.Network n  = ImportExport.loadNetwork("network2.txt", false);
		Set<Node> genes = NetworkFncs.findNLayerNeighbors(n.graph, 
				makeSet(n, new String[] {"a"}),new HashSet<Node>(), 1, false);
		assertEquals(4, genes.size()); 	
	}
	
	@Test
	public void nLayerNeighborTestUnDirected2() throws Exception {
		ImportExport.Network n = ImportExport.loadNetwork("network2.txt", false);
		Set<Node> genes = NetworkFncs.findNLayerNeighbors(n.graph, 
				makeSet(n, new String[] {"a"}),new HashSet<Node>(), 2, false);
		assertEquals(6, genes.size());
	}
		
	private Set<Node> makeSet(ImportExport.Network n, String[] arg){
		
		HashSet<Node> sig = new HashSet<Node>();
		for(String s : arg){
			sig.add(n.nodeMap.get(s));
		}
		return sig;
	}

}
