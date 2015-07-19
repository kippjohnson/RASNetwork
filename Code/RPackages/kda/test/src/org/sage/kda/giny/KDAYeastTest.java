package org.sage.kda.giny;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;
import org.sage.kda.giny.ImportExport.Network;

import giny.model.Node;

public class KDAYeastTest {
	
	private Network _network;
	private HashMap<String, Set<Node>> _modules = new HashMap<String,Set<Node>>();
	private String module = "HST02: Chr 2 560000";
	//String module = "HST08: Chr 12 680000";
	//private String module = "Chr 14 503000"
	
	
	@Before
	public void setup() throws Exception {
		_network = ImportExport.loadNetwork("./data/yeast_BN-full.pair",true);
		BufferedReader br = new BufferedReader(new FileReader("./data/yeast_eQTL-hotspots13.txt"));
		String line;
		while((line = br.readLine()) != null){
			line = line.trim();
			if(line.equals(""))
				continue;
			String[] arr = line.split("\t");
			Set<Node> set = _modules.get(arr[1]);
			if(set == null){
				set = new HashSet<Node>();
				_modules.put(arr[1], set);
			}
			Node node = _network.nodeMap.get(arr[0]);
			if(node == null){
				System.err.println("Warning: geneset node '" + arr[0] + "' not in network.");
			}else{
				set.add(node);
			}
		}
	}
	
	@Test
	public void _testPermutations(){
		testPermutations();
	}
	
	public double[] testPermutations(){
		return KeyDriverMain.runPermutations(_network.graph, _modules.get(module), 
				1, 6, true, 1000,true);
	}
	
	@Test
	public void _testKDA(){
		testKDA();
	}
	
	public double[] testKDA(){
		
		KeyDriverMain.KDAResult r = KeyDriverMain.run(_network.graph, _modules.get(module), 
				1, 6, true);
		double[] pvals = new double[r.keyDrivers.size()];
		for(int i = 0; i < pvals.length; ++i){
			pvals[i] = r.keyDrivers.get(i).pv;
		}
		return pvals;
	}
}
