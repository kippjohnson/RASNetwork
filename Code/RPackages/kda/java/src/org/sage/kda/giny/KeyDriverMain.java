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

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;

/**
 * Main entry point for running an analysis, either from command line or as API
 * @author justin.guinney@sagebase.org
 *
 */
public class KeyDriverMain {
	
	public static void main(String[] argv) throws Exception {
		Options options = new Options();
		options.addOption("n", true, "network file [required]");
		options.addOption("m", true, "module file [required]");
		options.addOption("o", true, "output file [optional]");
		options.addOption("p", true, "p-value threshold (default=.05)");
		options.addOption("r",true,"# of permutations (default=100)");
		options.addOption("u",false,"undirected network flag");
		options.addOption("k",true,"knn (default=[6, directed; 2, undirected])");
		options.addOption("e",true,"expansion (default=2)");
		options.addOption("v",false,"verbose");
		options.addOption("h",false,"show help");
		
		CommandLine line = new GnuParser().parse(options, argv);
		
		if(line.hasOption("h")){
			new HelpFormatter().printHelp( "KeyDriverMain", options );
			System.exit(0);
		}
		
		// defaults
		double pthreshold = .05;
		int knn = 6;
		int permutations = 100;
		int expansion = 2;
		
		String networkFile = line.getOptionValue("n");
		
		String moduleFile = line.getOptionValue("m");
		String outputFile = line.getOptionValue("o");
		
		boolean verbose = line.hasOption("v");
		boolean directed = !line.hasOption("u");
		
		if(line.hasOption("r")){
			permutations = Integer.parseInt(line.getOptionValue("r"));
		}
		if(line.hasOption("p")){
			pthreshold = Double.parseDouble(line.getOptionValue("p"));
		}
		if(line.hasOption("k")){
			knn = Integer.parseInt(line.getOptionValue("k"));
		}
		if(line.hasOption("e")){
			expansion = Integer.parseInt(line.getOptionValue("e"));
		}
		
		
		ImportExport.Network network = ImportExport.loadNetwork(networkFile, directed);
		Map<String,Set<Node>> modules = ImportExport.loadModules(moduleFile, network);
		
		Map<String,List<KeyDriver>> results = new  HashMap<String,List<KeyDriver>>();
		List<String> moduleNames = new ArrayList<String>(modules.keySet());
		Collections.sort(moduleNames);
		for(String s : moduleNames){
			if(verbose){ System.out.println("Computing module " + s); }
			Set<Node> module = modules.get(s);
			KDAResult result =  run(network.graph, module, 
					expansion, knn, directed, permutations);
			results.put(s, KeyDriver.pvalueFilter(result.keyDrivers,pthreshold));
		}
		
		// print out results to file or standard out
		Writer out = null;
		if(outputFile == null){
			out = new PrintWriter(System.out);
		}else{
			out = new FileWriter(outputFile);
		}
		if(verbose) { System.out.println("Writing results."); }
		ImportExport.output(out, results);
		out.flush();
		if(outputFile != null)
			out.close();
	}
	
	
	
	public static class KDAResult {
		public List<KeyDriver> keyDrivers;
		public int knnUsed;
	}
	
	static final RandomData _random = new RandomDataImpl(new MersenneTwister());
	
	/**
	 * 
	 * @param network
	 * @param module
	 * @param expansion
	 * @param knn negative value for knn will 
	 * @param directed
	 * @return
	 */
	public static KDAResult run(GraphPerspective network, Set<Node> module, 
			int expansion, int knn, boolean directed){
		KDAResult result = new KDAResult();
		
		GraphPerspective expandNet = NetworkFncs.findNLayerNeighbors_GraphPerspective(
				network, module, expansion, false);
		
		boolean searchKnn = knn < 0;
		knn = Math.abs(knn);
		
		if(searchKnn){
			double bestpvalue = 1.0;
			int bestIdx = 0;
			// over all knns, select keydriver with smallest p-value
			ArrayList<List<KeyDriver>> knnRuns = new ArrayList<List<KeyDriver>>();
			int idx = 0;
			for(int k = 1; k <= knn; ++k, ++idx){
				List<KeyDriver> run = KeyDriverFncs.keyDriverInSubnetwork(expandNet,
						module, k, directed);
				for(KeyDriver d : run){
					if(d.pv < bestpvalue){
						bestIdx = idx;
						bestpvalue = d.pv;
					}
				}
				knnRuns.add(run);
			}

			result.keyDrivers = knnRuns.get(bestIdx);
			result.knnUsed = bestIdx+1;
		}else{
			result.keyDrivers = KeyDriverFncs.keyDriverInSubnetwork(expandNet, module, knn, directed);
			result.knnUsed = knn;
		}
		return result;
	}
	
	public static KDAResult run(GraphPerspective network, Set<Node> module, 
			int expansion, int knn, boolean directed, int permutations){
		KDAResult result = run(network, module, expansion, knn, directed);
		if(permutations > 0){
			double[] perm_pvals = runPermutations(network, module, expansion, result.knnUsed, 
					directed, permutations,false);
			double[] obs_pvals = new double[result.keyDrivers.size()];
			for(int i = 0; i < result.keyDrivers.size(); ++i){
				obs_pvals[i] = result.keyDrivers.get(i).pv;
			}
			double[] fdr = MathUtils.fdr(obs_pvals, perm_pvals);
			double[] qval = MathUtils.qvals(obs_pvals, fdr);
			for(int  i = 0; i < qval.length; ++i){
				result.keyDrivers.get(i).qvalue = qval[i];
			}
		}
		return result;
	}
	
	public static double[] runPermutations(GraphPerspective network, Set<Node> module, 
			int expansion, int knn, boolean directed, int permutations, boolean verbose)
	{
		ArrayList<Double> all = new ArrayList<Double>();
		for(int i = 1; i <= permutations; ++i){
			if(verbose && i % 10 == 0){
				System.out.print(".");
			}
			if(verbose &&  i % 100 == 0){
				System.out.println();
			}
			Set<Node> permutedModule = permuteModule(network, module.size());
			KDAResult result = run(network, permutedModule, expansion, knn, directed);
			for(KeyDriver kd : result.keyDrivers){
				all.add(kd.pv);
			}
		}
		double[] pvals = new double[all.size()];
		for(int i = 0; i < pvals.length; ++i)
			pvals[i] = all.get(i);
		return pvals;
	}
	
	private static Set<Node> permuteModule(GraphPerspective network, int moduleSize){
		HashSet<Node> newModule = new HashSet<Node>(moduleSize);
		int[] arrayIndices = network.getNodeIndicesArray();
		int[] selectedIdxs = _random.nextPermutation(arrayIndices.length, moduleSize);
		for(int i = 0; i < moduleSize; ++i){
			newModule.add(network.getNode(arrayIndices[selectedIdxs[i]]));
		}
		return newModule;
	}
}
