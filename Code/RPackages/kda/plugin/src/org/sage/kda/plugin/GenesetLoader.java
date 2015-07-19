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
package org.sage.kda.plugin;

import giny.model.Node;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import cytoscape.CyNetwork;
import cytoscape.Cytoscape;

public class GenesetLoader {
	
	static public HashMap<String, Set<Node>> loadGenesetModules(CyNetwork network, String text) throws IOException {
		
		String[] lines = text.split("\n");
		for(int i = 0; i < lines.length; ++i){
			lines[i] = lines[i].trim();
		}
		return onLines(network, lines);
	}
	
	static public HashMap<String, Set<Node>> loadGenesetModules(CyNetwork network, File file) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		ArrayList<String> lines = new ArrayList<String>();
		String line;
		while((line = br.readLine()) != null){
			lines.add(line);
		}
		return onLines(network, lines.toArray(new String[lines.size()]));
	}
	
	static public HashMap<String, Set<Node>> onLines(CyNetwork network, String[] lines){
		HashMap<String,Set<Node>> modules = new HashMap<String,Set<Node>>();
		
		boolean multipleModules = false;
		for(int i = 0 ; i < lines.length; ++i){
		
			String[] arr = lines[i].trim().split("\t");
			if(i == 0 && arr.length > 1)
				multipleModules = true;

			Node n = Cytoscape.getCyNode(arr[0], false);
			if(n != null && network.containsNode(n)){
				String mname = multipleModules ? arr[1] : "module";
				Set<Node> gs = modules.get(mname);
				if(gs == null){
					gs = new HashSet<Node>();
					modules.put(mname, gs);
				}
				gs.add(n);
			}
		}
		return modules;
	}
}
	
	
	
	
