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

import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.sage.kda.giny.KeyDriver;

import cytoscape.CyNetwork;

public class KDARun {

	HashMap<String,KDAState> _network2State = new HashMap<String, KDAState>();
	HashMap<String,KDAState> _module2State = new HashMap<String,KDAState>();
	
	public void addKDAState(KDAState state){
		_network2State.put(state.network.getIdentifier(), state);
		_module2State.put(state.moduleName, state);
	}
	
	public void removeKDAState(KDAState state){
		_network2State.remove(state.network.getIdentifier());
		_module2State.remove(state.moduleName);
	}
	
	public KDAState getStateByNetworkId(String id){
		return _network2State.get(id);
	}
	public KDAState getStateByModuleName(String name){
		return _module2State.get(name);
	}
	
	static class KDAState {
		public KDAState(String moduleName, CyNetwork network, List<KeyDriver> drivers, Set<Node> moduleNodes){
			this.moduleName = moduleName;
			this.network = network;
			this.keyDrivers = drivers;
			this.moduleNodes = moduleNodes;
		}
		final CyNetwork network;
		final String moduleName;
		final List<KeyDriver> keyDrivers;
		final Set<Node> moduleNodes;
	}
	
	

}
