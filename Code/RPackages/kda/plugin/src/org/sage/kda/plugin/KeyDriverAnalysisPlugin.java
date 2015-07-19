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

import giny.model.GraphPerspective;
import giny.model.Node;

import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JDialog;
import javax.swing.JMenu;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import org.sage.kda.giny.KeyDriver;
import org.sage.kda.giny.KeyDriverMain;
import org.sage.kda.giny.KeyDriverMain.KDAResult;

import cytoscape.CyNetwork;
import cytoscape.Cytoscape;
import cytoscape.layout.CyLayoutAlgorithm;
import cytoscape.layout.CyLayouts;
import cytoscape.plugin.CytoscapePlugin;
import cytoscape.task.Task;
import cytoscape.task.TaskMonitor;
import cytoscape.task.ui.JTaskConfig;
import cytoscape.task.util.TaskManager;
import cytoscape.util.CytoscapeAction;
import cytoscape.view.CyNetworkView;
import cytoscape.visual.VisualStyle;

public class KeyDriverAnalysisPlugin extends CytoscapePlugin implements PropertyChangeListener, KDACallback{
		
	static final double THRESHOLD = .05; // todo, make input property
	
	private Settings _driverSettings = new Settings();
	private KDARun _currentRun;
	private DriverPanel _driverDialog;
	
	//public String 
	
	public KeyDriverAnalysisPlugin() {
		 KDAPluginAction kdaAction = new KDAPluginAction();
        
        //and add it to the menus
        //Cytoscape.getDesktop().getCyMenus().addAction(kdaAction);
        JMenu pluginMenu = Cytoscape.getDesktop().getCyMenus().getMenuBar().getMenu("Plugins");
        pluginMenu.add(kdaAction);
        Cytoscape.getDesktop().getSwingPropertyChangeSupport().addPropertyChangeListener(this);
        //        JMenu sageMenu = new JMenu("Key Driver Analysis");
//        pluginMenu.add(sageMenu);
//        sageMenu.add(kdaAction);
    }
	
	/**
	 * Listens for focus 
	 */
	public void propertyChange(PropertyChangeEvent e){
		if(_currentRun != null){
			if(e.getPropertyName().equals(Cytoscape.NETWORK_DESTROYED)){
				String networkId = (String)e.getNewValue();
				KDARun.KDAState state = _currentRun.getStateByNetworkId(networkId);
				KDAViz.cleanupRun(state);
				_currentRun.removeKDAState(state);
			}
		}
	}
	
	
	public void onKDA(final Settings settings){
		removeCurrent();
		Task task = new Task() {
			boolean halt = false;
			TaskMonitor taskMonitor;
			
			public void setTaskMonitor(TaskMonitor tm)
					throws IllegalThreadStateException {
				taskMonitor = tm;
			}
			
			public void run() {
				try{
					KDARun run = new KDARun();
					CyNetwork network = settings.network;
					
					Map<String,Set<Node>> modules = GenesetLoader.loadGenesetModules(
							network, settings.geneModuleText);
		
					final HashMap<String, List<KeyDriver>> allResults = new HashMap<String,List<KeyDriver>>();
					int noModules = modules.size();
					int count = 0;
					for(String moduleName : modules.keySet()){
						if(halt) break;
						
						taskMonitor.setStatus("Calculating module '" + moduleName + "'");
						CyNetwork old = Cytoscape.getNetwork(moduleName);
						if(old != null){	
							Cytoscape.destroyNetworkView(old);
							Cytoscape.destroyNetwork(old);
						}
						
						int knn = settings.knn;
						if(settings.knnAuto)
							knn = -1 * knn;
						
						KDAResult r = KeyDriverMain.run(network, 
								modules.get(moduleName), 
								settings.expansion, knn, 
								settings.directed, 
								settings.permutations);
						
						List<KeyDriver> sigKDs = KeyDriver.pvalueFilter(r.keyDrivers, settings.pval);
						if(sigKDs.size() > 0){
							GraphPerspective kdnetwork = sigKDs.get(0).graph;
							
							CyNetwork modNetwork = Cytoscape.createNetwork(
									kdnetwork.getNodeIndicesArray(),
									kdnetwork.getEdgeIndicesArray(), 
									moduleName, network);
							final KDARun.KDAState state = new KDARun.KDAState(
									moduleName, modNetwork, sigKDs, modules.get(moduleName));
							run.addKDAState(state);
							allResults.put(moduleName, sigKDs);
							try{
								SwingUtilities.invokeAndWait(new Runnable(){
									public void run(){ createNetworkView(state); }
								});
							}catch(Exception e){
								e.printStackTrace();
							}
						}
						taskMonitor.setPercentCompleted((int)((double)++count / noModules) * 100);
					}
					_currentRun = run;
					try{
						SwingUtilities.invokeAndWait(new Runnable(){
							public void run(){
								if(allResults.values().size() == 0){
									JOptionPane.showMessageDialog(null, "No key drivers found. Try changing the significance threshold.");
								}else{
									_driverDialog.postResults(allResults);
								}
							}
						});
					}catch(Exception e){
						e.printStackTrace();
					}
					
				}catch(IOException e){
					// TODO
					e.printStackTrace();
					taskMonitor.setException(e, e.getMessage());
				}
			}
			
			public void halt() {
				halt = true;
			}
			
			public String getTitle() {
				return "Key driver running...";
			}
		};
		
		JTaskConfig jTaskConfig = new JTaskConfig();
		jTaskConfig.setOwner(Cytoscape.getDesktop());
		jTaskConfig.displayCloseButton(true);

		jTaskConfig.displayCancelButton(true);

		jTaskConfig.displayStatus(true);
		jTaskConfig.setAutoDispose(true);

		// Execute Task in New Thread; pops open JTask Dialog Box.
		TaskManager.executeTask(task, jTaskConfig);
	}
	
	public void keyDriverSelect(String moduleName, KeyDriver driver) {
		//CyNetworkView view = createNetworkView(network);
		KDARun.KDAState state = null;
		if(_currentRun != null && (state = _currentRun.getStateByModuleName(moduleName)) != null){
			Cytoscape.setCurrentNetworkView(state.network.getIdentifier());
			Cytoscape.getDesktop().setFocus(state.network.getIdentifier());
			state.network.unselectAllNodes();
			state.network.setSelectedNodeState(driver.node, true);
		}
	}
	
	void removeCurrent(){
		if(_currentRun != null){
			KDAViz.cleanupRun(_currentRun);
			for(String networkId : _currentRun._network2State.keySet().toArray(new String[0])){
				CyNetwork network = Cytoscape.getNetwork(networkId);
				if(network != null){
					Cytoscape.destroyNetworkView(network);
					Cytoscape.destroyNetwork(network);
				}
			}
			_currentRun = null;
		}
	}
		

	CyNetworkView createNetworkView(KDARun.KDAState state){
		// = Cytoscape.getNetworkView(state.network.getIdentifier()); 
		VisualStyle style = KDAViz.setupAttributeVizForRun(state);
		
		CyNetworkView view= Cytoscape.createNetworkView(state.network);
		view.setVisualStyle(style.getName());
		view.applyVizmapper(style);
        CyLayoutAlgorithm layout = CyLayouts.getLayout("circular");
        if(layout == null)
        	layout = CyLayouts.getDefaultLayout();
        layout.doLayout(view);
	
		return view;
	}
	
	class KDAPluginAction extends CytoscapeAction   {
		
		public KDAPluginAction(){
			super("Key Driver Analysis");
		}
		
		public void actionPerformed(ActionEvent evt)  {
			if(_driverDialog == null){
				_driverDialog = new DriverPanel(_driverSettings, KeyDriverAnalysisPlugin.this);
				_driverDialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
				_driverDialog.setSize(500,400);
				_driverDialog.addWindowListener(new WindowAdapter() {
					public void windowClosed(WindowEvent evt){
						_driverDialog = null;
					}
				});
			}
            _driverDialog.setVisible(true);   
		}
	}

	
	
}