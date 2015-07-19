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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.DefaultListSelectionModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.ListSelectionModel;
import javax.swing.SpringLayout;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.DefaultTableCellRenderer;

import org.sage.kda.giny.ImportExport;
import org.sage.kda.giny.KeyDriver;

import cytoscape.CyNetwork;
import cytoscape.Cytoscape;


public class DriverPanel extends JDialog {

    /*--------------------------------------------------------------
    FIELDS.
    --------------------------------------------------------------*/

    private JTextArea 	_geneSetTextArea = new JTextArea();
    private JComboBox _networkCombo = new JComboBox();
    private JTable 		_resultsTbl = new JTable();
    private JFormattedTextField	_pvTextField = new JFormattedTextField(new DecimalFormat("0.000E0"));
    private JRadioButton _directedGraphRadio = new JRadioButton("Directed");
    private JRadioButton _undirectedGraphRadio = new JRadioButton("Undirected");
    private JComboBox _expansionCombo =new JComboBox(new Integer[] { 1,2,3});
    private JComboBox _knnCombo = new JComboBox(new Integer[] { 1,2,3,4,5,6,7,8,9,10});
    private JButton _exportButton = new JButton("Export Results");
    private JRadioButton _knnSearchRadio = new JRadioButton("Search");
    private JRadioButton _knnFixedRadio = new JRadioButton("Fixed");
    private JComboBox _permutationsCombo = new JComboBox(new Integer[]{ 0, 100, 500, 1000, 5000, 10000 });
     
    private Map<String,List<KeyDriver>> _results;
    
    private static File _currentDir;
    
    private Settings _settings;
    private KDACallback _kdaCallback;
    
    
    /**
     * This constructor creates the panel with its swing-components.
     */
    public DriverPanel(Settings settings, KDACallback kdaCallback) {
   
    	_settings = settings;
    	_kdaCallback = kdaCallback;
    	
    	JSplitPane main = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
    	getContentPane().add(main);
    	
    	//main.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
    	JPanel runPanel = new JPanel();
    	runPanel.setLayout(new BoxLayout(runPanel, BoxLayout.PAGE_AXIS));
    	
        //create border.
        runPanel.setBorder(BorderFactory.createTitledBorder
                (BorderFactory.createLineBorder(Color.black),
                        "Key Driver Analysis Settings",
                        0,
                        0,
                        new Font("Sage", Font.BOLD, 16),
                        Color.black));
        
        JButton loadModules = new JButton("Import gene modules");
        loadModules.addActionListener(new LoadModulesActionListener());
        loadModules.setAlignmentX(LEFT_ALIGNMENT);
        runPanel.add(loadModules);
        runPanel.add(Box.createRigidArea(new Dimension(0,5)));
        
        JScrollPane listScroller = new JScrollPane(_geneSetTextArea);
        listScroller.setPreferredSize(new Dimension(250, 300));
        listScroller.setAlignmentX(LEFT_ALIGNMENT);
        runPanel.add(listScroller);
        
        JPanel springPanel = new JPanel(new SpringLayout());
        
        springPanel.add(new JLabel("Network"));
        _networkCombo = new JComboBox(getNetworkIds());
        springPanel.add(_networkCombo);
        
        ButtonGroup group = new ButtonGroup();
        ActionListener radioListener = new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				int knn = e.getSource() == _directedGraphRadio
					? 6 : 2;
				_knnCombo.setSelectedItem(knn);
			}
        };
        group.add(_directedGraphRadio);
        group.add(_undirectedGraphRadio);
        _directedGraphRadio.addActionListener(radioListener);
        _undirectedGraphRadio.addActionListener(radioListener);
        JPanel graphTypePane = new JPanel();
        graphTypePane.setLayout(new BoxLayout(graphTypePane, BoxLayout.LINE_AXIS));
        graphTypePane.add(_directedGraphRadio);
        graphTypePane.add(_undirectedGraphRadio);
        graphTypePane.add(Box.createHorizontalGlue());
        graphTypePane.setAlignmentX(LEFT_ALIGNMENT);
        
        springPanel.add(new JLabel("Graph type:"));
        springPanel.add(graphTypePane);
        
        springPanel.add(new JLabel("Significance threshold"));
        springPanel.add(_pvTextField);
        _pvTextField.setColumns(5);

        
        _expansionCombo.setEditable(false);
        springPanel.add(new JLabel("Module network expansion"));
        springPanel.add(_expansionCombo);

        _knnCombo.setEditable(false);
        springPanel.add(new JLabel("Nearest neighbors"));
        JPanel knnPane = new JPanel();
        ButtonGroup knnGroup = new ButtonGroup();
        knnGroup.add(_knnFixedRadio);
        knnGroup.add(_knnSearchRadio);
        knnPane.setLayout(new BoxLayout(knnPane, BoxLayout.LINE_AXIS));
        knnPane.setAlignmentX(LEFT_ALIGNMENT);
        knnPane.add(_knnSearchRadio);
        knnPane.add(_knnFixedRadio);
        knnPane.add(_knnCombo);
        knnPane.add(Box.createHorizontalGlue());
        springPanel.add(knnPane);
        
        springPanel.add(new JLabel("Permutations"));
        springPanel.add(_permutationsCombo);
        
        springPanel.setAlignmentX(LEFT_ALIGNMENT);
        
        runPanel.add(springPanel);
        JPanel runButtonPanel = new JPanel();
        runButtonPanel.setLayout(new BoxLayout(runButtonPanel, BoxLayout.LINE_AXIS));
        JButton runButton = new JButton("Run");
        runButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) { onRun(); }
		});
        runButtonPanel.add(Box.createHorizontalGlue());
        runButtonPanel.add(runButton);
        runButtonPanel.setAlignmentX(LEFT_ALIGNMENT);
        runPanel.add(runButtonPanel);
        
        
        SpringUtilities.makeCompactGrid(springPanel,
                6, 2, //rows, cols
                6, 6,        //initX, initY
                6, 6);       //xPad, yPad
        
        JPanel resultsPanel = new JPanel(new BorderLayout());
        resultsPanel.setBorder(BorderFactory.createTitledBorder
                (BorderFactory.createLineBorder(Color.black),
                        "Key Driver Analysis Results",
                        0,
                        0,
                        new Font("Sage", Font.BOLD, 16),
                        Color.black));

        resultsPanel.add(new JScrollPane(_resultsTbl), BorderLayout.CENTER);
        
        
        _exportButton.setEnabled(false);
        _exportButton.addActionListener(new ActionListener(){
        	public void actionPerformed(ActionEvent e) {onExport(); }});
        
        JPanel exportPanel = new JPanel();
        exportPanel.setLayout(new BoxLayout(exportPanel, BoxLayout.LINE_AXIS));
        exportPanel.add(Box.createHorizontalGlue());
        exportPanel.add(_exportButton);
        resultsPanel.add(exportPanel, BorderLayout.SOUTH);
        
        main.add(runPanel);
        main.add(resultsPanel);
        main.setDividerLocation(.75);
        
        setup(settings);
    }
    
    private void onExport(){
    	if(_results == null)
    		return;
    	
    	JFileChooser fc = new JFileChooser(_currentDir);
    	int returnVal = fc.showSaveDialog(this);
    	if(returnVal == JFileChooser.APPROVE_OPTION){
    		try{    
	    		File f = fc.getSelectedFile();
	    		FileWriter fw = new FileWriter(f);
	    		ImportExport.output(fw, _results);
	    		fw.close();
    		}catch(IOException e){
    			e.printStackTrace();
    			// TODO 
    		}
    	}
    }
    
    private String[] getNetworkIds(){
    	ArrayList<String> networkIds = new ArrayList<String>();
    	for(CyNetwork network : Cytoscape.getNetworkSet())
    		networkIds.add(network.getTitle());
    	Collections.sort(networkIds);
    	return networkIds.toArray(new String[networkIds.size()]);
    }
    
    class LoadModulesActionListener implements ActionListener {

		public void actionPerformed(ActionEvent e) {
			try{
				JFileChooser fileChooser = new JFileChooser(_currentDir);
				if(JFileChooser.APPROVE_OPTION == fileChooser.showDialog(DriverPanel.this, "Select gene modules(s)")){
					File f = fileChooser.getSelectedFile();
					_currentDir = f;
					StringBuffer sb = new StringBuffer();
					char[] buf = new char[1024];
					BufferedReader br = new BufferedReader(new FileReader(f));
					int nchars;
					while((nchars = br.read(buf)) != -1){
						sb.append(buf, 0, nchars);
					}
					_geneSetTextArea.setText(sb.toString());
					
				}
			}catch(IOException ex){
				ex.printStackTrace(); // TODO
			}
		}
    }
    
    void setup(Settings settings){
    	
    	CyNetwork selectedNetwork = settings.network != null
    		? settings.network
    		: Cytoscape.getCurrentNetwork();
    	if(selectedNetwork != null)
	    	_networkCombo.setSelectedItem(selectedNetwork.getIdentifier());
    	
    	_pvTextField.setValue(new Double(settings.pval));
    	if(settings.directed)
    		_directedGraphRadio.setSelected(true);
    	else
    		_undirectedGraphRadio.setSelected(true);
    	_expansionCombo.setSelectedItem(_settings.expansion);
    	_knnCombo.setSelectedItem(_settings.knn);
    	if(settings.knnAuto)
    		_knnSearchRadio.setSelected(true);
    	else
    		_knnFixedRadio.setSelected(true);
        
    	_geneSetTextArea.setText(settings.geneModuleText);
    	_permutationsCombo.setSelectedItem(settings.permutations);
    }
    
    void onRun(){
    	_settings.network = Cytoscape.getNetwork(
    			(String)_networkCombo.getSelectedItem());
    	_settings.pval = (Double)_pvTextField.getValue();
    	_settings.directed = _directedGraphRadio.isSelected();
    	_settings.expansion = (Integer)_expansionCombo.getSelectedItem();
    	_settings.knn = (Integer)_knnCombo.getSelectedItem();
    	_settings.geneModuleText = _geneSetTextArea.getText();
    	_settings.knnAuto = _knnSearchRadio.isSelected();
    	_settings.permutations = (Integer)_permutationsCombo.getSelectedItem();
    	_kdaCallback.onKDA(_settings);
    }
    
    void postResults(Map<String,List<KeyDriver>> results){
    	_exportButton.setEnabled(true);
    	_results = results;
    	
    	ListSelectionModel model = new DefaultListSelectionModel();
    	model.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				int rowIdx =  _resultsTbl.getSelectedRow();
				KeyDriverTableModel model = (KeyDriverTableModel)_resultsTbl.getModel();
				Object[] obj = model.alldrivers.get(rowIdx);
				String mod = (String)obj[0];
				KeyDriver kd = (KeyDriver)obj[1];
				_kdaCallback.keyDriverSelect(mod, kd);
			}
		});
    	_resultsTbl.setSelectionModel(model);
    	_resultsTbl.setModel(new KeyDriverTableModel(results));
    	_resultsTbl.setDefaultRenderer(Double.class, new PValRenderer());
    	_resultsTbl.validate();
    }
    
    static class PValRenderer extends DefaultTableCellRenderer {
        DecimalFormat formatter;
        
        public void setValue(Object value) {
            if (formatter==null) {
                formatter = new DecimalFormat("0.##E0");
            }
            if(value == null) setText("");
            else if(((Double)value) == 1.0) setText("1.0");
            else setText(formatter.format(value));
        }
    }

    
    static class KeyDriverTableModel extends AbstractTableModel {
    	
		final static String[] colNames = new String[] { 
				"module", "keydriver","pval",
				"pval-adj","qval","hits",
				"downstream","type" };
		
    	ArrayList<String> moduleNames;
    	ArrayList<Object[]> alldrivers;
    	
    	
    	KeyDriverTableModel(Map<String,List<KeyDriver>> results){
    		moduleNames = new ArrayList<String>(results.keySet());
    		Collections.sort(moduleNames);
        	
        	alldrivers = new ArrayList<Object[]>();
        	for(String m : moduleNames){
        		List<KeyDriver> kds = results.get(m);
        		for(KeyDriver kd : kds){
        			alldrivers.add(new Object[]{ m, kd } );
        		}
        	}
    	}
    	
    	@Override
		public Class<?> getColumnClass(int columnIndex) {
    		if(columnIndex == 2 || columnIndex == 3 || columnIndex == 4)
    			return Double.class;
    		else
				return super.getColumnClass(columnIndex);
		}
    	
    	public String getColumnName(int col) {
	        return colNames[col];
	    }

		public Object getValueAt(int rowIndex, int columnIndex) {
			Object[] obj = alldrivers.get(rowIndex);
			String module = (String)obj[0];
			KeyDriver kd = (KeyDriver)obj[1];
			switch(columnIndex){
				case 0:
					return module;
				case 1:
					return Cytoscape.getNodeAttributes().getAttribute(
							kd.node.getIdentifier(), "canonicalName");
				case 2:
					return kd.pv;
				case 3:
					return kd.pv_bonferroni;
				case 4:
					return kd.qvalue;
				case 5:
					return kd.no_hits;
				case 6:
					return kd.downstreamGenes.size();
				case 7:
					return kd.islocal ? "local" : "global";
				default:
					throw new IllegalStateException("Bad column idx");
			}
		}
		
		public boolean isCellEditable(int row, int col){
			return false;
		}
		
		public int getRowCount() {
			return alldrivers.size();
		}
		
		public int getColumnCount() {
			return colNames.length;
		}
    }
}

