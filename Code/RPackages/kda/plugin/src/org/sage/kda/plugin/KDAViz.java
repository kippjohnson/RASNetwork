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

import java.awt.Color;

import org.sage.kda.giny.KeyDriver;

import cytoscape.Cytoscape;
import cytoscape.data.CyAttributes;
import cytoscape.visual.ArrowShape;
import cytoscape.visual.CalculatorCatalog;
import cytoscape.visual.EdgeAppearanceCalculator;
import cytoscape.visual.GlobalAppearanceCalculator;
import cytoscape.visual.LineStyle;
import cytoscape.visual.NodeAppearanceCalculator;
import cytoscape.visual.NodeShape;
import cytoscape.visual.VisualMappingManager;
import cytoscape.visual.VisualPropertyType;
import cytoscape.visual.VisualStyle;
import cytoscape.visual.calculators.BasicCalculator;
import cytoscape.visual.calculators.Calculator;
import cytoscape.visual.mappings.DiscreteMapping;
import cytoscape.visual.mappings.ObjectMapping;
import cytoscape.visual.mappings.PassThroughMapping;

public class KDAViz {
	
	public static final String KEYDRIVER_VIZ = "KDAViz:";
	
	/* public attributes */
	
	public static final String MODULE_ATTR = "kda:module:";
	public static final String PV_ATTR = "kda:pvalue:";
	public static final String KEYDRIVER_ATTR = "kda:keydriver:";
	
	public enum DRIVER_TYPE { global, local, module, none };
	
	static void cleanupRun(KDARun run){
		for(KDARun.KDAState state : run._module2State.values()){
			cleanupRun(state);
		}
	}
	
	static void cleanupRun(KDARun.KDAState state){
		CyAttributes nodeAttrs = Cytoscape.getNodeAttributes();
		for(String s : nodeAttrs.getAttributeNames()){
			if(s.endsWith(state.moduleName))
				nodeAttrs.deleteAttribute(s);
		}
		VisualMappingManager manager = Cytoscape.getVisualMappingManager();
        manager.getCalculatorCatalog().removeVisualStyle(makeAttribute(KEYDRIVER_VIZ, state));
	}
	
	private static String makeAttribute(String name, KDARun.KDAState state){
		return name + state.moduleName;
	}
	
	static VisualStyle setupAttributeVizForRun(KDARun.KDAState currentRun){
		//clearKDAAttributes();
		CyAttributes nodeAttrs = Cytoscape.getNodeAttributes();
		
		for(Node n : currentRun.moduleNodes){
			// mark module nodes; if key driver, will be overwritten below
			nodeAttrs.setAttribute(n.getIdentifier(), 
					makeAttribute(MODULE_ATTR, currentRun) , 
					Boolean.TRUE);
			nodeAttrs.setAttribute(n.getIdentifier(), 
					makeAttribute(KEYDRIVER_ATTR, currentRun), 
					DRIVER_TYPE.module.name());				
		}
		
		for(KeyDriver kd : currentRun.keyDrivers){
			nodeAttrs.setAttribute(kd.node.getIdentifier(), 
					makeAttribute(PV_ATTR, currentRun), kd.pv);
			nodeAttrs.setAttribute(kd.node.getIdentifier(), 
					makeAttribute(KEYDRIVER_ATTR,currentRun), 
						(kd.isboosted || kd.islocal) 
							? DRIVER_TYPE.local.name() 
							: DRIVER_TYPE.global.name());
			
		}	
		VisualStyle viz = createVisualStyle(currentRun);
		VisualMappingManager manager = Cytoscape.getVisualMappingManager();
        CalculatorCatalog catalog = manager.getCalculatorCatalog();
		catalog.addVisualStyle(viz);
		return viz;
	}
	
	static VisualStyle createVisualStyle(KDARun.KDAState run){
		
		double defaultNodeSize = 35.0;
		
	    //methods to access the node, edge, and global appearance calculators
	    NodeAppearanceCalculator nodeAppCalc = new NodeAppearanceCalculator();
	    EdgeAppearanceCalculator edgeAppCalc = new EdgeAppearanceCalculator();
	    GlobalAppearanceCalculator globalAppCalc = new GlobalAppearanceCalculator();
	
		globalAppCalc.setDefaultBackgroundColor(Color.white);
		nodeAppCalc.getDefaultAppearance().set(VisualPropertyType.NODE_FILL_COLOR, Color.LIGHT_GRAY);
		nodeAppCalc.getDefaultAppearance().set(VisualPropertyType.NODE_SHAPE, NodeShape.ELLIPSE);
		nodeAppCalc.getDefaultAppearance().set(VisualPropertyType.NODE_SIZE, defaultNodeSize);
		
		edgeAppCalc.getDefaultAppearance().set(VisualPropertyType.EDGE_TGTARROW_SHAPE, ArrowShape.ARROW);
		edgeAppCalc.getDefaultAppearance().set(VisualPropertyType.EDGE_TGTARROW_COLOR, Color.BLACK);
		
		
	    DiscreteMapping lineMapping = new DiscreteMapping(LineStyle.SOLID,
	    												makeAttribute(MODULE_ATTR,run),
	                                                    ObjectMapping.NODE_MAPPING);
	    lineMapping.putMapValue(Boolean.TRUE, LineStyle.LONG_DASH);
	    lineMapping.putMapValue(Boolean.FALSE, LineStyle.SOLID);
	    Calculator lineCalculator = new BasicCalculator("Line Calculator", 
	    		lineMapping, VisualPropertyType.NODE_LINE_STYLE);
	    nodeAppCalc.setCalculator(lineCalculator);
	
	    // COLOR
	    DiscreteMapping colorMapping = new DiscreteMapping(Color.LIGHT_GRAY, 
				makeAttribute(KEYDRIVER_ATTR,run),
	            ObjectMapping.NODE_MAPPING);
	    colorMapping.putMapValue(DRIVER_TYPE.global.name(), Color.GREEN);
	    colorMapping.putMapValue(DRIVER_TYPE.local.name(), Color.RED);
	    colorMapping.putMapValue(DRIVER_TYPE.module.name(), Color.BLUE);
	    BasicCalculator driverCalculator = new BasicCalculator("PV Calc", 
	    		colorMapping, VisualPropertyType.NODE_FILL_COLOR);
	    nodeAppCalc.setCalculator(driverCalculator);
	    
	    //SIZE
	    DiscreteMapping sizeMapping = new DiscreteMapping(defaultNodeSize, 
				makeAttribute(KEYDRIVER_ATTR,run),
	            ObjectMapping.NODE_MAPPING);
		sizeMapping.putMapValue(DRIVER_TYPE.global.name(), defaultNodeSize * 2.0);
		Calculator sizeCalculator = new BasicCalculator("Shape Calculator", 
				sizeMapping, VisualPropertyType.NODE_SIZE);
		nodeAppCalc.setCalculator(sizeCalculator);
	    
	//    ContinuousMapping colorMapping = new ContinuousMapping(defaultObj, ObjectMapping.NODE_MAPPING);
	//    colorMapping.setControllingAttributeName(PV_ATTR, Cytoscape.getCurrentNetwork(), false);
	//    Interpolator numToColor = new LinearNumberToColorInterpolator();
	//    BoundaryRangeValues bv1 = new BoundaryRangeValues(defaultObj, Color.RED, Color.RED);
	//    BoundaryRangeValues bv2 = new BoundaryRangeValues(Color.RED, Color.RED, defaultObj);
	//    colorMapping.setInterpolator(numToColor);
	//    colorMapping.addPoint(0.0, bv1);
	//    colorMapping.addPoint(0.05, bv2);
	//    
	//    BasicCalculator pvCalculator = new BasicCalculator("PV Calc", colorMapping, VisualPropertyType.NODE_FILL_COLOR);
	//    nodeAppCalc.setCalculator(pvCalculator);
	    
	    //disMapping.setControllingAttributeName(shapeAttribute, network, false);
	//    disMapping.putMapValue(new Integer(1), Color.red);
	//    Calculator colorCalculator = new BasicCalculator("Module Calculator",
	//                                                      disMapping,
	//                                                     VisualPropertyType.NODE_FILL_COLOR);
	//    nodeAppCalc.setCalculator(colorCalculator);
	//    
	//    nodeAppCalc.setCalculator(shapeCalculator);
	//    
	    PassThroughMapping pm = new PassThroughMapping(new String(), "canonicalName");
	    Calculator nlc = new BasicCalculator("Node Label Calculator",
	                                         pm, VisualPropertyType.NODE_LABEL);
	    nodeAppCalc.setCalculator(nlc);
	//
	//
	//
	//    DiscreteMapping discreteColorMapping = new DiscreteMapping(Color.WHITE, 
	//    											colorAttribute,
	//                                                ObjectMapping.NODE_MAPPING);
	//    
	//    discreteColorMapping.putAll(colorMapping);
	//    
	//    Calculator nodeColorCalculator = new BasicCalculator("R-mapped fill color",
	//                                                    discreteColorMapping,
	//                                                    VisualPropertyType.NODE_FILL_COLOR);
	//    nodeAppCalc.setCalculator(nodeColorCalculator);
	//    Calculator nodeBorderColorCalculator = new BasicCalculator("R-mapped edge color",
	//            discreteColorMapping,
	//            VisualPropertyType.NODE_BORDER_COLOR);
	//    nodeAppCalc.setCalculator(nodeBorderColorCalculator);
	//    
	//    nodeAppCalc.getDefaultAppearance().set(VisualPropertyType.NODE_LABEL_POSITION, 
	//    		new LabelPosition(LabelPosition.convert(LabelPosition.northName),
	//    						  LabelPosition.convert(LabelPosition.southName),
	//    						  LabelPosition.convert(LabelPosition.justifyCenterName),
	//    						  0.0, 0.0));
	//    
	//    edgeAppCalc.getDefaultAppearance().set(VisualPropertyType.EDGE_COLOR, Color.LIGHT_GRAY);
	//    
	    return new VisualStyle(makeAttribute(KEYDRIVER_VIZ,run), 
	    		nodeAppCalc, edgeAppCalc, globalAppCalc);
	}

}
