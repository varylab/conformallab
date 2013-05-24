package de.varylab.discreteconformal.plugin.visualizer;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import javax.swing.JCheckBox;
import javax.swing.JPanel;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.impl.StringAdapter;
import de.jtem.halfedgetools.adapter.type.Label;
import de.jtem.halfedgetools.plugin.VisualizerPlugin;
import de.varylab.discreteconformal.heds.CoVertex;

public class ThetaVisualizer extends VisualizerPlugin {

	private JPanel
		panel = new JPanel();
	private JCheckBox
		boundaryOnlyChecker = new JCheckBox("Boundary Only", true);
	
	
	public ThetaVisualizer() {
		panel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.weightx = GridBagConstraints.REMAINDER;
		c.fill = GridBagConstraints.BOTH;
		c.insets = new Insets(2, 2, 2, 2);
		panel.add(boundaryOnlyChecker);
	}
	
	
	@Label
	private class ThetaLabelAdapter extends StringAdapter {

		private NumberFormat
			nf = new DecimalFormat("0.00");
		
		public ThetaLabelAdapter() {
			super(true, false);
		}

		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return CoVertex.class.isAssignableFrom(nodeClass);
		}

		@Override
		public <
			V extends Vertex<V,E,F>, 
			E extends Edge<V,E,F>, 
			F extends Face<V,E,F>
		> String getV(V v, AdapterSet a) {
			if (!HalfEdgeUtils.isBoundaryVertex(v) && boundaryOnlyChecker.isSelected()) {
				return "";
			} else {
				CoVertex cv = CoVertex.class.cast(v);
				return nf.format(cv.getTheta() / Math.PI) + "pi";
			}
		}
		
		@Override
		public double getPriority() {
			return 1;
		}
		
	}
	
	@Override
	public JPanel getOptionPanel() {
		return panel;
	}
	
	
	@Override
	public AdapterSet getAdapters() {
		AdapterSet result = new AdapterSet();
		result.add(new ThetaLabelAdapter());
		return result;
	}
	
	
	@Override
	public String getName() {
		return "Theta Visualizer";
	}

}
