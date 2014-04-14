package de.varylab.discreteconformal.plugin;


import java.awt.Window;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Color;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.util.HyperellipticUtility;

public class EllipticImageGenerator extends AlgorithmPlugin {
	
	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Generator;
	}

	@Override
	public String getAlgorithmName() {
		return "Elliptic Image";
	}
	
	
	@Override
	public < 
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS h, AdapterSet a, HalfedgeInterface hif) {
		Window w = SwingUtilities.getWindowAncestor(hif.getShrinkPanel());
		String numString = JOptionPane.showInputDialog(w, "Number of extra points", 0);
		if (numString == null) return;
		int extraPoints = Integer.parseInt(numString);
		Selection sel = hif.getSelection();
		int[] branchIndices = new int[sel.getVertices().size()];
		int i = 0;
		for (V v : sel.getVertices(h)) {
			branchIndices[i++] = v.getIndex();
		}
		CoHDS hds = hif.get(new CoHDS());
		Set<CoEdge> glueSet = new HashSet<CoEdge>();
//		DiscreteEllipticUtility.generateEllipticImage(hds, extraPoints, glueSet, branchIndices);
		HyperellipticUtility.generateHyperellipticImage(hds, extraPoints, glueSet, branchIndices);
		PathVisualizer pathVisualizer = new PathVisualizer();
		for (CoEdge e : glueSet) {
			pathVisualizer.add(e);
			pathVisualizer.add(e.getOppositeEdge());
		}
		// show the result
		hif.addLayerAdapter(pathVisualizer, true);
		hif.set(hds);
	}
	
	
	@Color
	public static class PathVisualizer extends AbstractAdapter<double[]> {

		private Set<CoEdge>
			edges = new HashSet<CoEdge>();
		private double[]
		    pathColor = {1, 0, 0},
		    defaultColor = {1, 1, 1};
		
		public PathVisualizer() {
			super(double[].class, true, false);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return nodeClass == CoEdge.class;
		}

		@Override
		public double getPriority() {
			return 10;
		}
		
		@Override
		public <
			V extends de.jtem.halfedge.Vertex<V,E,F>, 
			E extends de.jtem.halfedge.Edge<V,E,F>, 
			F extends de.jtem.halfedge.Face<V,E,F>
		> double[] getE(E e, de.jtem.halfedgetools.adapter.AdapterSet a) {
			if (edges.contains(e)) {
				return pathColor;
			} else {
				return defaultColor;
			}
		};
		
		public void add(CoEdge edge) {
			edges.add(edge);
		}

	}

	
}

