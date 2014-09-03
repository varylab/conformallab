package de.varylab.discreteconformal.plugin;


import java.awt.EventQueue;
import java.awt.Window;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.SelectionInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.jrworkspace.plugin.Controller;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.util.HyperellipticUtility;

public class HyperellipticCurveGenerator extends AlgorithmPlugin {
	
	private final int
		SHEET_PATHS_CHANEL = 8723784;
	private int
		numExtraPoints = 0;
	
	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.Generator;
	}

	@Override
	public String getAlgorithmName() {
		return "Hyperelliptic Curve";
	}
	
	@Override
	public < 
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>,
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void execute(HDS h, AdapterSet a, final HalfedgeInterface hif) {
		Runnable r = new Runnable() {
			@Override
			public void run() {
				Window w = SwingUtilities.getWindowAncestor(hif.getShrinkPanel());
				String numString = JOptionPane.showInputDialog(w, "Number of extra points", numExtraPoints);
				if (numString == null) return;
				numExtraPoints = Integer.parseInt(numString);
			}
		};
		try {
			EventQueue.invokeAndWait(r);
		} catch (Exception e1) {
			throw new RuntimeException(e1);
		}
		
		Selection sel = hif.getSelection();
		int[] branchIndices = new int[sel.getVertices().size()];
		int i = 0;
		for (V v : sel.getVertices(h)) {
			branchIndices[i++] = v.getIndex();
		}
		CoHDS hds = hif.get(new CoHDS());
		Set<CoEdge> glueSet = new HashSet<CoEdge>();
		HyperellipticUtility.generateHyperellipticImage(hds, numExtraPoints, glueSet, branchIndices);
		hif.set(hds);
		Selection s = new Selection();
		s.addAll(glueSet, SHEET_PATHS_CHANEL);
		hif.addSelection(s);
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		SelectionInterface sif = c.getPlugin(SelectionInterface.class);
		sif.registerChannelName(SHEET_PATHS_CHANEL, "Sheet Paths");
	}
	
}

