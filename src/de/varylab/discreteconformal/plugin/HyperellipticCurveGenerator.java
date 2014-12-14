package de.varylab.discreteconformal.plugin;


import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import de.jreality.ui.LayoutFactory;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.SelectionInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmDialogPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.jrworkspace.plugin.Controller;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.util.DiscreteEllipticUtility;

public class HyperellipticCurveGenerator extends AlgorithmDialogPlugin {
	
	private final int
		SHEET_PATHS_CHANEL = 8723784;
	private JPanel
		options = new JPanel();
	private SpinnerNumberModel	
		numPointsModel = new SpinnerNumberModel(0, 0, 10000, 1);
	private JSpinner
		numPointsSpinner = new JSpinner(numPointsModel);
	private JCheckBox
		projectChecker = new JCheckBox("Project To Sphere", true),
		convexHullChecker = new JCheckBox("Project To Sphere", true);
	
	public HyperellipticCurveGenerator() {
		GridBagConstraints lc = LayoutFactory.createLeftConstraint();
		GridBagConstraints rc = LayoutFactory.createRightConstraint();
		options.setLayout(new GridBagLayout());
		options.add(new JLabel("Random Points"), lc);
		options.add(numPointsSpinner, rc);
		options.add(projectChecker, rc);
		options.add(convexHullChecker, rc);
	}
	
	@Override
	protected JPanel getDialogPanel() {
		return options;
	}
	
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
	> void executeAfterDialog(HDS h, AdapterSet a, final HalfedgeInterface hif) {
		int numPoints = numPointsModel.getNumber().intValue();
		boolean project = projectChecker.isSelected();
		boolean convecHull = convexHullChecker.isSelected();
		Selection sel = hif.getSelection();
		int[] branchIndices = new int[sel.getVertices().size()];
		int i = 0;
		for (V v : sel.getVertices(h)) {
			branchIndices[i++] = v.getIndex();
		}
		CoHDS hds = hif.get(new CoHDS());
		Set<CoEdge> glueSet = new HashSet<CoEdge>();
		DiscreteEllipticUtility.generateEllipticImage(hds, numPoints, project, convecHull, glueSet, branchIndices);
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

