package de.varylab.discreteconformal.plugin.algorithm;

import static de.jtem.halfedgetools.algorithm.triangulation.Triangulator.triangulateByCuttingCorners;
import static de.varylab.discreteconformal.uniformization.SurfaceCurveUtility.createIntersectionVertices;
import static de.varylab.discreteconformal.unwrapper.EuclideanLayout.doLayout;
import static de.varylab.discreteconformal.util.CuttingUtility.glueAll;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.LengthTex;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition2d;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmCategory;
import de.jtem.halfedgetools.plugin.algorithm.AlgorithmDialogPlugin;
import de.jtem.halfedgetools.selection.Selection;
import de.jtem.halfedgetools.selection.TypedSelection;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.PluginInfo;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.plugin.DiscreteConformalPlugin;
import de.varylab.discreteconformal.plugin.image.ImageHook;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class CutAndGlueConformalDomain extends AlgorithmDialogPlugin implements ActionListener {

	private DiscreteConformalPlugin
		dcp = null;
	private HalfedgeInterface
		hif = null;
	private JPanel
		options = new JPanel();
	private JRadioButton
		orthogonalRadio = new JRadioButton("Orthogonal to Boundary", true),
		fixedAngleRadio = new JRadioButton("Predefined Angle");
	private SpinnerNumberModel
		angleModel = new SpinnerNumberModel(0, 0, 360, 0.1);
	private JSpinner
		angleSpinner = new JSpinner(angleModel);
	
	public CutAndGlueConformalDomain() {
		options.setPreferredSize(new Dimension(350, 150));
		options.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.weightx = 1.0;
		c.fill = GridBagConstraints.BOTH;
		c.gridwidth = GridBagConstraints.REMAINDER;
		c.insets = new Insets(2, 2, 2, 2);
		options.add(orthogonalRadio, c);
		c.gridwidth = GridBagConstraints.RELATIVE;
		options.add(fixedAngleRadio, c);
		c.gridwidth = GridBagConstraints.REMAINDER;
		options.add(angleSpinner, c);
		
		orthogonalRadio.addActionListener(this);
		fixedAngleRadio.addActionListener(this);
		
		ButtonGroup modeGroup = new ButtonGroup();
		modeGroup.add(orthogonalRadio);
		modeGroup.add(fixedAngleRadio);
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>, 
		HDS extends HalfEdgeDataStructure<V, E, F>
	> void executeAfterDialog(HDS data, AdapterSet a, HalfedgeInterface hcp) throws Exception {
		if (!(data instanceof CoHDS)) {
			throw new RuntimeException("Can only work with conformal data structures");
		}
		CoHDS hds = (CoHDS) data;
		TypedSelection<CoVertex> sel = hif.getSelection().getVertices(hds);
		if (sel.isEmpty()) {
			throw new RuntimeException("Please select a boundary vertex to define the cut location and direction");
		}
		CoVertex v = sel.iterator().next();
		if (!HalfEdgeUtils.isBoundaryVertex(v)) {
			throw new RuntimeException("Please select a boundary vertex to define the cut location and direction");
		}
		double[][] segment = null;
		if (orthogonalRadio.isSelected()) {
			CoEdge be = HalfEdgeUtils.incomingBoundaryEdge(v);
			CoVertex v1 = be.getStartVertex(); 
			CoVertex v2 = be.getOppositeEdge().getPreviousEdge().getStartVertex();
			double[] pv = a.getD(TexturePosition2d.class, v);
			double[] pv1 = a.getD(TexturePosition2d.class, v1);
			double[] pv2 = a.getD(TexturePosition2d.class, v2);
			double[] vec1 = Rn.subtract(null, pv1, pv);
			double[] vec2 = Rn.subtract(null, pv2, pv);
			Rn.normalize(vec1, vec1);
			Rn.normalize(vec2, vec2);
			double[] difVec = Rn.subtract(null, Rn.add(null, pv, vec1), Rn.add(null, pv, vec2));
			segment = new double[][]{{pv[0], pv[1], 1}, {difVec[1] + pv[0], difVec[0] + pv[1], 1}};
		} 
		if (fixedAngleRadio.isSelected()) {
			double angle = Math.toRadians(angleModel.getNumber().doubleValue());
			double[] pv = a.getD(TexturePosition2d.class, v);
			double[] pv2 = {pv[0] + Math.cos(angle), pv[1] + Math.sin(angle)};
			segment = new double[][]{{pv[0], pv[1], 1}, {pv2[0], pv2[1], 1}};
		}

		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = dcp.getCurrentCutInfo();
		double snapTolerance = 1E-5;
		int signature = Pn.EUCLIDEAN;
		Set<CoVertex> newVertices = new HashSet<>();
		createIntersectionVertices(segment, false, hds, hds, cutInfo, snapTolerance, signature, newVertices);
		triangulateByCuttingCorners(hds, a);
		
		Map<CoEdge, Double> lMap = new HashMap<>();
		for (CoEdge edge : hds.getPositiveEdges()) {
			Double l = a.get(LengthTex.class, edge, Double.class);
			lMap.put(edge, l);
			lMap.put(edge.getOppositeEdge(), l);
		}
		Map<CoEdge, Double> alphaMap = new HashMap<>();
		for (CoEdge edge : hds.getEdges()) {
			if (edge.getLeftFace() == null) {
				alphaMap.put(edge, 0.0);
				continue;
			}
			double[] p = a.getD(TexturePosition2d.class, edge.getNextEdge().getTargetVertex());
			double[] p1 = a.getD(TexturePosition2d.class, edge.getTargetVertex());
			double[] p2 = a.getD(TexturePosition2d.class, edge.getStartVertex());
			double alpha = Rn.euclideanAngle(Rn.subtract(null, p1, p), Rn.subtract(null, p2, p));
			alphaMap.put(edge, alpha);
		}
		glueAll(hds, cutInfo);

		Selection pathSelection = new Selection();
		pathSelection.addAll(newVertices);
		LinkedList<CoEdge> newCut = DiscreteConformalPlugin.selectCutPath(hds, newVertices, pathSelection);
		CuttingInfo<CoVertex, CoEdge, CoFace> newCutInfo = new CuttingInfo<CoVertex, CoEdge, CoFace>();
		CuttingUtility.cutAlongPath(newCut, newCutInfo);
		for (CoEdge cutEdge : newCutInfo.edgeCutMap.keySet()) {
			if (!lMap.containsKey(cutEdge)) continue;
			double l = lMap.get(cutEdge);
			lMap.put(cutEdge.getOppositeEdge(), l);
		}
		dcp.setCutCurrentInfo(newCutInfo);
		
		doLayout(hds, lMap, alphaMap);
		
		hif.update();
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		dcp = c.getPlugin(DiscreteConformalPlugin.class);
		hif = c.getPlugin(HalfedgeInterface.class);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		updateStates();
	}
	private void updateStates() {
		angleSpinner.setEnabled(fixedAngleRadio.isSelected());
	}
	@Override
	protected JPanel getDialogPanel() {
		return options;
	}
	@Override
	public AlgorithmCategory getAlgorithmCategory() {
		return AlgorithmCategory.TextureRemeshing;
	}
	@Override
	public String getAlgorithmName() {
		return "Cut And Glue Texture Domain";
	}
	@Override
	public PluginInfo getPluginInfo() {
		PluginInfo info = super.getPluginInfo();
		info.icon = ImageHook.getIcon("cut.png");
		return info;
	}

}
