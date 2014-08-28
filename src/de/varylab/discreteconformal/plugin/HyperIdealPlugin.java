package de.varylab.discreteconformal.plugin;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static de.varylab.discreteconformal.plugin.TargetGeometry.Hyperbolic;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import javax.swing.JButton;

import org.junit.Assert;

import de.jreality.math.Pn;
import de.jreality.plugin.scene.SceneShrinkPanel;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jpetsc.Vec;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.functional.HyperIdealFunctionalTest;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.HyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperIdealApplication;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class HyperIdealPlugin extends SceneShrinkPanel implements ActionListener {

	private double[]
		lawsonData = {1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531};
	private JButton
		lawsonButton = new JButton("Lawson's Surface");
	private DiscreteConformalPlugin
		conformalPlugin = null;
	private ConformalDataPlugin
		conformalDataPlugin = null;
	
	public HyperIdealPlugin() {
		shrinkPanel.setTitle("Hyper-Ideal Uniformization");
		shrinkPanel.add(lawsonButton);
		lawsonButton.addActionListener(this);
	}
	
	@Override
	public void actionPerformed(ActionEvent ev) {
		Tao.Initialize();
		CoHDS hds = HyperIdealFunctionalTest.createLawsonsSurface();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		Vec u = new Vec(app.getDomainDimension());
		for (int i = 0; i < app.getDomainDimension(); i++) {
			u.setValue(i, lawsonData[i], INSERT_VALUES);
		}
		// initialize angles
		app.evaluateObjectiveAndGradient(u, null);

		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));
		
		Set<CoEdge> cutEdges = new LinkedHashSet<>();
//		cutEdges.add(hds.getEdge(0));
//		cutEdges.add(hds.getEdge(1));
//		cutEdges.add(hds.getEdge(3));
//		cutEdges.add(hds.getEdge(5));
//		cutEdges.add(hds.getEdge(6));
//		cutEdges.add(hds.getEdge(8));
//		cutEdges.add(hds.getEdge(11));
//		cutEdges.add(hds.getEdge(13));
//		cutEdges.add(hds.getEdge(14));
//		cutEdges.add(hds.getEdge(16));
//		cutEdges.add(hds.getEdge(19));
//		cutEdges.add(hds.getEdge(21));
//		cutEdges.add(hds.getEdge(22));
//		cutEdges.add(hds.getEdge(23));
		
		cutEdges.add(hds.getEdge(2));
		cutEdges.add(hds.getEdge(1));
		cutEdges.add(hds.getEdge(5));
		cutEdges.add(hds.getEdge(8));
		cutEdges.add(hds.getEdge(11));
		cutEdges.add(hds.getEdge(17));
		cutEdges.add(hds.getEdge(14));
		cutEdges.add(hds.getEdge(18));
		cutEdges.add(hds.getEdge(21));
		
		TopologyAlgorithms.flipEdge(hds.getEdge(27));
		TopologyAlgorithms.flipEdge(hds.getEdge(31));
		TopologyAlgorithms.flipEdge(hds.getEdge(35));
		
		
		CoVertex root = hds.getVertex(2);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		CuttingUtility.cutAtEdges(cutInfo, cutEdges);
		cutInfo.cutRoot = root;
		Assert.assertEquals(0, HalfEdgeUtils.getGenus(hds));
		
		hds.getVertex(0).P = new double[]{0,0,0,1};
		hds.getVertex(1).P = new double[]{1,0,0,1};
		hds.getVertex(2).P = new double[]{0,1,0,1};
		hds.getVertex(3).P = new double[]{1,1,0,1};
		hds.getVertex(4).P = new double[]{2,0,0,1};
		hds.getVertex(5).P = new double[]{2,1,0,1};
		hds.getVertex(6).P = new double[]{1,2,0,1};
		hds.getVertex(7).P = new double[]{2,2,0,1};
		hds.getVertex(8).P = new double[]{3,2,0,1};
		hds.getVertex(9).P = new double[]{3,1,0,1};
		hds.getVertex(10).P = new double[]{2,3,0,1};
		hds.getVertex(11).P = new double[]{3,3,0,1};
		hds.getVertex(12).P = new double[]{4,3,0,1};
//		hds.getVertex(13).P = new double[]{4,2,0,1};
		
		Map<CoEdge, Double> lMap = new LinkedHashMap<>();
		for (CoEdge e : hds.getPositiveEdges()) {
			if (e.getLeftFace() == null) {
				e = e.getOppositeEdge();
			}
			double l = app.getEdgeLength(e, u);
			lMap.put(e, l);
			lMap.put(e.getOppositeEdge(), l);
		}
		
		// write beta angles to alpha for the layout
		for (CoEdge e : hds.getEdges()) {
			e.setAlpha(e.getBeta());
		}
		
		HyperbolicLayout.doLayout(hds, root, lMap);
		for (CoEdge e : hds.getPositiveEdges()) {
			double[] s = e.getStartVertex().T;
			double[] t = e.getTargetVertex().T;
			double lExpected = lMap.get(e);
			double l = Pn.distanceBetween(s, t, Pn.HYPERBOLIC);
			Assert.assertEquals(lExpected, l, 1E-5);
		}
		
		conformalDataPlugin.addDiscreteMap("Uniformizing Map", hds, cutInfo);
		conformalPlugin.createUniformization(hds, Hyperbolic, cutInfo);
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		conformalPlugin = c.getPlugin(DiscreteConformalPlugin.class);
		conformalDataPlugin = c.getPlugin(ConformalDataPlugin.class);
	}
	
}
