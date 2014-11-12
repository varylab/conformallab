package de.varylab.discreteconformal.plugin;

import static de.jtem.jpetsc.InsertMode.INSERT_VALUES;
import static de.jtem.jpetsc.NormType.NORM_FROBENIUS;
import static de.varylab.discreteconformal.plugin.TargetGeometry.Hyperbolic;

import java.awt.GridLayout;
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
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.jpetsc.Vec;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jtao.Tao;
import de.varylab.discreteconformal.functional.HyperIdealGenerator;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.unwrapper.HyperbolicLayout;
import de.varylab.discreteconformal.unwrapper.numerics.CHyperIdealApplication;
import de.varylab.discreteconformal.util.CuttingUtility;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

public class HyperIdealPlugin extends SceneShrinkPanel implements ActionListener {

	private JButton
		lawsonSquareTiledButton = new JButton("Lawson's Square Tiled"),
		lawsonSquareTiledBranchButton = new JButton("Lawson's Square Tiled Branch");
	private DiscreteConformalPlugin
		conformalPlugin = null;
	private ConformalDataPlugin
		conformalDataPlugin = null;
	
	public HyperIdealPlugin() {
		shrinkPanel.setTitle("Hyper-Ideal Uniformization");
		shrinkPanel.setLayout(new GridLayout(2, 1));
		shrinkPanel.add(lawsonSquareTiledButton);
		shrinkPanel.add(lawsonSquareTiledBranchButton);
		lawsonSquareTiledButton.addActionListener(this);
		lawsonSquareTiledBranchButton.addActionListener(this);
	}
	
	static {
		Tao.Initialize();
	}
	
	@Length
	private class LawsonMetric extends AbstractAdapter<Double> {

		public LawsonMetric() {
			super(Double.class, true, false);
		}
		
		@Override
		public <
			V extends Vertex<V, E, F>, 
			E extends Edge<V, E, F>, 
			F extends Face<V, E, F>
		> Double getE(E e, AdapterSet a) {
			return e.getIndex() < 24 ? 1.0 : Math.sqrt(2.0);
		}
		
		@Override
		public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
			return true;
		}
		
	}
	
	@Override
	public void actionPerformed(ActionEvent ev) {
		if (lawsonSquareTiledButton == ev.getSource()) {
			uniformizeLawsonSquareTiled();
		}
		if (lawsonSquareTiledBranchButton == ev.getSource()) {
			uniformizeLawsonSquareTiledBranch();
		}
	}

	
	private void uniformizeLawsonSquareTiledBranch() {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiledWithBranchPoints();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		Vec u = new Vec(app.getDomainDimension());
		double[] lawsonData = {
			1.3287091215473186, 1.3287091235682074, 1.1828942571436754, 1.1828942582027029, 
			2.1314267018765376, 2.131426701401825, 2.1453936916905043, 2.2924316715616246, 
			2.2362567084524674, 2.2362567270890032, 2.292431670745553, 2.1453936858633007, 
			2.1314267006971437, 2.2672397560325046, 2.1314267038798467, 1.0416102536948757, 
			-0.4502780169195045, -0.4502780155597613, -0.041127823384474675, -0.04112782383611047, 
			-0.0899601865781942, -0.08996019016545298, -0.0308698979391615, -0.03086989627601265, 
			0.01091841034336131, 0.010918425260716142, -0.048171888154962574, -0.04817189956520202, 
			0.010918424084159396, 0.010918408554520824, -0.04817190541203753, -0.04817189322423887, 
			-0.0899601964256837, -0.08996019138046285, -0.030869894625751414, -0.030869899162547093, 
			-0.4502780123786227, -0.45027801548011975, -0.04112782420557621, -0.04112781875211804	
		};
		for (int i = 0; i < app.getDomainDimension(); i++) {
			u.setValue(i, lawsonData[i], INSERT_VALUES);
		}
		// initialize angles
		Vec G = new Vec(app.getDomainDimension());
		app.evaluateObjectiveAndGradient(u, G);
		Assert.assertEquals(0.0, G.norm(NORM_FROBENIUS), 1E-6);
		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));

		Set<CoEdge> cutEdges = new LinkedHashSet<>();
		// six quads around a vertex
		cutEdges.add(hds.getEdge(2));
		cutEdges.add(hds.getEdge(1));
		cutEdges.add(hds.getEdge(5));
		cutEdges.add(hds.getEdge(8));
		cutEdges.add(hds.getEdge(11));
		cutEdges.add(hds.getEdge(17));
		cutEdges.add(hds.getEdge(14));
		cutEdges.add(hds.getEdge(18));
		cutEdges.add(hds.getEdge(21));
		
		CoVertex root = hds.getVertex(2);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		CuttingUtility.cutAtEdges(cutInfo, cutEdges);
		cutInfo.cutRoot = root;
		Assert.assertEquals(0, HalfEdgeUtils.getGenus(hds));
		
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
			Assert.assertEquals(lExpected, l, 1E-7);
		}
		conformalDataPlugin.addHalfedgeMap("Uniformizing Map", hds, cutInfo);
		conformalPlugin.createUniformization(hds, Hyperbolic, cutInfo);
	}
	
	
	private void uniformizeLawsonSquareTiled() {
		CoHDS hds = HyperIdealGenerator.createLawsonSquareTiled();
		CHyperIdealApplication app = new CHyperIdealApplication(hds);
		Vec u = new Vec(app.getDomainDimension());
		double[] lawsonData = {
			1.1462158347805917, 1.1462158347805917, 1.1462158347805917, 1.1462158347805917, 
			1.7627471740392668, 1.7627471740392668, 1.7627471740392684, 1.7627471740392684, 
			1.7627471740392668, 1.7627471740392668, 1.7627471740392684, 1.7627471740392684, 
			1.7627471740392668, 1.7627471740392684, 1.7627471740392668, 1.7627471740392684, 
			2.6339157938497704, 2.6339157938497704, 2.6339157938497704, 2.6339157938497704, 
			2.6339157938497704, 2.6339157938497704
		};
		for (int i = 0; i < app.getDomainDimension(); i++) {
			u.setValue(i, lawsonData[i], INSERT_VALUES);
		}
		// initialize angles
		app.evaluateObjectiveAndGradient(u, null);

		Assert.assertEquals(2, HalfEdgeUtils.getGenus(hds));
		
		Set<CoEdge> cutEdges = new LinkedHashSet<>();
		// stair case of quads
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
		
		// six quads around a vertex
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

		conformalDataPlugin.addDiscreteMetric("Lawsons Squares", hds, new AdapterSet(new LawsonMetric()));
		
		CoVertex root = hds.getVertex(2);
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo = new CuttingInfo<>();
		CuttingUtility.cutAtEdges(cutInfo, cutEdges);
		cutInfo.cutRoot = root;
		Assert.assertEquals(0, HalfEdgeUtils.getGenus(hds));
		
		// stair case of quads
//		hds.getVertex(0).P = new double[]{0,0,0,1};
//		hds.getVertex(1).P = new double[]{1,0,0,1};
//		hds.getVertex(2).P = new double[]{0,1,0,1};
//		hds.getVertex(3).P = new double[]{1,1,0,1};
//		hds.getVertex(4).P = new double[]{2,0,0,1};
//		hds.getVertex(5).P = new double[]{2,1,0,1};
//		hds.getVertex(6).P = new double[]{1,2,0,1};
//		hds.getVertex(7).P = new double[]{2,2,0,1};
//		hds.getVertex(8).P = new double[]{3,2,0,1};
//		hds.getVertex(9).P = new double[]{3,1,0,1};
//		hds.getVertex(10).P = new double[]{2,3,0,1};
//		hds.getVertex(11).P = new double[]{3,3,0,1};
//		hds.getVertex(12).P = new double[]{4,3,0,1};
//		hds.getVertex(13).P = new double[]{4,2,0,1};
		
		// six quads around a vertex
		hds.getVertex(0).P = new double[]{1,0,0,1};
		hds.getVertex(1).P = new double[]{0,0,0,1};
		hds.getVertex(2).P = new double[]{1,1,0,1};
		hds.getVertex(3).P = new double[]{2,1,0,1};
		hds.getVertex(4).P = new double[]{2,0,0,1};
		hds.getVertex(5).P = new double[]{2,2,0,1};
		hds.getVertex(6).P = new double[]{1,0,0,1};
		hds.getVertex(7).P = new double[]{2,0,0,1};
		hds.getVertex(8).P = new double[]{2,1,0,1};
		hds.getVertex(9).P = new double[]{1,2,0,1};
		hds.getVertex(10).P = new double[]{2,2,0,1};
		hds.getVertex(11).P = new double[]{0,0,0,1};
		hds.getVertex(12).P = new double[]{0,1,0,1};
		
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
		
		conformalDataPlugin.addHalfedgeMap("Uniformizing Map", hds, cutInfo);
		conformalPlugin.createUniformization(hds, Hyperbolic, cutInfo);
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		conformalPlugin = c.getPlugin(DiscreteConformalPlugin.class);
		conformalDataPlugin = c.getPlugin(ConformalDataPlugin.class);
	}
	
}
