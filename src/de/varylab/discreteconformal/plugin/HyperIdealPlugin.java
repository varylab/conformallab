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
			1.3287091283689778, 1.328709144060893, 1.1828942665778752, 1.182894274711472, 
			2.1314267110516862, 2.1314266886793214, 2.145393699070382, 2.2924316728851686, 
			2.2362567280260284, 2.2362567650314595, 2.2924316732803294, 2.1453936986802735, 
			2.1314267097029105, 2.2672397657673327, 2.1314266893826184, 1.0416102421811257, 
			-0.4502780224399055, -0.45027802846560555, -0.0411278384808546, -0.0411278214456359, 
			-0.08996021141665035, -0.08996020110009696, -0.030869897541421235, -0.03086991142974196, 
			0.010918417011581107, 0.01091843405536372, -0.04817186915701951, -0.04817188236012954, 
			0.010918434042524154, 0.010918417025070244, -0.04817188210816582, -0.04817186915473929, 
			-0.08996020096657777, -0.0899602107928726, -0.03086991097650844, -0.030869902607411796, 
			-0.450278027759061, -0.450278021732768, -0.04112782131503768, -0.04112783768419421
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
			1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 1.1462158127870863, 
			1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 
			1.7627471360523435, 1.7627471360523435, 1.7627471360523428, 1.7627471360523428, 
			1.7627471360523435, 1.7627471360523428, 1.7627471360523435, 1.7627471360523428, 
			2.633915759978531, 2.633915759978531, 2.633915759978531, 2.633915759978531, 
			2.633915759978531, 2.633915759978531
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
