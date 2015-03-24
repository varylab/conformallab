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
			1.3169578891731424, 1.3169578927823975, 1.3169579011731716, 1.3169579020905602, 
			2.2924316583333115, 2.2924316392049833, 2.292431678989581, 2.292431643089526, 
			2.2924316941659173, 2.292431674031334, 2.2924316940930196, 2.2924316825837465, 
			2.292431682108128, 2.292431635825035, 2.292431663416895, 2.2924316882259688, 
			1.1959012101177978E-8, -2.0563476436963373E-8, -4.456976728393993E-8, -1.8429237951612702E-8, 
			-6.6205808102341355E-9, 1.5717234379207197E-8, -1.7293648067715866E-8, -4.483256556006885E-8, 
			-9.573468060377547E-10, -3.1141227023802106E-8, 1.0908622136947316E-8, 4.776674403327454E-8, 
			7.312334494445497E-9, 3.7181879656685486E-8, 2.14898248315642E-8, 2.8691886984531554E-9, 
			6.47831816268519E-9, -3.331867297470249E-9, 1.2954234834539E-8, 3.051358129457958E-8, 
			4.344545817244695E-9, 1.2223017352538969E-8, -1.6138776412832104E-8, -2.7556701966816714E-8
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
			Assert.assertEquals(lExpected, l, 1E-6);
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
