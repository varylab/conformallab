package de.varylab.discreteconformal.plugin;

import static de.jtem.halfedge.util.HalfEdgeUtils.getGenus;
import static de.jtem.halfedgetools.algorithm.triangulation.Delaunay.constructDelaunay;
import static de.varylab.discreteconformal.util.DiscreteRiemannUtility.getHolomorphicForms;
import static de.varylab.discreteconformal.util.LaplaceUtility.calculateCotanWeights;
import static javax.swing.JOptionPane.showMessageDialog;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

import javax.swing.JButton;

import de.jreality.plugin.JRViewer;
import de.jreality.plugin.basic.View;
import de.jreality.ui.LayoutFactory;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Color;
import de.jtem.halfedgetools.adapter.type.EdgeIndex;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.algorithm.triangulation.MappedLengthAdapter;
import de.jtem.halfedgetools.algorithm.triangulation.Triangulator;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.adapter.EuclideanLengthWeightAdapter;
import de.varylab.discreteconformal.adapter.MappedWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CanonicalBasisUtility;
import de.varylab.discreteconformal.util.DualityUtility;

public class DiscreteRiemannPlugin extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private JButton
		calculateButton = new JButton("Calculate Differentials"); 
	
	public DiscreteRiemannPlugin() {
		shrinkPanel.setLayout(new GridBagLayout());
		GridBagConstraints cr = LayoutFactory.createRightConstraint();
		shrinkPanel.add(calculateButton, cr);
		
		calculateButton.addActionListener(this);
	}
	
	public class HarmonicDifferentialAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

		private double[] 
		    dh = null;
		
		public HarmonicDifferentialAdapter(double[] dh) {
			super(null, CoEdge.class, null, Double.class, true, false);
			this.dh = dh;
		}
		
		@Override
		public Double getEdgeValue(CoEdge e, AdapterSet a) {
			Integer index = a.get(EdgeIndex.class, e, Integer.class);
			return dh[index];
		};
		
		@Override
		public String toString() {
			return "Harmonic Differential";
		}
		
	}
	
	public class HolomorphicDifferentialAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

		private Complex[] 
		    dh = null;
		private boolean
			showReal = false;
		private String
			name = "";
		
		public HolomorphicDifferentialAdapter(Complex[] dh, boolean showReal, String name) {
			super(null, CoEdge.class, null, Double.class, true, false);
			this.dh = dh;
			this.showReal = showReal;
			this.name = name;
		}
		
		@Override
		public Double getEdgeValue(CoEdge e, AdapterSet a) {
			Integer index = a.get(EdgeIndex.class, e, Integer.class);
			return showReal ? dh[index].re : dh[index].im;
		};
		
		@Override
		public String toString() {
			return "Holomorphic Differential " + name;
		}
		
	}
	
	@Color
	public class HolomorphicDifferentialColorAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

		private Complex[] 
		    dh = null;
		private boolean
			showReal = false;
		private String
			name = "";
		
		public HolomorphicDifferentialColorAdapter(Complex[] dh, boolean showReal, String name) {
			super(null, CoEdge.class, null, double[].class, true, false);
			this.dh = dh;
			this.showReal = showReal;
			this.name = name;
		}
		
		@Override
		public double[] getEdgeValue(CoEdge e, AdapterSet a) {
			Integer index = a.get(EdgeIndex.class, e, Integer.class);
			double val = showReal ? dh[index].re : dh[index].im;
			return new double[] {10*val, 0, 0};
		};
		
		@Override
		public String toString() {
			return "Holomorphic Differential " + name;
		}
		
	}
	
	@Color
	public class HarmonicDifferentialColor extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

		private double[] 
		    dh = null;
		
		public HarmonicDifferentialColor(double[] dh) {
			super(null, CoEdge.class, null, double[].class, true, false);
			this.dh = dh;
		}
		
		@Override
		public double[] getEdgeValue(CoEdge e, AdapterSet a) {
			Integer index = a.get(EdgeIndex.class, e, Integer.class);
			double val = Math.abs(dh[index]);
			return new double[] {10*val, 0, 0};
		};
		
		@Override
		public String toString() {
			return "Harmonic Differential";
		}
		
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		CoHDS S = hif.get(new CoHDS());
		int g = getGenus(S);
		if (g == 0) {
			showMessageDialog(
				shrinkPanel, 
				"No holomorphic differentials on genus 0 surfaces"
			);
			return;
		}
		AdapterSet a = hif.getAdapters();
		EuclideanLengthWeightAdapter wa = new EuclideanLengthWeightAdapter(null);
		
		// First make clear that we are working with a delaunay triangulation.
		List<CoEdge> newEdges = Triangulator.triangulateByCuttingCorners(S, hif.getAdapters());
		MappedLengthAdapter la = constructDelaunay(S, a);
		a.add(la);
		MappedWeightAdapter cotanWeights = calculateCotanWeights(S, a);
		a.add(cotanWeights);
		for (CoEdge edge : newEdges) {
			if (edge.isValid()) {
				TopologyAlgorithms.removeEdgeFill(edge);
			}
		}
		Complex[][] dhs = getHolomorphicForms(S, a, wa);
		
		int index = 0;
		for (Complex[] dh : dhs) {
			hif.addLayerAdapter(new HolomorphicDifferentialColorAdapter(dh, true, "dHRe" + index), false);
			hif.addLayerAdapter(new HolomorphicDifferentialColorAdapter(dh, false, "dHIm" + index++), false);
		}
		hif.update();
		
		// add introspection adapters
		index = 0;
		for (Complex[] dh : dhs) {
			hif.addLayerAdapter(new HolomorphicDifferentialAdapter(dh, true, "dHRe" + index), false);
			hif.addLayerAdapter(new HolomorphicDifferentialAdapter(dh, false, "dHIm" + index++), false);
		}
		
		CoVertex root = S.getVertex(0);
		List<List<CoEdge>> paths = CanonicalBasisUtility.getCanonicalHomologyBasis(root, a, wa);
		for (List<CoEdge> path : paths) {
			EdgeVectorAdapter eva = new EdgeVectorAdapter(path, "Homology Path " + path.size());
			hif.addLayerAdapter(eva, false);
		}
		
		List<List<CoEdge>> dualpaths = DualityUtility.getDualPaths(S,paths);
		for (List<CoEdge> path : dualpaths) {
			EdgeVectorAdapter eva = new EdgeVectorAdapter(path, "Dual Homology Path " + path.size());
			hif.addLayerAdapter(eva, false);
		}
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
	}
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}
	
	public static void main(String[] args) {
		JRViewer v = new JRViewer();
		v.addContentUI();
		v.addBasicUI();
		v.registerPlugin(new DiscreteRiemannPlugin());
		v.startup();
	}

}
