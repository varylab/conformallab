package de.varylab.discreteconformal.plugin;

import static de.jtem.halfedge.util.HalfEdgeUtils.getGenus;
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
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.varylab.discreteconformal.adapter.EuclideanLengthWeightAdapter;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CanonicalBasisUtility;
import de.varylab.discreteconformal.util.DiscreteHarmonicFormUtility;
import de.varylab.discreteconformal.util.DiscreteRiemannUtility;
import de.varylab.discreteconformal.util.DualityUtility;
import de.varylab.discreteconformal.util.ParallelPathUtility;

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
	
	private class HarmonicDifferentialAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

		private double[] 
		    dh = null;
		
		public HarmonicDifferentialAdapter(double[] dh) {
			super(null, CoEdge.class, null, Double.class, true, false);
			this.dh = dh;
		}
		
		@Override
		public Double getEdgeValue(CoEdge e, AdapterSet a) {
			int index = a.get(EdgeIndex.class, e, Integer.class);
			return dh[index];
		};
		
		@Override
		public String toString() {
			return "Harmonic Differential";
		}
		
	}
	
	private class HolomorphicDifferentialAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

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
			int index = a.get(EdgeIndex.class, e, Integer.class);
			return showReal ? dh[index].re : dh[index].im;
		};
		
		@Override
		public String toString() {
			return "Holomorphic Differential " + name;
		}
		
	}
	
	@Color
	private class HolomorphicDifferentialColorAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

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
			int index = a.get(EdgeIndex.class, e, Integer.class);
			double val = showReal ? dh[index].re : dh[index].im;
			return new double[] {10*val, 0, 0};
		};
		
		@Override
		public String toString() {
			return "Holomorphic Differential " + name;
		}
		
	}
	
	@Color
	private class HarmonicDifferentialColorAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, double[]> {

		private double[] 
		    dh = null;
		
		public HarmonicDifferentialColorAdapter(double[] dh) {
			super(null, CoEdge.class, null, double[].class, true, false);
			this.dh = dh;
		}
		
		@Override
		public double[] getEdgeValue(CoEdge e, AdapterSet a) {
			int index = a.get(EdgeIndex.class, e, Integer.class);
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
////		Complex[][] dhs = DiscreteRiemannUtility.getHolomorphicForms(S, a, wa);
//		double[][] dhs= DiscreteHarmonicFormUtility.getHarmonicFormsOnPrimalMesh(S, a, wa);
//		
//		int index = 0;
//		for (double[] dh : dhs) {
//			hif.addLayerAdapter(new HarmonicDifferentialColorAdapter(dh), false);
//		}
//		hif.update();
//		
//		// add introspection adapters
//		index = 0;
//		for (double[] dh : dhs) {
//			hif.addLayerAdapter(new HarmonicDifferentialAdapter(dh), false);
//			hif.addLayerAdapter(new HarmonicDifferentialAdapter(dh), false);
//		}
		
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
		
		List<List<CoEdge>> parallelpaths = ParallelPathUtility.getLeftParallelPaths(S,paths);
		for (List<CoEdge> path : parallelpaths) {
			EdgeVectorAdapter eva = new EdgeVectorAdapter(path, "Parallel Homology Path " + path.size());
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
