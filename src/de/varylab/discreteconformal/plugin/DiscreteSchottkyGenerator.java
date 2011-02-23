package de.varylab.discreteconformal.plugin;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jreality.plugin.basic.View;
import de.jreality.ui.LayoutFactory;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.java2dx.beans.Viewer2DWithInspector;
import de.jtem.java2dx.modelling.GraphicsModeller2D;
import de.jtem.java2dx.modelling.SimpleModeller2D;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.util.SurgeryUtility;

public class DiscreteSchottkyGenerator extends ShrinkPanelPlugin implements ActionListener, ChangeListener {

	private HalfedgeInterface
		hif = null;
	private SimpleModeller2D
		moddeller = new GraphicsModeller2D();
	private Viewer2DWithInspector
		viewer = moddeller.getViewer();
	private JScrollPane 
		protectorPane = new JScrollPane(viewer); 
	private JButton
		generateButton = new JButton("Generate Surface");
	private SpinnerNumberModel
		genusModel = new SpinnerNumberModel(2, 0, 20, 1);
	private JSpinner
		genusSpinner = new JSpinner(genusModel);
	
//	
//	private List<Complex>
//		fixPointsA = new LinkedList<Complex>(),
//		fixPointsB = new LinkedList<Complex>(),
//		muList = new LinkedList<Complex>();
		
	public DiscreteSchottkyGenerator() {
		GridBagConstraints c1 = LayoutFactory.createLeftConstraint();
		GridBagConstraints c2 = LayoutFactory.createRightConstraint();
		viewer.setPreferredSize(new Dimension(200, 200));
		viewer.setMinimumSize(viewer.getPreferredSize());
		shrinkPanel.setLayout(new GridBagLayout());
		shrinkPanel.add(new JLabel("Genus"), c1);
		shrinkPanel.add(genusSpinner, c2);
		shrinkPanel.add(protectorPane, c2);
		shrinkPanel.add(generateButton, c2);

		generateButton.addActionListener(this);
		genusSpinner.addChangeListener(this);
	}
	
	@Override
	public void actionPerformed(ActionEvent e) {
		Set<SchottkyPair> pairs = getSchottkyPairs();
		CoHDS hds = generate(hif.getAdapters(), pairs);
		hif.set(hds);
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		resetViewer();
	}
	
	
	private static class Circle {
		private boolean 
			orientation = true;
		private Complex 
			c = new Complex();
		private double 
			r = 1.0;
		
		public Circle(Complex c, double r, boolean o) {
			super();
			this.c = c;
			this.r = r;
			this.orientation = o;
		}
		
		public boolean isInside(Complex z, double tol) {
			Complex d = z.minus(c);
			if (orientation) {
				return d.abs() < r + tol;
			} else {
				return d.abs() > r - tol;
			}
		}
		
		public Complex getOrientedCenter() {
			if (orientation) {
				return c;
			} else { // invert
				// TODO fix this to be general
				return c.invert();
			}
		}
		
		@Override
		public String toString() {
			return "Circle: " + c + " r=" + r;
		}
	}
	
	
	private static class SchottkyPair {
		private Moebius 
			s = null;
		private Circle 
			source = null,
			target = null;
		public SchottkyPair(Moebius s, Circle source, Circle target) {
			super();
			this.s = s;
			this.source = source;
			this.target = target;
		}
	}
	
	
	private static double[] inverseStereographic(Complex Z) {
		if (Z.re == Double.POSITIVE_INFINITY) {
			return new double[] {0,0,1};
		}
		double X = Z.getRe();
		double Y = Z.getIm();
		double d = 1 + X*X + Y*Y;
		double x = 2*X / d;
		double y = 2*Y / d;
		double z = (X*X + Y*Y - 1) / d;
		return new double[] {x, y, z};
	}
	
	
	private static Complex stereographic(double[] p) {
		return new Complex(p[0] / (1 - p[2]), p[1] / (1 - p[2]));
	}
	
	
	private static Set<Circle> getAllCircles(Set<SchottkyPair> pairs) {
		Set<Circle> r = new HashSet<DiscreteSchottkyGenerator.Circle>();
		for (SchottkyPair p : pairs) {
			r.add(p.source);
			r.add(p.target);
		}
		return r;
	}
	
	
	private static Set<SchottkyPair> getSchottkyPairs() {
		Set<SchottkyPair> pairs = new HashSet<SchottkyPair>();
		// define transformation 1
		Complex A = new Complex(-0.2, 0.0);
		Complex B = new Complex(0.0, 0.0);
		Complex m = new Complex(0.3);
		Moebius s = new Moebius(A, B, m);
		Complex center = new Complex();
		Circle circle = new Circle(center, 1.0, false);
		Complex centerOfMappedCircle = new Complex();
		double sr = s.getRadiusOfMappedCircle(circle.c, circle.r, centerOfMappedCircle);
		Circle sCircle = new Circle(centerOfMappedCircle, sr, true);
		SchottkyPair pair = new SchottkyPair(s, circle, sCircle);
//		pairs.add(pair);
		
		// define transformation 2
		A = new Complex(0.5, 0.0);
		B = new Complex(-0.5, 0.0);
		m = new Complex(0.0, 0.2);
		s = new Moebius(A, B, m);
		center = new Complex(0.4, 0.0);
		circle = new Circle(center, 0.05, true);
		centerOfMappedCircle = new Complex();
		sr = s.getRadiusOfMappedCircle(circle.c, circle.r, centerOfMappedCircle);
		sCircle = new Circle(centerOfMappedCircle, sr, true);
		pair = new SchottkyPair(s, circle, sCircle);
		pairs.add(pair);
		return pairs;
	}
	
	
	
	private static CoHDS generate(AdapterSet a, Set<SchottkyPair> pairs) {
		CoHDS hds = new CoHDS();
		double zScale = 10.0; 
		int extraPointRes = 10;
		int circleRes = 10;
		
		Map<CoVertex, CoVertex> sMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, CoVertex> sInvMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, CoVertex> centerMap = new HashMap<CoVertex, CoVertex>();
		Map<CoFace, CoFace> sfMap = new HashMap<CoFace, CoFace>();
		Map<CoVertex, SchottkyPair> vertexPairMap = new HashMap<CoVertex, SchottkyPair>(); 
		Map<CoVertex, SchottkyPair> vertexPairInvMap = new HashMap<CoVertex, SchottkyPair>(); 
		Map<CoEdge, CoEdge> edgeMap = new HashMap<CoEdge, CoEdge>();
		Map<CoEdge, SchottkyPair> edgePairMap = new HashMap<CoEdge, SchottkyPair>();
		Map<CoEdge, SchottkyPair> edgePairInvMap = new HashMap<CoEdge, SchottkyPair>();
//		
//		// create extra grid points
//		Set<Complex> extraPoints = new HashSet<Complex>(); 
//		for (int i = 0; i < extraPointRes + 1; i++) {
//			for (int j = 0; j < extraPointRes + 1; j++) {
//				Complex z = new Complex(i / (0.5*extraPointRes) - 1, j / (0.5*extraPointRes) - 1);
//				extraPoints.add(z);
//			}
//		}
//		
//		// create circle density points
//		for (Circle c : getAllCircles(pairs)) {
//			for (int i = 1; i < 2; i++) {
//				int ringRes = circleRes * 3 / (i + 4);
//				for (int j = 0; j < ringRes; j++) {
//					double o = c.orientation ? 1.0 : -1.0;
//					double phi = 2*j*PI / ringRes;
//					double x = (c.r + o*c.r*i/10.0) * cos(phi) + c.c.re;
//					double y = (c.r + o*c.r*i/10.0) * sin(phi) + c.c.im;
//					Complex z = new Complex(x, y);
//					extraPoints.add(z);
//				}
//			}
//		}
//		
//		// exclude extra vertices from the circles
//		for (Complex ez : new HashSet<Complex>(extraPoints)) {
//			for (Circle c : getAllCircles(pairs)) {
//				if (c.isInside(ez, c.r*0.05)) {
//					extraPoints.remove(ez);
//				}
//			}
//		}
		
		// add the vertices on the source and target circles
		for (SchottkyPair p : pairs) {
			for (int i = 0; i < circleRes; i++) {
				Circle c = p.source;
				double phi = 2*i*PI / circleRes;
				double x = c.r * cos(phi) + c.c.re;
				double y = c.r * sin(phi) + c.c.im;
				Complex z = new Complex(x, y);
				Complex sz = p.s.applyTo(z);
				double[] zPos = new double[] {z.re, z.im, 0};
				double[] szPos = new double[] {sz.re, sz.im, 0};
				CoVertex v = hds.addNewVertex();
				CoVertex sv = hds.addNewVertex();
				a.set(Position.class, v, zPos);
				a.set(Position.class, sv, szPos);
				sMap.put(v, sv);
				sInvMap.put(sv, v);
				vertexPairMap.put(v, p);
				vertexPairInvMap.put(sv, p);
			}
			CoVertex center = hds.addNewVertex();
			CoVertex sCenter = hds.addNewVertex();
			double[] cPos = new double[] {p.source.c.re, p.source.c.im, 0};
			Complex mappedCenter = p.s.applyTo(p.source.c);
			double[] cPosMapped = new double[] {mappedCenter.re, mappedCenter.im, 0};
			a.set(Position.class, center, cPos);
			a.set(Position.class, sCenter, cPosMapped);
		}
//		// add the projected circle centers to act as handles
//		for (SchottkyPair p : pairs) {
//			CoVertex cv = hds.addNewVertex();
//			CoVertex scv = hds.addNewVertex();
//			centerMap.put(cv, scv);
//			Complex z = p.source.getOrientedCenter();
//			Complex sz = p.target.getOrientedCenter();
//			double[] zPos = new double[] {z.re, z.im, 0};
//			double[] szPos = new double[] {sz.re, sz.im, 0};
//			a.set(Position.class, cv, zPos);
//			a.set(Position.class, scv, szPos);
//		}
//		// add extra point vertices
//		for (Complex z : extraPoints) {
//			double[] zPos = new double[] {z.re, z.im, 0};
//			CoVertex v = hds.addNewVertex();
//			a.set(Position.class, v, zPos);	
//		}
//		
//		// scale and project 
//		for (CoVertex v : hds.getVertices()) {
//			double[] p = a.getD(Position3d.class, v);
//			Complex z = new Complex(p[0], p[1]);
//			z = z.times(zScale);
//			double[] pProjected = inverseStereographic(z);
//			a.set(Position.class, v, pProjected);
//		}
//		
//		// create triangulation
//		ConvexHull.convexHull(hds, a, 1E-8);
//		
//		// build edge identification opposites maps
//		for (CoVertex vt : sMap.keySet()) {
//			CoVertex svt = sMap.get(vt);
//			for (CoEdge e : incomingEdges(vt)) {
//				CoVertex vs = e.getStartVertex();
//				if (sMap.containsKey(vs)) {
//					CoVertex svs = sMap.get(vs);
//					for (CoEdge se : incomingEdges(svs)) {
//						if (se.getStartVertex() == svt) {
//							SchottkyPair p = vertexPairMap.get(vt);
//							edgeMap.put(e, se);
//							edgeMap.put(e.getOppositeEdge(), se.getOppositeEdge());
//							edgeMap.put(se, e);
//							edgeMap.put(se.getOppositeEdge(), e.getOppositeEdge());
//							edgePairMap.put(e, p);
//							edgePairMap.put(e.getOppositeEdge(), p);
//							edgePairInvMap.put(se, p);
//							edgePairInvMap.put(se.getOppositeEdge(), p);
//						}
//					}
//				}
//			}
//		}
//		// we know how many edges there are on the circles
//		assert edgeMap.size() == 2 * pairs.size() * circleRes * 2;
//		
//		// calculate length cross ratios
//		Map<CoEdge, Double> crMap = new HashMap<CoEdge, Double>();
//		for (CoEdge e : hds.getEdges()) {
//			CoVertex vi = e.getStartVertex();
//			CoVertex vj = e.getTargetVertex();
//			CoVertex vl = e.getNextEdge().getTargetVertex();
//			CoVertex vk = e.getOppositeEdge().getNextEdge().getTargetVertex();
//			Complex zi = stereographic(a.getD(Position3d.class, vi));
//			Complex zj = stereographic(a.getD(Position3d.class, vj));
//			Complex zl = stereographic(a.getD(Position3d.class, vl));
//			Complex zk = stereographic(a.getD(Position3d.class, vk));
//			
//			Moebius pullTransform = null;
//			boolean invertCrossRatio = false;
//			// check if we have to pull a vertex position
//			if (edgePairMap.containsKey(e)) {
//				assert centerMap.containsKey(vl) || centerMap.containsKey(vk);
//				// we are at a source circle edge
//				SchottkyPair pair = edgePairMap.get(e);
//				pullTransform = pair.s.invert();
//			} 
//			if (edgePairInvMap.containsKey(e)) {
//				assert centerMap.containsValue(vl) || centerMap.containsValue(vk);
//				// we are at a target circle edge
//				SchottkyPair pair = edgePairInvMap.get(e);
//				pullTransform = pair.s;
//				invertCrossRatio = true;
//			}
//			// pull one vertex from a different domain
//			if (pullTransform != null) {
//				CoEdge se = edgeMap.get(e);
//				CoVertex v1 = se.getOppositeEdge().getNextEdge().getTargetVertex();
//				CoVertex v2 = se.getNextEdge().getTargetVertex();
//				if (centerMap.containsKey(vl) || centerMap.containsValue(vl)) {
//					CoVertex mappedL = centerMap.containsKey(v1) || centerMap.containsValue(v1) ? v2 : v1;
//					assert !centerMap.containsKey(mappedL) && !centerMap.containsValue(mappedL);
//					zl = stereographic(a.getD(Position3d.class, mappedL));
//					zl = zl.times(1 / zScale);
//					zl = pullTransform.applyTo(zl);
//					zl = zl.times(zScale);
//				} else {
//					assert centerMap.containsKey(vk) || centerMap.containsValue(vk);
//					CoVertex mappedL = centerMap.containsKey(v1) || centerMap.containsValue(v1) ? v2 : v1;
//					assert !centerMap.containsKey(mappedL) && !centerMap.containsValue(mappedL);
//					zk = stereographic(a.getD(Position3d.class, mappedL));
//					zk = zk.times(1 / zScale);
//					zk = pullTransform.applyTo(zk);
//					zk = zk.times(zScale);
//				}
//			}
//			
//			double l_ik = zi.dist(zk); 
//			double l_kj = zk.dist(zj); 
//			double l_jl = zj.dist(zl); 
//			double l_li = zl.dist(zi); 
//			double cr = (l_kj * l_li) / (l_ik * l_jl);
//			if (invertCrossRatio) {
//				cr = 1 / cr;
//			}
//			// store the cross ratio
//			crMap.put(e, cr);
//		}
//
//		// remove circles from the triangulation
//		for (CoVertex v : centerMap.keySet()) {
//			CoVertex sv = centerMap.get(v);
//			CoFace f = TopologyAlgorithms.removeVertexFill(v);
//			CoFace sf = TopologyAlgorithms.removeVertexFill(sv);
//			sfMap.put(f, sf);
//		}
//		
//		// identify circles
//		try {
//			for (CoFace f : sfMap.keySet()) {
//				CoFace sf = sfMap.get(f);
//				CoVertex v = f.getBoundaryEdge().getTargetVertex();
//				CoVertex sv = sMap.get(v);
//				try {
//					SurgeryUtility.glueAlongFaces(f, sf, v, sv);
//				} catch (RuntimeException e) {
//					e.printStackTrace();
//				}
//			}
//		} catch (Exception e) {
//			System.err.println("Invalid Schottky configuration.");
//		}
//		
//		// check if identification was correct
//		for (CoEdge e : edgeMap.keySet()) {
//			if (!e.isValid()) continue;
//			CoEdge se = edgeMap.get(e);
//			assert e == se.getOppositeEdge() : "opposite";
//		}
//		
//		// assert mapped and opposite cross ratios are the same
//		for (CoEdge e : edgeMap.keySet()) {
//			if (!e.isValid()) continue;
//			CoEdge se = edgeMap.get(e);
//			double cr = crMap.get(e);
//			double scr = crMap.get(se);
//			assert cr - scr < 1E-8 : "Mapped cross ratio: " + cr + " - " + scr;
//		}
//		int count = 0;
//		for (CoEdge e : hds.getEdges()) {
//			if (!e.isValid()) continue;
//			CoEdge opp = e.getOppositeEdge();
//			double cr = crMap.get(e);
//			double scr = crMap.get(opp);
////			assert cr - scr < 1E-8 : "Opposite cross ratios: " + cr + " - " + scr;
//			if (cr - scr > 1E-8) count++;
//		}
//		System.out.println("wring opposite cross ratio counter " + count);
//		
		
		System.out.println("Generated surface of genus " + HalfEdgeUtils.getGenus(hds));
		return hds;
	}
	
	
	
	private void resetViewer() {
	}

	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
		hif.addGlobalAdapter(new CoPositionAdapter(), true);
		hif.addGlobalAdapter(new CoTexturePositionAdapter(true), true);
		resetViewer();
	}
	
	
	@Override
	public Class<? extends SideContainerPerspective> getPerspectivePluginClass() {
		return View.class;
	}

	
	
	public static void main(String[] args) {
		AdapterSet a = AdapterSet.createGenericAdapters();
		a.add(new CoPositionAdapter());
		a.add(new CoTexturePositionAdapter(true));
		Set<SchottkyPair> pairs = DiscreteSchottkyGenerator.getSchottkyPairs();
		DiscreteSchottkyGenerator.generate(a, pairs);
		
//		NativePathUtility.set("native");
//		JRViewer v = new JRViewer();
//		v.addContentUI();
//		v.addBasicUI();
//		v.registerPlugin(DiscreteSchottkyGenerator.class);
//		v.registerPlugin(new VertexRemoverPlugin());
//		v.registerPlugin(new RandomSphereGenerator());
//		v.startup();
	}
	
	
	
}
