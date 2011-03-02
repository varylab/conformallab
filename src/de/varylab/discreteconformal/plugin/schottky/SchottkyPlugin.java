package de.varylab.discreteconformal.plugin.schottky;

import static de.jtem.halfedge.util.HalfEdgeUtils.incomingEdges;
import static de.jtem.java2d.Annotation.SOUTHEAST;
import static java.awt.geom.AffineTransform.getTranslateInstance;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import java.awt.BasicStroke;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import de.jreality.math.Rn;
import de.jreality.plugin.JRViewer;
import de.jreality.plugin.basic.View;
import de.jreality.ui.LayoutFactory;
import de.jreality.util.NativePathUtility;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.algorithm.topology.TopologyAlgorithms;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.generator.RandomSphereGenerator;
import de.jtem.halfedgetools.plugin.algorithm.topology.VertexRemoverPlugin;
import de.jtem.java2d.Annotation;
import de.jtem.java2d.DragListener;
import de.jtem.java2d.SceneComponent;
import de.jtem.java2d.TransformedMouseEvent;
import de.jtem.java2dx.Ellipse2DDouble;
import de.jtem.java2dx.beans.Viewer2DWithInspector;
import de.jtem.java2dx.modelling.GraphicsModeller2D;
import de.jtem.java2dx.modelling.SimpleModeller2D;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import de.jtem.numericalMethods.geometry.meshGeneration.ruppert.Ruppert;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;
import de.varylab.discreteconformal.util.NodeIndexComparator;
import de.varylab.discreteconformal.util.SurgeryUtility;

public class SchottkyPlugin extends ShrinkPanelPlugin implements ActionListener, ChangeListener {

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
	private static SpinnerNumberModel
		genusModel = new SpinnerNumberModel(2, 0, 20, 1),
		stereographicScaleModel = new SpinnerNumberModel(7.0, 0.01, 100.0, 0.01),
		cirleResModel = new SpinnerNumberModel(20, 3, 1000, 1),
		delaunayAreaModel = new SpinnerNumberModel(0.05, 0.001, 10, 0.01);
	private JSpinner
		genusSpinner = new JSpinner(genusModel),
		delaunayAreaSpinner = new JSpinner(delaunayAreaModel),
		circleResSpinner = new JSpinner(cirleResModel),
		stereographicScaleSpinner = new JSpinner(stereographicScaleModel);
	private SceneComponent
		root = new SceneComponent();
	
	private static double 
		zScale = 7.0; 
	
	
	private class FixPointComponent extends SceneComponent {
		
		public FixPointComponent(double x, double y, int index, boolean isSource) {
			setTransform(AffineTransform.getTranslateInstance(x, y));
			setShape(new Ellipse2DDouble(-0.1, -0.1, 0.2, 0.2));
			setStroke(new BasicStroke(2));
			setOutlineDragEnabled(true);
			setAreaDragEnabled(true);
			setAnnotated(true);
			setAnnotationDragEnabled(false);
			String label = (isSource ? "A" : "B") + index;
			Annotation ann = new Annotation(label, 0, 0, SOUTHEAST);
			getAnnotations().add(ann);
			addDragListener(new DragListener() {
				double lastX, lastY;
				@Override
				public void dragStart(TransformedMouseEvent e) {
					lastX = e.getX();
					lastY = e.getY();
				}
				
				@Override
				public void dragEnd(TransformedMouseEvent e) {
				}
				
				@Override
				public void drag(TransformedMouseEvent e) {
					AffineTransform T = getTransform();
					double dx = e.getX() - lastX;
					double dy = e.getY() - lastY;
					T.concatenate(getTranslateInstance(dx, dy));
					setTransform(T);
					lastX = e.getX();
					lastY = e.getY();
					viewer.repaint();
				}
			});
		}

		@Override
		public AffineTransform getTransform() {
			if (super.getTransform() == null) {
				setTransform(new AffineTransform());
			}
			return super.getTransform();
		}
		
	}

	
	
	private class CircleComponent extends SceneComponent {
		
		public CircleComponent() {
			
		}
		
	}
	
		
	public SchottkyPlugin() {
		GridBagConstraints c1 = LayoutFactory.createLeftConstraint();
		GridBagConstraints c2 = LayoutFactory.createRightConstraint();
		viewer.setPreferredSize(new Dimension(200, 200));
		viewer.setMinimumSize(viewer.getPreferredSize());
		shrinkPanel.setLayout(new GridBagLayout());
		shrinkPanel.add(new JLabel("Genus"), c1);
		shrinkPanel.add(genusSpinner, c2);
		shrinkPanel.add(new JLabel("Max Triangle Area"), c1);
		shrinkPanel.add(delaunayAreaSpinner, c2);
		shrinkPanel.add(new JLabel("Stereographic Factor"), c1);
		shrinkPanel.add(stereographicScaleSpinner, c2);
		shrinkPanel.add(new JLabel("Circle Resolution"), c1);
		shrinkPanel.add(circleResSpinner, c2);
		shrinkPanel.add(protectorPane, c2);
		shrinkPanel.add(generateButton, c2);

		generateButton.addActionListener(this);
		genusSpinner.addChangeListener(this);
		
		moddeller.addScene(root);
		viewer.setScaleToolEnabled(true);
		viewer.setTranslateToolEnabled(true);
		viewer.setMenuToolEnabled(false);
	}
	
	
	@Override
	public void actionPerformed(ActionEvent e) {
		List<SchottkyGenerator> pairs = getSchottkyPairs();
		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
		CoHDS hds = generate(hif.getAdapters(), pairs, lMap);
		hif.set(hds);
		hif.addLayerAdapter(new SchottkyLengthAdapter(lMap), false);
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		resetViewer();
	}
	
	
	private static double[] inverseStereographic(Complex Z, double[]... result) {
		if (Z.re == Double.POSITIVE_INFINITY) {
			return new double[] {0,0,1};
		}
		double X = Z.getRe();
		double Y = Z.getIm();
		double d = 1 + X*X + Y*Y;
		double x = 2*X / d;
		double y = 2*Y / d;
		double z = (X*X + Y*Y - 1) / d;
		if (result.length != 0) {
			result[0][0] = x;
			result[0][1] = y;
			result[0][2] = z;
			return result[0];
		}
		return new double[] {x, y, z};
	}
	
	
	private static Complex stereographic(double[] p) {
		return new Complex(p[0] / (1 - p[2]), p[1] / (1 - p[2]));
	}
	
	
	private static List<SchottkyCircle> getAllCircles(List<SchottkyGenerator> pairs) {
		List<SchottkyCircle> r = new LinkedList<SchottkyCircle>();
		for (SchottkyGenerator p : pairs) {
			r.add(p.getSource());
			r.add(p.getTarget());
		}
		return r;
	}
	
	
	private static List<SchottkyGenerator> getSchottkyPairs() {
		List<SchottkyGenerator> pairs = new LinkedList<SchottkyGenerator>();
		// define transformation 1
		Complex A = new Complex(-1.5, -1.5);
		Complex B = new Complex(0.0, 0.0);
		Complex m = new Complex(0.01);
		Moebius s = new Moebius(A, B, m);
		Complex center = new Complex();
		SchottkyCircle circle = new SchottkyCircle(center, 1.0, false);
		Complex centerOfMappedCircle = new Complex();
		double sr = s.getRadiusOfMappedCircle(circle.getCenter(), circle.getRadius(), centerOfMappedCircle);
		SchottkyCircle sCircle = new SchottkyCircle(centerOfMappedCircle, sr, true);
		SchottkyGenerator pair = new SchottkyGenerator(s, circle, sCircle, "0");
		pairs.add(pair);
		
		// define transformation 2
		A = new Complex(0.23, 0.1);
		B = new Complex(-0.23, 0.1);
		m = new Complex(0.005);
		s = new Moebius(A, B, m);
		circle = new SchottkyCircle(A, 0.05, true);
		centerOfMappedCircle = new Complex();
		sr = s.getRadiusOfMappedCircle(circle.getCenter(), circle.getRadius(), centerOfMappedCircle);
		sCircle = new SchottkyCircle(centerOfMappedCircle, sr, true);
		pair = new SchottkyGenerator(s, circle, sCircle, "1");
		pairs.add(pair);
		
		// define transformation 3
		A = new Complex(0.1, 0.1);
		B = new Complex(0.1, -0.1);
		m = new Complex(0.01);
		s = new Moebius(A, B, m);
		circle = new SchottkyCircle(A, 0.01, true);
		centerOfMappedCircle = new Complex();
		sr = s.getRadiusOfMappedCircle(circle.getCenter(), circle.getRadius(), centerOfMappedCircle);
		sCircle = new SchottkyCircle(centerOfMappedCircle, sr, true);
		pair = new SchottkyGenerator(s, circle, sCircle, "2");
		pairs.add(pair);
		
		return pairs;
	}
	
	
	private static class StereographicRuppert extends Ruppert {

		private static final long serialVersionUID = 1L;
		private boolean inited = false;
		private double[] v1, v2, v3, vp1, vp2, vp3, vec1, vec2, cross;
		
		public StereographicRuppert(double[][] p) {
			super(p);
		}
		
		@Override
		protected double area(int f) {
			if (!inited) {
				v1 = new double[2];
				v2 = new double[2];
				v3 = new double[3];
				vp1 = new double[3];
				vp2 = new double[3];
				vp3 = new double[3];
				vec1 = new double[3];
				vec2 = new double[3];
				cross = new double[3];
				inited = true;
			}
			int v1Index = getIndex(f, 1);
			int v2Index = getIndex(f, 2);
			int v3Index = getIndex(f, 3);
			getPoint(v1Index, v1);
			getPoint(v2Index, v2);
			getPoint(v3Index, v3);
			Complex z1 = new Complex(v1[0], v1[1]);
			Complex z2 = new Complex(v2[0], v2[1]);
			Complex z3 = new Complex(v3[0], v3[1]);
			z1 = z1.times(zScale);
			z2 = z2.times(zScale);
			z3 = z3.times(zScale);
			inverseStereographic(z1, vp1);
			inverseStereographic(z2, vp2);
			inverseStereographic(z3, vp3);
			Rn.subtract(vec1, vp2, vp1);
			Rn.subtract(vec2, vp3, vp1);
			Rn.crossProduct(cross, vec1, vec2);
			return Rn.euclideanNorm(cross);
		}
		
	}
	
	
	
	private static void cutIdentificationHoles(CoHDS hds, List<ArrayList<CoVertex>> circles) {
		for (ArrayList<CoVertex> circle : circles) {
			int len = circle.size();
			for (int i = 0; i < len; i++) {
				CoVertex pre = circle.get((i + len - 1) % len);
				CoVertex act = circle.get((i + len + 0) % len);
				CoVertex pos = circle.get((i + len + 1) % len);
				List<CoEdge> star = incomingEdges(act);
				for (CoEdge e : star) {
					if (!e.isValid()) continue;
					if (e.getStartVertex() == pre) continue;
					if (e.getStartVertex() == pos) continue;
					if (circle.contains(e.getStartVertex())) {
						TopologyAlgorithms.removeEdge(e);
					}
				}
			}
		}
	}
	
	
	
	
	private static CoHDS generate(AdapterSet a, List<SchottkyGenerator> pairs, Map<CoEdge, Double> lMap) {
		zScale = stereographicScaleModel.getNumber().doubleValue();
		int circleRes = cirleResModel.getNumber().intValue();
		double ruppertArea = delaunayAreaModel.getNumber().doubleValue();
		
		CoHDS hds = new CoHDS();
		Map<CoVertex, CoVertex> sMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, CoVertex> sInvMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, SchottkyGenerator> vertexPairMap = new HashMap<CoVertex, SchottkyGenerator>(); 
		Map<CoVertex, SchottkyGenerator> vertexPairInvMap = new HashMap<CoVertex, SchottkyGenerator>(); 
		Map<CoEdge, CoEdge> edgeMap = new HashMap<CoEdge, CoEdge>();
		Map<CoEdge, SchottkyGenerator> edgePairMap = new HashMap<CoEdge, SchottkyGenerator>();
		Map<CoEdge, SchottkyGenerator> edgePairInvMap = new HashMap<CoEdge, SchottkyGenerator>();
		
		double[][] polygons = new double[pairs.size() * 2][circleRes * 2];
		List<ArrayList<CoVertex>> vertexCircles = new LinkedList<ArrayList<CoVertex>>();
		
		// add the vertices on the source and target circles
		int polygonIndex = 0;
		for (SchottkyGenerator p : pairs) {
			ArrayList<CoVertex> sourceCircle = new ArrayList<CoVertex>();
			ArrayList<CoVertex> targetCircle = new ArrayList<CoVertex>();
			vertexCircles.add(sourceCircle);
			vertexCircles.add(targetCircle);
			for (int i = 0; i < circleRes; i++) {
				SchottkyCircle c = p.getSource();
				double phi = 2*i*PI / circleRes;
				double x = c.getRadius() * cos(phi) + c.getCenter().re;
				double y = c.getRadius() * sin(phi) + c.getCenter().im;
				Complex z = new Complex(x, y);
				Complex sz = p.getS().applyTo(z);
				double[] zPos = new double[] {z.re, z.im, 0};
				double[] szPos = new double[] {sz.re, sz.im, 0};
				CoVertex v = hds.addNewVertex();
				CoVertex sv = hds.addNewVertex();
				sourceCircle.add(v);
				targetCircle.add(sv);
				a.set(Position.class, v, zPos);
				a.set(Position.class, sv, szPos);
				sMap.put(v, sv);
				sInvMap.put(sv, v);
				vertexPairMap.put(v, p);
				vertexPairInvMap.put(sv, p);
				polygons[polygonIndex][i * 2] = z.re;
				polygons[polygonIndex][i * 2 + 1] = z.im;
				polygons[polygonIndex + 1][i * 2] = sz.re;
				polygons[polygonIndex + 1][i * 2 + 1] = sz.im;
			}
			polygonIndex += 2;
		}
//		if (polygonIndex == 4) return hds;
		Ruppert ruppert = new StereographicRuppert(polygons);
		ruppert.setAreaConstraint(ruppertArea);
		ruppert.setMaximalNumberOfTriangles(100000);
		ruppert.refine();
		double[] verts = ruppert.getPoints();
		
		for (int i = 0; i < verts.length; i+=2) {
			Complex z = new Complex(verts[i], verts[i+1]);
			double[] zPos = new double[] {z.re, z.im, 0};
			CoVertex v = hds.addNewVertex();
			a.set(Position.class, v, zPos);
		}
		
		
		// exclude extra vertices from the circles
		for (CoVertex v : new HashSet<CoVertex>(hds.getVertices())) {
			if (sMap.containsKey(v) || sMap.containsValue(v)) {
				continue; // we are a circle vertex
			}
			for (SchottkyCircle c : getAllCircles(pairs)) {
				double[] vp = a.getD(Position3d.class, v);
				Complex vz = new Complex(vp[0], vp[1]);
				if (c.isInside(vz, c.getRadius() * 1E-4)) {
					hds.removeVertex(v);
				}
			}
		}
		
		// scale and project 
		Map<CoVertex, Complex> zMap = new HashMap<CoVertex, Complex>();
		for (CoVertex v : hds.getVertices()) {
			double[] p = a.getD(Position3d.class, v);
			Complex z = new Complex(p[0], p[1]);
			zMap.put(v, z); // store original positions
			z = z.times(zScale);
			double[] pProjected = inverseStereographic(z);
			a.set(Position.class, v, pProjected);
		}
		
		// create triangulation
		ConvexHull.convexHull(hds, a, 1E-8);
		cutIdentificationHoles(hds, vertexCircles);

		
		// build edge identification opposites maps
		for (CoVertex vt : sMap.keySet()) {
			CoVertex svt = sMap.get(vt);
			for (CoEdge e : incomingEdges(vt)) {
				CoVertex vs = e.getStartVertex();
				if (sMap.containsKey(vs)) {
					CoVertex svs = sMap.get(vs);
					for (CoEdge se : incomingEdges(svs)) {
						if (se.getStartVertex() == svt) {
							SchottkyGenerator p = vertexPairMap.get(vt);
							edgeMap.put(e, se);
							edgeMap.put(e.getOppositeEdge(), se.getOppositeEdge());
							edgeMap.put(se, e);
							edgeMap.put(se.getOppositeEdge(), e.getOppositeEdge());
							edgePairMap.put(e, p);
							edgePairMap.put(e.getOppositeEdge(), p);
							edgePairInvMap.put(se, p);
							edgePairInvMap.put(se.getOppositeEdge(), p);
						}
					}
				}
			}
		}
		// we know how many edges there are on the circles
		assert edgeMap.size() == 2 * pairs.size() * circleRes * 2 : "lengths of cycles";
		
		// calculate length cross ratios
		Map<CoEdge, Double> crMap = new HashMap<CoEdge, Double>();
		for (CoEdge e : hds.getEdges()) {
			CoVertex vi = e.getStartVertex();
			CoVertex vj = e.getTargetVertex();
			CoVertex vl = e.getNextEdge().getTargetVertex();
			CoVertex vk = e.getOppositeEdge().getNextEdge().getTargetVertex();
			Complex zi = stereographic(a.getD(Position3d.class, vi));
			Complex zj = stereographic(a.getD(Position3d.class, vj));
			Complex zl = stereographic(a.getD(Position3d.class, vl));
			Complex zk = stereographic(a.getD(Position3d.class, vk));
			
			Moebius pullTransform = null;
			// check if we have to pull a vertex position
			if (edgePairMap.containsKey(e)) {
				assert e.getLeftFace() == null || e.getRightFace() == null;
				// we are at a source circle edge
				SchottkyGenerator pair = edgePairMap.get(e);
				pullTransform = pair.getS().invert();
			} 
			if (edgePairInvMap.containsKey(e)) {
				assert e.getLeftFace() == null || e.getRightFace() == null;
				// we are at a target circle edge
				SchottkyGenerator pair = edgePairInvMap.get(e);
				pullTransform = pair.getS();
			}
			// pull one vertex from a different domain
			if (pullTransform != null) {
				CoEdge se = edgeMap.get(e);
				if (se.getLeftFace() == null) {
					CoVertex mappedL = se.getOppositeEdge().getNextEdge().getTargetVertex();
					zl = stereographic(a.getD(Position3d.class, mappedL));
					zl = zl.times(1 / zScale);
					zl = pullTransform.applyTo(zl);
					zl = zl.times(zScale);
				} else {
					CoVertex mappedL = se.getNextEdge().getTargetVertex();
					zk = stereographic(a.getD(Position3d.class, mappedL));
					zk = zk.times(1 / zScale);
					zk = pullTransform.applyTo(zk);
					zk = zk.times(zScale);
				}
			}
			double l_ik = zi.dist(zk); 
			double l_kj = zk.dist(zj); 
			double l_jl = zj.dist(zl); 
			double l_li = zl.dist(zi); 
			double cr = (l_kj * l_li) / (l_ik * l_jl);
			// store the cross ratio
			crMap.put(e, cr);
		}

		// identify circles
		try {
			Set<CoEdge> mappedEdgesSet = new TreeSet<CoEdge>(new NodeIndexComparator<CoEdge>());
			mappedEdgesSet.addAll(edgeMap.keySet());
			for (CoEdge e : mappedEdgesSet) {
				if (!e.isValid() || e.getLeftFace() != null) {
					continue;
				}
				CoEdge se = edgeMap.get(e);
				try {
					System.out.println("Identify " + e + " " + se);
					SurgeryUtility.glueAlongBoundaries(e, se);
				} catch (RuntimeException re) {
					re.printStackTrace();
				}
			}
		} catch (Exception e) {
			System.err.println("Invalid Schottky configuration.");
		}
		
		// check if identification was correct
		for (CoEdge e : edgeMap.keySet()) {
			if (!e.isValid()) continue;
			CoEdge se = edgeMap.get(e);
			assert e == se.getOppositeEdge() : "opposite";
		}
		// assert mapped and opposite cross ratios are the same
		for (CoEdge e : edgeMap.keySet()) {
			CoEdge se = edgeMap.get(e);
			double cr = crMap.get(e);
			double scr = crMap.get(se);
			assert abs(cr - scr) < 1E-8 : "Mapped cross ratio: " + cr + " - " + scr;
		}
		for (CoEdge e : hds.getEdges()) {
			CoEdge opp = e.getOppositeEdge();
			double cr = crMap.get(e);
			double scr = crMap.get(opp);
			assert abs(cr - scr) < 1E-8 : "Opposite cross ratios: " + cr + " - " + scr;
		}
		// check cross-ratio product around each vertex
		for (CoVertex v : hds.getVertices()) {
			double pr = 1.0;
			for (CoEdge e : incomingEdges(v)) {
				pr *= crMap.get(e);
			}
			assert abs(pr - 1.0) < 1E-8 : "Cross-ratio product around vertex not 1.0";
		}
		
		
		// define lengths from cross-ratios
		Map<CoEdge, Double> aMap = new HashMap<CoEdge, Double>();
		for (CoVertex v : hds.getVertices()) {
			double ai = 1.0;
			for (CoEdge e : incomingEdges(v)) {
				double q = crMap.get(e);
				ai *= q;
				aMap.put(e, ai);
			}
		}
		for (CoEdge e : hds.getEdges()) {
			double ai = aMap.get(e);
			double aj = aMap.get(e.getPreviousEdge());
			double l = sqrt(1 / (ai * aj));
			lMap.put(e, l);
		}
		
		// check opposite lengths
		for (CoEdge e : hds.getEdges()) {
			double l = lMap.get(e);
			double lo = lMap.get(e.getOppositeEdge());
			assert abs(l - lo) < 1E-8 : "lengths";
		}
		
		System.out.println("Generated surface of genus " + HalfEdgeUtils.getGenus(hds));
		return hds;
	}
	
	
	
	
	
	private void resetViewer() {
		root.removeAllChildren();
		for (SchottkyGenerator g : getSchottkyPairs()) {
			root.addChild(g.getEditor());
		}
		viewer.encompass(viewer.getViewport());
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
//		AdapterSet a = AdapterSet.createGenericAdapters();
//		a.add(new CoPositionAdapter());
//		a.add(new CoTexturePositionAdapter(true));
//		Map<CoEdge, Double> lMap = new HashMap<CoEdge, Double>();
//		List<SchottkyPair> pairs = DiscreteSchottkyGenerator.getSchottkyPairs();
//		DiscreteSchottkyGenerator.generate(a, pairs, lMap);
		
		NativePathUtility.set("native");
		JRViewer v = new JRViewer();
		v.setPropertiesFile("Schottky.jrw");
		v.setPropertiesResource(SchottkyPlugin.class, "Schottky.jrw");
		v.addContentUI();
		v.addBasicUI();
		v.registerPlugin(SchottkyPlugin.class);
		v.registerPlugin(new VertexRemoverPlugin());
		v.registerPlugin(new RandomSphereGenerator());
		v.startup();
	}
	
	
	
}
