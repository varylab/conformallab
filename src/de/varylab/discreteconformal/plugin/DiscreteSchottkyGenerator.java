package de.varylab.discreteconformal.plugin;

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
import de.jtem.java2dx.beans.Viewer2DWithInspector;
import de.jtem.java2dx.modelling.GraphicsModeller2D;
import de.jtem.java2dx.modelling.SimpleModeller2D;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
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
		generate();
	}
	
	@Override
	public void stateChanged(ChangeEvent e) {
		resetViewer();
	}
	
	
	private class Circle {
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
	
	
	private class SchottkyPair {
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
	
	
	private double[] inverseStereographic(Complex Z) {
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
	
	
//	private Complex stereographic(double[] p) {
//		return new Complex(p[0] / (1 - p[2]), p[1] / (1 - p[2]));
//	}
	
	
	private Set<Circle> getAllCircles(Set<SchottkyPair> pairs) {
		Set<Circle> r = new HashSet<DiscreteSchottkyGenerator.Circle>();
		for (SchottkyPair p : pairs) {
			r.add(p.source);
			r.add(p.target);
		}
		return r;
	}
	
	
	private void generate() {
		AdapterSet a = hif.getAdapters();
		CoHDS hds = new CoHDS();
		double zScale = 10.0; 
		int extraPointRes = 20;
		int circleRes = 50;
		
		Map<CoVertex, CoVertex> sMap = new HashMap<CoVertex, CoVertex>();
		Map<CoVertex, CoVertex> centerMap = new HashMap<CoVertex, CoVertex>();
		Map<CoFace, CoFace> sfMap = new HashMap<CoFace, CoFace>();
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
		pairs.add(pair);
		
		// define transformation 2
		A = new Complex(0.5, 0.0);
		B = new Complex(-0.5, 0.0);
		m = new Complex(1.1, 1.2);
		s = new Moebius(A, B, m);
		center = new Complex(-0.1, 0.0);
		circle = new Circle(center, 0.05, true);
		centerOfMappedCircle = new Complex();
		sr = s.getRadiusOfMappedCircle(circle.c, circle.r, centerOfMappedCircle);
		sCircle = new Circle(centerOfMappedCircle, sr, true);
		pair = new SchottkyPair(s, circle, sCircle);
		pairs.add(pair);
		
		
		// create extra grid points
		Set<Complex> extraPoints = new HashSet<Complex>(); 
		for (int i = 0; i < extraPointRes + 1; i++) {
			for (int j = 0; j < extraPointRes + 1; j++) {
				Complex z = new Complex(i / (0.5*extraPointRes) - 1, j / (0.5*extraPointRes) - 1);
				extraPoints.add(z);
			}
		}
		
		// create circle density points
		for (Circle c : getAllCircles(pairs)) {
			for (int i = 1; i < 10; i++) {
				int ringRes = circleRes * 3 / (i + 4);
				for (int j = 0; j < ringRes; j++) {
					double o = c.orientation ? 1.0 : -1.0;
					double phi = 2*j*PI / ringRes;
					double x = (c.r + o*c.r*i/10.0) * cos(phi) + c.c.re;
					double y = (c.r + o*c.r*i/10.0) * sin(phi) + c.c.im;
					Complex z = new Complex(x, y);
					extraPoints.add(z);
				}
			}
		}
		
		// exclude extra vertices from the circles
		for (Complex ez : new HashSet<Complex>(extraPoints)) {
			for (Circle c : getAllCircles(pairs)) {
				if (c.isInside(ez, c.r*0.05)) {
					extraPoints.remove(ez);
				}
			}
		}
		
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
			}
		}
		// add the projected circle centers to act as handles
		for (SchottkyPair p : pairs) {
			CoVertex cv = hds.addNewVertex();
			CoVertex scv = hds.addNewVertex();
			centerMap.put(cv, scv);
			Complex z = p.source.getOrientedCenter();
			Complex sz = p.target.getOrientedCenter();
			double[] zPos = new double[] {z.re, z.im, 0};
			double[] szPos = new double[] {sz.re, sz.im, 0};
			a.set(Position.class, cv, zPos);
			a.set(Position.class, scv, szPos);
		}
		// add extra point vertices
		for (Complex z : extraPoints) {
			double[] zPos = new double[] {z.re, z.im, 0};
			CoVertex v = hds.addNewVertex();
			a.set(Position.class, v, zPos);	
		}
		
		// scale and project 
		for (CoVertex v : hds.getVertices()) {
			double[] p = a.getD(Position3d.class, v);
			Complex z = new Complex(p[0], p[1]);
			z = z.times(zScale);
			double[] pProjected = inverseStereographic(z);
			a.set(Position.class, v, pProjected);
		}
		
		// create triangulation
		ConvexHull.convexHull(hds, a, 1E-8);
		
		// remove circles from the triangulation
		for (CoVertex v : centerMap.keySet()) {
			CoVertex sv = centerMap.get(v);
			CoFace f = TopologyAlgorithms.removeVertexFill(v);
			CoFace sf = TopologyAlgorithms.removeVertexFill(sv);
			sfMap.put(f, sf);
		}
		
		// identify circles
		try {
			for (CoFace f : sfMap.keySet()) {
				CoFace sf = sfMap.get(f);
				CoVertex v = f.getBoundaryEdge().getTargetVertex();
				CoVertex sv = sMap.get(v);
				try {
					SurgeryUtility.glueAlongFaces(f, sf, v, sv);
				} catch (RuntimeException e) {
					e.printStackTrace();
				}
			}
		} catch (Exception e) {
			System.err.println("Invalid Schottky configuration.");
		}
		
		System.out.println("Generated surface of genus " + HalfEdgeUtils.getGenus(hds));
		hif.set(hds);
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
		NativePathUtility.set("native");
		JRViewer v = new JRViewer();
		v.addContentUI();
		v.addBasicUI();
		DiscreteSchottkyGenerator dsg = new DiscreteSchottkyGenerator();
		v.registerPlugin(dsg);
		v.registerPlugin(new VertexRemoverPlugin());
		v.registerPlugin(new RandomSphereGenerator());
		v.startup();
	}
	
	
	
}
