package de.varylab.discreteconformal.plugin;

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.swing.BorderFactory;
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
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.generic.Position3d;
import de.jtem.halfedgetools.algorithm.computationalgeometry.ConvexHull;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.algorithm.generator.RandomSphereGenerator;
import de.jtem.java2dx.beans.Viewer2DWithInspector;
import de.jtem.java2dx.modelling.GraphicsModeller2D;
import de.jtem.java2dx.modelling.SimpleModeller2D;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.jtem.mfc.field.Complex;
import de.jtem.mfc.group.Moebius;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.CoPositionAdapter;
import de.varylab.discreteconformal.heds.adapter.CoTexturePositionAdapter;

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
	
	
	private Random
		rnd = new Random();
	
	
	private List<Complex>
		fixPointsA = new LinkedList<Complex>(),
		fixPointsB = new LinkedList<Complex>(),
		muList = new LinkedList<Complex>();
		
	
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
		boolean orientation = true;
		Complex c = new Complex();
		double r = 1.0;
		
		public Circle(Complex c, double r, boolean o) {
			super();
			this.c = c;
			this.r = r;
			this.orientation = o;
		}
		
		public boolean isInside(Complex z) {
			Complex d = z.minus(c);
			if (orientation) {
				return d.abs() < r;
			} else {
				return d.abs() > r;
			}
		}
	}
	
	
	private double[] inverseStereographic(Complex Z) {
		double X = Z.getRe();
		double Y = Z.getIm();
		double d = 1 + X*X + Y*Y;
		double x = 2*X / d;
		double y = 2*Y / d;
		double z = (X*X + Y*Y - 1) / d;
		return new double[] {x, y, z};
	}
	
	
	private double[] liftToR3(Complex z) {
		return new double[] {z.re, z.im, 0};
	}
	
	
	
	private Complex stereographic(double[] p) {
		return new Complex(p[0] / (1 - p[2]), p[1] / (1 - p[2]));
	}
	
	
	private void generate() {
		AdapterSet a = hif.getAdapters();
		CoHDS hds = new CoHDS();
		double zScale = 1.0; 
		int numExtraPoints = 200;
		int circleRes = 20;
		
		// create extra points
		Set<Complex> extraPoints = new HashSet<Complex>(); 
		for (int i = 0; i < numExtraPoints; i++) {
			double[] pos = {rnd.nextGaussian(), rnd.nextGaussian(), rnd.nextGaussian()};
			Rn.normalize(pos, pos);
			Complex z = stereographic(pos);
			extraPoints.add(z);
		}
		
		
		Complex A = new Complex(0, 0);
		Complex B = new Complex(0.2, 0.2);
		Complex m = new Complex(0.4);
		Moebius s = new Moebius(A, B, m);
		Circle c = new Circle(A, 1.0, false);
		Complex centerOfMappedCircle = new Complex();
		double sr = s.getRadiusOfMappedCircle(c.c, c.r, centerOfMappedCircle);
		Circle sc = new Circle(centerOfMappedCircle, sr, true);
		
		
		// exclude extra vertices from the circles
		for (Complex ez : new HashSet<Complex>(extraPoints)) {
			Complex scaledEz = ez.divide(zScale);
			if (c.isInside(scaledEz)) {
				extraPoints.remove(ez);
			}
			if (sc.isInside(scaledEz)) {
				extraPoints.remove(ez);
			}
		}
		
		// create the vertices on the source circle
		for (int i = 0; i < circleRes; i++) {
			double phi = 2*i*PI / circleRes;
			double x = c.r * cos(phi) + c.c.re;
			double y = c.r * sin(phi) + c.c.im;
			Complex z = new Complex(x, y);
			Complex sz = s.applyTo(z);
			double[] zPos = inverseStereographic(z.times(zScale));
			double[] szPos = inverseStereographic(sz.times(zScale));
			CoVertex v = hds.addNewVertex();
			CoVertex vs = hds.addNewVertex();
			a.set(Position.class, v, zPos);
			a.set(Position.class, vs, szPos);
		}
		
		for (Complex e : extraPoints) {
			double[] ePos = inverseStereographic(e);
			CoVertex v = hds.addNewVertex();
			a.set(Position.class, v, ePos);	
		}
		
		ConvexHull.convexHull(hds, a, 1E-8);
		
		hif.set(hds);
	}
	
	
	
	private void resetViewer() {
		int g = genusModel.getNumber().intValue();
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
		JRViewer v = new JRViewer();
		v.addContentUI();
		v.addBasicUI();
		DiscreteSchottkyGenerator dsg = new DiscreteSchottkyGenerator();
		v.registerPlugin(dsg);
		v.startup();
	}
	
	
	
}
