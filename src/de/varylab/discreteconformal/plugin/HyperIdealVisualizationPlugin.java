package de.varylab.discreteconformal.plugin;

import static java.lang.Math.cosh;
import static java.lang.Math.sinh;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Path2D;

import javax.swing.JCheckBox;

import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Radius;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.halfedgetools.plugin.texturespace.TextureSpacePlugin;
import de.jtem.java2d.SceneComponent;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.Plugin;
import de.jtem.jrworkspace.plugin.sidecontainer.widget.ShrinkPanel;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.adapter.HyperIdealRadiusAdapter;
import de.varylab.discreteconformal.uniformization.VisualizationUtility;

public class HyperIdealVisualizationPlugin extends Plugin implements TextureSpacePlugin, ActionListener {

	private HalfedgeInterface
		hif = null;
	private SceneComponent
		root = new SceneComponent();
	private ShrinkPanel
		ui = new ShrinkPanel("Hyperideal Uniformization");
	private JCheckBox
		vertexChecker = new JCheckBox("Vertex Circles", true),
		faceChecker = new JCheckBox("Face Circles", true);
	
	public HyperIdealVisualizationPlugin() {
		ui.setLayout(new GridLayout(2, 1));
		ui.add(vertexChecker);
		ui.add(faceChecker);
		vertexChecker.addActionListener(this);
		faceChecker.addActionListener(this);
	}

	@Override
	public SceneComponent getSceneComponent() {
		return root;
	}

	@Override
	public ShrinkPanel getOptionPanel() {
		return ui;
	}
	
	/**
	 * Takes center and radius in projective model of hyperbolic space
	 * and outputs euclidean center and radius for the circle in poincare model
	 * @param center
	 * @param radius
	 * @return circle[0-1], radius[2]
	 */
	protected static double[] getEuclideanCircleFromHyperbolic(double[] center, double radius) {
		double[] result = new double[3];
		double[] p1 = {sinh(radius), 0, 0, cosh(radius)};
		double[] p2 = {0.0, sinh(radius), 0, cosh(radius)};
		double[] p3 = {-sinh(radius), 0, 0, cosh(radius)};
		double[] c = {center[0], center[1], 0.0, center[3]};
		Matrix T = MatrixBuilder.hyperbolic().translate(new double[]{0,0,0,1}, center).getMatrix();
		T.transformVector(p1);
		T.transformVector(p2);
		T.transformVector(p3);
		p1[3] += 1.0; Pn.dehomogenize(p1, p1);
		p2[3] += 1.0; Pn.dehomogenize(p2, p2);
		p3[3] += 1.0; Pn.dehomogenize(p3, p3);
		c[3] += 1.0; Pn.dehomogenize(c, c);
		double[] ec = VisualizationUtility.getCircumCenter(p1, p2, p3);
		result[0] = ec[0];
		result[1] = ec[1];
		result[2] = Rn.euclideanDistance(ec, p1);
		return result;
	}

	protected void updateCircles(CoHDS hds) {
		root.removeAllChildren();
		AdapterSet a = hif.getActiveAdapters();
		boolean hasRadii = a.isAvailable(HyperIdealRadiusAdapter.class);
		if (!hasRadii) {
			root.fireAppearanceChange();
			return;
		}
		if (vertexChecker.isSelected()) {
			Path2D vertexCircles = new Path2D.Float();
			for (CoVertex v : hds.getVertices()) {
				Double r = a.get(Radius.class, v, Double.class);
				if (r < 1E-8) continue;
				double[] ec = getEuclideanCircleFromHyperbolic(v.T, r);
				Ellipse2D circle = new Ellipse2D.Double(ec[0] - ec[2], ec[1] - ec[2], 2*ec[2], 2*ec[2]);
				vertexCircles.append(circle, false);
			}
			SceneComponent vertexComponent = new SceneComponent();
			vertexComponent.setShape(vertexCircles);
			root.addChild(vertexComponent);
		}
		if (faceChecker.isSelected()) {
			Path2D faceCircles = new Path2D.Float();
			
			SceneComponent faceComponent = new SceneComponent();
			faceComponent.setShape(faceCircles);
			root.addChild(faceComponent);
		}
		root.fireAppearanceChange();
	}
	
	@Override
	public boolean getRenderOnTop() {
		return true;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		updateCircles(hif.get(new CoHDS()));
	}
	
	@Override
	public void install(Controller c) throws Exception {
		super.install(c);
		hif = c.getPlugin(HalfedgeInterface.class);
	}

}
