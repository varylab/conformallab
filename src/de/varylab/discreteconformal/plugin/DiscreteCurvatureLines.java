package de.varylab.discreteconformal.plugin;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.signum;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;

import de.jreality.math.Rn;
import de.jreality.plugin.basic.View;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.CurvatureField;
import de.jtem.halfedgetools.adapter.type.Normal;
import de.jtem.halfedgetools.adapter.type.generic.EdgeVector;
import de.jtem.halfedgetools.plugin.HalfedgeInterface;
import de.jtem.jrworkspace.plugin.Controller;
import de.jtem.jrworkspace.plugin.sidecontainer.SideContainerPerspective;
import de.jtem.jrworkspace.plugin.sidecontainer.template.ShrinkPanelPlugin;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;

public class DiscreteCurvatureLines extends ShrinkPanelPlugin implements ActionListener {

	private HalfedgeInterface
		hif = null;
	private JButton
		goButton = new JButton("Go");
	
	
	public DiscreteCurvatureLines() {
		shrinkPanel.setTitle("Curvature Line Parametrization");
		shrinkPanel.add(goButton);
		goButton.addActionListener(this);
	}
	
	private double normalizeAngle(double a) {
		a %= 2*PI;
		if (abs(a) > PI) {
			a -= signum(a) * 2*PI;
		}
		return a;
	}

	/**
	 * Calculate the angle at the target vertex of e
	 * @param e
	 * @return
	 */
	public static double getAngle(CoEdge e, AdapterSet a) {
		double[] v1 = a.get(EdgeVector.class, e.getOppositeEdge(), double[].class);
		double[] v2 = a.get(EdgeVector.class, e.getNextEdge(), double[].class);
		double[] cr = Rn.crossProduct(null, v1, v2);
		double x = Rn.innerProduct(v1, v2);
		double y = Rn.euclideanNorm(cr);
		return Math.atan2(y, x);
	}
	
	
	private void initAlphas(CoHDS hds, AdapterSet a) {
		for (CoFace f : hds.getFaces()) {
			double bSum = 0;
			for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
				double[] eVec = a.getD(EdgeVector.class, e.getPreviousEdge());
				double[] xVec = a.getD(CurvatureField.class, e.getPreviousEdge());
				double[] nVec = a.getD(Normal.class, e.getPreviousEdge());
				double[][] B = {nVec, eVec, xVec};
				double sign = Math.signum(Rn.determinant(B));
				double alphaP = normalizeAngle(Rn.euclideanAngle(eVec, xVec) * sign);
				if (alphaP < 0) {
					alphaP = alphaP + PI;
				} else {
					alphaP = alphaP - PI;
				}
				eVec = a.getD(EdgeVector.class, e);
				xVec = a.getD(CurvatureField.class, e);
				nVec = a.getD(Normal.class, e);
				B = new double[][] {nVec, eVec, xVec};
				sign = Math.signum(Rn.determinant(B));
				double alpha = normalizeAngle(Rn.euclideanAngle(eVec, xVec) * sign);
				
				double th = getAngle(e, a);
//				th *= surfaceOrientation;
				double gamma = normalizeAngle(th - alphaP + alpha);
				if(abs(gamma) > PI/2) { // flip x
//					System.out.println("flip at " + v.getIndex());
					alpha += PI;
					normalizeAngle(alpha);
//					Rn.times(xVec, -1, xVec); // flip vector
				}
				double a1 = alpha < 0 ? alpha + 2*PI : alpha;
				double a2 = alphaP < 0 ? alphaP + 2*PI : alphaP;
				double theta = a1 - a2;
				if (theta > 0) {
					theta = 2*PI - theta;
				}
				theta = abs(theta);
//				System.out.println("sum theta: " + th);
//				System.out.println("gamma: " + gamma);
//				System.out.println("a: " + alpha + "   ap: " + alphaP);
//				System.out.println("Theta at " + v.getIndex() + ": " + theta);
//				System.out.println("");
				bSum += Math.PI - theta;
				e.setAlpha(theta);

//				e = e.getNextEdge();
			}
		}
	}
	
	
	private void doLayout(CoHDS hds, AdapterSet a) {
		
	}
	
	
	@Override
	public void actionPerformed(ActionEvent ae) {
		AdapterSet a = hif.getAdapters();
		CoHDS hds = hif.get(new CoHDS());
		initAlphas(hds, a);
		doLayout(hds, a);
		
		for (CoFace f : hds.getFaces()) {
			double aSum = 0;
			for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
				aSum += e.getAlpha();
			}
			System.out.println(f + ": " + aSum);
		}
		
		for (CoVertex v : hds.getVertices()) {
			double aSum = 0;
			for (CoEdge e : HalfEdgeUtils.incomingEdges(v)) {
				aSum += e.getAlpha();
			}
			System.out.println(v + ": " + aSum);
		}
		
//		hif.update();
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

}
