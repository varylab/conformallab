package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jreality.math.Pn.normalize;
import static de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility.context;
import static java.awt.geom.Arc2D.OPEN;
import static java.lang.Math.cos;
import static java.lang.Math.signum;
import static java.lang.Math.sin;

import java.awt.Shape;
import java.awt.geom.Arc2D;
import java.awt.geom.Line2D;
import java.awt.geom.Path2D;
import java.math.BigDecimal;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import de.jreality.math.Matrix;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.discretegroup.core.DiscreteGroupElement;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.math.PnBig;
import de.varylab.discreteconformal.math.RnBig;

public class VisualizationUtility {
	
	private static BigDecimal 
		BIG_HALF = new BigDecimal(0.5);
	
	
	public static void createUniversalCover(
		FundamentalPolygon P,
		HyperbolicModel model,
		int maxDrawElements,
		double maxDrawDistance,
		boolean drawPolygon,
		boolean drawAxes,
		List<double[][]> axesSegments,
		List<double[][]> polygonSegments,
		Path2D axesPath,
		Path2D polygonPath,
		Path2D domainPath
	) {
		List<DiscreteGroupElement> G = P.createGoupElements(maxDrawElements, maxDrawDistance);
		boolean isFirst = true;
		for (DiscreteGroupElement s : G) {
			BigDecimal[] sBig = RnBig.toBig(null, s.getArray());
			createTransformedDomain(
				P, 
				sBig,
				drawPolygon, isFirst & drawAxes, 
				model, 
				axesSegments, polygonSegments, 
				axesPath, 
				polygonPath,
				isFirst ? domainPath : null
			);
			isFirst = false;
		}
	}
	
	public static void createTriangulation (
		CoHDS surface,
		FundamentalPolygon P,
		HyperbolicModel model,
		int maxDrawElements,
		double maxDrawDistance,
		Path2D triangulationPath,
		Path2D boundaryPath
	) {
		List<DiscreteGroupElement> G = P.createGoupElements(maxDrawElements, maxDrawDistance);
		for (DiscreteGroupElement g : G) {
			Matrix T = g.getMatrix();
			for (CoEdge e : surface.getPositiveEdges()) {
				double[] s = e.getStartVertex().T.clone();
				double[] t = e.getTargetVertex().T.clone();
				double[] m = Rn.linearCombination(null, 0.5, s, 0.5, t);
				T.transformVector(s);
				T.transformVector(t);
				T.transformVector(m);
				double[] p1 = dehomogenizeToModel(s, model);
				double[] p2 = dehomogenizeToModel(t, model);
				double[] p3 = dehomogenizeToModel(m, model);
				Shape arc = createArc(p1, p2, p3, model);
				triangulationPath.append(arc, false);
				if (HalfEdgeUtils.isBoundaryEdge(e)) {
					boundaryPath.append(arc, false);
				}
			}
		}
	}
	
		
	protected static void createTransformedDomain(
		FundamentalPolygon poly,
		BigDecimal[] T,
		boolean drawPolygon,
		boolean drawAxes,
		HyperbolicModel model,
		List<double[][]> axesSegments,
		List<double[][]> polygonSegments,
		Path2D axesPath,
		Path2D polygonPath,
		Path2D domainPath
	) {
		FundamentalPolygon rP = FundamentalPolygonUtility.copyPolygon(poly);
		for (FundamentalEdge ce : rP.getEdges()) {
			BigDecimal[] domainInv = RnBig.inverse(null, T, context);
			RnBig.times(ce.motionBig, ce.motionBig, T, context);
			RnBig.times(ce.motionBig, domainInv, ce.motionBig, context);
			RnBig.matrixTimesVector(ce.startPosition, T, ce.startPosition, context);
			PnBig.normalize(ce.startPosition, ce.startPosition, HYPERBOLIC, context);
		}
		if (drawAxes) {
			createPolygonAxes(rP, model, axesSegments, axesPath);
		}
		if (drawPolygon) {
			createPolygon(rP, model, polygonSegments, polygonPath, domainPath);
		}
	}
	
	
	protected static void createPolygonAxes(
		FundamentalPolygon poly, 
		HyperbolicModel model, 
		List<double[][]> segmentsOUT,
		Path2D axesPath
	) {
		Set<FundamentalEdge> axisDrawn = new HashSet<FundamentalEdge>();
		for (FundamentalEdge e : poly.getEdges()) {
			if (axisDrawn.contains(e) || axisDrawn.contains(e.partner)) {
				continue;
			}
			Shape axis = createAxis(e.motionBig, model, segmentsOUT);
			axisDrawn.add(e);
			if (axesPath != null) {
				axesPath.append(axis, false);
			}
		}
	}
	
	
	protected static Shape createAxis(
		BigDecimal[] Ta,		
		HyperbolicModel model, 
		List<double[][]> segmentsOUT
	) {
		double[][] Tad = new double[4][4];
		for (int i = 0; i < Ta.length; i++) {
			Tad[i/4][i%4] = Ta[i].doubleValue();
		}
		DenseMatrix T = new DenseMatrix(Tad);
		EVD evd = new EVD(4);
		try {
			evd.factor(T);
		} catch (NotConvergedException e1) {
			e1.printStackTrace();
			return new Line2D.Float();
		}
		DenseMatrix ev = evd.getRightEigenvectors();
		double[] evl = evd.getRealEigenvalues();
		int i1 = evl[0] > evl[1] ? (evl[0] > evl[2] ? 0 : 2) : (evl[1] > evl[2] ? 1 : 2);
		int i2 = evl[0] < evl[1] ? (evl[0] < evl[2] ? 0 : 2) : (evl[1] < evl[2] ? 1 : 2);
		double[] f1 = {ev.get(0, i1) / ev.get(3, i1), ev.get(1, i1) / ev.get(3, i1), 0, 1.0 + 1E-10};
		double[] f2 = {ev.get(0, i2) / ev.get(3, i2), ev.get(1, i2) / ev.get(3, i2), 0, 1.0 + 1E-10};
		double[] f3 = Rn.linearCombination(null, 0.5, f1, 0.5, f2);
		double[] p1 = dehomogenizeToModel(f1, model);
		double[] p2 = dehomogenizeToModel(f2, model);
		double[] p3 = dehomogenizeToModel(f3, model);
		if (segmentsOUT != null) {
			segmentsOUT.add(new double[][] {p1, p2});
		}
		return createArc(p1, p2, p3, model);
	}
	
	protected static void createPolygon(
		FundamentalPolygon poly, 
		HyperbolicModel model, 
		List<double[][]> segmentsOUT,
		Path2D polygonPath,
		Path2D domainPath
	) {
		for (FundamentalEdge e : poly.getEdges()) {
			BigDecimal[] p1a = e.startPosition; 
			BigDecimal[] p2a = e.nextEdge.startPosition;
			BigDecimal[] p3a = RnBig.linearCombination(null, BIG_HALF, p1a, BIG_HALF, p2a, FundamentalPolygonUtility.context);
			double[] p1ad = RnBig.toDouble(null, p1a);
			double[] p2ad = RnBig.toDouble(null, p2a);
			double[] p3ad = RnBig.toDouble(null, p3a);
			double[] p1 = dehomogenizeToModel(p1ad, model);
			double[] p2 = dehomogenizeToModel(p2ad, model);
			double[] p3 = dehomogenizeToModel(p3ad, model);
			if (segmentsOUT != null) {
				segmentsOUT.add(new double[][] {p1, p2});
			}
			Shape arc = createArc(p1, p2, p3, model);
			if (polygonPath != null) {
				polygonPath.append(arc, false);
			}
			if (domainPath != null) {
				domainPath.append(arc, true);
			}
		}
	}
	
	protected static Shape createArc(
		double[] p1,
		double[] p2,
		double[] p3,
		HyperbolicModel model
	) {
		Shape shape = null;
		try {
			if (model == HyperbolicModel.Klein) {
				shape = new Line2D.Float((float)p1[0], (float)p1[1], (float)p2[0], (float)p2[1]);
			} else {
				double[] center = getCircumCenter(p1, p2, p3);
				double[] vec1 = Rn.subtract(null, p1, center);
				double[] vec2 = Rn.subtract(null, p2, center);
				double angle = Rn.euclideanAngle(vec1, vec2);
				double startAngle = Math.atan2(-vec1[1], vec1[0]);
				angle *= -signum(vec2[0] * sin(startAngle) + vec2[1] * cos(startAngle));
				double degAngle = Math.toDegrees(angle);
				double degStartAngle = Math.toDegrees(startAngle);
				double radius = Rn.euclideanDistance(p1, center);
				double[] betaVec1 = Rn.subtract(null, p1, p3);
				double[] betaVec2 = Rn.subtract(null, p2, p3);
				double beta = Rn.euclideanAngle(betaVec1, betaVec2);
				if (beta < 3.1) {
					double cornerx = center[0] - radius;
					double cornery = center[1] - radius;
					double size = radius * 2;
					shape = new Arc2D.Float((float)cornerx, (float)cornery, (float)size, (float)size, (float)degStartAngle, (float)degAngle, OPEN);
				} else {
					shape = new Line2D.Float((float)p1[0], (float)p1[1], (float)p2[0], (float)p2[1]);
				}
			}
		} catch (Exception e) {
			shape = new Line2D.Float((float)p1[0], (float)p1[1], (float)p2[0], (float)p2[1]);
		}
		return shape;
	}
	
	
	/**
	 * Calculate the circum-center of a triangle in affine coordinates
	 * @return
	 */
	protected static double[] getCircumCenter(double[] A, double[] B, double[] C) {
		double a = Rn.euclideanDistance(B, C);
		double b = Rn.euclideanDistance(C, A);
		double c = Rn.euclideanDistance(A, B);
		double ca = a*a*(b*b + c*c - a*a);
		double cb = b*b*(c*c + a*a - b*b);
		double cc = c*c*(a*a + b*b - c*c);
		double l = ca + cb + cc;
		double[] r = {0, 0, 0};
		r[0] = (ca*A[0] + cb*B[0] + cc*C[0]) / l;
		r[1] = (ca*A[1] + cb*B[1] + cc*C[1]) / l;
		r[2] = (ca*A[2] + cb*B[2] + cc*C[2]) / l;
		return r;
	}
	
	
	protected static double[] dehomogenizeToModel(double[] p, HyperbolicModel model) {
		normalize(p, p, Pn.HYPERBOLIC);
		double[] result = null;
		switch (model) {
		case Klein:
			result = new double[] {p[0] / p[3], p[1] / p[3], 1};
			break;
		default:
		case Poincaré:
			result = new double[] {p[0] / (p[3] + 1), p[1] / (p[3] + 1), 1};
			break;
		case Halfplane:
			result = new double[] {p[1] / (p[3] - p[0]), 1 / (p[3] - p[0]), 1};
			break;
		}
		return result;
	}
	
}
