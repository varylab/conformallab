package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility.context;
import static java.awt.BasicStroke.CAP_SQUARE;
import static java.awt.BasicStroke.JOIN_ROUND;
import static java.awt.Color.BLACK;
import static java.awt.geom.Arc2D.OPEN;
import static java.awt.image.BufferedImage.TYPE_INT_ARGB;
import static java.lang.Math.cos;
import static java.lang.Math.signum;
import static java.lang.Math.sin;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.geom.Arc2D;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;

import de.jreality.math.Matrix;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.math.PnBig;
import de.varylab.discreteconformal.math.RnBig;

public class VisualizationUtility {
	
	private static BigDecimal 
		HALF = new BigDecimal(0.5);
	private static double 
		eps = 1E-4;

	
	public static Image drawUniversalCover(
		FundamentalPolygon poly, 
		int depth,
		HyperbolicModel model,
		int res,
		Color rootColor
	) {
		BufferedImage image = new BufferedImage(res, res, TYPE_INT_ARGB);
		Graphics2D g = (Graphics2D)image.getGraphics();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setColor(new Color(255, 255, 255, 0));
		g.fillRect(0, 0, res, res);
		
		List<BigDecimal[]> orbit = poly.getVertexPositions();
		g.setColor(BLACK);
		g.setStroke(new BasicStroke(2.0f));
		BigDecimal[] id = new BigDecimal[16];
		RnBig.setIdentityMatrix(id);
		drawPolygon(poly, orbit, id, g, res, 0, depth, model);
		
		g.scale(res / 2.0, -res / 2.0);
		g.translate(1.0, -1.0);
		
//		g.setColor(BLACK);
//		g.setStroke(new BasicStroke(5.0f/res));
//		drawCoverPolygon(poly, g, 0, depth+1, model);
		
		Stroke polygonStroke = new BasicStroke(15.0f/res);
		Color polygonColor = rootColor;
		Stroke axesStroke = new BasicStroke(10.0f/res, CAP_SQUARE, JOIN_ROUND, 1.0f, new float[] {0.02f, 0.02f}, 1.0f);
		Color axesColor = Color.ORANGE;
		drawPolygon(poly, model, g, polygonColor, polygonStroke, true, axesColor, axesStroke);
		
		g.setColor(BLACK);
		g.setStroke(new BasicStroke(15.0f/res));
		Arc2D arc = new Arc2D.Double(-1.0, -1.0, 2.0, 2.0, 0.0, 360, OPEN);
		g.draw(arc);
		return image;
	}
	
	
	private static boolean drawPolygon(
		FundamentalPolygon poly, 
		HyperbolicModel model, 
		Graphics2D g, 
		Color polygonColor,
		Stroke polygonStroke,
		boolean drawAxes,
		Color axesColor,
		Stroke axesStroke 
	) {
		boolean proceed = true;
		Set<FundamentalEdge> done = new HashSet<FundamentalEdge>();
		for (FundamentalEdge e : poly.getEdges()) {
			BigDecimal[] p1a = e.startPosition; 
			BigDecimal[] p2a = e.nextEdge.startPosition;
			BigDecimal[] p3a = RnBig.linearCombination(null, HALF, p1a, HALF, p2a, FundamentalPolygonUtility.context);
			PnBig.normalize(p3a, p3a, HYPERBOLIC, context);
			boolean drawArc = true;
			switch (model) {
			case Klein:
			case Poincaré:
				if (getDistToUnitCircle(p1a, FundamentalPolygonUtility.context) < eps && 
					getDistToUnitCircle(p2a, FundamentalPolygonUtility.context) < eps &&
					getDistToUnitCircle(p3a, FundamentalPolygonUtility.context) < eps) {
					proceed = drawArc = false;
				}
				break;
			case Halfplane:
				if (1 / (p1a[3].subtract(p1a[0], FundamentalPolygonUtility.context).doubleValue()) < eps ||
					1 / (p2a[3].subtract(p2a[0], FundamentalPolygonUtility.context).doubleValue()) < eps ||
					1 / (p3a[3].subtract(p3a[0], FundamentalPolygonUtility.context).doubleValue()) < eps) {
					proceed = drawArc = false;
				}
			}
			if (drawArc) {
				double[] p1ad = RnBig.toDouble(null, p1a);
				double[] p2ad = RnBig.toDouble(null, p2a);
				double[] p3ad = RnBig.toDouble(null, p3a);
				double[] p1, p2, p3;
				switch (model) {
					case Klein:
						p1 = new double[] {p1ad[0] / p1ad[3], p1ad[1] / p1ad[3], 0};
						p2 = new double[] {p2ad[0] / p2ad[3], p2ad[1] / p2ad[3], 0};
						p3 = new double[] {p3ad[0] / p3ad[3], p3ad[1] / p3ad[3], 0};		
						break;
					default:
					case Poincaré:
						p1 = new double[] {p1ad[0] / (p1ad[3] + 1), p1ad[1] / (p1ad[3] + 1), 0};
						p2 = new double[] {p2ad[0] / (p2ad[3] + 1), p2ad[1] / (p2ad[3] + 1), 0};
						p3 = new double[] {p3ad[0] / (p3ad[3] + 1), p3ad[1] / (p3ad[3] + 1), 0};
						break;
					case Halfplane:
						p1 = new double[] {p1ad[1] / (p1ad[3] - p1ad[0]), 1 / (p1ad[3] - p1ad[0]), 0};
						p2 = new double[] {p2ad[1] / (p2ad[3] - p2ad[0]), 1 / (p2ad[3] - p2ad[0]), 0};
						p3 = new double[] {p3ad[1] / (p3ad[3] - p3ad[0]), 1 / (p3ad[3] - p3ad[0]), 0};
						break;
				}
				g.setColor(polygonColor);
				g.setStroke(polygonStroke);
				drawArcNormalized(p1, p2, p3, g, model);
				if (done.contains(e) || done.contains(e.partner) || !drawAxes) {
					continue;
				}
				BigDecimal[] Ta = e.motionBig;
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
					continue;
				}
				DenseMatrix ev = evd.getRightEigenvectors();
				double[] evl = evd.getRealEigenvalues();
				int i1 = evl[0] > evl[1] ? (evl[0] > evl[2] ? 0 : 2) : (evl[1] > evl[2] ? 1 : 2);
				int i2 = evl[0] < evl[1] ? (evl[0] < evl[2] ? 0 : 2) : (evl[1] < evl[2] ? 1 : 2);
				double[] f1 = {ev.get(0, i1) / ev.get(3, i1), ev.get(1, i1) / ev.get(3, i1), 0, 1.0};
				double[] f2 = {ev.get(0, i2) / ev.get(3, i2), ev.get(1, i2) / ev.get(3, i2), 0, 1.0};
				double[] f3 = Rn.linearCombination(null, 0.5, f1, 0.5, f2);
				Pn.normalize(f3, f3, HYPERBOLIC);
				switch (model) {
					case Klein:
						p1 = new double[] {f1[0] / f1[3], f1[1] / f1[3], 0};
						p2 = new double[] {f2[0] / f2[3], f2[1] / f2[3], 0};
						p3 = new double[] {f3[0] / f3[3], f3[1] / f3[3], 0};		
						break;
					case Poincaré:
						p1 = new double[] {f1[0] / f1[3], f1[1] / f1[3], 0};
						p2 = new double[] {f2[0] / f2[3], f2[1] / f2[3], 0};
						p3 = new double[] {f3[0] / (f3[3] + 1), f3[1] / (f3[3] + 1), 0};
						break;
					case Halfplane:
						p1 = new double[] {f1[0] / f1[3], f1[1] / f1[3], 0};
						p2 = new double[] {f2[0] / f2[3], f2[1] / f2[3], 0};
						p3 = new double[] {f3[1] / (f3[3] - f3[0]), 1 / (f3[3] - f3[0]), 0};
						break;
				}
				g.setColor(axesColor);
				g.setStroke(axesStroke);
				drawArcNormalized(p1, p2, p3, g, model);
			}
			done.add(e);
		}
		return proceed;
	}
	
	
	private static void drawArcNormalized(
		double[] p1,
		double[] p2,
		double[] p3,
		Graphics2D g,
		HyperbolicModel model
	) {
		try {
			if (model == HyperbolicModel.Klein) {
				g.draw(new Line2D.Double(p1[0], p1[1], p2[0], p2[1]));
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
				Arc2D arc = new Arc2D.Double(center[0] - radius, center[1] - radius, radius * 2, radius * 2, degStartAngle, degAngle, OPEN);
				g.draw(arc);
			}
		} catch (Exception e) {
			g.draw(new Line2D.Double(p1[0], p1[1], p2[0], p2[1]));
		}
	}
	
	
	
//	private static void drawCoverPolygon(
//		FundamentalPolygon poly,
//		Graphics2D g,
//		int depth,
//		int maxDepth,
//		HyperbolicModel model
//	) {
//		boolean proceed = drawPolygon(poly, model, g);
//		if (depth + 1 > maxDepth) {
//			return;
//		}
//		if (proceed) {
//			for (FundamentalEdge e : poly.getEdges()) {
//				FundamentalPolygon rP = FundamentalPolygonUtility.copyPolygon(poly);
//				for (FundamentalEdge ce : rP.getEdges()) {
//					RnBig.times(ce.motionBig, e.motionBig, ce.motionBig, context);
//					RnBig.times(ce.motionBig, ce.motionBig, e.partner.motionBig, context);
//					ce.motion = Matrix.times(ce.motion, e.motion);
//					RnBig.matrixTimesVector(ce.startPosition, e.motionBig, ce.startPosition, context);
//					PnBig.normalize(ce.startPosition, ce.startPosition, HYPERBOLIC, context);
//				}
//				drawCoverPolygon(rP, g, depth + 1, maxDepth, model);
//			}
//		}
//	}
	
	
	
	private static void drawPolygon(
		FundamentalPolygon poly,
		List<BigDecimal[]> orbit,
		BigDecimal[] domain,
		Graphics2D g,
		int res,
		int depth,
		int maxDepth,
		HyperbolicModel model
	) {
		if (depth > maxDepth) {
			return;
		}
		double eps = 1E-4;
		for (FundamentalEdge fe : poly.getEdges()) {
			BigDecimal[] T = RnBig.times(null, domain, fe.motionBig, FundamentalPolygonUtility.context);
			boolean proceed = true;
			for (int i = 0; i < orbit.size(); i++) {
				BigDecimal[] p1a = orbit.get(i).clone();
				BigDecimal[] p2a = orbit.get((i + 1) % orbit.size()).clone();
				BigDecimal[] p3a = RnBig.linearCombination(null, HALF, p1a, HALF, p2a, FundamentalPolygonUtility.context);
				PnBig.normalize(p3a, p3a, HYPERBOLIC, FundamentalPolygonUtility.context);
				RnBig.matrixTimesVector(p1a, T, p1a, FundamentalPolygonUtility.context);
				RnBig.matrixTimesVector(p2a, T, p2a, FundamentalPolygonUtility.context);
				RnBig.matrixTimesVector(p3a, T, p3a, FundamentalPolygonUtility.context);
				PnBig.normalize(p1a, p1a, HYPERBOLIC, FundamentalPolygonUtility.context);
				PnBig.normalize(p2a, p2a, HYPERBOLIC, FundamentalPolygonUtility.context);
				PnBig.normalize(p3a, p3a, HYPERBOLIC, FundamentalPolygonUtility.context);
				boolean drawArc = true;
				switch (model) {
				case Klein:
				case Poincaré:
					if (getDistToUnitCircle(p1a, FundamentalPolygonUtility.context) < eps && 
						getDistToUnitCircle(p2a, FundamentalPolygonUtility.context) < eps &&
						getDistToUnitCircle(p3a, FundamentalPolygonUtility.context) < eps) {
						proceed = drawArc = false;
					}
					break;
				case Halfplane:
					if (1 / (p1a[3].subtract(p1a[0], FundamentalPolygonUtility.context).doubleValue()) < eps ||
						1 / (p2a[3].subtract(p2a[0], FundamentalPolygonUtility.context).doubleValue()) < eps ||
						1 / (p3a[3].subtract(p3a[0], FundamentalPolygonUtility.context).doubleValue()) < eps) {
						proceed = drawArc = false;
					}
				}
				if (drawArc) {
					double[] p1aDouble = RnBig.toDouble(null, p1a);
					double[] p2aDouble = RnBig.toDouble(null, p2a);
					double[] p3aDouble = RnBig.toDouble(null, p3a);
					drawArc(p1aDouble, p2aDouble, p3aDouble, g, res, model);
				}
			}
			if (proceed) {
				drawPolygon(poly, orbit, T, g, res, depth + 1, maxDepth, model);
			}
		}
	}
	
	private static double getDistToUnitCircle(BigDecimal[] p, MathContext c) {
		double x = p[0].divide(p[3].add(BigDecimal.ONE, FundamentalPolygonUtility.context), c).doubleValue();
		double y = p[1].divide(p[3].add(BigDecimal.ONE, FundamentalPolygonUtility.context), c).doubleValue();
		return 1 - Math.sqrt(x*x + y*y);
	}
	
	private static void drawArc(
		double[] p1a,
		double[] p2a,
		double[] p3a,
		Graphics2D g,
		int res,
		HyperbolicModel model
	) {
		double[] p1 = null;
		double[] p2 = null;
		double[] p3 = null;
		switch (model) {
			case Klein:
				p1 = new double[] {p1a[0] / p1a[3], p1a[1] / p1a[3], 0};
				p2 = new double[] {p2a[0] / p2a[3], p2a[1] / p2a[3], 0};
				p3 = new double[] {p3a[0] / p3a[3], p3a[1] / p3a[3], 0};		
				break;
				default:
			case Poincaré:
				p1 = new double[] {p1a[0] / (p1a[3] + 1), p1a[1] / (p1a[3] + 1), 0};
				p2 = new double[] {p2a[0] / (p2a[3] + 1), p2a[1] / (p2a[3] + 1), 0};
				p3 = new double[] {p3a[0] / (p3a[3] + 1), p3a[1] / (p3a[3] + 1), 0};
				break;
			case Halfplane:
				p1 = new double[] {p1a[1] / (p1a[3] - p1a[0]), 1 / (p1a[3] - p1a[0]), 0};
				p2 = new double[] {p2a[1] / (p2a[3] - p2a[0]), 1 / (p2a[3] - p2a[0]), 0};
				p3 = new double[] {p3a[1] / (p3a[3] - p3a[0]), 1 / (p3a[3] - p3a[0]), 0};
		}
		try {
			if (model == HyperbolicModel.Klein) {
				double sx = (p1[0] / 2 + 0.5) * res;
				double sy = (-p1[1] / 2 + 0.5) * res;
				double tx = (p2[0] / 2 + 0.5) * res;
				double ty = (-p2[1] / 2 + 0.5) * res;
				Line2D line = new Line2D.Double(sx, sy, tx, ty);
				g.draw(line);
			} else {
				double[] center = getCircumCenter(p1, p2, p3);
				double cRad = Rn.euclideanDistance(p1, center);
				double[] vec1 = Rn.subtract(null, p1, center);
				double[] vec2 = Rn.subtract(null, p2, center);
				double angle = Rn.euclideanAngle(vec1, vec2);
				double startAngle = Math.atan2(vec1[1], vec1[0]);
				angle *= signum(-vec2[0] * sin(startAngle) + vec2[1] * cos(startAngle));
				
				double degAngle = Math.toDegrees(angle);
				double degStartAngle = Math.toDegrees(startAngle);
				double centerX = (center[0] / 2 + 0.5) * res;
				double centerY = (-center[1] / 2 + 0.5) * res;
				double radius = (cRad / 2) * res;
				Arc2D arc = new Arc2D.Double(centerX - radius, centerY - radius, radius * 2, radius * 2, degStartAngle, degAngle, OPEN);
				g.draw(arc);
			}
		} catch (Exception e) {
			int p1x = (int)((p1[0] / 2 + 0.5) * res);
			int p1y = (int)((-p1[1] / 2 + 0.5) * res);
			int p2x = (int)((p2[0] / 2 + 0.5) * res);
			int p2y = (int)((-p2[1] / 2 + 0.5) * res);
			g.drawLine(p1x, p1y, p2x, p2y);
		}
		
	}
	
	
	/**
	 * Calculate the circum-center of a triangle in affine coordinates
	 * @param Ap
	 * @param Bp
	 * @param Cp
	 * @return
	 */
	private static double[] getCircumCenter(double[] A, double[] B, double[] C) {
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
	
	
	
	
}
