package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.varylab.discreteconformal.uniformization.FundamentalPolygonUtility.context;
import static java.awt.BasicStroke.CAP_SQUARE;
import static java.awt.BasicStroke.JOIN_ROUND;
import static java.awt.RenderingHints.KEY_ANTIALIASING;
import static java.awt.RenderingHints.VALUE_ANTIALIAS_ON;
import static java.awt.geom.Arc2D.OPEN;
import static java.lang.Math.cos;
import static java.lang.Math.signum;
import static java.lang.Math.sin;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashSet;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.scene.IndexedLineSet;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.math.PnBig;
import de.varylab.discreteconformal.math.RnBig;

public class VisualizationUtility {
	
	private static BigDecimal 
		HALF = new BigDecimal(0.5);
	private static double 
		eps = 1E-15;

	
	
	public static void drawTriangulation(
		CoHDS surface,
		HyperbolicModel model,
		Graphics2D g,
		int res,
		Color color
	) {
		float ls = res / 500f; // line scale
		g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
		g.scale(0.5, -0.5);
		g.translate(res, -res);
		g.setStroke(new BasicStroke(2 * ls));
		g.setColor(color);
		
		for (CoEdge e : surface.getPositiveEdges()) {
			double[] s = e.getStartVertex().T;
			double[] t = e.getTargetVertex().T;
			double[] m = Rn.linearCombination(null, 0.5, s, 0.5, t);
			Pn.normalize(m, m, Pn.HYPERBOLIC);
			double[] p1, p2, p3;
			switch (model) {
				case Klein:
					p1 = new double[] {s[0] / s[3], s[1] / s[3], 0};
					p2 = new double[] {t[0] / t[3], t[1] / t[3], 0};
					p3 = new double[] {m[0] / m[3], m[1] / m[3], 0};		
					break;
				default:
				case Poincaré:
					p1 = new double[] {s[0] / (s[3] + 1), s[1] / (s[3] + 1), 0};
					p2 = new double[] {t[0] / (t[3] + 1), t[1] / (t[3] + 1), 0};
					p3 = new double[] {m[0] / (m[3] + 1), m[1] / (m[3] + 1), 0};
					break;
				case Halfplane:
					p1 = new double[] {s[1] / (s[3] - s[0]), 1 / (s[3] - s[0]), 0};
					p2 = new double[] {t[1] / (t[3] - t[0]), 1 / (t[3] - t[0]), 0};
					p3 = new double[] {m[1] / (m[3] - m[0]), 1 / (m[3] - m[0]), 0};
					break;
			}
			drawArc(p1, p2, p3, g, model, res);
		}
		
		g.translate(-res, res);
		g.scale(2, -2);
	}
	
	
	public static IndexedLineSet createSurfaceCurves(
		FundamentalPolygon poly, 
		int depth,
		HyperbolicModel model,
		Graphics2D g,
		int res,
		Color rootColor
	) {
		
		
		
		return new IndexedLineSet();
	}
	
	
	
	
	public static void drawUniversalCoverImage(
		FundamentalPolygon poly, 
		int depth,
		HyperbolicModel model,
		Graphics2D g,
		int res,
		Color rootColor
	) {
		float ls = res / 500f; // line scale
		g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
		g.setColor(new Color(255, 255, 255, 0));
		g.fillRect(0, 0, res, res);
		g.scale(0.5, -0.5);
		g.translate(res, -res);

		if (model == HyperbolicModel.Klein || model == HyperbolicModel.Poincaré) {
			g.setColor(Color.BLACK);
			g.setStroke(new BasicStroke(4 * ls));
			Ellipse2D boundary = new Ellipse2D.Double(-res,-res, 2*res, 2*res);
			g.draw(boundary);
		}
		
//		g.setColor(Color.BLACK);
//		g.setStroke(new BasicStroke(1.5f * ls));
		g.setStroke(new BasicStroke(4 * ls));
		g.setColor(rootColor);
		drawUniversalCover(poly, depth, true, g, model, res);
		
		g.setStroke(new BasicStroke(4 * ls));
		g.setColor(rootColor);
		drawPolygon(poly, model, g, res);

		g.setStroke(new BasicStroke(3 * ls, CAP_SQUARE, JOIN_ROUND, 1.0f, new float[] {10 * ls, 10 * ls}, 1.0f));
		g.setColor(Color.BLUE);
		drawPolygonAxes(poly, model, g, res);
		
		g.translate(-res, res);
		g.scale(2, -2);
	}
	
	
	private static boolean drawPolygon(
		FundamentalPolygon poly, 
		HyperbolicModel model, 
		Graphics2D g,
		int resolution
	) {
		boolean proceed = true;
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
			if (!drawArc) {
				continue;
			}
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
			drawArc(p1, p2, p3, g, model, resolution);
		}
		return proceed;
	}
	
	
	private static void drawPolygonAxes(
		FundamentalPolygon poly, 
		HyperbolicModel model, 
		Graphics2D g,
		int resolution
	) {
		Set<FundamentalEdge> axisDrawn = new HashSet<FundamentalEdge>();
		for (FundamentalEdge e : poly.getEdges()) {
			if (axisDrawn.contains(e) || axisDrawn.contains(e.partner)) {
				continue;
			}
			drawAxis(e.motionBig, model, g, resolution);
			axisDrawn.add(e);
		}
	}
	
	
	private static void drawAxis(
		BigDecimal[] Ta,		
		HyperbolicModel model, 
		Graphics2D g,
		int resolution
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
			return;
		}
		DenseMatrix ev = evd.getRightEigenvectors();
		double[] evl = evd.getRealEigenvalues();
		int i1 = evl[0] > evl[1] ? (evl[0] > evl[2] ? 0 : 2) : (evl[1] > evl[2] ? 1 : 2);
		int i2 = evl[0] < evl[1] ? (evl[0] < evl[2] ? 0 : 2) : (evl[1] < evl[2] ? 1 : 2);
		double[] f1 = {ev.get(0, i1) / ev.get(3, i1), ev.get(1, i1) / ev.get(3, i1), 0, 1.0};
		double[] f2 = {ev.get(0, i2) / ev.get(3, i2), ev.get(1, i2) / ev.get(3, i2), 0, 1.0};
		double[] f3 = Rn.linearCombination(null, 0.5, f1, 0.5, f2);
		Pn.normalize(f3, f3, HYPERBOLIC);
		double[] p1, p2, p3;
		switch (model) {
			case Klein:
				p1 = new double[] {f1[0] / f1[3], f1[1] / f1[3], 0};
				p2 = new double[] {f2[0] / f2[3], f2[1] / f2[3], 0};
				p3 = new double[] {f3[0] / f3[3], f3[1] / f3[3], 0};		
				break;
			default:
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
		drawArc(p1, p2, p3, g, model, resolution);
	}
	
	
	private static void drawUniversalCover(
		FundamentalPolygon poly,
		int maxDepth,
		boolean drawAxes,
		Graphics2D g,
		HyperbolicModel model,
		int resolution
	) {
		BigDecimal[] id = new BigDecimal[16];
		RnBig.setIdentityMatrix(id);
		drawUniversalCoverR(poly, id, 0, maxDepth, drawAxes, g, model, resolution);
	}
		
		
	private static void drawUniversalCoverR(
		FundamentalPolygon poly,
		BigDecimal[] domain,
		int depth,
		int maxDepth,
		boolean drawAxes,
		Graphics2D g,
		HyperbolicModel model,
		int resolution
	) {
		FundamentalPolygon rP = FundamentalPolygonUtility.copyPolygon(poly);
		for (FundamentalEdge ce : rP.getEdges()) {
			BigDecimal[] domainInv = RnBig.inverse(null, domain, context);
			RnBig.times(ce.motionBig, ce.motionBig, domain, context);
			RnBig.times(ce.motionBig, domainInv, ce.motionBig, context);
			RnBig.matrixTimesVector(ce.startPosition, domain, ce.startPosition, context);
			PnBig.normalize(ce.startPosition, ce.startPosition, HYPERBOLIC, context);
		}
		if (drawAxes) {
			float ls = resolution / 500f;
			Stroke storedStroke = g.getStroke();
			Color storedColor = g.getColor();
//			g.setStroke(new BasicStroke(1.5f * ls, CAP_SQUARE, JOIN_ROUND, 1.0f, new float[] {10 * ls, 10 * ls}, 1.0f));
//			g.setColor(Color.BLUE);
			g.setStroke(new BasicStroke(3 * ls, CAP_SQUARE, JOIN_ROUND, 1.0f, new float[] {10 * ls, 10 * ls}, 1.0f));
			g.setColor(Color.BLUE);
			drawPolygonAxes(rP, model, g, resolution);
			g.setStroke(storedStroke);
			g.setColor(storedColor);
		}
		boolean proceed = drawPolygon(rP, model, g, resolution);
		if (!proceed || depth + 1 > maxDepth) {
			return;
		}
		for (FundamentalEdge e : poly.getEdges()) {
			BigDecimal[] newDomain = RnBig.times(null, domain, e.motionBig, context);
			drawUniversalCoverR(poly, newDomain, depth + 1, maxDepth, drawAxes, g, model, resolution);
		}
	}
	
	
	private static void drawArc(
		double[] p1,
		double[] p2,
		double[] p3,
		Graphics2D g,
		HyperbolicModel model,
		int resolution
	) {
		try {
			if (model == HyperbolicModel.Klein) {
				g.draw(new Line2D.Double(resolution * p1[0], resolution * p1[1], resolution * p2[0], resolution * p2[1]));
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
					double cornerx = resolution * (center[0] - radius);
					double cornery = resolution * (center[1] - radius);
					double size = resolution * radius * 2;
					Arc2D arc = new Arc2D.Double(cornerx, cornery, size, size, degStartAngle, degAngle, OPEN);
					g.draw(arc);
				} else {
					g.draw(new Line2D.Double(resolution * p1[0], resolution * p1[1], resolution * p2[0], resolution * p2[1]));
				}
			}
		} catch (Exception e) {
			g.draw(new Line2D.Double(resolution * p1[0], resolution * p1[1], resolution * p2[0], resolution * p2[1]));
		}
	}
	
	
	private static double getDistToUnitCircle(BigDecimal[] p, MathContext c) {
		double x = p[0].divide(p[3].add(BigDecimal.ONE, FundamentalPolygonUtility.context), c).doubleValue();
		double y = p[1].divide(p[3].add(BigDecimal.ONE, FundamentalPolygonUtility.context), c).doubleValue();
		return 1 - Math.sqrt(x*x + y*y);
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
