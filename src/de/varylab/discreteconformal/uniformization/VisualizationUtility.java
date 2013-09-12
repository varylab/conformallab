package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.jtem.halfedge.util.HalfEdgeUtils.isBoundaryEdge;
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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;
import no.uib.cipr.matrix.NotConvergedException;
import de.jreality.math.Matrix;
import de.jreality.math.P2;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.discretegroup.core.DiscreteGroup;
import de.jtem.discretegroup.core.DiscreteGroupConstraint;
import de.jtem.discretegroup.core.DiscreteGroupElement;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoHDS;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.heds.CustomVertexInfo;
import de.varylab.discreteconformal.math.PnBig;
import de.varylab.discreteconformal.math.RnBig;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;

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
		g.setStroke(new BasicStroke(1 * ls));
		g.setColor(color);
		
		for (CoEdge e : surface.getPositiveEdges()) {
			double[] s = e.getStartVertex().T;
			double[] t = e.getTargetVertex().T;
			double[] m = Rn.linearCombination(null, 0.5, s, 0.5, t);
			Pn.normalize(m, m, Pn.HYPERBOLIC);
			double[] p1, p2, p3;
			switch (model) {
				case Klein:
					p1 = new double[] {s[0] / s[3], s[1] / s[3], 1};
					p2 = new double[] {t[0] / t[3], t[1] / t[3], 1};
					p3 = new double[] {m[0] / m[3], m[1] / m[3], 1};		
					break;
				default:
				case Poincaré:
					p1 = new double[] {s[0] / (s[3] + 1), s[1] / (s[3] + 1), 1};
					p2 = new double[] {t[0] / (t[3] + 1), t[1] / (t[3] + 1), 1};
					p3 = new double[] {m[0] / (m[3] + 1), m[1] / (m[3] + 1), 1};
					break;
				case Halfplane:
					p1 = new double[] {s[1] / (s[3] - s[0]), 1 / (s[3] - s[0]), 1};
					p2 = new double[] {t[1] / (t[3] - t[0]), 1 / (t[3] - t[0]), 1};
					p3 = new double[] {m[1] / (m[3] - m[0]), 1 / (m[3] - m[0]), 1};
					break;
			}
			drawArc(p1, p2, p3, g, model, res);
		}
		
		g.translate(-res, res);
		g.scale(2, -2);
	}
	
	
	public static void drawUniversalCoverImage(
		FundamentalPolygon poly,
		boolean drawPolygon,
		boolean drawAxes,
		int maxDrawDepth,
		double maxDrawDistance,
		HyperbolicModel model,
		Graphics2D g,
		int res,
		Color polygonColor,
		Color axesColor
	) {
		float ls = res / 500f; // line scale
		g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
		g.scale(0.5, -0.5);
		g.translate(res, -res);

		g.setStroke(new BasicStroke(4 * ls));
		g.setColor(polygonColor);
		drawUniversalCover(poly, maxDrawDepth, maxDrawDistance, drawPolygon, drawAxes, g, model, res, polygonColor, axesColor, null, null);
		
		if (model == HyperbolicModel.Klein || model == HyperbolicModel.Poincaré) {
			Ellipse2D boundary = new Ellipse2D.Double(-res + ls, -res + ls, 2*res - 2*ls, 2*res - 2*ls);
			g.setColor(Color.BLACK);
			g.setStroke(new BasicStroke(4 * ls));
			g.draw(boundary);
		}
		
		g.translate(-res, res);
		g.scale(2, -2);
	}
	
	
	public static void drawDomainBackground(Graphics2D g, int res, HyperbolicModel model) {
		float ls = res / 500f; // line scale
		g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
		g.setColor(new Color(255, 255, 255, 0));
		g.fillRect(0, 0, res, res);
		g.scale(0.5, -0.5);
		g.translate(res, -res);
		
		if (model == HyperbolicModel.Klein || model == HyperbolicModel.Poincaré) {
			Ellipse2D boundary = new Ellipse2D.Double(-res + ls, -res + ls, 2*res - 2*ls, 2*res - 2*ls);
			g.setColor(Color.WHITE);
			g.fill(boundary);
		}
		g.translate(-res, res);
		g.scale(2, -2);
	}
	
	private static boolean drawPolygon(
		FundamentalPolygon poly, 
		HyperbolicModel model, 
		Graphics2D g,
		int resolution,
		List<double[][]> segmentsOUT
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
					p1 = new double[] {p1ad[0] / p1ad[3], p1ad[1] / p1ad[3], 1};
					p2 = new double[] {p2ad[0] / p2ad[3], p2ad[1] / p2ad[3], 1};
					p3 = new double[] {p3ad[0] / p3ad[3], p3ad[1] / p3ad[3], 1};
					if (segmentsOUT != null) {
						segmentsOUT.add(new double[][] {p1, p2});
					}
					break;
				default:
				case Poincaré:
					p1 = new double[] {p1ad[0] / (p1ad[3] + 1), p1ad[1] / (p1ad[3] + 1), 1};
					p2 = new double[] {p2ad[0] / (p2ad[3] + 1), p2ad[1] / (p2ad[3] + 1), 1};
					p3 = new double[] {p3ad[0] / (p3ad[3] + 1), p3ad[1] / (p3ad[3] + 1), 1};
					break;
				case Halfplane:
					p1 = new double[] {p1ad[1] / (p1ad[3] - p1ad[0]), 1 / (p1ad[3] - p1ad[0]), 1};
					p2 = new double[] {p2ad[1] / (p2ad[3] - p2ad[0]), 1 / (p2ad[3] - p2ad[0]), 1};
					p3 = new double[] {p3ad[1] / (p3ad[3] - p3ad[0]), 1 / (p3ad[3] - p3ad[0]), 1};
					break;
			}
			if (g != null) {
				drawArc(p1, p2, p3, g, model, resolution);
			}
		}
		return proceed;
	}
	
	
	protected static List<double[][]> getPolygonAxesSegments(FundamentalPolygon poly) {
		List<double[][]> result = new ArrayList<double[][]>();
		drawPolygonAxes(poly, HyperbolicModel.Klein, null, -1, result);
		return result;
	}
	
	
	private static void drawPolygonAxes(
		FundamentalPolygon poly, 
		HyperbolicModel model, 
		Graphics2D g,
		int resolution,
		List<double[][]> segmentsOUT
	) {
		Set<FundamentalEdge> axisDrawn = new HashSet<FundamentalEdge>();
		for (FundamentalEdge e : poly.getEdges()) {
			if (axisDrawn.contains(e) || axisDrawn.contains(e.partner)) {
				continue;
			}
			drawAxis(e.motionBig, model, g, resolution, segmentsOUT);
			axisDrawn.add(e);
		}
	}
	
	
	private static void drawAxis(
		BigDecimal[] Ta,		
		HyperbolicModel model, 
		Graphics2D g,
		int resolution,
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
				p1 = new double[] {f1[0] / f1[3], f1[1] / f1[3], 1};
				p2 = new double[] {f2[0] / f2[3], f2[1] / f2[3], 1};
				p3 = new double[] {f3[0] / f3[3], f3[1] / f3[3], 1};
				if (segmentsOUT != null) {
					segmentsOUT.add(new double[][] {p1, p2});
				}
				break;
			default:
			case Poincaré:
				p1 = new double[] {f1[0] / f1[3], f1[1] / f1[3], 1};
				p2 = new double[] {f2[0] / f2[3], f2[1] / f2[3], 1};
				p3 = new double[] {f3[0] / (f3[3] + 1), f3[1] / (f3[3] + 1), 1};
				break;
			case Halfplane:
				p1 = new double[] {f1[0] / f1[3], f1[1] / f1[3], 1};
				p2 = new double[] {f2[0] / f2[3], f2[1] / f2[3], 1};
				p3 = new double[] {f3[1] / (f3[3] - f3[0]), 1 / (f3[3] - f3[0]), 1};
				break;
		}
		if (g != null) {
			drawArc(p1, p2, p3, g, model, resolution);
		}
	}
	
	
	protected static void getUniversalCoverSegments(
		FundamentalPolygon poly,
		int maxDepth,
		double maxDrawDistance,
		boolean drawPolygon,
		boolean drawAxes,
		Color polygonColor,
		Color axesColor,
		List<double[][]> axesSegments,
		List<double[][]> polygonSegments
	) {
		drawUniversalCover(poly, maxDepth, maxDrawDistance, drawPolygon, drawAxes, null, HyperbolicModel.Klein, -1, polygonColor, axesColor, axesSegments, polygonSegments);
	}
	
	
	private static void drawUniversalCover(
		FundamentalPolygon poly,
		final int maxDrawElements,
		final double maxDrawDistance,
		boolean drawPolygon,
		boolean drawAxes,
		Graphics2D g,
		HyperbolicModel model,
		int resolution,
		Color polygonColor,
		Color axesColor,
		List<double[][]> axesSegments,
		List<double[][]> polygonSegments
	) {
		BigDecimal[] id = new BigDecimal[16];
		RnBig.setIdentityMatrix(id);
		DiscreteGroupConstraint constraint = new DiscreteGroupConstraint() {
			double[] zero = {0,0,0,1};
			@Override
			public void update() {
			}
			@Override
			public void setMaxNumberElements(int arg0) {
			}
			@Override
			public int getMaxNumberElements() {
				return maxDrawElements;
			}
			@Override
			public boolean acceptElement(DiscreteGroupElement s) {
				double[] sZero = s.getMatrix().multiplyVector(zero);
				double dist = Pn.distanceBetween(zero, sZero, Pn.HYPERBOLIC);
				return dist <= maxDrawDistance;
			}
		};
		DiscreteGroup G = poly.getDiscreteGroup();
		G.setConstraint(constraint);
		G.generateElements();
		boolean isFirst = true;
		for (DiscreteGroupElement s : G.getElementList()) {
			BigDecimal[] sBig = RnBig.toBig(null, s.getArray());
			drawUniversalCoverR(poly, sBig, 0, 0, drawPolygon, isFirst, g, model, resolution, polygonColor, axesColor, axesSegments, polygonSegments);
			isFirst = false;
		}
	}
		
		
	private static void drawUniversalCoverR(
		FundamentalPolygon poly,
		BigDecimal[] domain,
		int depth,
		int maxDepth,
		boolean drawPolygon,
		boolean drawAxes,
		Graphics2D g,
		HyperbolicModel model,
		int resolution,
		Color polygonColor,
		Color axesColor,
		List<double[][]> axesSegments,
		List<double[][]> polygonSegments
	) {
		FundamentalPolygon rP = FundamentalPolygonUtility.copyPolygon(poly);
		for (FundamentalEdge ce : rP.getEdges()) {
			BigDecimal[] domainInv = RnBig.inverse(null, domain, context);
			RnBig.times(ce.motionBig, ce.motionBig, domain, context);
			RnBig.times(ce.motionBig, domainInv, ce.motionBig, context);
			RnBig.matrixTimesVector(ce.startPosition, domain, ce.startPosition, context);
			PnBig.normalize(ce.startPosition, ce.startPosition, HYPERBOLIC, context);
		}
		if (drawAxes && depth == 0) {
			if (g != null) {
				float ls = resolution / 500f;
				Stroke storedStroke = g.getStroke();
				Color storedColor = g.getColor();
				g.setStroke(new BasicStroke(3 * ls, CAP_SQUARE, JOIN_ROUND, 1.0f, new float[] {10 * ls, 10 * ls}, 1.0f));
				g.setColor(axesColor);
				drawPolygonAxes(rP, model, g, resolution, axesSegments);
				g.setColor(storedColor);
				g.setStroke(storedStroke);
			} else {
				drawPolygonAxes(rP, model, null, resolution, axesSegments);
			}
		}
		boolean proceed = true;
		if (drawPolygon) {
			if (g != null) {
				g.setColor(polygonColor);
			}
			proceed = drawPolygon(rP, model, g, resolution, polygonSegments);
		}
		if (!proceed || depth + 1 > maxDepth) {
			return;
		}
		for (FundamentalEdge e : poly.getEdges()) {
			BigDecimal[] newDomain = RnBig.times(null, domain, e.motionBig, context);
			drawUniversalCoverR(poly, newDomain, depth + 1, maxDepth, drawPolygon, drawAxes, g, model, resolution, polygonColor, axesColor, axesSegments, polygonSegments);
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
	
	
	public static class ValenceComparator implements Comparator<CoVertex> {

		@Override
		public int compare(CoVertex o1, CoVertex o2) {
			int v1 = HalfEdgeUtils.incomingEdges(o1).size();
			int v2 = HalfEdgeUtils.incomingEdges(o2).size();
			return v1 - v2;
		}
		
	}
	
	
	public static Set<CoFace> reglueOutsideFaces(
		CoHDS hds, 
		int maxFaces, 
		FundamentalPolygon p, 
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo,
		int signature
	) {
		int counter = 0;
		Set<CoFace> reglued = new HashSet<CoFace>(); 
		Collection<CoVertex> bList = HalfEdgeUtils.boundaryVertices(hds);
		PriorityQueue<CoVertex> lowValenceFirst = new PriorityQueue<CoVertex>(bList.size(), new ValenceComparator());
		lowValenceFirst.addAll(bList);
		for (CoVertex v : lowValenceFirst) {
			if (!v.isValid()) continue;
			List<CoFace> faces = HalfEdgeUtils.facesIncidentWithVertex(v);
			for (CoFace f : faces) {
				if (counter >= maxFaces) {
					System.out.println(counter + " faces reglued");
					return reglued;
				}
				if (HalfEdgeUtils.isInteriorFace(f)) continue;
				if (!isFaceMovable(f)) continue;
				if (!isOutsideFundamentalPolygon(f, p)) continue;
				if (!isFaceMovedToFundamentalDomainByReglue(f, cutInfo, p, signature))continue;
				try {
					reglueFace(f, cutInfo, signature);
				} catch (AssertionError e) {
					System.err.println(e.getLocalizedMessage());
				}
				reglued.add(f);
				System.out.println("face " + f + " reglued.");
				counter++;
	//			assert HalfEdgeUtils.isValidSurface(hds, true) : "surface should be valid after face reglue";
			}
		}
		System.out.println(counter + " faces reglued");
		return reglued;
	}
	
	
	
	protected static boolean isFaceMovable(CoFace f) {
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
			if (HalfEdgeUtils.isBoundaryEdge(e)) continue;
			CoVertex v = e.getNextEdge().getTargetVertex();
			boolean ear = HalfEdgeUtils.incomingEdges(v).size() == 2;
			boolean sb = HalfEdgeUtils.isBoundaryVertex(e.getStartVertex());
			boolean tb = HalfEdgeUtils.isBoundaryVertex(e.getTargetVertex());
			if (sb && tb && !ear) return false;
		}
		return true;
	}
	
	
	
	protected static void reglueFace(
		CoFace f,
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo,
		int signature
	) {
		if (HalfEdgeUtils.isInteriorFace(f)) {
			throw new IllegalArgumentException("can only reglue boundary faces");
		}
		HalfEdgeDataStructure<CoVertex, CoEdge, CoFace> hds = f.getHalfEdgeDataStructure();
		CoEdge e1 = f.getBoundaryEdge();
		CoEdge e2 = e1.getNextEdge();
		CoEdge e3 = e2.getNextEdge();
		int numBoundaryEdges = 0;
		numBoundaryEdges += isBoundaryEdge(e1) ? 1 : 0;
		numBoundaryEdges += isBoundaryEdge(e2) ? 1 : 0;
		numBoundaryEdges += isBoundaryEdge(e3) ? 1 : 0;
		assert numBoundaryEdges < 3 : "a face must have at most two boundary edges";
		if (numBoundaryEdges == 1) {
			CoEdge e = null;
			for (CoEdge be : HalfEdgeUtils.boundaryEdges(f)) {
				if (be.getRightFace() == null) {
					e = be; break;
				}
			}
			assert e != null : "we should have found the boundary edge of face f";
			
			// treat source location
			CoEdge opp = e.getOppositeEdge();
			CoEdge oppNext = opp.getNextEdge();
			CoEdge oppPrev = opp.getPreviousEdge();
			CoEdge eNext = e.getNextEdge();
			CoEdge ePrev = e.getPreviousEdge();
			CoVertex sourceVertex = eNext.getTargetVertex();
			// vertex geometry
			double[] sourcePosT = sourceVertex.T;
			double[] sourcePosP = sourceVertex.P;
			Matrix A = new Matrix(createIsometryFromEdge(e, cutInfo, signature));
			double[] newPos = A.multiplyVector(sourcePosT);
			CoVertex newVertex = hds.addNewVertex();
			newVertex.T = newPos;
			newVertex.P = sourcePosP.clone();
			if (sourceVertex.info != null) {
				newVertex.info = new CustomVertexInfo(sourceVertex.info);
			}
			// relink boundary
			oppPrev.linkNextEdge(eNext);
			ePrev.linkNextEdge(oppNext);
			eNext.setLeftFace(null);
			ePrev.setLeftFace(null);
			// remove old boundary
			hds.removeEdge(e);
			hds.removeEdge(opp);
			// fix incoming edges
			oppPrev.setTargetVertex(eNext.getStartVertex());
			ePrev.setTargetVertex(oppNext.getStartVertex());
			
			// target location
			CoEdge pe = cutInfo.edgeCutMap.get(e);
			CoEdge peOpp = pe.getOppositeEdge();
			CoEdge peOppNext = peOpp.getNextEdge();
			CoEdge peOppPrev = peOpp.getPreviousEdge();
			CoEdge newEdge1 = hds.addNewEdge();
			CoEdge newEdge1Opp = hds.addNewEdge();
			CoEdge newEdge2 = hds.addNewEdge();			
			CoEdge newEdge2Opp = hds.addNewEdge();			
			// link boundary
			peOppPrev.linkNextEdge(newEdge1Opp);
			newEdge1Opp.linkNextEdge(newEdge2Opp);
			newEdge2Opp.linkNextEdge(peOppNext);
			peOpp.linkNextEdge(newEdge2);
			newEdge2.linkNextEdge(newEdge1);
			newEdge1.linkNextEdge(peOpp);
			newEdge2.setTargetVertex(newVertex);
			newEdge1Opp.setTargetVertex(newVertex);
			newEdge1.setTargetVertex(pe.getTargetVertex());
			newEdge2Opp.setTargetVertex(pe.getStartVertex());
			newEdge1.setLeftFace(f);
			newEdge2.setLeftFace(f);
			peOpp.setLeftFace(f);
			newEdge1.linkOppositeEdge(newEdge1Opp);
			newEdge2.linkOppositeEdge(newEdge2Opp);
			
			// fix cut info
			cutInfo.edgeCutMap.put(ePrev, newEdge1Opp);
			cutInfo.edgeCutMap.put(newEdge1Opp, ePrev);
			cutInfo.edgeCutMap.put(eNext, newEdge2Opp);
			cutInfo.edgeCutMap.put(newEdge2Opp, eNext);
			cutInfo.edgeCutMap.put(ePrev.getOppositeEdge(), newEdge1);
			cutInfo.edgeCutMap.put(newEdge1, ePrev.getOppositeEdge());
			cutInfo.edgeCutMap.put(eNext.getOppositeEdge(), newEdge2);
			cutInfo.edgeCutMap.put(newEdge2, eNext.getOppositeEdge());
			cutInfo.edgeCutMap.remove(e);
			cutInfo.edgeCutMap.remove(opp);
			cutInfo.edgeCutMap.remove(pe);
			cutInfo.edgeCutMap.remove(peOpp);
			cutInfo.vertexCopyMap.put(sourceVertex, newVertex);
		} else {
			CoEdge e = null;
			for (CoEdge be : HalfEdgeUtils.boundaryEdges(f)) {
				if (be.getRightFace() != null) {
					e = be; break;
				}
			}
			assert e != null : "we should have found the non-boundary edge of face f";
			CoEdge b1 = e.getNextEdge();
			CoEdge b2 = e.getPreviousEdge();
			CoEdge pb1 = cutInfo.edgeCutMap.get(b1);
			CoEdge pb2 = cutInfo.edgeCutMap.get(b2);
			if (pb1.getOppositeEdge().getNextEdge() != pb2.getOppositeEdge()) {
				System.out.println("no continuous edge identification at face " + f);
				return;
			}
			CoEdge b1Opp = b1.getOppositeEdge();
			CoEdge b2Opp = b2.getOppositeEdge();
			CoEdge b1OppNext = b1Opp.getNextEdge();
			CoEdge b2OppPrev = b2Opp.getPreviousEdge();
			CoVertex v = b1.getTargetVertex();
			CoVertex v1 = b1.getStartVertex();
			CoVertex v2 = e.getStartVertex();
			// relink source boundary
			b2OppPrev.linkNextEdge(e);
			e.linkNextEdge(b1OppNext);
			e.setLeftFace(null);
			// delete old boundary
			hds.removeEdge(b1);
			hds.removeEdge(b2);
			hds.removeEdge(b1Opp);
			hds.removeEdge(b2Opp);
			hds.removeVertex(v);
			// fix incoming edges
			b2OppPrev.setTargetVertex(v2);
			e.setTargetVertex(v1);
			
			// target location
			CoEdge newEdge = hds.addNewEdge();
			CoEdge newEdgeOpp = hds.addNewEdge();
			CoEdge pb1Opp = pb1.getOppositeEdge();
			CoEdge pb2Opp = pb2.getOppositeEdge();
			CoEdge pb1OppPrev = pb1Opp.getPreviousEdge();
			CoEdge pb2OppNext = pb2Opp.getNextEdge();
			// link new boundary
			pb1OppPrev.linkNextEdge(newEdgeOpp);
			newEdgeOpp.linkNextEdge(pb2OppNext);
			// relink face
			newEdge.linkOppositeEdge(newEdgeOpp);
			pb2Opp.linkNextEdge(newEdge);
			newEdge.linkNextEdge(pb1Opp);
			newEdge.setLeftFace(f);
			pb1Opp.setLeftFace(f);
			pb2Opp.setLeftFace(f);
			newEdge.setTargetVertex(pb1.getTargetVertex());
			newEdgeOpp.setTargetVertex(pb2.getStartVertex());
			
			// fix cut info
			cutInfo.edgeCutMap.put(e, newEdgeOpp);
			cutInfo.edgeCutMap.put(newEdgeOpp, e);
			cutInfo.edgeCutMap.put(e.getOppositeEdge(), newEdge);
			cutInfo.edgeCutMap.put(newEdge, e.getOppositeEdge());
			cutInfo.edgeCutMap.remove(b1);
			cutInfo.edgeCutMap.remove(b2);
			cutInfo.edgeCutMap.remove(b1Opp);
			cutInfo.edgeCutMap.remove(b1Opp);
			cutInfo.vertexCopyMap.remove(v);
		}
	}
	

	protected static boolean isOutsideFundamentalPolygon(CoFace f, FundamentalPolygon p) {
		for (CoVertex v : HalfEdgeUtils.boundaryVertices(f)) {
			if (isInsideFundamentalPolygon(v, p)) return false;
		}
		return true;
	}
	
	protected static boolean isInsideFundamentalPolygon(CoVertex v, FundamentalPolygon p) {
		double[] vt = P2.projectP3ToP2(null, v.T); 
		for (FundamentalEdge e : p.getEdges()) {
			double[] s = P2.projectP3ToP2(null, RnBig.toDouble(null, e.startPosition));
			double[] t = P2.projectP3ToP2(null, RnBig.toDouble(null, e.nextEdge.startPosition));
			double[] line = P2.lineFromPoints(null, s, t);
			double dot = Rn.innerProduct(line, vt);
			if (dot < 0) return false;
		}
		return true;
	}
	
	protected static boolean isFaceMovedToFundamentalDomainByReglue(
		CoFace f, 
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, 
		FundamentalPolygon p, 
		int signature
	) {
		for (CoEdge e : HalfEdgeUtils.boundaryEdges(f)) {
			CoEdge pe = cutInfo.edgeCutMap.get(e);
			if (pe != null) {
				if (!isInsideFundamentalPolygon(pe.getStartVertex(), p)) return false;
				if (!isInsideFundamentalPolygon(pe.getTargetVertex(), p)) return false;
			}
		}
		return true;
	}
	
	
	protected static double[] createIsometryFromEdge(CoEdge e, CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo, int signature) {
		if (!HalfEdgeUtils.isBoundaryEdge(e)) {
			throw new IllegalArgumentException("No boundary edge");
		}
		CoEdge ePartner = cutInfo.edgeCutMap.get(e);
		assert ePartner != null : "every boundary edge has to be in the cut info";
		double[] s1 = P2.projectP3ToP2(null, e.getStartVertex().T);
		double[] t1 = P2.projectP3ToP2(null, e.getTargetVertex().T);
		double[] s2 = P2.projectP3ToP2(null, ePartner.getStartVertex().T);
		double[] t2 = P2.projectP3ToP2(null, ePartner.getTargetVertex().T);
		double[] a = P2.makeDirectIsometryFromFrames(null, s1, t1, t2, s2, signature);
		return P2.imbedMatrixP2InP3(null, a);
	}
	
}
