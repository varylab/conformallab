package de.varylab.discreteconformal.uniformization;

import static de.jreality.math.Pn.HYPERBOLIC;
import static de.varylab.discreteconformal.uniformization.UniformizationUtility.context;
import static java.awt.Color.BLACK;
import static java.awt.Color.WHITE;
import static java.awt.geom.Arc2D.OPEN;
import static java.awt.image.BufferedImage.TYPE_INT_ARGB;
import static java.lang.Math.cos;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.geom.Arc2D;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.List;

import de.jreality.geometry.IndexedLineSetFactory;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.scene.SceneGraphComponent;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.math.PnBig;
import de.varylab.discreteconformal.math.RnBig;

public class FundamentalDomainUtility {

	
	
	public static void createFundamentalPolygon(
		FundamentalPolygon fundamentalPolygon,
		double[] pRoot,
		SceneGraphComponent c, 
		int resolution, 
		HyperbolicModel model
	) {
		final IndexedLineSetFactory ilsf = new IndexedLineSetFactory();
		int n = fundamentalPolygon.getLength();
		BigDecimal[][] verts = null;
		int[][] edges = null;
		switch (model) {
			default:
			case Klein:
				ilsf.setVertexCount(n);
				ilsf.setEdgeCount(n);	
				verts = new BigDecimal[n][];
				edges = new int[n][];
				break;
			case Poincaré:
			case Halfplane:
				ilsf.setVertexCount(n * resolution);
				ilsf.setEdgeCount(n * resolution);
				verts = new BigDecimal[n * resolution][];
				edges = new int[n * resolution][];	
				break;
		}
		
		double[] root = new double[] {pRoot[0], pRoot[1], 0.0, pRoot[3]};
		List<BigDecimal[]> orbit = fundamentalPolygon.getOrbit(root);
		for (int i = 0; i < n; i++) {
			BigDecimal[] pos = orbit.get(i);
			BigDecimal[] posNext = orbit.get((i + 1) % n);
			switch (model) {
				default:
				case Klein:
					verts[i] = new BigDecimal[] {pos[0], pos[1], ZERO, pos[3]};
					edges[i] = new int[] {i, (i + 1) % (n)};
					break;
				case Poincaré:
					for (int j = 0; j < resolution; j++) {
						int index = i * resolution + j;
						double t = j / (resolution - 1.0);
						BigDecimal tBig = new BigDecimal(t);
						BigDecimal[] p = RnBig.linearCombination(null, tBig, posNext, ONE.subtract(tBig), pos, context);
						PnBig.normalize(p, p, Pn.HYPERBOLIC, context);
						verts[index] = new BigDecimal[] {p[0], p[1], ZERO, p[3].add(ONE)};
						edges[index] = new int[] {index, (index + 1) % (n * resolution)};
					}
					break;
				case Halfplane:
					for (int j = 0; j < resolution; j++) {
						int index = i * resolution + j;
						double t = j / (resolution - 1.0);
						BigDecimal tBig = new BigDecimal(t);
						BigDecimal[] p = RnBig.linearCombination(null, tBig, posNext, ONE.subtract(tBig), pos, context);
						PnBig.normalize(p, p, Pn.HYPERBOLIC, context);
						verts[index] = new BigDecimal[] {p[1], ONE, ZERO, p[3].subtract(p[0], context)};
						edges[index] = new int[] {index, (index + 1) % (n * resolution)};
					}
					break;
			}
		}
		ilsf.setVertexCoordinates(RnBig.toDouble(null, verts));
		ilsf.setEdgeIndices(edges);
		ilsf.update();
		c.setGeometry(ilsf.getGeometry());
	}
	
	
	
	
	
	public static Image createCoverTexture(
		double[] root,
		FundamentalPolygon poly, 
		int depth,
		HyperbolicModel model,
		int res
	) {
		BufferedImage image = new BufferedImage(res, res, TYPE_INT_ARGB);
		Graphics2D g = (Graphics2D)image.getGraphics();
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		g.setColor(new Color(0, 0, 0, 0));
		g.fillRect(0, 0, res, res);
		g.setColor(WHITE);
		g.setStroke(new BasicStroke(2.0f));
		
		switch (model) {
			case Klein:
			case Poincaré:
				g.fillArc(0, 0, res, res, 0, 360);
				break;
			case Halfplane:
				g.fillRect(0, 0, res, res);
				break;
		}
		
		List<BigDecimal[]> orbit = poly.getOrbit(root);
		g.setColor(BLACK);
		BigDecimal[] id = new BigDecimal[16];
		RnBig.setIdentityMatrix(id);
		drawPolygon(poly, orbit, id, g, res, 0, depth, model);
//		g.setColor(RED);
//		RnBig.setIdentityMatrix(id);
//		drawPolygon(poly, orbit, id, g, res, 0, 0, model);
//		try {
//			ImageIO.write(image, "png", new File("test.png"));
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		return image;
	}
	
	
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
		BigDecimal HALF = new BigDecimal(0.5);
		for (FundamentalEdge fe : poly.edgeList) {
			BigDecimal[] T = RnBig.times(null, domain, fe.motionBig, context);
			boolean proceed = true;
			for (int i = 0; i < orbit.size(); i++) {
				BigDecimal[] p1a = orbit.get(i).clone();
				BigDecimal[] p2a = orbit.get((i + 1) % orbit.size()).clone();
				BigDecimal[] p3a = RnBig.linearCombination(null, HALF, p1a, HALF, p2a, context);
				PnBig.normalize(p3a, p3a, HYPERBOLIC, context);
				RnBig.matrixTimesVector(p1a, T, p1a, context);
				RnBig.matrixTimesVector(p2a, T, p2a, context);
				RnBig.matrixTimesVector(p3a, T, p3a, context);
				PnBig.normalize(p1a, p1a, HYPERBOLIC, context);
				PnBig.normalize(p2a, p2a, HYPERBOLIC, context);
				PnBig.normalize(p3a, p3a, HYPERBOLIC, context);
				boolean drawArc = true;
				switch (model) {
				case Klein:
				case Poincaré:
					if (getDistToUnitCircle(p1a, context) < eps && 
						getDistToUnitCircle(p2a, context) < eps &&
						getDistToUnitCircle(p3a, context) < eps) {
						proceed = drawArc = false;
					}
					break;
				case Halfplane:
					if (1 / (p1a[3].subtract(p1a[0], context).doubleValue()) < eps ||
						1 / (p2a[3].subtract(p2a[0], context).doubleValue()) < eps ||
						1 / (p3a[3].subtract(p3a[0], context).doubleValue()) < eps) {
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
		double x = p[0].divide(p[3].add(BigDecimal.ONE, context), c).doubleValue();
		double y = p[1].divide(p[3].add(BigDecimal.ONE, context), c).doubleValue();
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
	public static double[] getCircumCenter(double[] A, double[] B, double[] C) {
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
