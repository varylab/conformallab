package de.varylab.discreteconformal.util;

import static de.jreality.math.Pn.HYPERBOLIC;
import static java.awt.Color.BLACK;
import static java.awt.Color.WHITE;
import static java.awt.geom.Arc2D.OPEN;
import static java.awt.image.BufferedImage.TYPE_INT_ARGB;
import static java.lang.Math.cos;
import static java.lang.Math.signum;
import static java.lang.Math.sin;
import geom3d.Circle;
import geom3d.Point;
import geom3d.Triangle;
import geom3d.Vector;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.RenderingHints;
import java.awt.geom.Arc2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;

import javax.imageio.ImageIO;

import de.jreality.geometry.IndexedLineSetFactory;
import de.jreality.math.Matrix;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jreality.scene.SceneGraphComponent;
import de.varylab.discreteconformal.adapter.HyperbolicModel;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;
import de.varylab.discreteconformal.util.CuttingUtility.CuttingInfo;
import de.varylab.discreteconformal.util.UniformizationUtility.FundamentalEdge;
import de.varylab.discreteconformal.util.UniformizationUtility.FundamentalPolygon;

public class FundamentalDomainUtility {

	
	
	public static void createFundamentalPolygon(
		FundamentalPolygon fundamentalPolygon,
		CuttingInfo<CoVertex, CoEdge, CoFace> cutInfo,
		SceneGraphComponent c, 
		int resolution, 
		HyperbolicModel model
	) {
		if (fundamentalPolygon == null || cutInfo == null) {
			return;
		}
		final IndexedLineSetFactory ilsf = new IndexedLineSetFactory();
		int n = fundamentalPolygon.getLength();
		double[][] verts = null;
		int[][] edges = null;
		switch (model) {
			default:
			case Klein:
				ilsf.setVertexCount(n);
				ilsf.setEdgeCount(n);	
				verts = new double[n][];
				edges = new int[n][];
				break;
			case Poincaré:
			case Halfplane:
				ilsf.setVertexCount(n * resolution);
				ilsf.setEdgeCount(n * resolution);
				verts = new double[n * resolution][];
				edges = new int[n * resolution][];	
				break;
		}
		Point pRoot = cutInfo.cutRoot.getTextureCoord();
		
		double[] root = new double[] {pRoot.x(), pRoot.y(), 0.0, pRoot.z()};
		List<double[]> orbit = fundamentalPolygon.getOrbit(root);
		for (int i = 0; i < n; i++) {
			double[] pos = orbit.get(i);
			double[] posNext = orbit.get((i + 1) % n);
			switch (model) {
				default:
				case Klein:
					verts[i] = new double[] {pos[0], pos[1], 0.0, pos[3]};
					edges[i] = new int[] {i, (i + 1) % (n)};
					break;
				case Poincaré:
					for (int j = 0; j < resolution; j++) {
						int index = i * resolution + j;
						double t = j / (resolution - 1.0);
						double[] p = Rn.linearCombination(null, t, posNext, 1-t, pos);
						Pn.normalize(p, p, Pn.HYPERBOLIC);
						verts[index] = new double[] {p[0], p[1], 0.0, p[3] + 1};
						edges[index] = new int[] {index, (index + 1) % (n * resolution)};
					}
					break;
				case Halfplane:
					for (int j = 0; j < resolution; j++) {
						int index = i * resolution + j;
						double t = j / (resolution - 1.0);
						double[] p = Rn.linearCombination(null, t, posNext, 1-t, pos);
						Pn.normalize(p, p, Pn.HYPERBOLIC);
						verts[index] = new double[] {p[1], 1.0, 0.0, (p[3] - p[0])};
						edges[index] = new int[] {index, (index + 1) % (n * resolution)};
					}
					break;
			}
		}
		ilsf.setVertexCoordinates(verts);
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
		
		List<double[]> orbit = poly.getOrbit(root);
		g.setColor(BLACK);
		drawPolygon(poly, orbit, new Matrix(), g, res, 0, depth, model);
		try {
			ImageIO.write(image, "png", new File("test.png"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		return image;
	}
	
	
	private static void drawPolygon(
		FundamentalPolygon poly,
		List<double[]> orbit,
		Matrix domain,
		Graphics2D g,
		int res,
		int depth,
		int maxDepth,
		HyperbolicModel model
	) {
		if (depth > maxDepth) {
			return;
		}
		for (FundamentalEdge fe : poly.edgeList) {
			Matrix T = Matrix.times(domain, fe.motion);
			boolean proceed = true;
			for (int i = 0; i < orbit.size(); i++) {
				double[] p1a = orbit.get(i).clone();
				double[] p2a = orbit.get((i + 1) % orbit.size()).clone();
				double[] p3a = Rn.linearCombination(null, 0.5, p1a, 0.5, p2a);
				Pn.normalize(p3a, p3a, HYPERBOLIC);
				T.transformVector(p1a);
				T.transformVector(p2a);
				T.transformVector(p3a);
				Pn.normalize(p1a, p1a, HYPERBOLIC);
				Pn.normalize(p2a, p2a, HYPERBOLIC);
				Pn.normalize(p3a, p3a, HYPERBOLIC);
				boolean drawArc = true;
				switch (model) {
				case Klein:
				case Poincaré:
					if (getDistToUnitCircle(p1a) < 0.001 && 
						getDistToUnitCircle(p2a) < 0.001 &&
						getDistToUnitCircle(p3a) < 0.001) {
						proceed = drawArc = false;
					}
					break;
				case Halfplane:
					if (1 / (p1a[3] - p1a[0]) < 0.0001 ||
						1 / (p2a[3] - p2a[0]) < 0.0001 ||
						1 / (p3a[3] - p3a[0]) < 0.0001) {
						proceed = drawArc = false;
					}
				}
				if (drawArc) {
					drawArc(p1a, p2a, p3a, g, res, model);
				}
			}
			if (proceed) {
				drawPolygon(poly, orbit, T, g, res, depth + 1, maxDepth, model);
			}
		}
	}
	
	private static double getDistToUnitCircle(double[] p) {
		double x = p[0] / (p[3] + 1);
		double y = p[1] / (p[3] + 1);
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
		Point p1 = null;
		Point p2 = null;
		Point p3 = null;
		switch (model) {
			case Klein:
				p1 = new Point(p1a[0] / p1a[3], p1a[1] / p1a[3], 0);
				p2 = new Point(p2a[0] / p2a[3], p2a[1] / p2a[3], 0);
				p3 = new Point(p3a[0] / p3a[3], p3a[1] / p3a[3], 0);		
				break;
				default:
			case Poincaré:
				p1 = new Point(p1a[0] / (p1a[3] + 1), p1a[1] / (p1a[3] + 1), 0);
				p2 = new Point(p2a[0] / (p2a[3] + 1), p2a[1] / (p2a[3] + 1), 0);
				p3 = new Point(p3a[0] / (p3a[3] + 1), p3a[1] / (p3a[3] + 1), 0);
				break;
			case Halfplane:
				p1 = new Point(p1a[1] / (p1a[3] - p1a[0]), 1 / (p1a[3] - p1a[0]), 0);
				p2 = new Point(p2a[1] / (p2a[3] - p2a[0]), 1 / (p2a[3] - p2a[0]), 0);
				p3 = new Point(p3a[1] / (p3a[3] - p3a[0]), 1 / (p3a[3] - p3a[0]), 0);
		}

		Triangle t = new Triangle(p1, p2, p3);
		try {
			Circle c = t.getCircumCircle();
			Point center = c.getCenter();
			Vector vec1 = center.vectorTo(p1);
			Vector vec2 = center.vectorTo(p2);
			double angle = vec1.getAngle(vec2);
			double startAngle = Math.atan2(vec1.y(), vec1.x());
			angle *= signum(-vec2.x()*sin(startAngle) + vec2.y()*cos(startAngle));
			
			double degAngle = Math.toDegrees(angle);
			double degStartAngle = Math.toDegrees(startAngle);
			double centerX = (center.x() / 2 + 0.5) * res;
			double centerY = (-center.y() / 2 + 0.5) * res;
			double radius = (c.getRadius() / 2) * res;
			Arc2D arc = new Arc2D.Double(centerX - radius, centerY - radius, radius * 2, radius * 2, degStartAngle, degAngle, OPEN);
			g.draw(arc);
			
		} catch (Exception e) {
			int p1x = (int)((p1.x() / 2 + 0.5) * res);
			int p1y = (int)((-p1.y() / 2 + 0.5) * res);
			int p2x = (int)((p2.x() / 2 + 0.5) * res);
			int p2y = (int)((-p2.y() / 2 + 0.5) * res);
			g.drawLine(p1x, p1y, p2x, p2y);
		}
		
	}
	
	
	
	
}