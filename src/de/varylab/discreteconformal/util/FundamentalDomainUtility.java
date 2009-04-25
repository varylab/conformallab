package de.varylab.discreteconformal.util;

import static de.jreality.math.Pn.HYPERBOLIC;
import static java.awt.Color.BLACK;
import static java.awt.Color.WHITE;
import static java.awt.RenderingHints.KEY_ANTIALIASING;
import static java.awt.RenderingHints.KEY_RENDERING;
import static java.awt.RenderingHints.VALUE_ANTIALIAS_ON;
import static java.awt.RenderingHints.VALUE_RENDER_QUALITY;
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
import java.awt.geom.Arc2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.List;

import javax.imageio.ImageIO;

import de.jreality.math.Matrix;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.varylab.discreteconformal.util.UniformizationUtility.FundamentalEdge;
import de.varylab.discreteconformal.util.UniformizationUtility.FundamentalPolygon;

public class FundamentalDomainUtility {

	public static Image createCoverTexture(
		double[] root,
		FundamentalPolygon poly, 
		int depth
	) {
		int res = 1000;
		BufferedImage image = new BufferedImage(res, res, TYPE_INT_ARGB);
		Graphics2D g = (Graphics2D)image.getGraphics();
		g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
		g.setRenderingHint(KEY_RENDERING, VALUE_RENDER_QUALITY);
		g.setColor(new Color(0, 0, 0, 0));
		g.fillRect(0, 0, res, res);
		g.setColor(WHITE);
		g.fillArc(0, 0, res, res, 0, 360);
		g.setRenderingHint(KEY_ANTIALIASING, VALUE_ANTIALIAS_ON);
		List<double[]> orbit = poly.getOrbit(root);
		drawArc(poly, orbit, new Matrix(), g, res, 0, depth);
		try {
			ImageIO.write(image, "png", new File("test.png"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		return image;
	}
	
	
	private static void drawArc(
		FundamentalPolygon poly,
		List<double[]> orbit,
		Matrix domain,
		Graphics2D g,
		int res,
		int depth,
		int maxDepth
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
				if (getDistToUnitCircle(p1a) < 0.001 && 
					getDistToUnitCircle(p2a) < 0.001 &&
					getDistToUnitCircle(p3a) < 0.001) {
					proceed = drawArc = false;
				}
				if (drawArc) {
					drawArc(p1a, p2a, p3a, g, res);
				}
			}
			if (proceed) {
				drawArc(poly, orbit, T, g, res, depth + 1, maxDepth);
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
		int res
	) {
		Point p1 = new Point(p1a[0] / (p1a[3] + 1), p1a[1] / (p1a[3] + 1), 0);
		Point p2 = new Point(p2a[0] / (p2a[3] + 1), p2a[1] / (p2a[3] + 1), 0);
		Point p3 = new Point(p3a[0] / (p3a[3] + 1), p3a[1] / (p3a[3] + 1), 0);
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
			double radius = (c.getRadius() / 2) * 1000;
			g.setColor(BLACK);
			g.setStroke(new BasicStroke(3.0f));
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
