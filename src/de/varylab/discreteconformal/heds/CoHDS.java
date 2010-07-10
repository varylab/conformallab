package de.varylab.discreteconformal.heds;

import static java.lang.Double.MAX_VALUE;
import geom3d.Point;
import de.jtem.halfedge.HalfEdgeDataStructure;

public class CoHDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace> {

	private boolean
		texCoordinatesValid = false;
	private Point
		normalizeCenter = new Point(0, 0, 0);
	private double
		normalizeFactor = 1.0;
	
	public CoHDS() {
		super(CoVertex.class, CoEdge.class, CoFace.class);
	}
	
	
	public boolean isTexCoordinatesValid() {
		return texCoordinatesValid;
	}

	public void setTexCoordinatesValid(boolean texCoordinatesValid) {
		this.texCoordinatesValid = texCoordinatesValid;
	}
	
	
	public void normalizeCoordinates() {
		double[][] bounds = new double[][]{{MAX_VALUE, -MAX_VALUE},{MAX_VALUE, -MAX_VALUE},{MAX_VALUE, -MAX_VALUE}};
		for (CoVertex v : getVertices()) {
			Point p = v.getPosition();
			bounds[0][0] = p.x() < bounds[0][0] ? p.x() : bounds[0][0];
			bounds[0][1] = p.x() > bounds[0][1] ? p.x() : bounds[0][1];
			bounds[1][0] = p.y() < bounds[1][0] ? p.y() : bounds[1][0];
			bounds[1][1] = p.y() > bounds[1][1] ? p.y() : bounds[1][1];
			bounds[2][0] = p.z() < bounds[2][0] ? p.z() : bounds[2][0];
			bounds[2][1] = p.z() > bounds[2][1] ? p.z() : bounds[2][1];
		}
		double xExtend = Math.abs(bounds[0][1] - bounds[0][0]);
		double yExtend = Math.abs(bounds[1][1] - bounds[1][0]);
		double zExtend = Math.abs(bounds[2][1] - bounds[2][0]);
		double max = Math.max(xExtend, Math.max(yExtend, zExtend));
		double xCenter = (bounds[0][1] + bounds[0][0]) / 2 / max;
		double yCenter = (bounds[1][1] + bounds[1][0]) / 2 / max;
		double zCenter = (bounds[2][1] + bounds[2][0]) / 2 / max;
		normalizeCenter = new Point(-xCenter, -yCenter, -zCenter);
		normalizeFactor = 1 / max;
		for (CoVertex v : getVertices()) {
			Point p = v.getPosition();
			p.times(normalizeFactor);
			p.move(normalizeCenter);
		}
	}
	
	
	public void revertNormalization() {
		normalizeCenter.times(-1);
		for (CoVertex v : getVertices()) {
			Point p = v.getPosition();
			p.move(normalizeCenter);
			p.times(1 / normalizeFactor);
		}
		normalizeCenter = new Point(0, 0, 0);
		normalizeFactor = 1.0;
	}
	
}
