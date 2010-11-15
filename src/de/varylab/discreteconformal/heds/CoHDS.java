package de.varylab.discreteconformal.heds;

import static java.lang.Double.MAX_VALUE;
import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.HalfEdgeDataStructure;

public class CoHDS extends HalfEdgeDataStructure<CoVertex, CoEdge, CoFace> {

	private boolean
		texCoordinatesValid = false;
	private double[]
		normalizeCenter = {0,0,0,1};
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
			double[] p = v.P;
			Pn.dehomogenize(p, p);
			bounds[0][0] = p[0] < bounds[0][0] ? p[0] : bounds[0][0];
			bounds[0][1] = p[0] > bounds[0][1] ? p[0] : bounds[0][1];
			bounds[1][0] = p[1] < bounds[1][0] ? p[1] : bounds[1][0];
			bounds[1][1] = p[1] > bounds[1][1] ? p[1] : bounds[1][1];
			bounds[2][0] = p[2] < bounds[2][0] ? p[2] : bounds[2][0];
			bounds[2][1] = p[2] > bounds[2][1] ? p[2] : bounds[2][1];
		}
		double xExtend = Math.abs(bounds[0][1] - bounds[0][0]);
		double yExtend = Math.abs(bounds[1][1] - bounds[1][0]);
		double zExtend = Math.abs(bounds[2][1] - bounds[2][0]);
		double max = Math.max(xExtend, Math.max(yExtend, zExtend));
		double xCenter = (bounds[0][1] + bounds[0][0]) / 2 / max;
		double yCenter = (bounds[1][1] + bounds[1][0]) / 2 / max;
		double zCenter = (bounds[2][1] + bounds[2][0]) / 2 / max;
		normalizeCenter = new double[] {-xCenter, -yCenter, -zCenter};
		normalizeFactor = 1 / max;
		for (CoVertex v : getVertices()) {
			double[] p = v.P;
			Rn.times(p, normalizeFactor, p);
			Rn.add(p, normalizeCenter, p);
			p[3] = 1.0;
		}
	}
	
	
	public void revertNormalization() {
		Rn.times(normalizeCenter, -1, normalizeCenter);
		for (CoVertex v : getVertices()) {
			double[] p = v.P;
			Pn.dehomogenize(p, p);
			Rn.add(p, normalizeCenter, p);
			p[3] = normalizeFactor;
		}
		normalizeCenter = new double[] {0,0,0,1};
		normalizeFactor = 1.0;
	}
	
}
