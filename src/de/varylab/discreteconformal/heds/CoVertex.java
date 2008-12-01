package de.varylab.discreteconformal.heds;

import geom3d.Point;
import geom3d.Vector;
import de.jtem.halfedge.functional.conformal.node.ConformalVertex;
import de.varylab.discreteconformal.heds.bsp.HasPosition;

public class CoVertex extends ConformalVertex<CoVertex, CoEdge, CoFace> implements HasPosition {

	private Point
	    P = new Point(),
	    T = new Point();
	private Vector
		N = new Vector();

	public Point getPosition() {
		return P;
	}

	public void setPosition(Point p) {
		P.set(p);
	}

	public Point getTextureCoord() {
		return T;
	}
	
	public void setTextureCoord(Point t) {
		T.set(t);
	}
	
	public Vector getNormal() {
		return N;
	}

	public void setNormal(Vector normal) {
		N.set(normal);
	}

}
