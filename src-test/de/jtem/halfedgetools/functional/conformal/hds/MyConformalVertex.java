package de.jtem.halfedgetools.functional.conformal.hds;

import javax.vecmath.Point3d;

import de.jtem.halfedgetools.functional.HasPosition;
import de.varylab.discreteconformal.functional.node.ConformalVertex;


public class MyConformalVertex extends ConformalVertex<MyConformalVertex, MyConformalEdge, MyConformalFace> implements HasPosition {

	private Point3d
		pos = new Point3d();
	
	public void setPosition(Point3d p) {
		pos.set(p);
	}
	
	public Point3d getPosition() {
		return pos;
	}
	
}
