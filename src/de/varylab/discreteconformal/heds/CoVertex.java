package de.varylab.discreteconformal.heds;

import geom3d.Point;
import geom3d.Vector;
import de.varylab.discreteconformal.functional.node.ConformalVertex;
import de.varylab.discreteconformal.heds.bsp.HasBspPos;

public class CoVertex extends ConformalVertex<CoVertex, CoEdge, CoFace> implements HasBspPos {

	private Point
	    P = new Point(),
	    T = new Point();
	private Vector
		N = new Vector();
	private CustomVertexInfo
		info = null;

	public Point getPosition() {
		return P;
	}
	public void setPosition(Point p) {
		P.set(p);
	}
	
	@Override
	public Point getBspPos() {
		return T;
	}
	@Override
	public void setBspPos(Point p) {
		T.set(p);
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
	
	public CustomVertexInfo getCustomInfo() {
		return info;
	}
	public void setCustomInfo(CustomVertexInfo info) {
		this.info = info;
	}
	
	@Override
	public void copyData(CoVertex v) {
		super.copyData(v);
		setPosition(v.getPosition());
		setTextureCoord(v.getTextureCoord());
		setNormal(v.getNormal());
		setCustomInfo(v.getCustomInfo());
	}
	
}
