package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import geom3d.Point;
import geom3d.Vector;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.heds.bsp.HasPosition;

public class CVertex extends Vertex<CVertex, CEdge, CFace> implements HasPosition {

	private Point
	    P = new Point(),
	    T = new Point();
	private Vector
		N = new Vector();
	private double 
		theta = 2 * PI;
	private Integer
		solverIndex = -1;


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

	public double getTheta() {
		return theta;
	}

	public void setTheta(double theta) {
		this.theta = theta;
	}

	public Integer getSolverIndex() {
		return solverIndex;
	}

	public void setSolverIndex(Integer solverIndex) {
		this.solverIndex = solverIndex;
	}

	
}
