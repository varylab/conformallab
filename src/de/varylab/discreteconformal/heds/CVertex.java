package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import geom3d.Point;
import geom3d.Vector;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.heds.decoration.HasNormal;
import de.varylab.discreteconformal.heds.decoration.HasPosition;

public class CVertex extends Vertex<CVertex, CEdge, CFace> implements HasPosition, HasNormal{

	private Point
	    P = new Point();
	private Vector
		N = new Vector();
	private double 
		u = 1.0,
		theta = 2 * PI;

	
	public double getU() {
		return u;
	}

	public void setU(double u) {
		this.u = u;
	}

	public Point getPosition() {
		return P;
	}

	public void setPosition(Point p) {
		P.set(P);
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


}
