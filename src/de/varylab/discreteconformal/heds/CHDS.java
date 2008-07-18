package de.varylab.discreteconformal.heds;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.log;
import geom3d.Point;
import geom3d.Vector;
import de.jtem.halfedge.HalfEdgeDataStructure;
import de.varylab.discreteconformal.math.Lob;

public class CHDS extends HalfEdgeDataStructure<CVertex, CEdge, CFace> {

	public CHDS() {
		super(CVertex.class, CEdge.class, CFace.class);
	}

	
	public double angleAt(Point a, Point b, Point c) {
		Vector ab = a.vectorTo(b);
		Vector ac = a.vectorTo(c);
		Vector cr = ab.cross(ac);
		double angle = Math.atan2(Math.sqrt(cr.dot(cr)), ab.dot(ac));
		if (!(0 <= angle && angle <= PI))
			throw new IllegalArgumentException("Illegal triangle in angleAt()");
		return angle;
	}
	
	
	public void prepareData(double[] theta) {
		for (CEdge e : getPositiveEdges()) {
			double l = e.getLength();
			e.setLambda(2.0 * log(l));
			e.getOppositeEdge().setLambda(e.getLambda());
		}
		for (CVertex v : getVertices()) {
			v.setTheta(theta[v.getIndex()]);
		}
		for (CFace f : getFaces()) {
			CEdge e1 = f.getBoundaryEdge();
			CEdge e2 = e1.getNextEdge();
			CEdge e3 = e1.getPreviousEdge();
			Point p1 = e1.getTargetVertex().getPosition();
			Point p2 = e2.getTargetVertex().getPosition();
			Point p3 = e3.getTargetVertex().getPosition();
			double a1 = angleAt(p2, p1, p3);
			double a2 = angleAt(p3, p2, p1);
			double a3 = angleAt(p1, p3, p2);
			double E1 = a1*e1.getLambda() + a2*e2.getLambda() + a3*e3.getLambda();
			double E2 = Lob.valueAt(a1) + Lob.valueAt(a2) + Lob.valueAt(a3);
			f.setEnergy(E1 + E2);
		}
	}
	
	
	public void writeUState(double[] u) {
		for (CVertex v : getVertices())
			v.setU(u[v.getIndex()]);
	}
	
	
	
	
	private double triangleEnergyAndAngles(CFace f, double[] a123) {
		CEdge e1 = f.getBoundaryEdge();
		CEdge e2 = e1.getNextEdge();
		CEdge e3 = e1.getPreviousEdge();
		CVertex v1 = e1.getTargetVertex();
		CVertex v2 = e2.getTargetVertex();
		CVertex v3 = e3.getTargetVertex();
		a123[0] = 0.0;
		a123[1] = 0.0;
		a123[2] = 0.0;
		double u1 = v1.getU();
		double u2 = v2.getU();
		double u3 = v3.getU();
		double umean = (u1+u2+u3)/3;
		double x12 = e2.getLambda() + (u1+u2 - 2*umean);
		double x23 = e3.getLambda() + (u2+u3 - 2*umean);
		double x31 = e1.getLambda() + (u3+u1 - 2*umean);
		double l12 = exp(x12);
		double l23 = exp(x23);
		double l31 = exp(x31);
		double t31 = +l12+l23-l31;
		double t23 = +l12-l23+l31;
		double t12 = -l12+l23+l31;
		if (t31 > 0 && t23 > 0 && t12 > 0) {
			double l123 = l12 + l23 + l31;
			double denom = Math.sqrt(t12 * t23 * t31 * l123);
			a123[0] = 2 * Math.atan2(t12 * t31, denom);
			a123[1] = 2 * Math.atan2(t23 * t12, denom);
			a123[2] = 2 * Math.atan2(t31 * t23, denom);
		} else if (t31 <= 0) {
			a123[1] = 0.0;
		} else if (t23 <= 0) {
			a123[0] = 0.0;
		} else if (t12 <= 0) {
			a123[2] = 0.0;
		}
		double E1 = a123[0]*x23 + a123[1]*x31 + a123[2]*x12;
		double E2 = Lob.valueAt(a123[0]) + Lob.valueAt(a123[1]) + Lob.valueAt(a123[2]);
		double E3 = - PI * umean - f.getEnergy();
		return E1 + E2 + E3; 
	}
	
	
	
	
	public void conformalEnergy(double[] u, double[] E, double[] G, double[][] H) {
		writeUState(u);
		// Vertex Energy
		if (E != null) 
			E[0] = 0.0;
		for (CVertex v : getVertices()) {
			if (E != null)
				E[0] += v.getTheta() * v.getU();
			if (G != null)
				G[v.getIndex()] = v.getTheta();
		}
		if (E != null)
			E[0] *= 0.5;
		// Face Energy
		double[] a123 = {0, 0, 0};
		for (CFace t : getFaces()) {
			double e = triangleEnergyAndAngles(t, a123);
			if (E != null)
				E[0] += e;
//			if (G != null)
//				G
			//:TODO Keep going
		}
		// Gradient
		
		
	}
	
	
	
	
	
}
