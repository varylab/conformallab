package de.varylab.discreteconformal.heds;

import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryEdges;
import static java.lang.Math.PI;
import static java.lang.Math.exp;
import geom3d.Point;
import geom3d.Triangle;
import geom3d.Vector;

import java.util.List;

import de.jtem.halfedge.Face;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.math.Lob;

public class CFace extends Face<CVertex, CEdge, CFace> {

    private Vector
    	normal = null;
    private Double
    	curvature = null;
    
    public double getFaceEneergy() {
    	CEdge e1 = getBoundaryEdge();
    	CEdge e2 = e1.getNextEdge();
    	CEdge e3 = e2.getNextEdge();
    	
    	double u1 = e2.getTargetVertex().getU();
    	double u2 = e3.getTargetVertex().getU();
    	double u3 = e1.getTargetVertex().getU();
    	
    	double la1 = e1.getLambda();
    	double la2 = e2.getLambda();
    	double la3 = e3.getLambda();
    	
    	double l1 = exp(la1 / 2);
    	double l2 = exp(la2 / 2);
    	double l3 = exp(la3 / 2);
    	
    	double a1 = 2 * Math.atan2((l2 + l3 - l1) * (l1 + l2 + l3), (l1 + l2 - l3) * (l1 + l3 - l2));
    	double a2 = 2 * Math.atan2((l1 + l3 - l2) * (l1 + l2 + l3), (l2 + l3 - l1) * (l1 + l2 - l3));
    	double a3 = 2 * Math.atan2((l1 + l2 - l3) * (l1 + l2 + l3), (l2 + l3 - l1) * (l1 + l3 - l2));
    		
		double f1 = 0.5 * (la1*a1 + la2*a2 + la3*a3);
		double f2 = Lob.valueAt(a1) + Lob.valueAt(a2) + Lob.valueAt(a3);
		
		return f1 + f2 - 0.5 * PI * (u1 + u2 + u3);
    }

    
    
    public double getCurvature() {
    	if (curvature == null) {
    		List<CEdge> b = boundaryEdges(this);
    		curvature = 0.0;
    		for (CEdge e : b)
    			curvature += e.getAngle();
    		curvature /= b.size();
		}
		return curvature;
    }
    
    
    public Triangle toTriangle() {
    	List<CEdge> b = HalfEdgeUtils.boundaryEdges(this);
    	if (b.size() != 3)
    		throw new RuntimeException("No triangle in toTriangle()");
    	return new Triangle(b.get(0).getTargetVertex().getPosition(), b.get(1).getTargetVertex().getPosition(), b.get(2).getTargetVertex().getPosition());
    }
    
    
    public Vector getNormal() {
    	if (normal == null) {
    		List<CEdge> boundary = HalfEdgeUtils.boundaryEdges(this);
    		Point a = boundary.get(0).getTargetVertex().getPosition();
    		Point b = boundary.get(1).getTargetVertex().getPosition();
    		Point c = boundary.get(2).getTargetVertex().getPosition();
    		Vector v1 = a.vectorTo(b);
    		Vector v2 = a.vectorTo(c);
    		normal = v1.getNormal(v2);
		}
    	return normal;
	}
}
