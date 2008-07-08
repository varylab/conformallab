package de.varylab.discreteconformal.heds;

import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryEdges;
import geom3d.Point;
import geom3d.Triangle;
import geom3d.Vector;

import java.util.List;

import de.jtem.halfedge.Face;
import de.jtem.halfedge.util.HalfEdgeUtils;

public class CFace extends Face<CVertex, CEdge, CFace> {

    private Vector
    	normal = null;
    private Double
    	curvature = null;
    
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
