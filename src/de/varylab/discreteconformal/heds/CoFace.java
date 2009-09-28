package de.varylab.discreteconformal.heds;

import static de.jtem.halfedge.util.HalfEdgeUtils.boundaryEdges;
import geom3d.Point;
import geom3d.Triangle;
import geom3d.Vector;

import java.util.List;

import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.functional.conformal.node.ConformalFace;

public class CoFace extends ConformalFace<CoVertex, CoEdge, CoFace> {

    private Vector
    	normal = null;
    private Double
    	curvature = null;
    
    public double getCurvature() {
    	if (curvature == null) {
    		List<CoEdge> b = boundaryEdges(this);
    		curvature = 0.0;
    		for (CoEdge e : b)
    			curvature += e.getCurvature();
    		curvature /= b.size();
		}
		return curvature;
    }
    
    
    public Triangle toTriangle() {
    	List<CoEdge> b = HalfEdgeUtils.boundaryEdges(this);
    	if (b.size() != 3)
    		throw new RuntimeException("No triangle in toTriangle()");
    	return new Triangle(b.get(0).getTargetVertex().getPosition(), b.get(1).getTargetVertex().getPosition(), b.get(2).getTargetVertex().getPosition());
    }
    
    
    public Vector getNormal() {
    	if (normal == null) {
    		List<CoEdge> boundary = HalfEdgeUtils.boundaryEdges(this);
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
