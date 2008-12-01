package de.varylab.discreteconformal.heds;

import geom3d.Point;
import geom3d.Vector;
import de.varylab.discreteconformal.heds.util.MeshUtility;

public class CoEdge extends de.jtem.halfedge.functional.conformal.node.ConformalEdge<CoVertex, CoEdge, CoFace> {

    private Double
    	curvature = null;
    
    public double getCurvature() {
    	if (curvature == null) {
			curvature = MeshUtility.getCurvature(this);
		}
		return curvature;
    }

    
    public double getLength(){
		Point start = getStartVertex().getPosition();
		Point target = getTargetVertex().getPosition();
		return start.distanceTo(target);
	}
    
    public Vector getVector(){
    	Point start = getStartVertex().getPosition();
		Point target = getTargetVertex().getPosition();
		return start.vectorTo(target);
    	
    }
	
}
