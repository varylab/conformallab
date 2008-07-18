package de.varylab.discreteconformal.heds;

import geom3d.Point;
import geom3d.Vector;
import de.jtem.halfedge.Edge;
import de.varylab.discreteconformal.heds.util.MeshUtility;

public class CEdge extends Edge<CVertex, CEdge, CFace> {


    private static final long
        serialVersionUID = 1L;

    private Double
    	curvature = null;
    private double
    	lambda = 1.0;
    
    
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


	public double getLambda() {
		return lambda;
	}
	
	public void setLambda(double lambda) {
		this.lambda = lambda;
	}


}
