package de.varylab.discreteconformal.heds;

import geom3d.Point;
import geom3d.Vector;
import de.jreality.math.Matrix;
import de.jreality.math.MatrixBuilder;
import de.jtem.halfedge.Edge;

public class CEdge extends Edge<CVertex, CEdge, CFace> {


    private static final long
        serialVersionUID = 1L;

    private Double
    	curvature = null;
    
    
    public double getAngle() {
    	if (curvature == null) {
			curvature = signedAngle();
		}
		return curvature;
    }
    
    
	public double signedAngle() {
		CFace lf = getLeftFace();
		CFace rf = getRightFace();
		if (lf == null || rf == null)
			return 0;
		
		return curvatureSign()*lf.getNormal().getAngle(rf.getNormal());
	}
	
	/*
	 * 
	 * @param e, an MEdge
	 * @return -1,0,1 the sign of the angle between the left and the right face.
	 * 			negative is concave, positive if convex
	 *          
	 */
	private double curvatureSign(){
		Matrix m = MatrixBuilder.euclidean().getMatrix();
		for (int i = 0; i < 3; i++) {
			m.setEntry(i, 0, getVector().get(i));
			m.setEntry(i, 1, getLeftFace().getNormal().get(i));
			m.setEntry(i, 2, getRightFace().getNormal().get(i));
		}
		double det = m.getDeterminant() ;
		if(Math.abs(det)< 1E-7)
			return 0;
		else if (det<0)
			return -1;
		else 
			return 1;
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
