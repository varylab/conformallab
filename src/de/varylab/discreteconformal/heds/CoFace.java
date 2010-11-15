package de.varylab.discreteconformal.heds;

import java.util.List;

import de.jreality.math.Pn;
import de.jreality.math.Rn;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.varylab.discreteconformal.functional.node.ConformalFace;

public class CoFace extends ConformalFace<CoVertex, CoEdge, CoFace> {

    private double[]
    	normal = null;
    
    public double[] getNormal() {
    	if (normal == null) {
    		List<CoEdge> boundary = HalfEdgeUtils.boundaryEdges(this);
    		double[] a = boundary.get(0).getTargetVertex().P;
    		double[] b = boundary.get(1).getTargetVertex().P;
    		double[] c = boundary.get(2).getTargetVertex().P;
    		Pn.dehomogenize(a, a);
    		Pn.dehomogenize(b, b);
    		Pn.dehomogenize(c, c);
    		double[] v1 = Rn.subtract(null, b, a);
    		double[] v2 = Rn.subtract(null, c, a);
    		normal = Rn.crossProduct(null, v1, v2);
		}
    	return normal;
	}


}
