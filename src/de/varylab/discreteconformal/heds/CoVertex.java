package de.varylab.discreteconformal.heds;

import de.varylab.discreteconformal.functional.node.ConformalVertex;

public class CoVertex extends ConformalVertex<CoVertex, CoEdge, CoFace> {

	public double[]
	    P = {0,0,0,1},
	    T = {0,0,0,1};
	public CustomVertexInfo
		info = null;
	
	@Override
	public void copyData(CoVertex v) {
		super.copyData(v);
		P = v.P.clone();
		T = v.T.clone();
		if (v.info != null) {
			info = new CustomVertexInfo(v.info);
		}
	}
	
}
