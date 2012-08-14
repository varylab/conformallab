package de.varylab.discreteconformal.heds;

import de.jreality.math.Pn;
import de.varylab.discreteconformal.functional.node.ConformalEdge;

public class CoEdge extends ConformalEdge<CoVertex, CoEdge, CoFace> {

	public CustomEdgeInfo
		info = null;
	
	public double getLength() {
		double[] s = getStartVertex().P;
		double[] t = getTargetVertex().P;
		return Pn.distanceBetween(s, t, Pn.EUCLIDEAN);
	}
	
	public double getTexLength() {
		double[] s = getStartVertex().T;
		double[] t = getTargetVertex().T;
		return Pn.distanceBetween(s, t, Pn.EUCLIDEAN);
	}
	
	@Override
	public void copyData(CoEdge e) {
		super.copyData(e);
		if (e.info != null) {
			info = new CustomEdgeInfo(e.info);
		}
	}
	
}
