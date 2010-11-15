package de.varylab.discreteconformal.heds;

import de.jreality.math.Pn;
import de.varylab.discreteconformal.functional.node.ConformalEdge;

public class CoEdge extends ConformalEdge<CoVertex, CoEdge, CoFace> {

	public double getLength() {
		double[] s = getStartVertex().P;
		double[] t = getTargetVertex().P;
		return Pn.distanceBetween(s, t, Pn.EUCLIDEAN);
	}
	
}
