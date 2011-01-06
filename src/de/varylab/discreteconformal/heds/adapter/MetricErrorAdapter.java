package de.varylab.discreteconformal.heds.adapter;

import java.util.HashMap;
import java.util.Map;

import de.jreality.math.Pn;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public class MetricErrorAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private Map<CoEdge, Double>
		lengthMap = new HashMap<CoEdge, Double>();
	private int
		signature = Pn.EUCLIDEAN;
	
	public MetricErrorAdapter() {
		super(null, CoEdge.class, null, Double.class, true, false);
	}

	public void setLengthMap(Map<CoEdge, Double> lengthMap) {
		this.lengthMap = lengthMap;
	}
	
	public void setSignature(int signature) {
		this.signature = signature;
	}
	
	@Override
	public Double getEdgeValue(CoEdge e, AdapterSet a) {
		double d = lengthMap.get(e) == null ? 0.0 : lengthMap.get(e);
		double[] s = e.getStartVertex().T;
		double[] t = e.getTargetVertex().T;
		double rd = Pn.distanceBetween(s, t, signature);
		return Math.abs(d - rd);
	}
	
}
