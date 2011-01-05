package de.varylab.discreteconformal.heds.adapter;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

public class MetricEdgeLengthAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Number> {

	private Map<CoEdge, Double>
		lengthMap = new HashMap<CoEdge, Double>();
	
	public MetricEdgeLengthAdapter() {
		super(null, CoEdge.class, null, Number.class, true, false);
	}

	public void setLengthMap(Map<CoEdge, Double> lengthMap) {
		this.lengthMap = lengthMap;
	}
	
	@Override
	public Number getEdgeValue(CoEdge e, AdapterSet a) {
		Double d = lengthMap.get(e);
		if (d == null) {
			d = 0.0;
		}
		return d;
	}
	
}
