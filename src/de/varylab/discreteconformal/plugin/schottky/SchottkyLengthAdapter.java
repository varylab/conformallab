package de.varylab.discreteconformal.plugin.schottky;

import java.util.Map;

import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Length
public class SchottkyLengthAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private Map<CoEdge, Double>
		lMap = null;
	
	public SchottkyLengthAdapter(Map<CoEdge, Double> map) {
		super(null, CoEdge.class, null, Double.class, true, false);
		lMap = map;
	}
	
	@Override
	public Double getEdgeValue(CoEdge e, AdapterSet a) {
		return lMap.get(e);
	}	
	
	@Override
	public double getPriority() {
		return 10;
	}
	
}