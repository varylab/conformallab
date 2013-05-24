package de.varylab.discreteconformal.heds.adapter;

import static de.jreality.math.Pn.EUCLIDEAN;

import java.util.HashMap;
import java.util.Map;

import de.jreality.math.Pn;
import de.jtem.halfedge.Node;
import de.jtem.halfedgetools.adapter.AbstractTypedAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoFace;
import de.varylab.discreteconformal.heds.CoVertex;

@Length
public class MappedEdgeLengthAdapter extends AbstractTypedAdapter<CoVertex, CoEdge, CoFace, Double> {

	private Map<CoEdge, Double>
		lengthMap = new HashMap<CoEdge, Double>();
	private double 
		priority = 0;
	
	public MappedEdgeLengthAdapter(double priority) {
		super(null, CoEdge.class, null, Double.class, true, true);
		this.priority = priority;
	}
	
	public MappedEdgeLengthAdapter(Map<CoEdge, Double> lengthMap, double priority) {
		this(priority);
		this.lengthMap = lengthMap;
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return CoEdge.class.isAssignableFrom(nodeClass);
	}
	
	@Override
	public Double getEdgeValue(CoEdge e, AdapterSet a) {
		if (lengthMap.containsKey(e)) {
			return lengthMap.get(e);
		} else {
			CoVertex s = e.getStartVertex();
			CoVertex t = e.getTargetVertex();
			return Pn.distanceBetween(s.P, t.P, EUCLIDEAN);
		}
	}
	
	@Override
	public void setEdgeValue(CoEdge e, Double value, AdapterSet a) {
		lengthMap.put(e, value);
	}
	
	
	@Override
	public double getPriority() {
		return priority;
	}

}
