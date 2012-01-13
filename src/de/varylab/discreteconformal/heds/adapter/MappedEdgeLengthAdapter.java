package de.varylab.discreteconformal.heds.adapter;

import static de.jreality.math.Pn.EUCLIDEAN;

import java.util.HashMap;
import java.util.Map;

import de.jreality.math.Pn;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.varylab.discreteconformal.heds.CoEdge;
import de.varylab.discreteconformal.heds.CoVertex;

@Length
public class MappedEdgeLengthAdapter extends AbstractAdapter<Double> {

	private Map<CoEdge, Double>
		lengthMap = new HashMap<CoEdge, Double>();
	private double 
		priority = 0;
	
	public MappedEdgeLengthAdapter(double priority) {
		super(Double.class, true, true);
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
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getE(E e, AdapterSet a) {
		if (lengthMap.containsKey(e)) {
			return lengthMap.get(e);
		} else {
			CoVertex s = (CoVertex)e.getStartVertex();
			CoVertex t = (CoVertex)e.getTargetVertex();
			return Pn.distanceBetween(s.P, t.P, EUCLIDEAN);
		}
	}	
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void setE(E e, Double value, AdapterSet a) {
		lengthMap.put((CoEdge)e, value);
	}
	
	@Override
	public double getPriority() {
		return priority;
	}

}
