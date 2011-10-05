package de.varylab.discreteconformal.adapter;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Weight;

@Weight
public class MappedWeightAdapter extends AbstractAdapter<Double> {

	private Map<Edge<?, ?, ?>, Double>
		wMap = new HashMap<Edge<?,?,?>, Double>();
	
	public MappedWeightAdapter() {
		super(Double.class, true, true);
	}
	public MappedWeightAdapter(Map<Edge<?,?,?>, Double> map) {
		super(Double.class, true, true);
		this.wMap = map;
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getE(E e, AdapterSet a) {
		if (wMap.containsKey(e)) {
			return wMap.get(e);
		} else {
			return 0.0;
		}
	}	
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void setE(E e, Double value, AdapterSet a) {
		wMap.put(e, value);
		wMap.put(e.getOppositeEdge(), value);
	}
	
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Edge.class.isAssignableFrom(nodeClass);
	}
	
	@Override
	public double getPriority() {
		return 10;
	}
	
}