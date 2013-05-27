package de.varylab.discreteconformal.heds.adapter;

import java.util.HashMap;
import java.util.Map;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.varylab.discreteconformal.heds.adapter.types.LengthCrossRatio;

@LengthCrossRatio
public class LengthCrossRatioAdapter extends AbstractAdapter<Double> {

	private Map<? extends Edge<?,?,?>, Double>
		lcrMap = new HashMap<Edge<?,?,?>, Double>(); 
	
	/**
	 * Uses the lengths defined by the Length adapter
	 */
	public LengthCrossRatioAdapter() {
		super(Double.class, true, false);
	}
	
	/**
	 * Uses the length cross ratios defines in the given map. If an edge
	 * is not contained the respective length adapter is used
	 * @param lcrMap
	 */
	public LengthCrossRatioAdapter(Map<? extends Edge<?,?,?>, Double> lcrMap) {
		this();
		this.lcrMap = lcrMap;
	}
	
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getE(E e, AdapterSet a) {
		if (lcrMap.containsKey(e)) {
			return lcrMap.get(e);
		}
		Double l_kj = a.get(Length.class, e.getOppositeEdge().getPreviousEdge(), Double.class);
		Double l_li = a.get(Length.class, e.getPreviousEdge(), Double.class);
		Double l_ik = a.get(Length.class, e.getOppositeEdge().getNextEdge(), Double.class);
		Double l_jl = a.get(Length.class, e.getNextEdge(), Double.class);
		return (l_kj * l_li) / (l_ik * l_jl);
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Edge.class.isAssignableFrom(nodeClass);
	}

}
