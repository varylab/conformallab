package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedge.util.HalfEdgeUtils;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility;

public class IndexMedialGraphAdapter extends AbstractAdapter<Double> {
	
	private double digits = 2;
	
	public IndexMedialGraphAdapter() {
		super(Double.class, true, false);
	}

	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getV(V v, AdapterSet a) {
		if(HalfEdgeUtils.isBoundaryVertex(v)) {
			return null;
		}
		return Math.round(digits*QuasiisothermicUtility.alphaRotationFromAdapterSet(v, a) / (2.0 * Math.PI))/digits;
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getF(F f, AdapterSet a) {
		return Math.round(digits*QuasiisothermicUtility.alphaRotationFromAdapterSet(f, a) / (2.0 * Math.PI))/digits;
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Vertex.class.isAssignableFrom(nodeClass) || Face.class.isAssignableFrom(nodeClass);
	}
	
	@Override
	public String toString() {
		return "Index Medial Graph";
	}
	
}