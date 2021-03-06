package de.varylab.discreteconformal.adapter;

import static de.varylab.discreteconformal.unwrapper.quasiisothermic.QuasiisothermicUtility.alphaOrientationFromAdapterSet;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;

public class TriangleOrientationAdapter extends AbstractAdapter<Double> {
	
	public TriangleOrientationAdapter() {
		super(Double.class, true, false);
	}

	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getF(F f, AdapterSet a) {
		try {
			return (Double)alphaOrientationFromAdapterSet(f, a);
		} catch (Exception e) {
			return -1.0;
		}
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Face.class.isAssignableFrom(nodeClass);
	}
	
	@Override
	public String toString() {
		return "Flipped Triangles";
	}
	
}