package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.TexturePosition;
import de.jtem.halfedgetools.adapter.type.generic.TexturePosition4d;
import de.varylab.discreteconformal.heds.CoVertex;

@TexturePosition
@TexturePosition4d
public class CoDirectTextureAdapter extends AbstractAdapter<double[]> {

	public CoDirectTextureAdapter() {
		super(double[].class, true, false);
	}

	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return CoVertex.class.isAssignableFrom(nodeClass);
	}
	
	@Override
	public double getPriority() {
		return Double.MAX_VALUE;
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double[] getV(V v, AdapterSet a) {
		CoVertex cv = (CoVertex)v;
		return cv.T;
	}
}