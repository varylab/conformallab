package de.varylab.discreteconformal.heds.adapter;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Position;
import de.jtem.halfedgetools.adapter.type.TexturePosition;

@TexturePosition
@Position
public class CoTextureDomainPositionAdapter extends AbstractAdapter<double[]> {
	
	public CoTextureDomainPositionAdapter() {
		super(double[].class, true, true);
	}
	
	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return Vertex.class.isAssignableFrom(nodeClass);
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double[] getV(V v, AdapterSet a) {
		a.setPriorityBound(getPriority());
		double[] tp = a.getD(TexturePosition.class, v);
		a.removePriorityBound();
		return tp;
	}
	
	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> void setV(V v, double[] coords, AdapterSet a) {
		a.setPriorityBound(getPriority());
		a.set(TexturePosition.class, v, coords);
		a.removePriorityBound();
	}
	@Override
	public double getPriority() {
		return 1000.0;
	}
	
}