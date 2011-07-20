package de.varylab.discreteconformal.adapter;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.jtem.halfedgetools.adapter.AbstractAdapter;
import de.jtem.halfedgetools.adapter.AdapterSet;
import de.jtem.halfedgetools.adapter.type.Length;
import de.jtem.halfedgetools.adapter.type.Weight;

@Weight
public class CotanWeightAdapter extends AbstractAdapter<Double> {

	public CotanWeightAdapter() {
		super(Double.class, true, false);
	}

	@Override
	public <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Double getE(E e, AdapterSet adapters) {
		if (!adapters.contains(Length.class, e.getClass(), Double.class)) {
			throw new RuntimeException(
					"Need adapter for length of edges to calculate cotan weights.");
		}

		double leftcotanalpha = getLeftWeight(e, adapters);
		double rightcotanalpha = getLeftWeight(e.getOppositeEdge(), adapters);

		return .5*(leftcotanalpha+rightcotanalpha);
	}

	private <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> double getLeftWeight(E e, AdapterSet adapters) {
	
		double a = adapters.get(Length.class, e, Double.class);
		double b = adapters
				.get(Length.class, e.getNextEdge(), Double.class);
		double c = adapters.get(Length.class,
				e.getNextEdge().getNextEdge(), Double.class);

		if (!e.getNextEdge().getNextEdge().getNextEdge().equals(e)) {
			throw new RuntimeException("Face is not a triangle.");
		}

		double cosalpha = (b * b + c * c - a * a) / (2 * b * c);
		double cotanalpha = cosalpha / (Math.sqrt(1 - cosalpha * cosalpha));
		return cotanalpha;
	};

	@Override
	public <
		N extends Node<?, ?, ?>
	> boolean canAccept(Class<N> nodeClass) {
		return Edge.class.isAssignableFrom(nodeClass);
	}

}