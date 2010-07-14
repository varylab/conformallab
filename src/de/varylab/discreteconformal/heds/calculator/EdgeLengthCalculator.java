package de.varylab.discreteconformal.heds.calculator;

import geom3d.Point;
import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Node;
import de.jtem.halfedge.Vertex;
import de.varylab.discreteconformal.heds.CoEdge;

public class EdgeLengthCalculator extends
		de.varylab.discreteconformal.util.Delaunay.EdgeLengthCalculator {

	@Override
	public <N extends Node<?, ?, ?>> boolean canAccept(Class<N> nodeClass) {
		return nodeClass.isAssignableFrom(CoEdge.class);
	}

	@Override
	public <
		V extends Vertex<V, E, F>, 
		E extends Edge<V, E, F>, 
		F extends Face<V, E, F>
	> double getLength(E e) {
		CoEdge edge = (CoEdge)e;
		Point s = edge.getStartVertex().getPosition();
		Point t = edge.getTargetVertex().getPosition();
		return s.distanceTo(t);
	}

}
