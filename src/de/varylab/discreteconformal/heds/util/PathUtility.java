package de.varylab.discreteconformal.heds.util;

import java.util.HashSet;
import java.util.Set;

import de.jtem.halfedge.Edge;
import de.jtem.halfedge.Face;
import de.jtem.halfedge.Vertex;

public class PathUtility {

	public static <
		V extends Vertex<V, E, F>,
		E extends Edge<V, E, F>,
		F extends Face<V, E, F>
	> Set<V> getVerticesOnPath(Set<E> path) {
		HashSet<V> result = new HashSet<V>();
		for (E e : path) {
			result.add(e.getStartVertex());
			result.add(e.getTargetVertex());
		}
		return result;
	}
	
	
}
